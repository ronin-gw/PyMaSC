"""Unified calculation handler using strategy pattern.

This module provides a unified handler that replaces both CCCalcHandler
and BACalcHandler using the strategy pattern. It delegates algorithm-specific
logic to strategies while maintaining a single, consistent interface.

Key components:
- UnifiedCalcHandler: Main handler that uses strategies for algorithm selection
- Integrates with factory pattern for calculator and worker creation
- Maintains backward compatibility while simplifying the architecture
"""
import logging
from multiprocessing import Queue, Lock
from pathlib import Path
from typing import Optional, List, Dict, Any, Union

import pysam

from PyMaSC.core.strategy import CalculationContext
from PyMaSC.core.factory import CalculatorFactory, WorkerFactory
from PyMaSC.core.configuration import ConfigurationService
from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig,
    AlgorithmType, ExecutionMode, WorkerConfig
)
from PyMaSC.core.interfaces import CrossCorrelationCalculator

from PyMaSC.utils.progress import ProgressBar, ProgressHook, MultiLineProgressManager
from PyMaSC.utils.calc import filter_chroms, exec_worker_pool
from PyMaSC.core.readlen import estimate_readlen
from PyMaSC.core.progress_adapter import ProgressBarAdapter, ProgressManager, get_progress_manager

logger = logging.getLogger(__name__)


class InputUnseekable(Exception):
    """Exception raised when input stream cannot be rewound."""
    pass


class NothingToCalc(Exception):
    """Exception raised when no chromosomes match filtering criteria."""
    pass


class UnifiedCalcHandler:
    """Unified cross-correlation calculation handler.
    
    This handler replaces both CCCalcHandler and BACalcHandler by using
    the strategy pattern for algorithm-specific logic. It maintains the
    same public interface as the original handlers while delegating
    algorithm selection to strategies.
    
    The handler coordinates:
    - Algorithm strategy selection
    - BAM file processing and validation
    - Calculator and worker creation via factories
    - Result collection and aggregation
    - Progress reporting
    
    Attributes:
        path: Input BAM file path
        config: Calculation configuration
        execution_config: Execution mode configuration
        mappability_config: Optional mappability configuration
        references: List of target chromosome names
        lengths: List of chromosome lengths
        read_len: Estimated or specified read length
        calc_context: Calculation context with current strategy
    """
    
    def __init__(self, 
                 path: Union[str, Path],
                 config: CalculationConfig,
                 execution_config: Optional[ExecutionConfig] = None,
                 mappability_config: Optional[MappabilityConfig] = None):
        """Initialize unified handler with configuration.
        
        Args:
            path: Path to input BAM file
            config: Calculation configuration
            execution_config: Execution configuration (defaults to single process)
            mappability_config: Optional mappability configuration
            
        Raises:
            ValueError: If BAM file has no sequences defined
            NothingToCalc: If no chromosomes match filtering criteria
        """
        self.path = str(path)
        self.config = config
        self.execution_config = execution_config or ExecutionConfig()
        self.mappability_config = mappability_config
        
        # Initialize calculation context with appropriate strategy
        self.calc_context = CalculationContext.create_from_config(config)
        
        # Open and validate BAM file
        try:
            self.align_file = pysam.AlignmentFile(self.path)
        except ValueError:
            logger.error("File has no sequences defined.")
            raise
        
        # Filter chromosomes if specified
        chromfilter = getattr(config, 'chromfilter', None)
        target_references = filter_chroms(self.align_file.references, chromfilter)
        if not target_references:
            logger.error("There is no targeted chromosomes.")
            raise NothingToCalc
        
        # Build reference lists
        self.references = []
        self.lengths = []
        for reference, length in zip(self.align_file.references, self.align_file.lengths):
            if reference in target_references:
                self.references.append(reference)
                self.lengths.append(length)
        
        # Update config with filtered references
        self.config.references = self.references
        self.config.lengths = self.lengths
        
        # Check for index in multiprocess mode
        if (self.execution_config.mode == ExecutionMode.MULTI_PROCESS and 
            not self.align_file.has_index()):
            logger.error("Need indexed alignment file for multi-processing. "
                        "Calculation will be executed by single process.")
            self.execution_config.mode = ExecutionMode.SINGLE_PROCESS
            self.execution_config.worker_count = 1
        elif self.execution_config.mode == ExecutionMode.MULTI_PROCESS:
            self.align_file.close()
        
        # Initialize result storage
        self.ref2forward_sum = {}
        self.ref2reverse_sum = {}
        self.ref2ccbins = {}
        
        self.mappable_ref2forward_sum = {}
        self.mappable_ref2reverse_sum = {}
        self.mappable_ref2ccbins = {}
        self.ref2mappable_len = {}
        
        # Read length will be set later
        self.read_len = None
        
        # Mappability handler reference
        self.mappability_handler = None
        
        # Legacy attributes for backward compatibility
        self._esttype = getattr(config, 'esttype', 'mean')
        
        # Observer pattern support
        self.use_observer_progress = True
        self._progress_observers = []
        self._progress_manager = None
    
    def set_readlen(self, readlen: Optional[int] = None):
        """Set or estimate read length.
        
        Args:
            readlen: Explicit read length to use, or None to estimate
            
        Raises:
            InputUnseekable: If estimation needed but input is unseekable
            ValueError: If read length exceeds maximum shift distance
        """
        if readlen:
            self.read_len = readlen
        elif self.path == '-':
            logger.error("Cannot execute read length checking for unseekable input.")
            raise InputUnseekable
        else:
            logger.info(f"Check read length... : {self.path}")
            # Use config's estimation type if available
            esttype = getattr(self.config, 'esttype', 'mean')
            self.read_len = estimate_readlen(
                self.path, esttype, self.config.mapq_criteria
            )
        
        if self.read_len > self.config.max_shift:
            logger.error(f"Read length ({self.read_len}) seems to be longer than "
                        f"shift size ({self.config.max_shift}).")
            raise ValueError
        
        # Update config with read length
        self.config.read_length = self.read_len
    
    def set_mappability_handler(self, mappability_handler: Any):
        """Set mappability handler and validate chromosome sizes.
        
        Args:
            mappability_handler: MappabilityHandler instance
        """
        self.mappability_handler = mappability_handler
        bw_chromsizes = self.mappability_handler.chromsizes
        
        for i, reference in enumerate(self.references):
            if reference not in bw_chromsizes:
                logger.debug(f"mappability for '{reference}' not found")
                continue
            self._compare_refsize(reference, bw_chromsizes[reference])
    
    def _compare_refsize(self, reference: str, bw_chr_size: int):
        """Compare and validate chromosome sizes between BAM and BigWig."""
        i = self.references.index(reference)
        bam_chr_size = self.lengths[i]
        if bw_chr_size != bam_chr_size:
            logger.warning(f"'{reference}' reference length mismatch: "
                          f"SAM/BAM -> {bam_chr_size:,}, BigWig -> {bw_chr_size:,}")
            if bam_chr_size < bw_chr_size:
                logger.warning(f"Use longer length '{bw_chr_size:d}' for "
                              f"'{reference}' anyway")
                self.lengths[i] = bw_chr_size
                self.config.lengths[i] = bw_chr_size
    
    def run_calcuration(self):
        """Execute cross-correlation calculation workflow.
        
        Note: Method name kept as 'calcuration' for backward compatibility
        """
        if self.execution_config.mode == ExecutionMode.MULTI_PROCESS:
            # Disable single-line progress for multiprocess
            ProgressBar.global_switch = False
            ProgressHook.global_switch = True
            self._run_multiprocess_calculation()
        else:
            self._run_singleprocess_calculation()
    
    def _run_singleprocess_calculation(self):
        """Execute calculation in single process mode."""
        # Create calculator using strategy
        calculator = self.calc_context.create_calculator(
            self.config,
            self.mappability_config
        )
        
        # Create progress bar with observer support if enabled
        use_observer = getattr(self, 'use_observer_progress', True)
        if use_observer:
            progress = ProgressBarAdapter()
            # Get or create progress manager
            progress_manager = getattr(self, '_progress_manager', None) or get_progress_manager()
            # Attach any configured observers
            if hasattr(self, '_progress_observers'):
                for observer in self._progress_observers:
                    progress._subject.attach(observer)
        else:
            progress = ProgressBar()
        
        ref2genomelen = dict(zip(self.references, self.lengths))
        current_chr = None
        
        # Reopen file if needed
        if not hasattr(self, 'align_file') or self.align_file.closed:
            self.align_file = pysam.AlignmentFile(self.path)
        
        for read in self.align_file:
            chrom = read.reference_name
            if chrom is None or chrom not in self.references:
                continue
            
            # Update progress for new chromosome
            if chrom != current_chr:
                if current_chr is not None:
                    progress.clean()
                current_chr = chrom
                progress.set(chrom, ref2genomelen[chrom])
            
            progress.update(read.reference_start)
            
            # Skip reads based on criteria
            if (read.is_read2 or read.mapping_quality < self.config.mapq_criteria or
                read.is_unmapped or read.is_duplicate):
                continue
            
            # Feed read to calculator
            if read.is_reverse:
                calculator.feed_reverse_read(
                    chrom, read.reference_start + 1, read.infer_query_length()
                )
            else:
                calculator.feed_forward_read(
                    chrom, read.reference_start + 1, read.infer_query_length()
                )
        
        progress.clean()
        self.align_file.close()
        
        # Finalize calculation
        calculator.finishup_calculation()
        
        # Collect results
        self._collect_calculator_results(calculator)
        self._calc_unsolved_mappability()
    
    def _run_multiprocess_calculation(self):
        """Execute calculation using multiple processes."""
        self._order_queue = Queue()
        self._report_queue = Queue()
        self._logger_lock = Lock()
        
        # Create worker configuration
        worker_config = WorkerConfig(
            calculation_config=self.config,
            mappability_config=self.mappability_config,
            progress_enabled=True,
            logger_lock=self._logger_lock
        )
        
        # Create workers using factory
        workers = []
        num_workers = min(self.execution_config.worker_count, len(self.references))
        for _ in range(num_workers):
            worker = WorkerFactory.create_worker(
                worker_config,
                self._order_queue,
                self._report_queue,
                self._logger_lock,
                self.path
            )
            workers.append(worker)
        
        # Run calculation with workers
        self._run_calculation_with_workers(workers)
        self._calc_unsolved_mappability()
    
    def _run_calculation_with_workers(self, workers: List[Any]):
        """Coordinate worker processes and collect results."""
        _chrom2finished = {c: False for c in self.references}
        progress = MultiLineProgressManager()
        
        with exec_worker_pool(workers, self.references, self._order_queue):
            while True:
                chrom, obj = self._report_queue.get()
                if chrom is None:  # Progress update
                    chrom, body = obj
                    with self._logger_lock:
                        progress.update(chrom, body)
                else:
                    # Result received
                    mappable_len, cc_stats, masc_stats = obj
                    self._receive_results(chrom, mappable_len, cc_stats, masc_stats)
                    
                    _chrom2finished[chrom] = True
                    if all(_chrom2finished.values()):
                        break
                    
                    with self._logger_lock:
                        progress.erase(chrom)
        
        progress.clean()
    
    def _receive_results(self, chrom: str, mappable_len: Optional[int],
                        cc_stats: tuple, masc_stats: tuple):
        """Process and store results from worker processes."""
        f_sum, r_sum, ccbins = cc_stats
        mf_sum, mr_sum, mccbins = masc_stats
        
        if mappable_len is not None and self.mappability_handler:
            self.mappability_handler.chrom2mappable_len[chrom] = mappable_len
            self.mappability_handler.chrom2is_called[chrom] = True
        
        if None not in (f_sum, r_sum, ccbins):
            self.ref2forward_sum[chrom] = f_sum
            self.ref2reverse_sum[chrom] = r_sum
            self.ref2ccbins[chrom] = ccbins
        
        if None not in (mappable_len, mf_sum, mr_sum, mccbins):
            self.mappable_ref2forward_sum[chrom] = mf_sum
            self.mappable_ref2reverse_sum[chrom] = mr_sum
            self.mappable_ref2ccbins[chrom] = mccbins
    
    def _collect_calculator_results(self, calculator: CrossCorrelationCalculator):
        """Collect results from calculator instance."""
        # Basic NCC results
        if not self.config.skip_ncc:
            self.ref2forward_sum = calculator.ref2forward_sum
            self.ref2reverse_sum = calculator.ref2reverse_sum
            self.ref2ccbins = calculator.ref2ccbins
        
        # MSCC/BitArray results if available
        if hasattr(calculator, 'ref2mappable_forward_sum'):
            adapter = calculator  # Should be CalculatorAdapter
            if adapter.ref2mappable_forward_sum:
                self.mappable_ref2forward_sum = adapter.ref2mappable_forward_sum
                self.mappable_ref2reverse_sum = adapter.ref2mappable_reverse_sum
                self.mappable_ref2ccbins = adapter.ref2mascbins
        
        # Mappable length handling for BitArray
        if hasattr(calculator, '_calculator') and hasattr(calculator._calculator, 'ref2mappable_len'):
            mappable_len = {k: v for k, v in calculator._calculator.ref2mappable_len.items() 
                           if v is not None}
            if mappable_len and self.mappability_handler:
                self.ref2mappable_len = mappable_len
                self.mappability_handler.chrom2mappable_len = mappable_len
                self.mappability_handler.mappable_len = list(map(sum, zip(*mappable_len.values())))
    
    def _calc_unsolved_mappability(self):
        """Complete any remaining mappability calculations."""
        if self.mappability_handler is not None:
            if not self.mappability_handler.is_called:
                self.mappability_handler.is_called = all(
                    self.mappability_handler.chrom2is_called.values()
                )
                self.mappability_handler.calc_mappability()
            self.ref2mappable_len = self.mappability_handler.chrom2mappable_len
    
    # Properties for backward compatibility
    @property
    def skip_ncc(self) -> bool:
        """Whether NCC calculation is skipped."""
        return self.config.skip_ncc
    
    @property
    def max_shift(self) -> int:
        """Maximum shift distance."""
        return self.config.max_shift
    
    @property
    def mapq_criteria(self) -> int:
        """Minimum mapping quality."""
        return self.config.mapq_criteria
    
    @property
    def nworker(self) -> int:
        """Number of worker processes."""
        return self.execution_config.worker_count
    
    @property
    def esttype(self) -> str:
        """Read length estimation type."""
        return self._esttype
    
    @esttype.setter
    def esttype(self, value: str):
        """Set read length estimation type."""
        self._esttype = value
    
    # Progress observer methods
    def attach_progress_observer(self, observer: 'ProgressObserver') -> None:
        """Attach a progress observer.
        
        Args:
            observer: Progress observer to attach
        """
        if observer not in self._progress_observers:
            self._progress_observers.append(observer)
    
    def detach_progress_observer(self, observer: 'ProgressObserver') -> None:
        """Detach a progress observer.
        
        Args:
            observer: Progress observer to detach
        """
        if observer in self._progress_observers:
            self._progress_observers.remove(observer)
    
    def set_progress_manager(self, manager: ProgressManager) -> None:
        """Set a custom progress manager.
        
        Args:
            manager: Progress manager to use
        """
        self._progress_manager = manager
    
    def enable_observer_progress(self, enabled: bool = True) -> None:
        """Enable or disable observer-based progress reporting.
        
        Args:
            enabled: Whether to use observer pattern for progress
        """
        self.use_observer_progress = enabled