"""Unified calculation handler using factory pattern.

This module provides a handler that replaces both CCCalcHandler
and BACalcHandler by directly using the factory pattern for calculator creation.
This eliminates unnecessary abstraction layers while maintaining a single,
consistent interface.

Key components:
- CalcHandler: Main handler that uses factory for calculator creation
- Direct integration with CalculatorFactory for optimal performance
- Maintains backward compatibility while simplifying the architecture
"""
from __future__ import annotations

import logging
import os
from multiprocessing import Queue, Lock
from typing import Optional, List, Tuple, Union, Any

import pysam

from PyMaSC.core.exceptions import InputUnseekable, NothingToCalc
from PyMaSC.core.factory import CalculatorFactory
from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig,
    ExecutionMode, WorkerConfig
)
from PyMaSC.core.bam_processor import BAMFileProcessor, BAMValidationError
from PyMaSC.core.result import (
    ChromResult, GenomeWideResult
)
from PyMaSC.core.result import aggregate_results
from PyMaSC.core.readlen import estimate_readlen
from PyMaSC.core.worker import CalcWorker
from PyMaSC.utils.calc import exec_worker_pool
from PyMaSC.utils.progress import ProgressBar, ProgressHook, MultiLineProgressManager
from PyMaSC.handler.read import create_read_processor

logger = logging.getLogger(__name__)


class CalcHandler:
    """Cross-correlation calculation handler.

    This handler replaces both CCCalcHandler and BACalcHandler by using
    the factory pattern for direct calculator creation. It maintains the
    same public interface as the original handlers while providing
    optimized performance through reduced abstraction layers.

    The handler coordinates:
    - Direct calculator creation via CalculatorFactory
    - BAM file processing and validation (via BAMFileProcessor)
    - Calculator and worker creation via factories
    - Result collection and aggregation
    - Progress reporting (via ProgressCoordinator)

    Attributes:
        path: Input BAM file path
        config: Calculation configuration
        execution_config: Execution mode configuration
        mappability_config: Optional mappability configuration
        references: List of target chromosome names
        lengths: List of chromosome lengths
        read_len: Estimated or specified read length
        bam_processor: BAM file processor for file operations
        bam_metadata: BAM file metadata
    """

    def __init__(self,
                 path: os.PathLike[str],
                 config: CalculationConfig,
                 execution_config: ExecutionConfig,
                 mappability_config: Optional[MappabilityConfig] = None):
        """Initialize handler with configuration.

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
        self.execution_config = execution_config
        self.mappability_config = mappability_config

        # Initialize BAM file processor
        self.bam_processor = BAMFileProcessor(self.path)

        # Validate BAM file and extract metadata
        try:
            chromfilter = getattr(config, 'chromfilter', None)
            self.bam_metadata = self.bam_processor.validate_and_open(chromfilter)
        except BAMValidationError as e:
            if "no sequences defined" in str(e):
                logger.error("File has no sequences defined.")
                raise ValueError("File has no sequences defined.")
            elif "No chromosomes match" in str(e):
                logger.error("There is no targeted chromosomes.")
                raise NothingToCalc
            elif "no reference sequences" in str(e):
                logger.error("There is no targeted chromosomes.")
                raise NothingToCalc
            else:
                raise

        # Extract filtered references and lengths
        self.references = self.bam_metadata.filtered_references
        self.lengths = self.bam_metadata.filtered_lengths

        # Update config with filtered references
        self.config.references = self.references
        self.config.lengths = self.lengths

        # Check for index in multiprocess mode
        if (self.execution_config.mode is ExecutionMode.MULTI_PROCESS and
                not self.bam_processor.check_multiprocess_compatibility()):
            self.execution_config.mode = ExecutionMode.SINGLE_PROCESS
            self.execution_config.worker_count = 1
        elif self.execution_config.mode is ExecutionMode.MULTI_PROCESS:
            self.bam_processor.close()

        # For backward compatibility, maintain align_file access
        self.align_file: Optional[pysam.AlignmentFile] = None

        # Read length will be set later
        self.read_len: Optional[int] = None

        # Mappability handler reference
        self.mappability_handler: Any = None

        # Legacy attributes for backward compatibility
        self._esttype = getattr(config, 'esttype', 'mean')

    def set_or_estimate_readlen(self, readlen: Optional[int] = None) -> None:
        """Set or estimate read length.

        Args:
            readlen: Explicit read length to use, or None to estimate

        Raises:
            InputUnseekable: If estimation needed but input is unseekable
            ValueError: If read length exceeds maximum shift distance
        """
        if readlen is not None:
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

        if self.read_len is not None and self.read_len > self.config.max_shift:
            logger.error(f"Read length ({self.read_len}) seems to be longer than "
                         f"shift size ({self.config.max_shift}).")
            raise ValueError

        # Update config with read length
        self.config.read_length = self.read_len

    def set_mappability_handler(self, mappability_handler: Any) -> None:
        """Set mappability handler and validate chromosome sizes.

        Args:
            mappability_handler: MappabilityHandler instance
        """
        self.mappability_handler = mappability_handler
        if self.mappability_handler is not None:
            bw_chromsizes = self.mappability_handler.chromsizes
            # Use BAM processor to validate and reconcile chromosome sizes
            updated_sizes = self.bam_processor.validate_chromosome_sizes(bw_chromsizes)

            # Update lengths based on validation results
            for i, reference in enumerate(self.references):
                if reference in updated_sizes:
                    new_length = updated_sizes[reference]
                    if new_length != self.lengths[i]:
                        self.lengths[i] = new_length
                        self.config.lengths[i] = new_length

    def run_calculation(self) -> GenomeWideResult:
        """Execute cross-correlation calculation workflow."""
        if self.execution_config.mode is ExecutionMode.MULTI_PROCESS:
            # Configure progress for multiprocess mode
            ProgressBar.global_switch = False
            ProgressHook.global_switch = True
            return self._run_multiprocess_calculation()
        else:
            # Single process mode uses default progress settings
            return self._run_singleprocess_calculation()

    def _run_singleprocess_calculation(self) -> GenomeWideResult:
        """Execute calculation in single process mode."""
        # Create calculator using factory
        calculator = CalculatorFactory.create_calculator(
            self.config.target,
            self.config.implementation,
            self.config,
            self.mappability_config
        )

        #
        read_processor = create_read_processor(calculator, self.config.mapq_criteria)
        progress_bar = ProgressBar()

        # Setup progress tracking for single process
        ref2genomelen = dict(zip(self.references, self.lengths))

        # Reopen file if needed
        self.align_file = self.bam_processor.reopen()

        processed_bases = 0
        current_chrom = None

        for read in self.align_file:
            chrom = read.reference_name

            if chrom is None or chrom not in self.references:
                continue

            processed_bases = max(processed_bases, read.reference_start)
            if chrom != current_chrom:
                current_chrom = chrom
                progress_bar.set(chrom, ref2genomelen[chrom])
                progress_bar.clean()

            read_processor.process_read(read)
            progress_bar.update(read.reference_start)

        # progress_bar.clean()
        self.bam_processor.close()

        # Finalize calculation
        calculator.finishup_calculation()
        self._calc_unsolved_mappability()
        return calculator.get_whole_result()

    def _run_multiprocess_calculation(self) -> GenomeWideResult:
        """Execute calculation using multiple processes."""

        # Queue for distributing tasks to workers
        self._order_queue: Queue = Queue()
        # IMPORTANT: `self._report_queue` receives both results and progress updates.
        # - results are collected for aggregation
        # If the 1st item is a string, it shoud be a chromosome and the 2nd item must be a `WorkerResult`.
        # Otherwise, it's a progress report from a `ProgressHook` and that should be avoid here.
        self._report_queue: Queue[Tuple[Optional[str], Union[ChromResult, Tuple[str, int]]]] = Queue()
        self._logger_lock = Lock()

        # Create worker configuration
        # Handle None mappability_config by providing a default or skipping
        if self.mappability_config is not None:
            worker_config = WorkerConfig(
                calculation_config=self.config,
                mappability_config=self.mappability_config,
                progress_enabled=True,
                logger_lock=self._logger_lock
            )
        else:
            # Create a minimal config when mappability is not needed
            default_mappability = MappabilityConfig()  # Use default config
            worker_config = WorkerConfig(
                calculation_config=self.config,
                mappability_config=default_mappability,
                progress_enabled=True,
                logger_lock=self._logger_lock
            )

        # Create workers using factory
        workers = []
        num_workers = min(self.execution_config.worker_count, len(self.references))
        for _ in range(num_workers):
            worker = CalcWorker(
                worker_config,
                self._order_queue,
                self._report_queue,
                self._logger_lock,
                self.path
            )
            workers.append(worker)

        # Run calculation with workers
        result = self._run_calculation_with_workers(workers)
        self._calc_unsolved_mappability()

        return result

    def _run_calculation_with_workers(self, workers: List[Any]) -> GenomeWideResult:
        """Coordinate worker processes and collect results."""
        _chrom2finished = {c: False for c in self.references}
        _worker_results = {}  # Collect results for new aggregation system
        progress = MultiLineProgressManager()

        with exec_worker_pool(workers, self.references, self._order_queue):
            while True:
                chrom, obj = self._report_queue.get()

                #
                if chrom == '__ERROR__':
                    raise RuntimeError(f"Worker error:\n{obj}")

                elif chrom is None:
                    # This is a progress update from ProgressHook
                    assert not isinstance(obj, ChromResult)
                    chrom_text, body = obj  # 'body' is the progress bar string
                    with self._logger_lock:
                        progress.update(chrom_text, body)
                    continue

                # Result received
                if obj is not None:
                    assert isinstance(obj, ChromResult)
                    # Store for aggregation system
                    _worker_results[chrom] = obj

                    _chrom2finished[chrom] = True

                    # Erase progress for completed chromosome (only when result is received)
                    with self._logger_lock:
                        progress.erase(chrom)

                    if all(_chrom2finished.values()):
                        break

        # Clean up progress display
        progress.clean()

        # Apply new aggregation system for worker results
        return aggregate_results(_worker_results)

    def _calc_unsolved_mappability(self) -> None:
        """Complete any remaining mappability calculations."""
        if self.mappability_handler is not None:
            if not self.mappability_handler.is_called:
                self.mappability_handler.is_called = all(
                    self.mappability_handler.chrom2is_called.values()
                )
                self.mappability_handler.calc_mappability()
