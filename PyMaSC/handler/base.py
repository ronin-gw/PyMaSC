"""Base handler class with common calculation workflow functionality.

This module provides the BaseCalcHandler class that extracts common functionality
from all calculation handlers, eliminating duplication and providing a consistent
interface for calculation workflows.

Key features:
- Common BAM file initialization and validation
- Unified chromosome filtering and reference setup
- Standardized result aggregation patterns
- Shared progress management logic
- Common multiprocessing workflow coordination
"""
import logging
from abc import ABC, abstractmethod
from multiprocessing import Queue, Lock
from typing import Dict, List, Tuple, Optional, Any

from pysam import AlignmentFile

from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.progress import ProgressBar
from PyMaSC.utils.calc import filter_chroms

logger = logging.getLogger(__name__)


class NothingToCalc(Exception):
    """Exception raised when no chromosomes are available for calculation."""
    pass


class BaseCalcHandler(ABC):
    """Base class for all calculation handlers.
    
    Provides common functionality for PyMaSC calculation handlers including
    BAM file handling, chromosome filtering, multiprocessing coordination,
    and result aggregation. Subclasses implement algorithm-specific logic.
    
    This eliminates 85-105 lines of duplicate code across handlers by
    centralizing common initialization and workflow patterns.
    
    Attributes:
        path: Path to BAM file
        mapq_criteria: Minimum mapping quality threshold
        max_shift: Maximum shift distance for correlation
        references: List of chromosome names
        lengths: List of chromosome lengths
        mappability_handler: Optional mappability handler
        nworker: Number of worker processes
        skip_ncc: Whether to skip NCC calculation
    """
    
    def __init__(self, 
                 path: str,
                 mapq_criteria: int,
                 max_shift: int,
                 nworker: int = 1,
                 mappability_handler: Optional[MappabilityHandler] = None,
                 skip_ncc: bool = False,
                 **kwargs):
        """Initialize base calculation handler.
        
        Args:
            path: Path to BAM file
            mapq_criteria: Minimum mapping quality threshold
            max_shift: Maximum shift distance
            nworker: Number of worker processes
            mappability_handler: Optional mappability handler
            skip_ncc: Whether to skip NCC calculation
            **kwargs: Additional algorithm-specific parameters
        """
        self.path = path
        self.mapq_criteria = mapq_criteria
        self.max_shift = max_shift
        self.nworker = nworker
        self.mappability_handler = mappability_handler
        self.skip_ncc = skip_ncc
        
        # Initialize common attributes
        self._order_queue = None
        self._report_queue = None
        self._logger_lock = None
        
        # Result storage
        self.ref2forward_sum: Dict[str, int] = {}
        self.ref2reverse_sum: Dict[str, int] = {}
        self.ref2ccbins: Dict[str, Any] = {}
        self.ref2mappable_len: Dict[str, int] = {}
        self.mappable_ref2forward_sum: Dict[str, int] = {}
        self.mappable_ref2reverse_sum: Dict[str, int] = {}
        self.mappable_ref2ccbins: Dict[str, Any] = {}
        
        # Initialize BAM file and references
        self._initialize_bam_file()
        self._initialize_references(**kwargs)
    
    def _initialize_bam_file(self) -> None:
        """Initialize and validate BAM file.
        
        Common BAM file initialization logic extracted from all handlers.
        Handles file opening, sequence validation, and index checking.
        
        Raises:
            RuntimeError: If BAM file has no sequences defined
            IOError: If BAM file cannot be opened or is not indexed for multiprocessing
        """
        # Open BAM file and validate
        self.align_file = AlignmentFile(self.path)
        
        # Check for sequences
        if not self.align_file.references:
            self.align_file.close()
            logger.error("File has no sequences defined.")
            raise RuntimeError("BAM file has no sequences defined")
        
        # Check index requirement for multiprocessing
        if self.nworker > 1:
            try:
                # Test if file is indexed by attempting to fetch from first reference
                list(self.align_file.fetch(self.align_file.references[0], 1, 2))
            except (OSError, ValueError) as e:
                self.align_file.close()
                logger.error(f"BAM file must be indexed for multiprocessing: {e}")
                raise IOError("BAM file must be indexed for multiprocessing") from e
    
    def _initialize_references(self, **kwargs) -> None:
        """Initialize chromosome references and filtering.
        
        Common reference initialization logic extracted from all handlers.
        Applies chromosome filtering and validates that chromosomes exist.
        
        Args:
            **kwargs: May contain 'chroms' for chromosome filtering
            
        Raises:
            NothingToCalc: If no chromosomes remain after filtering
        """
        # Get all references from BAM file
        all_references = list(self.align_file.references)
        all_lengths = list(self.align_file.lengths)
        
        # Apply chromosome filtering if specified
        chroms = kwargs.get('chroms')
        if chroms:
            filtered_refs, filtered_lengths = filter_chroms(
                all_references, all_lengths, chroms
            )
        else:
            filtered_refs, filtered_lengths = all_references, all_lengths
        
        # Validate that we have chromosomes to process
        if not filtered_refs:
            self.align_file.close()
            raise NothingToCalc("No chromosomes available for calculation after filtering")
        
        self.references = filtered_refs
        self.lengths = filtered_lengths
        
        logger.info(f"Processing {len(self.references)} chromosomes: {', '.join(self.references)}")
    
    def _setup_progress_management(self) -> None:
        """Setup progress bar management for calculations.
        
        Common progress setup logic extracted from all handlers.
        Configures global progress bar behavior for single vs multiprocess execution.
        """
        if self.nworker == 1:
            # Single process - enable progress bar
            ProgressBar.global_switch = True
        else:
            # Multiprocess - disable progress bar to avoid conflicts
            ProgressBar.global_switch = False
    
    def _create_worker_queues(self) -> None:
        """Create multiprocessing queues for worker communication.
        
        Common queue creation logic extracted from all handlers.
        """
        self._order_queue = Queue()
        self._report_queue = Queue()
        self._logger_lock = Lock()
    
    def _distribute_work(self, workers: List[Any]) -> None:
        """Distribute chromosomes to workers and collect results.
        
        Common multiprocessing workflow extracted from all handlers.
        Sends work orders and collects results from worker processes.
        
        Args:
            workers: List of worker process instances
        """
        # Start all workers
        for worker in workers:
            worker.start()
        
        # Send work orders
        for chrom in self.references:
            self._order_queue.put(chrom)
        
        # Send termination signals
        for _ in workers:
            self._order_queue.put(None)
        
        # Collect results
        self._receive_results(len(self.references))
        
        # Wait for workers to complete
        for worker in workers:
            worker.join()
    
    def _receive_results(self, expected_results: int) -> None:
        """Receive and aggregate results from worker processes.
        
        Common result aggregation logic extracted from all handlers.
        Handles the standard result format and stores in appropriate dictionaries.
        
        Args:
            expected_results: Number of results to expect
        """
        for _ in range(expected_results):
            chrom, (mappable_len, ncc_stats, mscc_stats) = self._report_queue.get()
            
            # Store mappable length if available
            if mappable_len is not None:
                self.ref2mappable_len[chrom] = mappable_len
            
            # Store NCC results if available
            if ncc_stats[0] is not None:  # Check if NCC data exists
                f_sum, r_sum, ccbins = ncc_stats
                self.ref2forward_sum[chrom] = f_sum
                self.ref2reverse_sum[chrom] = r_sum
                self.ref2ccbins[chrom] = ccbins
            
            # Store MSCC results if available
            if mscc_stats[0] is not None:  # Check if MSCC data exists
                f_sum, r_sum, ccbins = mscc_stats
                self.mappable_ref2forward_sum[chrom] = f_sum
                self.mappable_ref2reverse_sum[chrom] = r_sum
                self.mappable_ref2ccbins[chrom] = ccbins
    
    def _cleanup(self) -> None:
        """Clean up resources after calculation.
        
        Common cleanup logic for all handlers.
        """
        if hasattr(self, 'align_file') and self.align_file:
            self.align_file.close()
    
    def run(self) -> None:
        """Execute the calculation workflow.
        
        Main entry point that orchestrates the calculation using template method pattern.
        Subclasses implement algorithm-specific logic through abstract methods.
        """
        try:
            self._setup_progress_management()
            
            if self.nworker == 1:
                self._run_singleprocess_calculation()
            else:
                self._run_multiprocess_calculation()
                
        finally:
            self._cleanup()
    
    @abstractmethod
    def _run_singleprocess_calculation(self) -> None:
        """Run calculation in single process mode.
        
        Subclasses must implement algorithm-specific single process logic.
        """
        pass
    
    @abstractmethod
    def _run_multiprocess_calculation(self) -> None:
        """Run calculation in multiprocess mode.
        
        Subclasses must implement algorithm-specific multiprocess logic.
        """
        pass
    
    # Backward compatibility properties
    @property
    def align_file(self):
        """Access to alignment file for backward compatibility."""
        return self._align_file
    
    @align_file.setter
    def align_file(self, value):
        """Set alignment file."""
        self._align_file = value
    
    @property
    def read_len(self) -> Optional[int]:
        """Read length for backward compatibility."""
        return getattr(self, '_read_len', None)
    
    @read_len.setter
    def read_len(self, value: int):
        """Set read length."""
        self._read_len = value