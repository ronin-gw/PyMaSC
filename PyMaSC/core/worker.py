"""Unified worker architecture for multiprocessing calculations.

This module provides a unified architecture for worker processes that
eliminates code duplication and provides a consistent interface for
all calculation algorithms.

Key components:
- BaseWorker: Abstract base class with common worker functionality
- ReadProcessor: Protocol for read processing logic
- UnifiedWorker: Single worker implementation using strategy pattern
- Result collectors for different calculation types
"""
import logging
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from multiprocessing import Process, Queue, Lock
from typing import Optional, Dict, Any, Protocol, runtime_checkable

from pysam import AlignmentFile, AlignedSegment

from PyMaSC.core.interfaces import CrossCorrelationCalculator
from PyMaSC.core.models import WorkerConfig, CalculationConfig, AlgorithmType
from PyMaSC.utils.progress import ProgressHook
from PyMaSC.utils.read_processing import ReadFilter, ReadProcessor as ReadProcessorUtil, create_standard_filter

logger = logging.getLogger(__name__)


@runtime_checkable
class ReadProcessor(Protocol):
    """Protocol for read processing logic.

    Defines the interface for components that process individual reads
    from BAM files. This allows different algorithms to implement their
    own read handling logic while maintaining a consistent interface.
    """

    def should_skip_read(self, read: AlignedSegment) -> bool:
        """Check if a read should be skipped.

        Args:
            read: Read to check

        Returns:
            True if read should be skipped
        """
        ...

    def process_read(self, read: AlignedSegment) -> None:
        """Process a single read.

        Args:
            read: Read to process
        """
        ...


@dataclass
class WorkerResult:
    """Container for worker calculation results.

    Provides a standardized format for returning results from workers
    to the main process, supporting both NCC and MSCC results.
    """
    chromosome: str
    mappable_length: Optional[int] = None
    ncc_forward_sum: Optional[int] = None
    ncc_reverse_sum: Optional[int] = None
    ncc_bins: Optional[Any] = None
    mscc_forward_sum: Optional[int] = None
    mscc_reverse_sum: Optional[int] = None
    mscc_bins: Optional[Any] = None
    metadata: Optional[Dict[str, Any]] = None

    def to_legacy_format(self) -> tuple:
        """Convert to legacy result format for backward compatibility.

        Returns:
            Tuple in legacy format: (chrom, (mappable_len, ncc_stats, mscc_stats))
        """
        ncc_stats = (self.ncc_forward_sum, self.ncc_reverse_sum, self.ncc_bins)
        mscc_stats = (self.mscc_forward_sum, self.mscc_reverse_sum, self.mscc_bins)
        return (self.chromosome, (self.mappable_length, ncc_stats, mscc_stats))


class BaseWorker(Process, ABC):
    """Abstract base class for all calculation workers.

    Provides common functionality for multiprocessing workers including
    process management, communication, and basic workflow. Subclasses
    need only implement the specific calculation logic.

    This design eliminates code duplication and provides a consistent
    interface for all worker types.
    """

    def __init__(self,
                 worker_config: WorkerConfig,
                 order_queue: Queue,
                 report_queue: Queue,
                 logger_lock: Lock,
                 bam_path: str):
        """Initialize base worker.

        Args:
            worker_config: Configuration for this worker
            order_queue: Queue for receiving work orders
            report_queue: Queue for sending results
            logger_lock: Lock for thread-safe logging
            bam_path: Path to BAM file to process
        """
        super().__init__()

        self.config = worker_config
        self.order_queue = order_queue
        self.report_queue = report_queue
        self.logger_lock = logger_lock
        self.bam_path = bam_path

        # Create progress reporter
        self._progress = ProgressHook(report_queue) if worker_config.progress_enabled else None
        self._current_chrom = None

        # Reference information
        calc_config = worker_config.calculation_config
        self._ref2genomelen = dict(zip(calc_config.references, calc_config.lengths))

    def run(self) -> None:
        """Main worker process loop."""
        with self.logger_lock:
            logger.debug(f"{self.name}: Starting worker, PID {os.getpid()}")

        try:
            # Initialize worker-specific components
            self._initialize()

            with AlignmentFile(self.bam_path) as alignfile:
                while True:
                    # Get next chromosome to process
                    chrom = self.order_queue.get()
                    if chrom is None:
                        break

                    # Process chromosome
                    self._process_chromosome(alignfile, chrom)

                    # Report results
                    result = self._collect_results(chrom)
                    self._report_result(result)

        except Exception as e:
            with self.logger_lock:
                logger.error(f"{self.name}: Error in worker: {e}")
            raise

        finally:
            with self.logger_lock:
                logger.debug(f"{self.name}: Shutting down worker")
            self._cleanup()

    def _process_chromosome(self, alignfile: AlignmentFile, chrom: str) -> None:
        """Process all reads from a chromosome.

        Args:
            alignfile: Open BAM file
            chrom: Chromosome to process
        """
        # Update progress tracking
        if self._progress and chrom != self._current_chrom:
            self._current_chrom = chrom
            self._progress.set(chrom, self._ref2genomelen[chrom])

        with self.logger_lock:
            logger.debug(f"{self.name}: Processing {chrom}")

        # Get read processor
        processor = self._get_read_processor()

        # Process reads
        for read in alignfile.fetch(chrom):
            if self._progress:
                self._progress.update(read.reference_start)

            # Skip if needed
            if processor.should_skip_read(read):
                continue

            # Process read
            try:
                processor.process_read(read)
            except Exception as e:
                logger.error(f"Error processing read: {e}")
                raise

    def _report_result(self, result: WorkerResult) -> None:
        """Report result to main process.

        Args:
            result: Result to report
        """
        # Convert to legacy format for backward compatibility
        self.report_queue.put(result.to_legacy_format())

    @abstractmethod
    def _initialize(self) -> None:
        """Initialize worker-specific components."""
        pass

    @abstractmethod
    def _get_read_processor(self) -> ReadProcessor:
        """Get the read processor for this worker.

        Returns:
            ReadProcessor instance
        """
        pass

    @abstractmethod
    def _collect_results(self, chrom: str) -> WorkerResult:
        """Collect results for a chromosome.

        Args:
            chrom: Chromosome to collect results for

        Returns:
            WorkerResult containing calculation results
        """
        pass

    @abstractmethod
    def _cleanup(self) -> None:
        """Clean up worker resources."""
        pass


class StandardReadProcessor:
    """Standard read processor implementation.

    Provides the common read filtering and processing logic used by
    most PyMaSC calculations. Now uses the centralized read processing
    utilities to eliminate duplicate filtering logic.
    """

    def __init__(self,
                 calculator: CrossCorrelationCalculator,
                 mapq_criteria: int):
        """Initialize read processor.

        Args:
            calculator: Calculator to feed reads to
            mapq_criteria: Minimum mapping quality threshold
        """
        self.calculator = calculator
        self.mapq_criteria = mapq_criteria
        # Use centralized read processing utilities
        self._filter = create_standard_filter(mapq_criteria)
        self._processor = ReadProcessorUtil(calculator, self._filter)

    def should_skip_read(self, read: AlignedSegment) -> bool:
        """Check if read should be skipped based on quality criteria.

        Uses centralized filtering logic to ensure consistency.

        Args:
            read: Read to check

        Returns:
            True if read should be skipped
        """
        return self._filter.should_skip_read(read)

    def process_read(self, read: AlignedSegment) -> None:
        """Process a single read by feeding to calculator.

        Uses centralized processing logic to eliminate duplication.

        Args:
            read: Read to process
        """
        self._processor.process_read(read)


class UnifiedWorker(BaseWorker):
    """Unified worker implementation using strategy pattern.

    This single worker class replaces all the specific worker classes
    by using the strategy pattern to delegate algorithm-specific logic
    to the appropriate calculator and processor implementations.

    This eliminates code duplication and provides a consistent interface
    for all calculation types while maintaining flexibility.
    """

    def __init__(self,
                 worker_config: WorkerConfig,
                 order_queue: Queue,
                 report_queue: Queue,
                 logger_lock: Lock,
                 bam_path: str,
                 calculator_factory: Optional[Any] = None):
        """Initialize unified worker.

        Args:
            worker_config: Worker configuration
            order_queue: Queue for receiving work orders
            report_queue: Queue for sending results
            logger_lock: Lock for thread-safe logging
            bam_path: Path to BAM file
            calculator_factory: Optional calculator factory (for testing)
        """
        super().__init__(worker_config, order_queue, report_queue,
                        logger_lock, bam_path)

        self._calculator_factory = calculator_factory
        self._calculator = None
        self._processor = None
        self._secondary_calculator = None  # For dual NCC/MSCC mode
        self._bwfeeder = None

    def _initialize(self) -> None:
        """Initialize calculators based on configuration."""
        from PyMaSC.core.factory import CalculatorFactory

        factory = self._calculator_factory or CalculatorFactory
        calc_config = self.config.calculation_config
        map_config = self.config.mappability_config

        # Create primary calculator
        self._calculator = factory.create_calculator(
            calc_config.algorithm,
            calc_config,
            map_config,
            self.logger_lock,
            self._progress
        )

        # For dual NCC/MSCC mode, create secondary NCC calculator
        if (calc_config.algorithm == AlgorithmType.MSCC and
            not calc_config.skip_ncc):
            # Create NCC calculator
            ncc_config = CalculationConfig(
                algorithm=AlgorithmType.NAIVE_CC,
                max_shift=calc_config.max_shift,
                mapq_criteria=calc_config.mapq_criteria,
                references=calc_config.references,
                lengths=calc_config.lengths,
                read_length=calc_config.read_length
            )
            self._secondary_calculator = factory.create_calculator(
                AlgorithmType.NAIVE_CC,
                ncc_config,
                None,  # No mappability for NCC
                self.logger_lock,
                None   # No separate progress for secondary
            )

        # Initialize BigWig feeder if needed
        if map_config and map_config.is_enabled():
            from PyMaSC.reader.bigwig import BigWigReader
            self._bwfeeder = BigWigReader(str(map_config.mappability_path))

    def _get_read_processor(self) -> ReadProcessor:
        """Get appropriate read processor."""
        if self._processor is None:
            # Create processor based on configuration
            calc_config = self.config.calculation_config

            if self._secondary_calculator:
                # Dual processor for NCC+MSCC
                self._processor = DualReadProcessor(
                    self._calculator,
                    self._secondary_calculator,
                    calc_config.mapq_criteria
                )
            else:
                # Standard single processor
                self._processor = StandardReadProcessor(
                    self._calculator,
                    calc_config.mapq_criteria
                )

        return self._processor

    def _collect_results(self, chrom: str) -> WorkerResult:
        """Collect results from calculators."""
        result = WorkerResult(chromosome=chrom)

        # Flush calculators
        if hasattr(self._calculator, 'flush'):
            self._calculator.flush()
        if self._secondary_calculator and hasattr(self._secondary_calculator, 'flush'):
            self._secondary_calculator.flush()
        if self._bwfeeder and hasattr(self._bwfeeder, 'flush'):
            self._bwfeeder.flush()

        # Get mappable length if available
        if self._bwfeeder and hasattr(self._bwfeeder, 'chrom2mappable_len'):
            result.mappable_length = self._bwfeeder.chrom2mappable_len.get(chrom)
        elif hasattr(self._calculator, 'ref2mappable_len'):
            result.mappable_length = self._calculator.ref2mappable_len.get(chrom)

        # Collect NCC results
        if self._secondary_calculator:
            # From secondary NCC calculator
            result.ncc_forward_sum = self._secondary_calculator.ref2forward_sum.get(chrom)
            result.ncc_reverse_sum = self._secondary_calculator.ref2reverse_sum.get(chrom)
            result.ncc_bins = self._secondary_calculator.ref2ccbins.get(chrom)
        elif not self.config.calculation_config.skip_ncc:
            # From primary calculator if it includes NCC
            result.ncc_forward_sum = self._calculator.ref2forward_sum.get(chrom)
            result.ncc_reverse_sum = self._calculator.ref2reverse_sum.get(chrom)
            result.ncc_bins = self._calculator.ref2ccbins.get(chrom)

        # Collect MSCC results
        if hasattr(self._calculator, 'ref2mappable_forward_sum'):
            mappable_forward = self._calculator.ref2mappable_forward_sum
            mappable_reverse = self._calculator.ref2mappable_reverse_sum
            mappable_bins = self._calculator.ref2mascbins

            if mappable_forward is not None:
                result.mscc_forward_sum = mappable_forward.get(chrom)
                result.mscc_reverse_sum = mappable_reverse.get(chrom)
                result.mscc_bins = mappable_bins.get(chrom)

        if (result.mscc_forward_sum is None and
            self.config.calculation_config.algorithm == AlgorithmType.MSCC):
            # For pure MSCC mode using standard ref2* attributes
            result.mscc_forward_sum = self._calculator.ref2forward_sum.get(chrom)
            result.mscc_reverse_sum = self._calculator.ref2reverse_sum.get(chrom)
            result.mscc_bins = self._calculator.ref2ccbins.get(chrom)

        return result

    def _cleanup(self) -> None:
        """Clean up resources."""
        if self._bwfeeder:
            self._bwfeeder.close()

        # Call finishup on calculators if available
        if hasattr(self._calculator, 'finishup_calculation'):
            self._calculator.finishup_calculation()
        if self._secondary_calculator and hasattr(self._secondary_calculator, 'finishup_calculation'):
            self._secondary_calculator.finishup_calculation()


class DualReadProcessor:
    """Read processor that feeds reads to two calculators.

    Used for simultaneous NCC and MSCC calculations, feeding each
    read to both calculators to compute both standard and
    mappability-corrected results. Now uses centralized utilities.
    """

    def __init__(self,
                 primary_calculator: CrossCorrelationCalculator,
                 secondary_calculator: CrossCorrelationCalculator,
                 mapq_criteria: int):
        """Initialize dual processor.

        Args:
            primary_calculator: Primary calculator (usually MSCC)
            secondary_calculator: Secondary calculator (usually NCC)
            mapq_criteria: Minimum mapping quality threshold
        """
        self.primary_calculator = primary_calculator
        self.secondary_calculator = secondary_calculator
        self.mapq_criteria = mapq_criteria
        # Use centralized read processing utilities
        from PyMaSC.utils.read_processing import create_dual_processor
        self._dual_processor = create_dual_processor(
            primary_calculator,
            secondary_calculator,
            mapq_criteria
        )

    def should_skip_read(self, read: AlignedSegment) -> bool:
        """Check if read should be skipped using centralized logic."""
        return self._dual_processor.filter.should_skip_read(read)

    def process_read(self, read: AlignedSegment) -> None:
        """Process read in both calculators using centralized logic."""
        self._dual_processor.process_read(read)


# Backward compatibility aliases
def create_legacy_worker(worker_type: str, *args, **kwargs) -> BaseWorker:
    """Create worker using legacy interface.

    This function provides backward compatibility for code expecting
    the old worker classes.

    Args:
        worker_type: Type of worker ('ncc', 'mscc', 'ncc_mscc', 'bitarray')
        *args, **kwargs: Legacy arguments

    Returns:
        UnifiedWorker configured appropriately
    """
    # Extract common arguments
    order_queue = args[0]
    report_queue = args[1]
    logger_lock = args[2]
    bam_path = args[3]
    mapq_criteria = args[4]
    max_shift = args[5]
    references = args[6]
    lengths = args[7]

    # Create appropriate configuration
    if worker_type == 'bitarray':
        from PyMaSC.core.models import AlgorithmType
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=max_shift,
            mapq_criteria=mapq_criteria,
            references=references,
            lengths=lengths,
            read_length=kwargs.get('read_len', 50),
            skip_ncc=kwargs.get('skip_ncc', False)
        )
    else:
        # Standard NCC/MSCC workers
        from PyMaSC.core.models import AlgorithmType
        if worker_type == 'ncc':
            algorithm = AlgorithmType.NAIVE_CC
            skip_ncc = False
        elif worker_type == 'mscc':
            algorithm = AlgorithmType.MSCC
            skip_ncc = True
        else:  # ncc_mscc
            algorithm = AlgorithmType.MSCC
            skip_ncc = False

        calc_config = CalculationConfig(
            algorithm=algorithm,
            max_shift=max_shift,
            mapq_criteria=mapq_criteria,
            references=references,
            lengths=lengths,
            read_length=kwargs.get('read_len', 50),
            skip_ncc=skip_ncc
        )

    # Create mappability config if needed
    map_config = None
    if 'mappable_path' in kwargs:
        from PyMaSC.core.models import MappabilityConfig
        map_config = MappabilityConfig(
            mappability_path=kwargs['mappable_path']
        )

    # Create worker config
    worker_config = WorkerConfig(
        calculation_config=calc_config,
        mappability_config=map_config,
        progress_enabled=True,
        logger_lock=logger_lock
    )

    return UnifiedWorker(
        worker_config,
        order_queue,
        report_queue,
        logger_lock,
        bam_path
    )