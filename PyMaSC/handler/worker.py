"""Worker architecture for multiprocessing calculations.

Provides worker processes for cross-correlation calculations.

Key components:
- BaseWorker: Abstract base class with common worker functionality
- ReadProcessor: Protocol for read processing logic
- CalcWorker: Worker implementation using strategy pattern
- Result collectors for different calculation types
"""
import logging
import os
import traceback
import sys
from abc import ABC, abstractmethod
from multiprocessing import Process, Queue
from multiprocessing.synchronize import Lock
from typing import Optional, Any

from PyMaSC.reader.bam import BAMFileProcessor, AlignmentLike
from PyMaSC.interfaces.config import PyMaSCConfig, CalculationTarget
from PyMaSC.interfaces.calculator import CrossCorrelationCalculator
from .read import ReadProcessor, create_read_processor
from .factory import create_calculator
from PyMaSC.result import EmptyResult, EmptyNCCResult, EmptyMSCCResult, EmptyBothChromResult
from PyMaSC.utils.progress import ProgressHook

logger = logging.getLogger(__name__)


class BaseWorker(Process, ABC):
    """Abstract base class for calculation workers.

    Provides common functionality for multiprocessing workers including
    process management and communication.
    """

    def __init__(self,
                 config: PyMaSCConfig,
                 order_queue: Queue,
                 report_queue: Queue,
                 logger_lock: Lock,
                 bam_path: str):
        """Initialize base worker.

        Args:
            config: PyMaSC configuration
            order_queue: Queue for receiving chromosome processing orders
            report_queue: Queue for sending results back to main process
            logger_lock: Lock for thread-safe logging
            bam_path: Path to BAM file for read access
        """
        super().__init__()

        self.config = config
        self.order_queue = order_queue
        self.report_queue = report_queue
        self.logger_lock = logger_lock
        self.bam_path = bam_path

        # Create progress reporter
        self._progress = ProgressHook(report_queue)
        self._current_chrom: Optional[str] = None

        # Reference information
        self._ref2genomelen = config.ref2lengths

    def run(self) -> None:
        """Main worker process loop."""
        try:
            # Initialize worker-specific components
            self._initialize()

            with BAMFileProcessor(self.bam_path) as alignfile:
                while True:
                    # Get next chromosome to process
                    chrom = self.order_queue.get()
                    if chrom is None:
                        break

                    # Process chromosome
                    self._process_chromosome(alignfile, chrom)

                    # Report results
                    self._report_result(chrom)

        except KeyboardInterrupt:
            if 0 < logger.level <= logging.DEBUG:
                raise

        except Exception as e:
            with self.logger_lock:
                logger.error(f"{self.name}: Error in worker: {e}")

            tb = traceback.format_exc()
            self.report_queue.put(('__ERROR__', tb))
            sys.stderr.write(tb)
            sys.stderr.flush()
            os._exit(1)

        finally:
            with self.logger_lock:
                logger.debug(f"{self.name}: Shutting down worker")
            self._cleanup()

    def _process_chromosome(self, alignfile: AlignmentLike, chrom: str) -> None:
        """Process all reads from a chromosome.

        Args:
            alignfile: Open BAM file
            chrom: Chromosome to process
        """
        # Update progress tracking
        if self._progress and chrom != self._current_chrom:
            self._current_chrom = chrom
            self._progress.enable_bar()
            self._progress.reset_progress("<1II1>" * 12, '>', '<', chrom, self._ref2genomelen[chrom])

        with self.logger_lock:
            logger.debug(f"{self.name}: Processing {chrom}")

        # Get read processor
        processor = self._get_read_processor()

        # Process reads
        for read in alignfile.fetch(chrom):
            # Update progress (always call if progress exists)
            if self._progress:
                self._progress.update(read.reference_start)

            # Process read
            processor.process_read(read)

    def _report_result(self, chrom: str) -> None:
        """Report result to main process.

        Args:
            result: Result to report
        """
        pass

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
    def _cleanup(self) -> None:
        """Clean up worker resources."""
        pass


class CalcWorker(BaseWorker):
    """Unified worker implementation using strategy pattern.

    This single worker class replaces all the specific worker classes
    by using the strategy pattern to delegate algorithm-specific logic
    to the appropriate calculator and processor implementations.

    This eliminates code duplication and provides a consistent interface
    for all calculation types while maintaining flexibility.
    """

    def __init__(self,
                 config: PyMaSCConfig,
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
        super().__init__(config, order_queue, report_queue,
                         logger_lock, bam_path)

        self._calculator_factory = calculator_factory
        self._calculator: Optional[CrossCorrelationCalculator] = None
        self._processor: Optional[ReadProcessor] = None
        self._bwfeeder: Optional[Any] = None

    def _initialize(self) -> None:
        """Initialize calculators based on configuration."""
        self._calculator, self._bwfeeder = create_calculator(
            self.config,
            self.logger_lock,
            self._progress
        )

    def _get_read_processor(self) -> ReadProcessor:
        """Get appropriate read processor."""
        if self._processor is None:
            assert self._calculator is not None
            self._processor = create_read_processor(
                self._calculator,
                self.config.mapq_criteria
            )

        return self._processor

    def _report_result(self, chrom: str) -> None:
        """Report result to main process, creating empty result if no data processed.

        This method ensures that all chromosomes produce a result, even if they
        contain no reads. Empty results maintain genome length consistency
        between single-process and parallel processing modes.
        """
        assert self._calculator is not None

        self._calculator.flush(chrom)

        try:
            result = self._calculator.get_result(chrom)
        except KeyError:
            # No data was processed for this chromosome
            # Create appropriate empty result to maintain genome length consistency
            result = self._create_empty_result(chrom)
        self.report_queue.put((chrom, result))

    def _create_empty_result(self, chrom: str) -> EmptyResult:
        """Create empty result for chromosomes with no data.

        Args:
            chrom: Chromosome name

        Returns:
            Appropriate empty result based on calculation target
        """
        genome_length = self._ref2genomelen[chrom]
        max_shift = self.config.max_shift
        read_len = self.config.read_length

        # read_length is required to create meaningful empty results
        if read_len is None:
            raise ValueError(f"read_length is required to create empty result for chromosome {chrom}")

        # Create appropriate empty result based on calculation target
        if self.config.target is CalculationTarget.NCC:
            return EmptyNCCResult.create_empty(genome_length, max_shift, read_len)
        elif self.config.target is CalculationTarget.MSCC:
            return EmptyMSCCResult.create_empty(genome_length, max_shift, read_len)
        elif self.config.target is CalculationTarget.BOTH:
            return EmptyBothChromResult.create_empty(genome_length, max_shift, read_len)
        else:
            raise ValueError(f"Unknown calculation target: {self.config.target}")

    def _cleanup(self) -> None:
        """Clean up resources."""
        if self._bwfeeder:
            self._bwfeeder.close()
