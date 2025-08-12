"""Calculation handler using factory pattern.

Provides the main handler for cross-correlation calculations using
factory-based calculator creation.
"""
from __future__ import annotations

import logging
import os
from multiprocessing import Queue, Lock
from typing import Optional, List, Tuple, Union, Any

from PyMaSC.core.exceptions import InputUnseekable, NothingToCalc
from PyMaSC.interfaces.config import PyMaSCConfig
from PyMaSC.reader.bam import BAMFileProcessor, BAMNoReadsError, BAMNoTargetChroms
from PyMaSC.core.readlen import estimate_readlen
from .read import create_read_processor
from .mappability import MappabilityHandler
from .factory import create_calculator
from .worker import CalcWorker
from PyMaSC.result import ChromResult, GenomeWideResult, aggregate_results
from PyMaSC.utils.calc import exec_worker_pool
from PyMaSC.utils.progress import ProgressBar, ProgressHook, MultiLineProgressManager

logger = logging.getLogger(__name__)


class CalcHandler:
    """Main calculation handler for cross-correlation analysis."""

    def __init__(self, path: os.PathLike[str], config: PyMaSCConfig) -> None:
        """Initialize calculation handler.

        Args:
            path: Path to BAM file
            config: Configuration for calculation parameters
        """
        self.path = str(path)
        self.config = config

        # Initialize BAM file processor
        self.bam_processor = BAMFileProcessor(self.path)

        # Validate BAM file and extract metadata
        try:
            references, lengths = self.bam_processor.apply_chromfilter(self.config.chromfilter)
        except BAMNoReadsError:
            raise ValueError("File has no sequences defined.")
        except BAMNoTargetChroms:
            raise NothingToCalc

        # Update config with filtered references
        self.config.ref2lengths = dict(zip(references, lengths))

        # Check for index in multiprocess mode
        if self.config.multiprocess and not self.bam_processor.check_multiprocess_compatibility():
            self.config.nproc = 1
        elif self.config.multiprocess:
            self.bam_processor.close()

        # Mappability handler reference
        self.mappability_handler: Optional[MappabilityHandler] = None

    @property
    def read_len(self) -> Optional[int]:
        """Get the read length from the configuration."""
        return self.config.read_length

    @read_len.setter
    def read_len(self, value: int) -> None:
        """Set the read length in the configuration."""
        self.config.read_length = value

    def estimate_readlen(self) -> int:
        """Set or estimate read length.

        Raises:
            InputUnseekable: If estimation needed but input is unseekable
            ValueError: If read length exceeds maximum shift distance
        """
        if self.path == '-':
            logger.error("Cannot execute read length checking for unseekable input.")
            raise InputUnseekable

        logger.info(f"Check read length... : {self.path}")
        # Use config-based estimation for cleaner interface
        read_len = estimate_readlen(
            path=self.path,
            esttype=self.config.esttype.value,
            mapq_criteria=self.config.mapq_criteria
        )

        if read_len > self.config.max_shift:
            logger.error(f"Read length ({read_len}) seems to be longer than "
                         f"shift size ({self.config.max_shift}).")
            raise ValueError

        return read_len

    def set_mappability_handler(self, mappability_handler: MappabilityHandler) -> None:
        """Set mappability handler and validate chromosome sizes.

        Args:
            mappability_handler: MappabilityHandler instance
        """
        self.mappability_handler = mappability_handler

        bw_chromsizes = self.mappability_handler.chromsizes
        # Use BAM processor to validate and reconcile chromosome sizes
        updated_sizes = self.bam_processor.validate_chromosome_sizes(bw_chromsizes)

        # Update lengths based on validation results
        for chrom, length in updated_sizes.items():
            if chrom in self.config.ref2lengths:
                self.config.ref2lengths[chrom] = length

    def run_calculation(self) -> GenomeWideResult:
        """Execute cross-correlation calculation workflow."""
        if self.config.multiprocess:
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
        calculator, _ = create_calculator(self.config)

        #
        read_processor = create_read_processor(calculator, self.config.mapq_criteria)
        progress_bar = ProgressBar()

        processed_bases = 0
        current_chrom = None

        for read in self.bam_processor:
            chrom = read.reference_name

            if chrom is None or chrom not in self.config.references:
                continue

            processed_bases = max(processed_bases, read.reference_start)
            if chrom != current_chrom:
                current_chrom = chrom
                progress_bar.set(chrom, self.config.ref2lengths[chrom])
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
        self._report_queue: Queue[Tuple[Optional[str], Union[ChromResult, Tuple[str, str]]]] = Queue()
        self._logger_lock = Lock()

        # Create workers using factory
        workers = []
        num_workers = min(self.config.nproc, len(self.config.references))
        for _ in range(num_workers):
            worker = CalcWorker(
                self.config,
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
        _chrom2finished = {c: False for c in self.config.references}
        _worker_results = {}  # Collect results for new aggregation system
        progress = MultiLineProgressManager()

        with exec_worker_pool(workers, self.config.references, self._order_queue):
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
