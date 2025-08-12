"""Mappability management and calculation for PyMaSC.

Provides mappability calculations with parallel processing, caching,
and statistics persistence. Mappability data is essential for MSCC
calculations to correct for genomic regions with variable mappability.

Key components:
- MappabilityHandler: Main interface for mappability operations
- MappabilityCalcWorker: Worker process for parallel calculation
- NumpyEncoder: JSON serialization for numpy data types
- Exception classes for error handling

The module handles BigWig file access, parallel computation across
chromosomes, and persistent storage of mappability statistics to
avoid redundant calculations.
"""
from __future__ import annotations

import logging
import os
import json
from multiprocessing import Queue, Process, Lock
from pathlib import Path
from typing import Any, Dict, Optional, Union

import sys
if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

import numpy as np

from PyMaSC.interfaces.config import MappabilityConfig
from PyMaSC.core.mappability import MappableLengthCalculator
from PyMaSC.utils.progress import ProgressHook, MultiLineProgressManager
from PyMaSC.utils.output import prepare_outdir
from PyMaSC.utils.calc import exec_worker_pool

logger = logging.getLogger(__name__)


class BWIOError(IOError):
    """Exception raised for BigWig file I/O errors.

    Indicates problems accessing or reading BigWig mappability files,
    including file not found, permission denied, or file corruption.
    """
    pass


class JSONIOError(IOError):
    """Exception raised for JSON statistics file I/O errors.

    Indicates problems with mappability statistics file operations,
    including read/write failures or directory access issues.
    """
    pass


class NeedUpdate(Exception):
    """Exception indicating mappability statistics need updating.

    Raised when existing statistics are incompatible with current
    analysis parameters and must be recalculated.
    """
    pass


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder for numpy data types.

    Extends the standard JSON encoder to handle numpy numeric types
    that are commonly used in mappability calculations. Converts
    numpy integers and floats to standard Python types for JSON
    serialization.
    """

    def default(self, obj: Any) -> Union[int, float, Any]:
        """Convert numpy types to JSON-serializable Python types.

        Args:
            obj: Object to serialize

        Returns:
            JSON-serializable representation of the object

        Raises:
            TypeError: If object type is not supported
        """
        if isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        else:
            return super(NumpyEncoder, self).default(obj)


class MappabilityHandler(MappableLengthCalculator):
    """Mappability calculation and management.

    Handles mappability analysis including parallel calculation,
    statistics caching, and file management. Extends MappableLengthCalculator
    with persistence and workflow management capabilities.

    Manages BigWig file access, calculates mappability statistics
    across chromosomes using multiple workers, and persists results
    to JSON files for future reuse.

    Attributes:
        path: Path to BigWig mappability file
        max_shift: Maximum shift distance for analysis
        map_path: Path for statistics cache file
        nworker: Number of parallel worker processes
        chrom2mappable_len: Per-chromosome mappability arrays
        mappable_len: Genome-wide mappability array
        need_save_stats: Whether statistics need to be saved
    """

    @staticmethod
    def calc_mappable_len_required_shift_size(readlen: int, max_shift: int) -> int:
        """Calculate required shift size for mappability analysis.

        Determines the appropriate shift size based on read length
        and maximum shift parameters to ensure complete coverage
        of the analysis window.

        Args:
            readlen: Read length in base pairs
            max_shift: Maximum shift distance requested

        Returns:
            Required shift size for mappability calculation
        """
        return max_shift - readlen + 1 if max_shift > 2*readlen - 1 else readlen

    @classmethod
    def from_config(cls, config: MappabilityConfig) -> Self:
        """Create MappabilityHandler from configuration.

        Args:
            config: Configuration containing mappability parameters

        Returns:
            MappabilityHandler instance configured with provided parameters
        """
        return cls(
            path=config.mappability_path,
            max_shift=config.max_shift,
            readlen=config.read_length,
            map_path=config.mappability_stats_path,
            nworker=config.nproc
        )

    def __init__(
        self,
        path: os.PathLike[str],
        max_shift: int = 0,
        readlen: int = 0,
        map_path: Optional[os.PathLike[str]] = None,
        nworker: int = 1
    ) -> None:
        """Initialize mappability handler.

        Sets up mappability analysis with specified parameters,
        validates file access, and determines whether calculation
        or cached statistics loading is required.

        Args:
            path: Path to BigWig mappability file
            max_shift: Maximum shift distance for cross-correlation
            readlen: Read length in base pairs
            map_path: Path for statistics cache file (optional)
            nworker: Number of parallel worker processes

        Raises:
            BWIOError: If BigWig file is not accessible
            JSONIOError: If statistics file operations fail
        """
        max_shift = self.calc_mappable_len_required_shift_size(readlen, max_shift)
        self.nworker = nworker

        # Convert path to string for file access and parent class compatibility
        path_str = os.fspath(path)

        if not os.access(path_str, os.R_OK):
            reason = "file is unreadable." if Path(path).is_file() else "no such file."
            logger.critical("Failed to open '{}': {}".format(path, reason))
            raise BWIOError

        super(MappabilityHandler, self).__init__(path_str, max_shift)
        self.close()
        self._progress.disable_bar()
        self.need_save_stats = True

        if map_path:
            self.map_path = map_path
        else:
            path_obj = Path(path)
            # Create: /path/to/file.bw -> /path/to/file_mappability.json
            stem_with_suffix = path_obj.with_suffix("").name + "_mappability"
            self.map_path = path_obj.parent / f"{stem_with_suffix}.json"

        if not Path(self.map_path).exists():
            self._check_saving_directory_is_writable()
            logger.info("Calcurate mappable length with max shift size {}.".format(max_shift))
        elif not Path(self.map_path).is_file():
            logger.critical("Specified path is not file: '{}'".format(self.map_path))
            raise JSONIOError
        elif not os.access(self.map_path, os.R_OK):
            logger.error("Failed to read '{}'".format(self.map_path))
        else:
            self._try_load_mappability_stats()
            if self.need_save_stats:
                self._check_stats_is_overwritable()
                logger.info("Calcurate mappable length with max shift size {}.".format(max_shift))
            else:
                logger.info("Use mappability stats read from '{}'".format(self.map_path))

    def _check_saving_directory_is_writable(self) -> None:
        dirname = str(Path(self.map_path).parent)
        dirname = dirname if dirname else '.'
        if not prepare_outdir(dirname, logger):
            raise JSONIOError

    def _try_load_mappability_stats(self) -> None:
        try:
            stats = self._read_mappability_stats()
        except IOError as e:
            logger.error("Failed to read '{}'".format(self.map_path))
            logger.error("[Errno {}] {}".format(e.errno, str(e)))
        except (TypeError, OverflowError, ValueError, KeyError, IndexError):
            logger.error("Failed to load json file: '{}'".format(self.map_path))
        except NeedUpdate:
            pass
        else:
            self._load_mappability_stats(stats)

    def _read_mappability_stats(self) -> Dict[str, Any]:
        with open(self.map_path) as f:
            stats = json.load(f)

        for k in ("max_shift", "__whole__", "references"):
            if k not in stats:
                logger.error("Mandatory key '{}' not found.".format(k))
                raise KeyError(k)

        if stats["max_shift"] < self.max_shift:
            logger.info("Specified shift length longer than former analysis. The stats will be updated.")
            raise NeedUpdate

        if stats["max_shift"] != len(stats["__whole__"]) - 1:
            logger.error("Max shift length for whole genome unmatched.")
            raise IndexError

        for ref in self.chromsizes:
            if ref not in stats["references"]:
                logger.error("Reference '{}' not found.".format(ref))
                raise KeyError(ref)

            if stats["max_shift"] != len(stats["references"][ref]) - 1:
                logger.error("Max shift length for '{}' unmatched.".format(ref))
                raise IndexError

        return stats

    def _load_mappability_stats(self, stats: Dict[str, Any]) -> None:
        self.mappable_len = stats["__whole__"][:self.max_shift + 1]
        self.chrom2mappable_len = {ref: b[:self.max_shift + 1] for ref, b in stats["references"].items()}
        self.chrom2is_called = {ref: True for ref in self.chromsizes}
        self.is_called = True
        self.need_save_stats = False

    def _check_stats_is_overwritable(self) -> None:
        if not os.access(self.map_path, os.W_OK):
            logger.critical("Failed to overwrite '{}'".format(self.map_path))
            raise JSONIOError
        else:
            logger.warning("Existing file '{}' will be overwritten.".format(self.map_path))

    def save_mappability_stats(self) -> None:
        """Save mappability statistics to JSON cache file.

        Persists calculated mappability statistics to a JSON file
        for future reuse, avoiding redundant calculations. Uses
        NumpyEncoder to handle numpy data types.

        The statistics include per-chromosome and genome-wide
        mappability arrays for all shift distances up to max_shift.

        Note:
            Only saves if need_save_stats flag is True
        """
        if not self.need_save_stats:
            return logger.info("Mappability stats updating is not required.")

        logger.info("Save mappable length to '{}'".format(self.map_path))
        try:
            with open(self.map_path, 'w') as f:
                json.dump({
                    "max_shift": self.max_shift,
                    "__whole__": self.mappable_len,
                    "references": self.chrom2mappable_len
                }, f, indent=4, sort_keys=True, cls=NumpyEncoder)
        except IOError as e:
            logger.error("Faild to output: {}\n[Errno {}] {}".format(
                         e.filename, e.errno, str(e)))

        self.need_save_stats = False

    def calc_mappability(self, _chrom: Optional[str] = None) -> None:
        """Calculate mappability statistics using parallel workers.

        Orchestrates parallel calculation of mappability statistics
        across chromosomes using multiple worker processes. Manages
        progress reporting and result collection from workers.

        The method only processes chromosomes that haven't been
        calculated yet, enabling incremental updates and resumption
        of interrupted calculations.

        Note:
            Uses multiprocessing for efficient parallel computation
            across available CPU cores
        """
        target_chroms = [c for c, b in self.chrom2is_called.items() if b is False]
        if not target_chroms:
            return self._sumup_mappability()

        order_queue: Queue = Queue()
        report_queue: Queue = Queue()
        logger_lock = Lock()
        progress = MultiLineProgressManager()

        workers = [MappabilityCalcWorker(Path(self.path), self.max_shift, order_queue, report_queue, logger_lock)
                   for _ in range(min(self.nworker, len(target_chroms)))]

        with exec_worker_pool(workers, target_chroms, order_queue):
            while not self.is_called:
                chrom, obj = report_queue.get()
                if chrom is None:  # update progress
                    chrom, body = obj
                    with logger_lock:
                        progress.update(chrom, body)
                else:
                    length = obj
                    self.chrom2mappable_len[chrom] = tuple(length)
                    self.chrom2is_called[chrom] = True
                    if all(self.chrom2is_called.values()):
                        self.is_called = True
                    with logger_lock:
                        progress.erase(chrom)

        progress.clean()
        self._sumup_mappability()

    def _sumup_mappability(self) -> None:
        for length in self.chrom2mappable_len.values():
            for i in range(self.max_shift + 1):
                self.mappable_len[i] += length[i]


class MappabilityCalcWorker(Process):
    """Worker process for parallel mappability calculation.

    Implements a multiprocessing worker that calculates mappability
    statistics for assigned chromosomes. Each worker operates
    independently on different chromosomes to enable parallel
    computation across the genome.

    The worker receives chromosome assignments through a queue,
    processes each chromosome using MappableLengthCalculator,
    and reports results back to the main process.

    Attributes:
        path: Path to BigWig mappability file
        max_shift: Maximum shift distance for calculation
        order_queue: Queue for receiving chromosome assignments
        report_queue: Queue for sending results and progress updates
        logger_lock: Lock for thread-safe logging
    """

    def __init__(self, path: os.PathLike[str], max_shift: int, order_queue: Queue, report_queue: Queue, logger_lock: Any) -> None:
        """Initialize mappability calculation worker.

        Args:
            path: Path to BigWig mappability file
            max_shift: Maximum shift distance for calculation
            order_queue: Queue for receiving work assignments
            report_queue: Queue for sending results
            logger_lock: Lock for thread-safe logging operations
        """
        super(MappabilityCalcWorker, self).__init__()

        self.path = os.fspath(path)
        self.max_shift = max_shift
        self.order_queue = order_queue
        self.report_queue = report_queue
        self.logger_lock = logger_lock

    def run(self) -> None:
        """Execute worker process main loop.

        Continuously processes chromosome assignments from the order
        queue until a termination signal (None) is received. For each
        chromosome, calculates mappability statistics and reports
        results back through the report queue.

        The worker maintains thread-safe logging and handles proper
        cleanup of resources upon completion.
        """
        super().run()

        # Create calculator instance in the worker process to avoid pickle issues
        calculator = MappableLengthCalculator(self.path, self.max_shift, self.logger_lock)
        calculator._progress.disable_bar()
        calculator._progress = ProgressHook(self.report_queue)

        try:
            while True:
                chrom = self.order_queue.get()
                if chrom is None:
                    break

                with self.logger_lock:
                    logger.debug("{}: Process {}...".format(self.name, chrom))
                calculator.calc_mappability(chrom)
                self.report_queue.put((chrom, calculator.chrom2mappable_len[chrom]))
        finally:
            with self.logger_lock:
                logger.debug("{}: Goodbye.".format(self.name))
            calculator.close()
