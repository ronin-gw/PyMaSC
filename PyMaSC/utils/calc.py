"""Calculation utilities for PyMaSC analysis.

Provides common calculation functions including filtering, smoothing,
and multiprocessing utilities for cross-correlation analysis.
"""
import fnmatch
from itertools import groupby, chain
import logging
from typing import Any, Iterable, List, Sequence, Optional, Set, Tuple, Union
from functools import wraps
from multiprocessing import Process
from multiprocessing.queues import Queue

import numpy as np
import numpy.typing as npt

from scipy.stats import norm
from typing import Callable, TypeVar
from typing_extensions import ParamSpec

logger = logging.getLogger(__name__)


def moving_avr_filter(arr: np.ndarray, window: int) -> np.ndarray:
    """Apply moving average filter to array.

    Computes a moving average filter over the input array using a specified
    window size. Handles edge effects by using partial windows at array boundaries.

    Args:
        arr: Input array to filter
        window: Window size for moving average (must be positive)

    Returns:
        Filtered array of same length as input with smoothed values

    Note:
        Edge values use progressively larger windows to avoid boundary artifacts
    """
    f = np.repeat(1, window) / float(window)
    avr = np.correlate(arr, f, mode="same")
    h_w = window // 2
    for i in range(h_w):
        avr[i] = np.average(arr[0:(h_w + i)])
        avr[-(i + 1)] = np.average(arr[-(h_w + i):])
    return avr


def filter_chroms(chroms: Union[List[str], Set[str], Iterable[str]], filters: Optional[List[Tuple[bool, List[str]]]]) -> Set[str]:
    """Filter chromosome list using Unix-style patterns.

    Applies include/exclude patterns to filter a list of chromosome names.
    Supports Unix shell-style wildcards (*, ?, []) for flexible pattern matching.

    Args:
        chroms: List of chromosome names to filter
        filters: List of (include_flag, patterns) tuples where include_flag
                is True for include patterns and False for exclude patterns

    Returns:
        Set of chromosome names that match the filtering criteria

    Note:
        Filters are applied in order, with exclude patterns removing chromosomes
        from the current set and include patterns adding them back
    """
    if filters is None:
        logger.debug("There is no chromosome filters.")
        return set(chroms)

    logger.debug("Filtering chromosome lists: " + repr(chroms))
    chroms = set(chroms)
    include_chroms = set()

    for i, (to_include, values) in enumerate(groupby(filters, key=lambda f: f[0])):
        patterns = set(chain(*(f[1] for f in values)))
        logger.debug("{} filters: {}".format("Include" if to_include else "Exclude", repr(patterns)))

        filtered_chroms = set.union(*(set(fnmatch.filter(chroms, p)) for p in patterns))
        logger.debug("Matches: " + repr(filtered_chroms))

        if not to_include:
            include_chroms |= chroms - filtered_chroms
        chroms = filtered_chroms

        logger.debug("current include chroms: " + repr(include_chroms))

    if to_include:
        include_chroms |= chroms

    logger.debug("Include chroms: " + repr(include_chroms))
    return include_chroms


class exec_worker_pool(object):
    """Context manager for multiprocessing worker pool execution.

    Manages the lifecycle of worker processes including startup, task distribution,
    and cleanup. Provides a clean context manager interface for coordinating
    parallel processing across multiple chromosomes.

    Attributes:
        workers: List of worker process instances
        tasks: List of tasks to distribute to workers
        task_queue: Queue for distributing tasks to workers
    """
    def __init__(self, workers: Sequence[Process], tasks: Sequence[str], task_queue: Queue) -> None:
        """Initialize worker pool manager.

        Args:
            workers: List of worker process instances
            tasks: List of tasks to distribute to workers
            task_queue: Queue object for task distribution
        """
        self.workers = workers
        self.tasks = tasks
        self.task_queue = task_queue

    def __enter__(self) -> None:
        """Start all workers and distribute tasks.

        Starts all worker processes, distributes tasks to the queue,
        and sends termination signals for clean shutdown.
        """
        for w in self.workers:
            w.start()
        for t in self.tasks:
            self.task_queue.put(t)
        for _ in range(len(self.workers)):
            self.task_queue.put(None)

    def __exit__(self, type: Optional[type], value: Optional[Exception], traceback: Optional[Any]) -> None:
        """Clean up worker processes.

        Terminates any remaining worker processes to ensure clean shutdown.

        Args:
            type: Exception type (if any)
            value: Exception value (if any)
            traceback: Exception traceback (if any)
        """
        for w in self.workers:
            if w.is_alive():
                w.terminate()
            w.join()


P = ParamSpec('P')
T = TypeVar('T')


def npcalc_with_logging_warn(func: Callable[P, T]) -> Callable[P, T]:
    """Decorator for handling numpy floating point errors gracefully.

    This decorator wraps numerical calculation functions to handle common
    floating point errors that can occur during cross-correlation calculations,
    such as division by zero or invalid operations.
    """
    @wraps(func)
    def _inner(*args: P.args, **kwargs: P.kwargs) -> T:
        try:
            with np.errstate(divide="raise", invalid="raise"):
                return func(*args, **kwargs)
        except (FloatingPointError, ZeroDivisionError) as e:
            logger.debug("catch numpy warning: " + repr(e))
            logger.debug("continue anyway.")
            with np.errstate(divide="ignore", invalid="ignore"):
                return func(*args, **kwargs)
    return _inner


def merge_correlations(
    genome_lengths: npt.NDArray[np.int64],
    correlation_arrays: List[npt.NDArray[np.float64]],
    read_length: int,
    confidence_interval: float = 0.99
) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Merge cross-correlation results using Fisher z-transformation.

    This method replicates the exact behavior of CCResult._merge_cc(),
    combining per-chromosome correlations into genome-wide profiles with
    confidence intervals.

    Mathematical Process:
    1. Fisher z-transformation: z = arctanh(r) for each correlation
    2. Weighted averaging: z_avg = Σ(w_i * z_i) / Σ(w_i) where w_i = n_i - 3
    3. Confidence intervals: CI = tanh(z_avg ± z_α/2 * SE)
    4. Inverse transformation: r = tanh(z_avg) back to correlation

    Args:
        genome_lengths: Array of genome lengths for each chromosome
        correlation_arrays: List of correlation arrays, one per chromosome
        read_length: Read length for position-dependent weighting

    Returns:
        Tuple of (merged_correlations, lower_confidence, upper_confidence)
    """
    ns = genome_lengths

    merged_r = []
    interval_upper = []
    interval_lower = []

    for i, __ccs in enumerate(zip(*correlation_arrays)):
        _ccs: Tuple[np.float64, ...] = __ccs
        # Filter out NaN values
        nans = np.isnan(_ccs)
        ccs = np.array(_ccs)[~nans]

        # Calculate weights (effective sample size - 3)
        if len(ns.shape) == 1:
            _ns = ns[~nans] - 3
        else:
            _ns = ns[~nans, abs(read_length - i)] - 3

        # Fisher z-transformation
        zs = np.arctanh(ccs)

        # Filter out infinite values
        infs = np.isinf(zs)
        zs = zs[~infs]
        _ns = _ns[~infs]

        # Weighted average in z-space
        avr_z = np.average(zs, weights=_ns)

        # Calculate confidence interval
        z_interval = norm.ppf(1 - (1 - confidence_interval) / 2) * np.sqrt(1 / np.sum(_ns))
        z_interval_upper = avr_z + z_interval
        z_interval_lower = avr_z - z_interval

        # Transform back to correlation space
        merged_r.append(np.tanh(avr_z))
        interval_upper.append(np.tanh(z_interval_upper))
        interval_lower.append(np.tanh(z_interval_lower))

    return (
        np.array(merged_r, dtype=np.float64),
        np.array(interval_lower, dtype=np.float64),
        np.array(interval_upper, dtype=np.float64)
    )
