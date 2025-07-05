"""Enhanced read length estimation with observer pattern support.

This module provides an enhanced version of read length estimation that
supports the new observer pattern for progress reporting while maintaining
full backward compatibility.
"""
import logging
from typing import Optional, Dict, Any, List

from .readlen import estimate_readlen as _original_estimate_readlen
from ..utils.progress import ReadCountProgressBar
from .progress_adapter import ProgressManager, get_progress_manager
from .observer import ProgressObserver

logger = logging.getLogger(__name__)


class ReadLengthEstimator:
    """Enhanced read length estimator with observer pattern support.

    This class wraps the original read length estimation functionality
    and adds support for progress observers, allowing flexible progress
    monitoring during read length estimation.
    """

    def __init__(self, 
                 use_observer: bool = True,
                 progress_manager: Optional[ProgressManager] = None):
        """Initialize read length estimator.

        Args:
            use_observer: Whether to use observer pattern for progress
            progress_manager: Optional progress manager to use
        """
        self.use_observer = use_observer
        self._progress_manager = progress_manager or (
            get_progress_manager() if use_observer else None
        )
        self._observers: List[ProgressObserver] = []

    def attach_observer(self, observer: ProgressObserver) -> None:
        """Attach a progress observer.

        Args:
            observer: Observer to attach
        """
        if observer not in self._observers:
            self._observers.append(observer)

    def detach_observer(self, observer: ProgressObserver) -> None:
        """Detach a progress observer.

        Args:
            observer: Observer to detach
        """
        if observer in self._observers:
            self._observers.remove(observer)

    def estimate(self, 
                 path: str, 
                 esttype: str, 
                 mapq_criteria: int,
                 **kwargs) -> int:
        """Estimate read length with optional observer notifications.

        Args:
            path: Path to BAM file
            esttype: Estimation method
            mapq_criteria: Minimum mapping quality
            **kwargs: Additional arguments for observers

        Returns:
            Estimated read length
        """
        # If not using observers, just call original function
        if not self.use_observer:
            return _original_estimate_readlen(path, esttype, mapq_criteria)

        # Notify observers about start
        if self._progress_manager:
            self._progress_manager.notify_start(
                source="read_length_estimation",
                message=f"Estimating read length from {path}",
                esttype=esttype,
                mapq_criteria=mapq_criteria,
                **kwargs
            )

        for observer in self._observers:
            if hasattr(observer, 'notify_start'):
                observer.notify_start("read_length_estimation", path=path)

        try:
            # Perform estimation
            result = _original_estimate_readlen(path, esttype, mapq_criteria)

            # Notify success
            if self._progress_manager:
                self._progress_manager.notify_complete(
                    source="read_length_estimation",
                    message=f"Read length estimated: {result}bp",
                    result=result
                )

            for observer in self._observers:
                if hasattr(observer, 'notify_complete'):
                    observer.notify_complete("read_length_estimation", result=result)

            return result

        except Exception as e:
            # Notify failure
            if self._progress_manager:
                self._progress_manager.notify_failure(
                    source="read_length_estimation",
                    message=f"Failed to estimate read length: {e}",
                    error=str(e)
                )

            for observer in self._observers:
                if hasattr(observer, 'notify_failure'):
                    observer.notify_failure("read_length_estimation", error=str(e))

            raise


def estimate_readlen_enhanced(path: str, 
                            esttype: str, 
                            mapq_criteria: int,
                            use_observer: bool = True,
                            progress_manager: Optional[ProgressManager] = None,
                            observers: Optional[List[ProgressObserver]] = None,
                            **kwargs) -> int:
    """Enhanced read length estimation with observer pattern support.

    This function provides a convenient interface for read length estimation
    with optional progress observer support. It maintains full backward
    compatibility while adding new capabilities.

    Args:
        path: Path to BAM file
        esttype: Estimation method ('mean', 'median', 'mode', 'min', 'max')
        mapq_criteria: Minimum mapping quality threshold
        use_observer: Whether to use observer pattern (default: True)
        progress_manager: Optional progress manager
        observers: Optional list of observers to attach
        **kwargs: Additional arguments passed to observers

    Returns:
        Estimated read length in base pairs

    Example:
        # Basic usage (backward compatible)
        >>> readlen = estimate_readlen_enhanced(
        ...     "input.bam", "mean", 20, use_observer=False
        ... )

        # With observer
        >>> from PyMaSC.core.observer import FileProgressObserver
        >>> observer = FileProgressObserver("progress.log")
        >>> readlen = estimate_readlen_enhanced(
        ...     "input.bam", "mean", 20, observers=[observer]
        ... )
    """
    estimator = ReadLengthEstimator(
        use_observer=use_observer,
        progress_manager=progress_manager
    )

    if observers:
        for observer in observers:
            estimator.attach_observer(observer)

    return estimator.estimate(path, esttype, mapq_criteria, **kwargs)


# Monkey-patch support for gradual migration
def enable_readlen_observers(enabled: bool = True) -> None:
    """Enable or disable observer pattern for all read length estimation.

    This function allows global configuration of observer support for
    read length estimation, useful for gradual migration.

    Args:
        enabled: Whether to enable observer pattern
    """
    import PyMaSC.core.readlen
    if enabled:
        logger.info("Enabling observer pattern for read length estimation")
        # Store original if not already stored
        if not hasattr(PyMaSC.core.readlen, '_original_estimate_readlen'):
            PyMaSC.core.readlen._original_estimate_readlen = (
                PyMaSC.core.readlen.estimate_readlen
            )

        # Replace with enhanced version
        def _wrapper(path, esttype, mapq_criteria):
            return estimate_readlen_enhanced(
                path, esttype, mapq_criteria, 
                use_observer=True
            )

        PyMaSC.core.readlen.estimate_readlen = _wrapper
    else:
        logger.info("Disabling observer pattern for read length estimation")
        # Restore original if available
        if hasattr(PyMaSC.core.readlen, '_original_estimate_readlen'):
            PyMaSC.core.readlen.estimate_readlen = (
                PyMaSC.core.readlen._original_estimate_readlen
            )