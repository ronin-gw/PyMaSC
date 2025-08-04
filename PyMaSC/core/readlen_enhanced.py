"""Enhanced read length estimation - Observer functionality disabled.

This module previously provided observer pattern support for read length
estimation. The observer functionality has been removed, and this module
now serves as a simple wrapper around the original read length estimation.
"""
import logging
from typing import Optional, Any

from .readlen import estimate_readlen as _original_estimate_readlen

logger = logging.getLogger(__name__)


class ReadLengthEstimator:
    """Read length estimator - Observer functionality disabled.

    This class now serves as a simple wrapper around the original
    read length estimation functionality. Observer pattern support
    has been removed.
    """

    def __init__(self,
                 use_observer: bool = True,
                 progress_manager: Optional[Any] = None):
        """Initialize read length estimator.

        Args:
            use_observer: Ignored - observer functionality disabled
            progress_manager: Ignored - observer functionality disabled
        """
        # Observer functionality disabled
        pass

    def attach_observer(self, observer: Any) -> None:
        """Attach a progress observer - DISABLED.
        
        Args:
            observer: Ignored
        """
        # Observer functionality disabled
        pass

    def detach_observer(self, observer: Any) -> None:
        """Detach a progress observer - DISABLED.
        
        Args:
            observer: Ignored
        """
        # Observer functionality disabled
        pass

    def estimate(self,
                 path: str,
                 esttype: str,
                 mapq_criteria: int,
                 **kwargs: Any) -> int:
        """Estimate read length - simplified version.

        Args:
            path: Path to BAM file
            esttype: Estimation method
            mapq_criteria: Minimum mapping quality
            **kwargs: Ignored

        Returns:
            Estimated read length
        """
        # Always use original function, observer support removed
        return _original_estimate_readlen(path, esttype, mapq_criteria)


# Convenience function for backward compatibility
def estimate_readlen_enhanced(path: str,
                            esttype: str,
                            mapq_criteria: int,
                            use_observer: bool = True,
                            **kwargs: Any) -> int:
    """Enhanced read length estimation - simplified version.
    
    Observer functionality has been removed. This function now
    directly calls the original estimate_readlen function.
    
    Args:
        path: Path to BAM file
        esttype: Estimation method
        mapq_criteria: Minimum mapping quality
        use_observer: Ignored
        **kwargs: Ignored
    
    Returns:
        Estimated read length
    """
    return _original_estimate_readlen(path, esttype, mapq_criteria)