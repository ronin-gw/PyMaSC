"""Enhanced statistics calculator for cross-correlation analysis.

This module extends the original CCStats functionality with better
organization, cleaner interfaces, and improved error handling.

Key improvements:
- Separation of calculation logic from initialization complexity
- Better type safety and validation
- Modular calculation methods
- Enhanced error handling with meaningful messages
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Union
import logging

import numpy as np

from PyMaSC.core.result_container import NCCData, MSCCData
from PyMaSC.utils.calc import moving_avr_filter
from PyMaSC.utils.stats_utils import CrossCorrelationMetrics

logger = logging.getLogger(__name__)


@dataclass
class StatisticsConfig:
    """Configuration for statistics calculations.
    
    Attributes:
        min_calc_width: Width for minimum value calculation
        mv_avr_filter_len: Moving average filter length
        filter_mask_len: Mask length around read length
        output_warnings: Whether to output warning messages
        do_fragment_estimation: Whether to perform fragment length estimation
        expected_fragment_length: User-specified expected fragment length
    """
    min_calc_width: int = 10
    mv_avr_filter_len: int = 15
    filter_mask_len: int = 5
    output_warnings: bool = True
    do_fragment_estimation: bool = False
    expected_fragment_length: Optional[int] = None


@dataclass
class StatisticsResult:
    """Container for calculated statistics.
    
    Attributes:
        cc_min: Minimum cross-correlation value
        cc_at_read_length: Cross-correlation at read length
        cc_at_fragment_length: Cross-correlation at fragment length
        nsc: Normalized Strand Coefficient
        rsc: Relative Strand Coefficient
        estimated_fragment_length: Estimated fragment length
        fragment_peak_width: Full Width at Half Maximum
        vsn: Variance Stabilizing Normalization factor
    """
    cc_min: float
    cc_at_read_length: float
    cc_at_fragment_length: Optional[float] = None
    nsc: Optional[float] = None
    rsc: Optional[float] = None
    estimated_fragment_length: Optional[int] = None
    fragment_peak_width: Optional[Union[int, bool]] = None
    vsn: Optional[float] = None


class StatisticsCalculator:
    """Enhanced cross-correlation statistics calculator.
    
    Provides comprehensive statistical analysis of cross-correlation data
    with improved organization and error handling compared to the original
    CCStats implementation.
    
    This class maintains compatibility with existing code while providing
    a cleaner interface for statistical calculations.
    """
    
    def __init__(self, 
                 data: Union[NCCData, MSCCData],
                 read_length: int,
                 config: Optional[StatisticsConfig] = None):
        """Initialize statistics calculator.
        
        Args:
            data: Cross-correlation data (NCC or MSCC)
            read_length: Read length in base pairs
            config: Statistics calculation configuration
        """
        self.data = data
        self.read_length = read_length
        self.config = config or StatisticsConfig()
        
        # Computed values
        self._smoothed_cc: Optional[np.ndarray] = None
        self._result: Optional[StatisticsResult] = None
        
        # Validation
        self._validate_inputs()
    
    def _validate_inputs(self) -> None:
        """Validate input data and configuration."""
        if self.read_length <= 0:
            raise ValueError("Read length must be positive")
        
        if isinstance(self.data, NCCData):
            if len(self.data.ccbins) == 0:
                raise ValueError("Cross-correlation data cannot be empty")
        elif isinstance(self.data, MSCCData):
            if len(self.data.ccbins) == 0:
                raise ValueError("MSCC data cannot be empty")
        else:
            raise TypeError("Data must be NCCData or MSCCData")
    
    def calculate_statistics(self) -> StatisticsResult:
        """Calculate comprehensive cross-correlation statistics.
        
        Returns:
            StatisticsResult with all calculated metrics
        """
        if self._result is not None:
            return self._result
        
        # Get cross-correlation array
        cc = self.data.ccbins
        
        # Calculate basic statistics
        cc_min = self._calculate_minimum(cc)
        cc_at_read_length = self._calculate_cc_at_read_length(cc)
        
        # Initialize result
        result = StatisticsResult(
            cc_min=cc_min,
            cc_at_read_length=cc_at_read_length
        )
        
        # Calculate fragment-based metrics if fragment length is known
        if self.config.expected_fragment_length:
            fragment_metrics = self._calculate_fragment_metrics(
                cc, self.config.expected_fragment_length
            )
            result.cc_at_fragment_length = fragment_metrics[0]
            result.nsc = fragment_metrics[1]
            result.rsc = fragment_metrics[2]
            
            # Calculate peak width
            try:
                result.fragment_peak_width = self._calculate_peak_width(
                    self.config.expected_fragment_length
                )
            except Exception as e:
                if self.config.output_warnings:
                    logger.warning(f"Failed to calculate peak width: {e}")
                result.fragment_peak_width = None
        
        # Estimate fragment length if requested
        if self.config.do_fragment_estimation:
            estimated_length = self._estimate_fragment_length(cc)
            result.estimated_fragment_length = estimated_length
            
            if estimated_length:
                # Calculate metrics for estimated length
                est_metrics = self._calculate_fragment_metrics(cc, estimated_length)
                # Store as estimated values (could extend result for this)
        
        # Calculate VSN if we have fragment metrics and width
        if (result.cc_at_fragment_length is not None and 
            result.fragment_peak_width is not None and
            isinstance(result.fragment_peak_width, int)):
            result.vsn = self._calculate_vsn(
                result.cc_at_fragment_length,
                result.fragment_peak_width
            )
        
        self._result = result
        return result
    
    def _calculate_minimum(self, cc: np.ndarray) -> float:
        """Calculate minimum cross-correlation value.
        
        Uses the median of the last min_calc_width values as a robust
        estimate of the background correlation level.
        
        Args:
            cc: Cross-correlation array
            
        Returns:
            Minimum correlation value
        """
        width = min(self.config.min_calc_width, len(cc))
        if width == 0:
            return 0.0
        
        # Use median of last values as minimum estimate
        tail_values = cc[-width:]
        cc_min = float(np.median(tail_values))
        
        # Warning check for potential issues
        if (self.config.output_warnings and len(cc) > 10):
            near_zero_median = np.median(cc[:10])
            if near_zero_median < cc_min:
                logger.warning(
                    "Detected minimum seems larger than beginning values. "
                    "Consider increasing shift size (--max-shift)."
                )
        
        return cc_min
    
    def _calculate_cc_at_read_length(self, cc: np.ndarray) -> float:
        """Calculate cross-correlation at read length position.
        
        Args:
            cc: Cross-correlation array
            
        Returns:
            Cross-correlation value at read length
        """
        if self.read_length > len(cc):
            return 0.0
        return float(cc[self.read_length - 1])
    
    def _calculate_fragment_metrics(self, 
                                  cc: np.ndarray, 
                                  fragment_length: int) -> Tuple[float, float, float]:
        """Calculate fragment-based metrics (NSC, RSC).
        
        Args:
            cc: Cross-correlation array
            fragment_length: Fragment length position
            
        Returns:
            Tuple of (cc_at_fragment, nsc, rsc)
        """
        if fragment_length > len(cc):
            fragment_length = len(cc)
        
        cc_at_fragment = float(cc[fragment_length - 1])
        cc_min = self._result.cc_min if self._result else self._calculate_minimum(cc)
        cc_at_read = self._result.cc_at_read_length if self._result else self._calculate_cc_at_read_length(cc)
        
        # Calculate NSC and RSC with error handling
        try:
            nsc = cc_at_fragment / cc_min if cc_min > 0 else float('inf')
            
            denominator = cc_at_read - cc_min
            if denominator > 0:
                rsc = (cc_at_fragment - cc_min) / denominator
            else:
                rsc = float('inf')
                
        except (ZeroDivisionError, FloatingPointError):
            if self.config.output_warnings:
                logger.warning("Numerical issue in NSC/RSC calculation")
            nsc = rsc = float('inf')
        
        return cc_at_fragment, nsc, rsc
    
    def _estimate_fragment_length(self, cc: np.ndarray) -> Optional[int]:
        """Estimate fragment length from cross-correlation peak.
        
        Args:
            cc: Cross-correlation array
            
        Returns:
            Estimated fragment length or None if estimation fails
        """
        # Apply smoothing
        smoothed = moving_avr_filter(cc, self.config.mv_avr_filter_len)
        self._smoothed_cc = smoothed
        
        # Find initial peak
        estimated_length = int(np.argmax(smoothed) + 1)
        
        # Check if peak is too close to read length
        distance_to_read_length = abs(estimated_length - self.read_length)
        
        if distance_to_read_length <= self.config.filter_mask_len:
            if self.config.output_warnings:
                logger.warning(
                    f"Estimated fragment length ({estimated_length}) is close to "
                    f"read length ({self.read_length}). Applying mask..."
                )
            
            # Apply mask around read length
            masked_cc = smoothed.copy()
            mask_start = max(0, self.read_length - 1 - self.config.filter_mask_len)
            mask_end = min(len(masked_cc), self.read_length + self.config.filter_mask_len)
            
            masked_cc[mask_start:mask_end] = -float("inf")
            estimated_length = int(np.argmax(masked_cc) + 1)
            
            # Check if we're still at the mask boundary
            if estimated_length - 1 in (mask_start - 1, mask_end):
                if self.config.output_warnings:
                    logger.error(
                        "Estimated fragment length is still at mask boundary. "
                        "Check output plots for quality assessment."
                    )
        
        # Final validation
        if (self.config.output_warnings and 
            abs(estimated_length - self.read_length) <= 5):
            logger.error(
                "Estimated fragment length is very close to read length! "
                "Please check output plots."
            )
        
        return estimated_length
    
    def _calculate_peak_width(self, fragment_length: int) -> Union[int, bool]:
        """Calculate Full Width at Half Maximum (FWHM) of the peak.
        
        Args:
            fragment_length: Fragment length position
            
        Returns:
            Peak width in base pairs, or False if calculation fails
        """
        if self._smoothed_cc is None:
            self._smoothed_cc = moving_avr_filter(
                self.data.ccbins, self.config.mv_avr_filter_len
            )
        
        cc = self._smoothed_cc
        max_idx = fragment_length - 1
        
        if max_idx >= len(cc) or max_idx < 0:
            raise ValueError(f"Fragment length {fragment_length} out of bounds")
        
        cc_max = cc[max_idx]
        cc_min = self._result.cc_min if self._result else self._calculate_minimum(self.data.ccbins)
        
        if cc_max <= cc_min:
            raise ValueError("Peak value is not greater than minimum")
        
        target = cc_min + (cc_max - cc_min) / 2
        
        # Find forward half-width
        forward_shift = 0
        forward_failed = False
        while max_idx + forward_shift < len(cc) and cc[max_idx + forward_shift] > target:
            forward_shift += 1
        
        if max_idx + forward_shift >= len(cc):
            if self.config.output_warnings:
                logger.warning(
                    "Failed to calculate forward half-width. "
                    "Consider increasing shift size."
                )
            forward_failed = True
            forward_shift -= 1
        
        # Find backward half-width
        backward_shift = 0
        backward_failed = False
        while max_idx - backward_shift >= 0 and cc[max_idx - backward_shift] > target:
            backward_shift += 1
        
        if max_idx - backward_shift < 0:
            if self.config.output_warnings:
                logger.warning("Failed to calculate backward half-width.")
            backward_failed = True
            backward_shift -= 1
        
        # Determine final width
        if forward_failed and backward_failed:
            if self.config.output_warnings:
                logger.error("Failed to calculate full width at half maximum.")
            return False
        elif forward_failed:
            if self.config.output_warnings:
                logger.warning("Using doubled backward half-width")
            return backward_shift * 2 + 1
        elif backward_failed:
            if self.config.output_warnings:
                logger.warning("Using doubled forward half-width")
            return forward_shift * 2 + 1
        else:
            return backward_shift + forward_shift + 1
    
    def _calculate_vsn(self, cc_at_fragment: float, peak_width: int) -> float:
        """Calculate Variance Stabilizing Normalization factor.
        
        Args:
            cc_at_fragment: Cross-correlation at fragment length
            peak_width: Full width at half maximum
            
        Returns:
            VSN factor
        """
        if isinstance(self.data, NCCData):
            total_reads = self.data.forward_sum + self.data.reverse_sum
        else:  # MSCCData
            # For MSCC, use sum of arrays
            total_reads = np.sum(self.data.forward_sum) + np.sum(self.data.reverse_sum)
        
        if total_reads == 0:
            return 0.0
        
        return 2 * cc_at_fragment * peak_width / total_reads
    
    @property
    def smoothed_correlation(self) -> Optional[np.ndarray]:
        """Get smoothed cross-correlation data if available."""
        return self._smoothed_cc


# Compatibility function for legacy code
def create_legacy_ccstats(data: Union[NCCData, MSCCData],
                         read_length: int,
                         min_calc_width: int = 10,
                         mv_avr_filter_len: int = 15,
                         filter_mask_len: int = 5,
                         output_warnings: bool = True,
                         do_fragment_estimation: bool = False,
                         estimated_fragment_length: Optional[int] = None,
                         expected_fragment_length: Optional[int] = None) -> StatisticsResult:
    """Create statistics using legacy parameter interface.
    
    This function provides backward compatibility with the original
    CCStats initialization interface.
    
    Args:
        data: Cross-correlation data
        read_length: Read length in base pairs
        min_calc_width: Width for minimum calculation
        mv_avr_filter_len: Moving average filter length
        filter_mask_len: Mask length around read length
        output_warnings: Whether to output warnings
        do_fragment_estimation: Whether to estimate fragment length
        estimated_fragment_length: Pre-computed estimated length
        expected_fragment_length: User-specified expected length
        
    Returns:
        Calculated statistics result
    """
    config = StatisticsConfig(
        min_calc_width=min_calc_width,
        mv_avr_filter_len=mv_avr_filter_len,
        filter_mask_len=filter_mask_len,
        output_warnings=output_warnings,
        do_fragment_estimation=do_fragment_estimation,
        expected_fragment_length=expected_fragment_length
    )
    
    calculator = StatisticsCalculator(data, read_length, config)
    return calculator.calculate_statistics()