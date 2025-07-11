"""Unified statistics module for PyMaSC cross-correlation analysis.

This module consolidates the statistical calculation functionality from CCStats
and PyMaSCStats into a unified, maintainable architecture. It provides:

1. BaseCorrelationStats: Common statistical calculations
2. NCCStats: Naive cross-correlation statistics
3. MSCCStats: Mappability-sensitive cross-correlation statistics
4. UnifiedStats: Combined NCC and MSCC analysis

Key improvements:
- Eliminates code duplication between CCStats and PyMaSCStats
- Provides cleaner separation of concerns
- Maintains backward compatibility
- Reduces overall codebase size

Mathematical background:
- NSC (Normalized Strand Coefficient): Signal-to-noise ratio
- RSC (Relative Strand Coefficient): Enrichment relative to phantom peak
- Fragment length estimation using moving average filtering
- FWHM (Full Width at Half Maximum) for peak characterization
"""

import logging
from abc import ABC, abstractmethod
from functools import wraps
from typing import Union, Tuple, Optional, Dict, Any

import numpy as np
from scipy.stats import norm

from PyMaSC.utils.calc import moving_avr_filter

logger = logging.getLogger(__name__)

# Constants
NEAR_READLEN_ERR_CRITERION = 5
NEAR_ZERO_MIN_CALC_LEN = 10


def npcalc_with_logging_warn(func):
    """Decorator for handling numpy floating point errors gracefully."""
    @wraps(func)
    def _inner(*args, **kwargs):
        try:
            with np.errstate(divide="raise", invalid="raise"):
                return func(*args, **kwargs)
        except FloatingPointError as e:
            logger.debug("catch numpy warning: " + repr(e))
            logger.debug("continue anyway.")
            with np.errstate(divide="ignore", invalid="ignore"):
                return func(*args, **kwargs)
    return _inner


class BaseCorrelationStats(ABC):
    """Base class for cross-correlation statistical calculations.
    
    This class provides common functionality for both NCC and MSCC statistics,
    including background calculation, quality metrics, and fragment length estimation.
    
    Attributes:
        cc (np.ndarray): Cross-correlation values
        genomelen (int): Genome length for normalization
        forward_sum (int): Total forward reads
        reverse_sum (int): Total reverse reads
        read_len (int): Read length
        cc_min (float): Background correlation level
        ccrl (float): Correlation at read length (phantom peak)
        ccfl (float): Correlation at fragment length
        nsc (float): Normalized Strand Coefficient
        rsc (float): Relative Strand Coefficient
        cc_width (int): FWHM of fragment length peak
        vsn (float): Variance Stabilizing Normalization
    """
    
    def __init__(self, cc: np.ndarray, genomelen: int, forward_sum: int, reverse_sum: int,
                 read_len: int, min_calc_width: int, mv_avr_filter_len: int,
                 filter_mask_len: int, output_warnings: bool = True,
                 estimated_library_len: Optional[int] = None,
                 expected_library_len: Optional[int] = None):
        """Initialize base correlation statistics.
        
        Args:
            cc: Cross-correlation array
            genomelen: Genome length for normalization
            forward_sum: Total forward reads
            reverse_sum: Total reverse reads
            read_len: Read length
            min_calc_width: Width for background calculation
            mv_avr_filter_len: Moving average filter length
            filter_mask_len: Masking window around read length
            output_warnings: Whether to output warnings
            estimated_library_len: Estimated fragment length
            expected_library_len: Expected fragment length
        """
        # Ensure cc is a numpy array
        self.cc = np.array(cc, dtype=np.float_) if cc is not None else None
        self.genomelen = genomelen
        self.forward_sum = forward_sum
        self.reverse_sum = reverse_sum
        self.read_len = read_len
        self.min_calc_width = min_calc_width
        self.mv_avr_filter_len = mv_avr_filter_len
        self.filter_mask_len = filter_mask_len
        self.output_warnings = output_warnings
        self.est_lib_len = estimated_library_len
        self.library_len = expected_library_len
        
        # Set up logging functions
        if self.output_warnings:
            self.error = logger.error
            self.warning = logger.warning
        else:
            self.error = self.warning = logger.debug
        
        # Calculate basic statistics
        self._calc_cc_min()
        self._calc_cc_rl()
        
        # Initialize quality metrics
        self.ccfl = self.nsc = self.rsc = None
        self.est_ccfl = self.est_nsc = self.est_rsc = None
        self.cc_width = self.est_cc_width = None
        self.vsn = self.est_vsn = None
        self.avr_cc = None
        
        # Calculate metrics for expected library length
        if self.library_len:
            self.ccfl, self.nsc, self.rsc = self._calc_lib_metrics(self.library_len)
            try:
                self.cc_width = self._get_FWHM(self.library_len)
                self.vsn = self._calc_vsn(self.ccfl, self.cc_width)
            except AssertionError as e:
                self.warning(e)
        
        # Calculate metrics for estimated library length
        if self.est_lib_len:
            self.est_ccfl, self.est_nsc, self.est_rsc = self._calc_lib_metrics(self.est_lib_len)
            try:
                self.est_cc_width = self._get_FWHM(self.est_lib_len)
                self.est_vsn = self._calc_vsn(self.est_ccfl, self.est_cc_width)
            except AssertionError as e:
                self.warning(e)
    
    def _calc_cc_min(self) -> None:
        """Calculate minimum cross-correlation value for background estimation."""
        self.cc_min = np.sort(self.cc[-self.min_calc_width:])[
            min(self.min_calc_width, self.cc.size) // 2
        ]
        if (np.median(self.cc[:NEAR_ZERO_MIN_CALC_LEN]) < self.cc_min and 
            self.output_warnings):
            logger.warning(
                "Detected minimum coefficient seems to be larger than beginning part minimum. "
                "Consider increasing shift size (-d/--max-shift)."
            )
    
    def _calc_cc_rl(self) -> None:
        """Calculate cross-correlation at read length (phantom peak)."""
        if self.read_len > self.cc.size:
            self.ccrl = 0
        else:
            self.ccrl = self.cc[self.read_len - 1]
    
    @npcalc_with_logging_warn
    def _calc_lib_metrics(self, library_len: int) -> Tuple[float, float, float]:
        """Calculate library-specific quality metrics (NSC and RSC)."""
        # Add bounds checking to prevent IndexError
        if library_len < 1 or library_len > len(self.cc):
            if self.output_warnings:
                self.warning(f"Library length {library_len} is out of bounds for correlation array of size {len(self.cc)}")
            return 0.0, 0.0, 0.0
        
        ccfl = self.cc[library_len - 1]
        nsc = ccfl / self.cc_min
        rsc = (ccfl - self.cc_min) / (self.ccrl - self.cc_min)
        return ccfl, nsc, rsc
    
    def _calc_vsn(self, ccfl: float, cc_width: int) -> float:
        """Calculate Variance Stabilizing Normalization value."""
        if ccfl is None or cc_width is None:
            return None
        return 2 * ccfl * cc_width / (self.forward_sum + self.reverse_sum)
    
    def _get_FWHM(self, library_len: int) -> Union[int, bool]:
        """Calculate Full Width at Half Maximum (FWHM) of fragment length peak."""
        if self.avr_cc is None:
            self.avr_cc = moving_avr_filter(self.cc, self.mv_avr_filter_len)
        
        max_i = library_len - 1
        assert max_i >= 0  # Original implementation check
        cc_max = self.avr_cc[max_i - 1]  # Original implementation uses max_i - 1  
        assert cc_max > self.cc_min  # Original implementation check
        
        target = self.cc_min + (cc_max - self.cc_min) / 2
        
        # Forward direction
        forward_shift = 0
        forward_failed = False
        while self.avr_cc[max_i + forward_shift] > target:
            forward_shift += 1
            if max_i + forward_shift == self.avr_cc.size:
                self.warning(
                    "Failed to calc the half width at half maximum in the forward "
                    "side of the peak. Consider increasing shift size (-d/--max-shift)."
                )
                forward_failed = True
                forward_shift -= 1
                break
        
        # Backward direction
        backward_shift = 0
        backward_failed = False
        while self.avr_cc[max_i - backward_shift] > target:
            backward_shift += 1
            if max_i < backward_shift:
                self.warning(
                    "Failed to calc the half width at half maximum in the backward side of the peak."
                )
                backward_failed = True
                backward_shift -= 1
                break
        
        logger.debug((forward_shift, backward_shift))
        
        if forward_failed and backward_failed:
            self.error("Failed to calcurate the full width at half maximum.")
            return False
        elif forward_failed:
            self.warning("Use twice width of the half width at half maximum in the backward side")
            return backward_shift * 2 + 1
        elif backward_failed:
            self.warning("Use twice width of the half width at half maximum in the forward side")
            return forward_shift * 2 + 1
        else:
            return backward_shift + forward_shift + 1
    
    def estimate_fragment_length(self) -> None:
        """Estimate fragment length from cross-correlation peak."""
        self.avr_cc = moving_avr_filter(self.cc, self.mv_avr_filter_len)
        self.est_lib_len = np.argmax(self.avr_cc) + 1
        
        need_warning = False
        
        # Handle phantom peak masking
        if (self.filter_mask_len and 
            abs(self.est_lib_len - self.read_len) <= self.filter_mask_len):
            self.warning("Estimated library length is close to the read length.")
            self.warning(
                f"Trying to masking around the read length +/- {self.filter_mask_len}bp..."
            )
            
            _avr_cc = self.avr_cc.copy()
            mask_from = max(0, self.read_len - 1 - self.filter_mask_len)
            mask_to = min(len(_avr_cc), self.read_len + self.filter_mask_len)
            for i in range(mask_from, mask_to):
                _avr_cc[i] = -float("inf")
            self.est_lib_len = np.argmax(_avr_cc) + 1
            
            if self.est_lib_len - 1 in (mask_from - 1, mask_to):
                need_warning = True
        
        elif (self.output_warnings and 
              abs(self.est_lib_len - self.read_len) <= NEAR_READLEN_ERR_CRITERION):
            need_warning = True
        
        if self.output_warnings and need_warning:
            logger.error(
                "Estimated library length is close to the read length! "
                "Please check output plots."
            )
        
        # Update estimated metrics
        if self.est_lib_len:
            self.est_ccfl, self.est_nsc, self.est_rsc = self._calc_lib_metrics(self.est_lib_len)
            try:
                self.est_cc_width = self._get_FWHM(self.est_lib_len)
                self.est_vsn = self._calc_vsn(self.est_ccfl, self.est_cc_width)
            except AssertionError as e:
                self.warning(e)


class NCCStats(BaseCorrelationStats):
    """Naive Cross-Correlation statistics calculator.
    
    This class handles standard cross-correlation analysis without
    mappability correction. It extends BaseCorrelationStats with
    NCC-specific functionality.
    """
    
    def __init__(self, cc: np.ndarray, genomelen: int, forward_sum: int, reverse_sum: int,
                 read_len: int, min_calc_width: int, mv_avr_filter_len: int,
                 filter_mask_len: int, output_warnings: bool = True,
                 do_llestimation: bool = False, estimated_library_len: Optional[int] = None,
                 expected_library_len: Optional[int] = None):
        """Initialize NCC statistics calculator.
        
        Args:
            cc: Cross-correlation array
            genomelen: Genome length for normalization
            forward_sum: Total forward reads
            reverse_sum: Total reverse reads
            read_len: Read length
            min_calc_width: Width for background calculation
            mv_avr_filter_len: Moving average filter length
            filter_mask_len: Masking window around read length
            output_warnings: Whether to output warnings
            do_llestimation: Whether to perform library length estimation
            estimated_library_len: Pre-estimated fragment length
            expected_library_len: Expected fragment length
        """
        super().__init__(
            cc, genomelen, forward_sum, reverse_sum, read_len,
            min_calc_width, mv_avr_filter_len, filter_mask_len,
            output_warnings, estimated_library_len, expected_library_len
        )
        
        # Perform fragment length estimation if requested
        if do_llestimation:
            self.estimate_fragment_length()


class MSCCStats(BaseCorrelationStats):
    """Mappability-Sensitive Cross-Correlation statistics calculator.
    
    This class handles cross-correlation analysis with mappability correction.
    It extends BaseCorrelationStats with MSCC-specific functionality.
    """
    
    def __init__(self, cc: np.ndarray, mappable_len: np.ndarray, 
                 forward_sum: np.ndarray, reverse_sum: np.ndarray,
                 read_len: int, min_calc_width: int, mv_avr_filter_len: int,
                 filter_mask_len: int, output_warnings: bool = True,
                 do_llestimation: bool = False, estimated_library_len: Optional[int] = None,
                 expected_library_len: Optional[int] = None):
        """Initialize MSCC statistics calculator.
        
        Args:
            cc: Cross-correlation array
            mappable_len: Mappable genome length array
            forward_sum: Mappable forward reads array
            reverse_sum: Mappable reverse reads array
            read_len: Read length
            min_calc_width: Width for background calculation
            mv_avr_filter_len: Moving average filter length
            filter_mask_len: Masking window around read length
            output_warnings: Whether to output warnings
            do_llestimation: Whether to perform library length estimation
            estimated_library_len: Pre-estimated fragment length
            expected_library_len: Expected fragment length
        """
        # Convert arrays to numpy arrays and ensure proper types
        mappable_len = np.array(mappable_len, dtype=np.float_) if mappable_len is not None else None
        forward_sum = np.array(forward_sum, dtype=np.float_) if forward_sum is not None else None
        reverse_sum = np.array(reverse_sum, dtype=np.float_) if reverse_sum is not None else None
        
        # For MSCC, use mappable_len as the effective genome length
        # Use the first value for initialization, full array is stored separately
        effective_genomelen = mappable_len[0] if mappable_len is not None and len(mappable_len) > 0 else 1
        effective_forward_sum = forward_sum[0] if forward_sum is not None and len(forward_sum) > 0 else 0
        effective_reverse_sum = reverse_sum[0] if reverse_sum is not None and len(reverse_sum) > 0 else 0
        
        super().__init__(
            cc, int(effective_genomelen), int(effective_forward_sum), int(effective_reverse_sum),
            read_len, min_calc_width, mv_avr_filter_len, filter_mask_len,
            output_warnings, estimated_library_len, expected_library_len
        )
        
        # Store MSCC-specific arrays
        self.mappable_len = mappable_len
        self.mappable_forward_sum = forward_sum
        self.mappable_reverse_sum = reverse_sum
        
        # Override scalar attributes with arrays for backward compatibility
        # The output module expects these to be arrays indexed by shift distance
        self.genomelen = mappable_len  # Array of mappable lengths
        self.forward_sum = forward_sum  # Array of mappable forward sums
        self.reverse_sum = reverse_sum  # Array of mappable reverse sums
        
        # Perform fragment length estimation if requested
        if do_llestimation:
            self.estimate_fragment_length()


class UnifiedStats:
    """Unified statistics calculator for both NCC and MSCC analysis.
    
    This class combines NCC and MSCC statistics into a single interface,
    providing compatibility with the existing PyMaSCStats API while
    reducing code duplication.
    """
    
    def __init__(self, read_len: int, mv_avr_filter_len: int = 15,
                 expected_library_len: Optional[int] = None,
                 # NCC parameters
                 genomelen: Optional[int] = None,
                 forward_sum: Optional[int] = None,
                 reverse_sum: Optional[int] = None,
                 ccbins: Optional[np.ndarray] = None,
                 cc: Optional[np.ndarray] = None,
                 # MSCC parameters
                 mappable_len: Optional[np.ndarray] = None,
                 mappable_forward_sum: Optional[np.ndarray] = None,
                 mappable_reverse_sum: Optional[np.ndarray] = None,
                 mappable_ccbins: Optional[np.ndarray] = None,
                 masc: Optional[np.ndarray] = None,
                 # Configuration
                 filter_mask_len: int = NEAR_READLEN_ERR_CRITERION,
                 warning: bool = False,
                 min_calc_width: Optional[int] = None):
        """Initialize unified statistics calculator.
        
        Args:
            read_len: Read length
            mv_avr_filter_len: Moving average filter length
            expected_library_len: Expected fragment length
            genomelen: Genome length for NCC
            forward_sum: Total forward reads for NCC
            reverse_sum: Total reverse reads for NCC
            ccbins: Cross-correlation bins for NCC
            cc: Pre-computed NCC array
            mappable_len: Mappable genome length array
            mappable_forward_sum: Mappable forward reads array
            mappable_reverse_sum: Mappable reverse reads array
            mappable_ccbins: Mappable cross-correlation bins
            masc: Pre-computed MSCC array
            filter_mask_len: Masking window around read length
            warning: Whether to output warnings
            min_calc_width: Width for background calculation
        """
        self.read_len = read_len
        self.mv_avr_filter_len = mv_avr_filter_len
        self.library_len = expected_library_len
        self.filter_mask_len = filter_mask_len
        self.output_warnings = warning
        self.min_calc_width = min_calc_width or 10
        
        # Initialize statistics objects
        self.cc = None
        self.masc = None
        self.est_lib_len = None
        
        # Determine what calculations are possible
        self.calc_ncc = self._can_calc_ncc(genomelen, forward_sum, reverse_sum, ccbins, cc)
        self.calc_masc = self._can_calc_masc(mappable_len, mappable_forward_sum, 
                                           mappable_reverse_sum, mappable_ccbins, masc)
        
        if not self.calc_ncc and not self.calc_masc:
            return
        
        # Calculate from bins if provided
        if ccbins is not None or mappable_ccbins is not None:
            self._calc_from_bins(genomelen, forward_sum, reverse_sum, ccbins,
                               mappable_len, mappable_forward_sum, mappable_reverse_sum,
                               mappable_ccbins)
        else:
            # Calculate from pre-computed arrays
            self._calc_from_arrays(genomelen, forward_sum, reverse_sum, cc,
                                 mappable_len, mappable_forward_sum, mappable_reverse_sum, masc)
    
    def _can_calc_ncc(self, genomelen, forward_sum, reverse_sum, ccbins, cc):
        """Check if NCC calculation is possible."""
        if ccbins is not None:
            return all(x is not None for x in (forward_sum, reverse_sum, ccbins, genomelen))
        else:
            return all(x is not None for x in (forward_sum, reverse_sum, genomelen, cc))
    
    def _can_calc_masc(self, mappable_len, mappable_forward_sum, mappable_reverse_sum, 
                      mappable_ccbins, masc):
        """Check if MSCC calculation is possible."""
        if mappable_ccbins is not None:
            return all(x is not None for x in (mappable_forward_sum, mappable_reverse_sum,
                                             mappable_ccbins, mappable_len))
        else:
            return all(x is not None for x in (mappable_forward_sum, mappable_reverse_sum,
                                             mappable_len, masc))
    
    def _calc_from_bins(self, genomelen, forward_sum, reverse_sum, ccbins,
                       mappable_len, mappable_forward_sum, mappable_reverse_sum, mappable_ccbins):
        """Calculate statistics from raw count bins."""
        # Determine maximum shift distance
        if self.calc_ncc:
            ncc_max_shift = len(ccbins) - 1
        else:
            ncc_max_shift = None
        
        if self.calc_masc:
            masc_max_shift = min(map(len, (mappable_forward_sum, mappable_reverse_sum, 
                                         mappable_ccbins))) - 1
        else:
            masc_max_shift = None
        
        if ncc_max_shift and masc_max_shift:
            max_shift = min(ncc_max_shift, masc_max_shift)
        elif ncc_max_shift:
            max_shift = ncc_max_shift
        elif masc_max_shift:
            max_shift = masc_max_shift
        else:
            max_shift = 0
        
        self.max_shift = max_shift
        
        # Calculate correlation arrays
        cc_array = self._calc_naive_cc(genomelen, forward_sum, reverse_sum, ccbins, max_shift)
        masc_array = self._calc_masc_cc(mappable_len, mappable_forward_sum, mappable_reverse_sum,
                                       mappable_ccbins, max_shift)
        
        # Create statistics objects
        self._create_stats_objects(genomelen, forward_sum, reverse_sum, cc_array,
                                 mappable_len, mappable_forward_sum, mappable_reverse_sum, masc_array)
    
    def _calc_from_arrays(self, genomelen, forward_sum, reverse_sum, cc,
                         mappable_len, mappable_forward_sum, mappable_reverse_sum, masc):
        """Calculate statistics from pre-computed correlation arrays."""
        if self.calc_ncc and self.calc_masc:
            assert len(cc) == len(masc), "NCC and MASC arrays must have same length"
            self.max_shift = min(len(cc), len(masc)) - 1
        elif self.calc_ncc:
            self.max_shift = len(cc) - 1
        elif self.calc_masc:
            self.max_shift = len(masc) - 1
        
        # Create statistics objects
        self._create_stats_objects(genomelen, forward_sum, reverse_sum, cc,
                                 mappable_len, mappable_forward_sum, mappable_reverse_sum, masc)
    
    def _create_stats_objects(self, genomelen, forward_sum, reverse_sum, cc,
                            mappable_len, mappable_forward_sum, mappable_reverse_sum, masc):
        """Create NCC and MSCC statistics objects."""
        # Create MSCC first to estimate fragment length
        if self.calc_masc:
            self.masc = MSCCStats(
                cc=masc,
                mappable_len=mappable_len,
                forward_sum=mappable_forward_sum,
                reverse_sum=mappable_reverse_sum,
                read_len=self.read_len,
                min_calc_width=self.min_calc_width,
                mv_avr_filter_len=self.mv_avr_filter_len,
                filter_mask_len=self.filter_mask_len,
                output_warnings=self.output_warnings,
                do_llestimation=True,
                expected_library_len=self.library_len
            )
            self.est_lib_len = self.masc.est_lib_len
        
        # Create NCC with estimated fragment length from MSCC
        if self.calc_ncc:
            self.cc = NCCStats(
                cc=cc,
                genomelen=genomelen,
                forward_sum=forward_sum,
                reverse_sum=reverse_sum,
                read_len=self.read_len,
                min_calc_width=self.min_calc_width,
                mv_avr_filter_len=self.mv_avr_filter_len,
                filter_mask_len=self.filter_mask_len,
                output_warnings=self.output_warnings,
                do_llestimation=not self.calc_masc,  # Only estimate if no MSCC
                estimated_library_len=self.est_lib_len,
                expected_library_len=self.library_len
            )
            
            # If no MSCC, use NCC fragment length estimation
            if not self.calc_masc:
                self.est_lib_len = self.cc.est_lib_len
    
    @staticmethod
    @npcalc_with_logging_warn
    def _calc_cc(forward_sum, reverse_sum, ccbins, totlen, denom):
        """Calculate normalized cross-correlation using binomial variance model."""
        forward_mean = forward_sum / totlen
        reverse_mean = reverse_sum / totlen
        
        forward_var = forward_mean * (1 - forward_mean)
        reverse_var = reverse_mean * (1 - reverse_mean)
        
        sum_prod = forward_mean * reverse_mean
        var_geomean = (forward_var * reverse_var) ** 0.5
        return (ccbins / denom - sum_prod) / var_geomean
    
    def _calc_naive_cc(self, genomelen, forward_sum, reverse_sum, ccbins, max_shift):
        """Calculate standard normalized cross-correlation (NCC)."""
        if not self.calc_ncc:
            return None
        
        denom = genomelen - np.array(range(max_shift + 1), dtype=np.float_)
        return self._calc_cc(
            float(forward_sum),
            float(reverse_sum),
            ccbins[:max_shift + 1],
            genomelen,
            denom
        )
    
    @npcalc_with_logging_warn
    def _calc_masc_cc(self, mappable_len, mappable_forward_sum, mappable_reverse_sum,
                     mappable_ccbins, max_shift):
        """Calculate mappability-sensitive cross-correlation (MSCC)."""
        if not self.calc_masc:
            return None
        
        totlen = np.array(mappable_len, dtype=np.float_)
        totlen = np.concatenate((
            totlen[:self.read_len][::-1], totlen[1:]
        ))[:max_shift + 1]
        
        return self._calc_cc(
            np.array(mappable_forward_sum[:max_shift + 1], dtype=np.float_),
            np.array(mappable_reverse_sum[:max_shift + 1], dtype=np.float_),
            mappable_ccbins[:max_shift + 1],
            totlen,
            totlen
        )