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
from typing import Union, Tuple, Optional, Dict, Any, List
from dataclasses import dataclass

import numpy as np
from scipy.stats import norm

from PyMaSC.utils.calc import moving_avr_filter
from PyMaSC.core.interfaces.result import GenomeWideResult, NCCGenomeWideResult, MSCCGenomeWideResult, BothGenomeWideResult

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


@dataclass
class StatisticsResult:
    """New statistics result class for clean architecture without compatibility layers.

    This class replaces CCResult and uses only the new statistics system from
    statistics.py, eliminating legacy compatibility code entirely.

    Attributes:
        chromosome_stats: Per-chromosome UnifiedStats objects
        genome_wide_stats: Genome-wide aggregated UnifiedStats
        references: List of chromosome names
        read_length: Read length in base pairs
        skip_ncc: Whether NCC calculation was skipped
        has_mappability: Whether mappability data is available
        nsc: Normalized Strand Coefficient (signal-to-noise ratio)
        rsc: Relative Strand Coefficient (fragment-to-phantom peak ratio)
        estimated_fragment_length: Estimated fragment length from correlation peak
    """
    chromosome_stats: Dict[str, UnifiedStats]
    genome_wide_stats: UnifiedStats
    references: List[str]
    read_length: int
    skip_ncc: bool
    has_mappability: bool
    nsc: float
    rsc: float
    estimated_fragment_length: Optional[int]

    @classmethod
    def from_handler(
        cls,
        handler: Any,  # UnifiedCalcHandler
        mv_avr_filter_len: int = 15,
        expected_library_len: Optional[int] = None,
        filter_mask_len: int = 5,
        min_calc_width: int = 10
    ) -> 'StatisticsResult':
        """Create StatisticsResult from UnifiedCalcHandler using new statistics system.

        This method bypasses all legacy compatibility layers and uses only the
        new container-based data access and UnifiedStats for calculations.

        Args:
            handler: UnifiedCalcHandler with completed calculations
            mv_avr_filter_len: Moving average filter length for peak detection
            expected_library_len: Expected fragment length (optional)
            filter_mask_len: Masking window around read length for phantom peak
            min_calc_width: Minimum width for background correlation calculation

        Returns:
            StatisticsResult with complete statistical analysis
        """
        # Get container-based results (no legacy attributes)
        aggregation_result = handler.get_aggregation_result()
        if aggregation_result is None:
            raise ValueError("No calculation results available from handler")

        container = aggregation_result.container
        if container is None:
            raise ValueError("No container results available")

        # Calculate per-chromosome statistics using UnifiedStats
        chromosome_stats = {}
        for ref in handler.references:
            chrom_result = container.chromosome_results.get(ref)
            if chrom_result is None:
                continue

            # Extract NCC data (raw ccbins, same as CCResult._init_from_handler)
            ncc_data = chrom_result.ncc_data
            genomelen = chrom_result.length if chrom_result else 1
            forward_sum = ncc_data.forward_sum if ncc_data else 0
            reverse_sum = ncc_data.reverse_sum if ncc_data else 0
            ccbins = ncc_data.ccbins if ncc_data else None

            # Extract MSCC data (raw ccbins, same as CCResult._init_from_handler)
            mscc_data = chrom_result.mscc_data
            mappable_len = mscc_data.mappable_len if mscc_data else None
            mappable_forward_sum = mscc_data.forward_sum if mscc_data else None
            mappable_reverse_sum = mscc_data.reverse_sum if mscc_data else None
            mappable_ccbins = mscc_data.ccbins if mscc_data else None

            # Create UnifiedStats for this chromosome (same parameters as CCResult._init_from_handler)
            chromosome_stats[ref] = UnifiedStats(
                read_len=handler.read_len,
                mv_avr_filter_len=mv_avr_filter_len,
                expected_library_len=expected_library_len,
                genomelen=genomelen,
                forward_sum=forward_sum,
                reverse_sum=reverse_sum,
                ccbins=ccbins,  # Raw ccbins, UnifiedStats will normalize internally
                mappable_len=mappable_len,
                mappable_forward_sum=mappable_forward_sum,
                mappable_reverse_sum=mappable_reverse_sum,
                mappable_ccbins=mappable_ccbins,  # Raw mappable ccbins
                filter_mask_len=filter_mask_len,
                min_calc_width=min_calc_width
            )

        # Calculate genome-wide statistics
        genome_wide_stats = StatisticsResult._calculate_genome_wide_stats(
            container, handler, mv_avr_filter_len,
            expected_library_len, filter_mask_len, min_calc_width
        )

        # Extract quality metrics
        nsc = 0.0
        rsc = 0.0
        estimated_fragment_length = None

        if genome_wide_stats.cc:
            nsc = getattr(genome_wide_stats.cc, 'nsc', 0.0) or 0.0
            rsc = getattr(genome_wide_stats.cc, 'rsc', 0.0) or 0.0

        estimated_fragment_length = getattr(genome_wide_stats, 'est_lib_len', None)

        return cls(
            chromosome_stats=chromosome_stats,
            genome_wide_stats=genome_wide_stats,
            references=handler.references,
            read_length=handler.read_len,
            skip_ncc=container.skip_ncc,
            has_mappability=container.has_mappability,
            nsc=nsc,
            rsc=rsc,
            estimated_fragment_length=estimated_fragment_length
        )

    @staticmethod
    def _calculate_genome_wide_stats(
        container: Any,  # ResultContainer
        handler: Any,    # UnifiedCalcHandler
        mv_avr_filter_len: int,
        expected_library_len: Optional[int],
        filter_mask_len: int,
        min_calc_width: int
    ) -> UnifiedStats:
        """Calculate genome-wide statistics from container data.

        Args:
            container: ResultContainer with chromosome results
            handler: UnifiedCalcHandler for metadata
            mv_avr_filter_len: Moving average filter length
            expected_library_len: Expected fragment length
            filter_mask_len: Phantom peak masking window
            min_calc_width: Background calculation width

        Returns:
            UnifiedStats object with genome-wide statistics
        """
        # Get aggregated data from container
        all_ncc_data = container.get_ncc_data()
        all_mscc_data = container.get_mscc_data()

        # Calculate totals
        total_genomelen = sum(handler.lengths)
        total_forward = sum(data.forward_sum for data in all_ncc_data.values()) if all_ncc_data else 0
        total_reverse = sum(data.reverse_sum for data in all_ncc_data.values()) if all_ncc_data else 0

        # Create per-chromosome UnifiedStats to calculate correlation arrays (same as CCResult)
        ref2stats = {}
        for ref in handler.references:
            chrom_result = container.chromosome_results.get(ref)
            if chrom_result is None:
                continue

            # Extract data for this chromosome (raw ccbins, same as CCResult)
            ncc_data = chrom_result.ncc_data
            mscc_data = chrom_result.mscc_data

            ref2stats[ref] = UnifiedStats(
                read_len=handler.read_len,
                mv_avr_filter_len=mv_avr_filter_len,
                expected_library_len=expected_library_len,
                genomelen=chrom_result.length,
                forward_sum=ncc_data.forward_sum if ncc_data else 0,
                reverse_sum=ncc_data.reverse_sum if ncc_data else 0,
                ccbins=ncc_data.ccbins if ncc_data else None,
                mappable_len=mscc_data.mappable_len if mscc_data else None,
                mappable_forward_sum=mscc_data.forward_sum if mscc_data else None,
                mappable_reverse_sum=mscc_data.reverse_sum if mscc_data else None,
                mappable_ccbins=mscc_data.ccbins if mscc_data else None,
                filter_mask_len=filter_mask_len,
                min_calc_width=min_calc_width
            )

        # Aggregate correlation arrays using Fisher z-transformation (same as CCResult lines 158-186)
        merged_cc = None
        if not container.skip_ncc:
            # Extract computed correlation arrays from UnifiedStats objects (same as CCResult)
            ncc_data = [(handler.lengths[i], ref2stats[ref].cc.cc)
                       for i, ref in enumerate(handler.references)
                       if ref in ref2stats and ref2stats[ref].cc is not None]
            if ncc_data:
                genome_lengths, correlation_arrays = zip(*ncc_data)
                from PyMaSC.core.constants import MERGED_CC_CONFIDENCE_INTERVAL
                from PyMaSC.utils.stats_utils import CrossCorrelationMerger
                merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
                merged_cc, _, _ = merger.merge_correlations(
                    np.array(genome_lengths), correlation_arrays, handler.read_len
                )

        merged_mscc = None
        if container.has_mappability:
            # Extract computed MSCC arrays from UnifiedStats objects (same as CCResult)
            # Need to get mappable lengths correctly
            mscc_data = []
            for ref in handler.references:
                if ref in ref2stats and ref2stats[ref].masc is not None:
                    chrom_result = container.chromosome_results.get(ref)
                    if chrom_result and chrom_result.mscc_data and chrom_result.mscc_data.mappable_len is not None:
                        # Use the final mappable length value (at read_len position)
                        mappable_len_value = chrom_result.mscc_data.mappable_len[handler.read_len - 1] if len(chrom_result.mscc_data.mappable_len) > handler.read_len - 1 else chrom_result.mscc_data.mappable_len[-1]
                        mscc_data.append((mappable_len_value, ref2stats[ref].masc.cc))

            if mscc_data:
                mappable_lengths, correlation_arrays = zip(*mscc_data)
                merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
                merged_mscc, _, _ = merger.merge_correlations(
                    np.array(mappable_lengths), correlation_arrays, handler.read_len
                )

        # Aggregate mappability data
        total_mappable_len = None
        total_mappable_forward = None
        total_mappable_reverse = None

        if all_mscc_data:
            from PyMaSC.utils.stats_utils import ArrayAggregator
            mappable_lens = [data.mappable_len for data in all_mscc_data.values() if data.mappable_len is not None]
            mappable_forwards = [data.forward_sum for data in all_mscc_data.values() if data.forward_sum is not None]
            mappable_reverses = [data.reverse_sum for data in all_mscc_data.values() if data.reverse_sum is not None]

            if mappable_lens:
                total_mappable_len = ArrayAggregator.safe_array_sum(mappable_lens)
            if mappable_forwards:
                total_mappable_forward = ArrayAggregator.safe_array_sum(mappable_forwards)
            if mappable_reverses:
                total_mappable_reverse = ArrayAggregator.safe_array_sum(mappable_reverses)

        return UnifiedStats(
            read_len=handler.read_len,
            mv_avr_filter_len=mv_avr_filter_len,
            expected_library_len=expected_library_len,
            genomelen=total_genomelen,
            forward_sum=total_forward,
            reverse_sum=total_reverse,
            cc=merged_cc,  # Pre-computed correlation array from Fisher z-transformation
            mappable_len=total_mappable_len,
            mappable_forward_sum=total_mappable_forward,
            mappable_reverse_sum=total_mappable_reverse,
            masc=merged_mscc,  # Pre-computed MSCC array from Fisher z-transformation
            filter_mask_len=filter_mask_len,
            warning=True,
            min_calc_width=min_calc_width
        )

    @classmethod
    def from_file_data(
        cls,
        mv_avr_filter_len: int,
        chi2_pval: float,
        filter_mask_len: int,
        min_calc_width: int,
        expected_library_len: Optional[int] = None,
        read_len: Optional[int] = None,
        references: Optional[List[str]] = None,
        ref2genomelen: Optional[Dict[str, int]] = None,
        ref2forward_sum: Optional[Dict[str, int]] = None,
        ref2reverse_sum: Optional[Dict[str, int]] = None,
        ref2cc: Optional[Dict[str, Any]] = None,
        ref2mappable_len: Optional[Dict[str, int]] = None,
        mappable_ref2forward_sum: Optional[Dict[str, Any]] = None,
        mappable_ref2reverse_sum: Optional[Dict[str, Any]] = None,
        ref2masc: Optional[Dict[str, Any]] = None
    ) -> 'StatisticsResult':
        """Create StatisticsResult from file data (for pymasc-plot).

        This method creates a StatisticsResult from data loaded from files,
        replicating the functionality of CCResult.from_file_data_with_builder()
        for the new StatisticsResult system.

        Args:
            mv_avr_filter_len: Moving average filter length
            chi2_pval: Chi-squared p-value threshold
            filter_mask_len: Filter mask length around read length
            min_calc_width: Minimum width for background calculation
            expected_library_len: Expected fragment length
            read_len: Read length in base pairs
            references: List of chromosome names
            ref2genomelen: Genome lengths by chromosome
            ref2forward_sum: Forward read counts by chromosome
            ref2reverse_sum: Reverse read counts by chromosome
            ref2cc: Cross-correlation data by chromosome
            ref2mappable_len: Mappable lengths by chromosome
            mappable_ref2forward_sum: Mappable forward counts by chromosome
            mappable_ref2reverse_sum: Mappable reverse counts by chromosome
            ref2masc: MSCC data by chromosome

        Returns:
            StatisticsResult instance constructed from file data
        """
        if not read_len:
            raise ValueError("Read length is required")
        if not references:
            raise ValueError("References list is required")
        if not ref2genomelen:
            raise ValueError("Genome lengths are required")

        # Determine skip_ncc and has_mappability
        skip_ncc = not (ref2forward_sum and ref2reverse_sum and ref2cc)
        has_mappability = bool(ref2mappable_len and ref2masc)

        if skip_ncc and not has_mappability:
            raise ValueError("Either NCC or MSCC data must be provided")

        # Create per-chromosome stats
        chromosome_stats: Dict[str, UnifiedStats] = {}

        for ref in references:
            genomelen = ref2genomelen[ref]

            # NCC data
            forward_sum = ref2forward_sum.get(ref) if ref2forward_sum else None
            reverse_sum = ref2reverse_sum.get(ref) if ref2reverse_sum else None
            cc = ref2cc.get(ref) if ref2cc else None

            # MSCC data
            mappable_len = ref2mappable_len.get(ref) if ref2mappable_len else None
            mappable_forward_sum = mappable_ref2forward_sum.get(ref) if mappable_ref2forward_sum else None
            mappable_reverse_sum = mappable_ref2reverse_sum.get(ref) if mappable_ref2reverse_sum else None
            masc = ref2masc.get(ref) if ref2masc else None

            # Create UnifiedStats for this chromosome
            chromosome_stats[ref] = UnifiedStats(
                read_len=read_len,
                mv_avr_filter_len=mv_avr_filter_len,
                expected_library_len=expected_library_len,
                genomelen=genomelen,
                forward_sum=forward_sum or 0,
                reverse_sum=reverse_sum or 0,
                cc=cc,
                mappable_len=mappable_len,
                mappable_forward_sum=mappable_forward_sum,
                mappable_reverse_sum=mappable_reverse_sum,
                masc=masc,
                filter_mask_len=filter_mask_len,
                min_calc_width=min_calc_width
            )

        # Calculate genome-wide statistics using Fisher z-transformation
        from PyMaSC.utils.stats_utils import CrossCorrelationMerger, ArrayAggregator

        # Aggregate NCC statistics
        merged_cc = None
        if not skip_ncc:
            merger = CrossCorrelationMerger(confidence_interval=0.95)
            ncc_data = [(ref2genomelen[ref], chromosome_stats[ref].cc.cc)
                       for ref in references if chromosome_stats[ref].cc and chromosome_stats[ref].cc.cc is not None]
            if ncc_data:
                genome_lengths, correlation_arrays = zip(*ncc_data)
                merged_cc, _, _ = merger.merge_correlations(
                    np.array(genome_lengths), correlation_arrays, read_len
                )

        # Aggregate MSCC statistics
        merged_mscc = None
        if has_mappability:
            merger = CrossCorrelationMerger(confidence_interval=0.95)
            mscc_data = [(ref2mappable_len[ref], chromosome_stats[ref].masc.cc)
                        for ref in references if chromosome_stats[ref].masc and chromosome_stats[ref].masc.cc is not None]
            if mscc_data:
                mappable_lengths, correlation_arrays = zip(*mscc_data)
                merged_mscc, _, _ = merger.merge_correlations(
                    np.array(mappable_lengths), correlation_arrays, read_len
                )

        # Calculate totals for genome-wide stats
        total_genomelen = sum(ref2genomelen.values())
        total_forward = sum(ref2forward_sum.values()) if ref2forward_sum else 0
        total_reverse = sum(ref2reverse_sum.values()) if ref2reverse_sum else 0

        total_mappable_len = None
        total_mappable_forward = None
        total_mappable_reverse = None

        if has_mappability:
            total_mappable_len = ArrayAggregator.safe_array_sum(ref2mappable_len.values())
            if mappable_ref2forward_sum:
                total_mappable_forward = ArrayAggregator.safe_array_sum(mappable_ref2forward_sum.values())
            if mappable_ref2reverse_sum:
                total_mappable_reverse = ArrayAggregator.safe_array_sum(mappable_ref2reverse_sum.values())

        # Create genome-wide UnifiedStats
        genome_wide_stats = UnifiedStats(
            read_len=read_len,
            mv_avr_filter_len=mv_avr_filter_len,
            expected_library_len=expected_library_len,
            genomelen=total_genomelen,
            forward_sum=total_forward,
            reverse_sum=total_reverse,
            cc=merged_cc,
            mappable_len=total_mappable_len,
            mappable_forward_sum=total_mappable_forward,
            mappable_reverse_sum=total_mappable_reverse,
            masc=merged_mscc,
            filter_mask_len=filter_mask_len,
            warning=True,
            min_calc_width=min_calc_width
        )

        # Extract quality metrics
        nsc = 0.0
        rsc = 0.0
        estimated_fragment_length = None

        if genome_wide_stats.cc:
            nsc = getattr(genome_wide_stats.cc, 'nsc', 0.0) or 0.0
            rsc = getattr(genome_wide_stats.cc, 'rsc', 0.0) or 0.0

        estimated_fragment_length = getattr(genome_wide_stats, 'est_lib_len', None)

        return cls(
            chromosome_stats=chromosome_stats,
            genome_wide_stats=genome_wide_stats,
            references=references,
            read_length=read_len,
            skip_ncc=skip_ncc,
            has_mappability=has_mappability,
            nsc=nsc,
            rsc=rsc,
            estimated_fragment_length=estimated_fragment_length
        )


class GenomeWideStatisticsResult(StatisticsResult):
    """StatisticsResult for GenomeWideResult input data.

    This class extends StatisticsResult to handle GenomeWideResult objects
    directly, providing the same statistical analysis capabilities while
    working with the new result container system.
    """

    @classmethod
    def from_genome_wide_result(
        cls,
        genome_wide_result: GenomeWideResult,
        references: List[str],
        lengths: List[int],
        read_length: int,
        mv_avr_filter_len: int = 15,
        expected_library_len: Optional[int] = None,
        filter_mask_len: int = 5,
        min_calc_width: int = 10
    ) -> 'GenomeWideStatisticsResult':
        """Create GenomeWideStatisticsResult from GenomeWideResult.

        Args:
            genome_wide_result: GenomeWideResult object containing calculation results
            references: List of chromosome names
            lengths: List of chromosome lengths
            read_length: Read length in base pairs
            mv_avr_filter_len: Moving average filter length for peak detection
            expected_library_len: Expected fragment length (optional)
            filter_mask_len: Masking window around read length for phantom peak
            min_calc_width: Minimum width for background correlation calculation

        Returns:
            GenomeWideStatisticsResult with complete statistical analysis
        """
        # Convert GenomeWideResult to statistics containers
        chromosome_stats = cls._convert_to_chromosome_stats(
            genome_wide_result, references, lengths, read_length,
            mv_avr_filter_len, expected_library_len, filter_mask_len, min_calc_width
        )

        # Calculate genome-wide statistics
        genome_wide_stats = cls._calculate_genome_wide_stats_from_containers(
            genome_wide_result, references, lengths, read_length,
            mv_avr_filter_len, expected_library_len, filter_mask_len, min_calc_width
        )

        # Extract quality metrics
        nsc = 0.0
        rsc = 0.0
        estimated_fragment_length = None

        if genome_wide_stats.cc:
            nsc = getattr(genome_wide_stats.cc, 'nsc', 0.0) or 0.0
            rsc = getattr(genome_wide_stats.cc, 'rsc', 0.0) or 0.0

        estimated_fragment_length = getattr(genome_wide_stats, 'est_lib_len', None)

        # Determine capabilities
        skip_ncc = not cls._has_ncc_data(genome_wide_result)
        has_mappability = cls._has_mappability_data(genome_wide_result)

        return cls(
            chromosome_stats=chromosome_stats,
            genome_wide_stats=genome_wide_stats,
            references=references,
            read_length=read_length,
            skip_ncc=skip_ncc,
            has_mappability=has_mappability,
            nsc=nsc,
            rsc=rsc,
            estimated_fragment_length=estimated_fragment_length
        )

    @staticmethod
    def _has_ncc_data(genome_wide_result: GenomeWideResult) -> bool:
        """Check if GenomeWideResult contains NCC data."""
        return isinstance(genome_wide_result, (NCCGenomeWideResult, BothGenomeWideResult))

    @staticmethod
    def _has_mappability_data(genome_wide_result: GenomeWideResult) -> bool:
        """Check if GenomeWideResult contains mappability data."""
        return isinstance(genome_wide_result, (MSCCGenomeWideResult, BothGenomeWideResult))

    @staticmethod
    def _convert_to_chromosome_stats(
        genome_wide_result: GenomeWideResult,
        references: List[str],
        lengths: List[int],
        read_length: int,
        mv_avr_filter_len: int,
        expected_library_len: Optional[int],
        filter_mask_len: int,
        min_calc_width: int
    ) -> Dict[str, UnifiedStats]:
        """Convert GenomeWideResult to per-chromosome UnifiedStats."""
        chromosome_stats = {}

        for i, ref in enumerate(references):
            genomelen = lengths[i]

            # Extract NCC data
            ncc_forward_sum = None
            ncc_reverse_sum = None
            ccbins = None

            if isinstance(genome_wide_result, (NCCGenomeWideResult, BothGenomeWideResult)):
                ncc_result = genome_wide_result.chroms.get(ref)
                if ncc_result:
                    ncc_forward_sum = ncc_result.forward_sum
                    ncc_reverse_sum = ncc_result.reverse_sum
                    ccbins = ncc_result.ccbins

            # Extract MSCC data
            mappable_len = None
            mappable_forward_sum = None
            mappable_reverse_sum = None
            mappable_ccbins = None

            if isinstance(genome_wide_result, MSCCGenomeWideResult):
                mscc_result = genome_wide_result.chroms.get(ref)
                if mscc_result:
                    mappable_len = mscc_result.mappable_len
                    mappable_forward_sum = mscc_result.forward_sum
                    mappable_reverse_sum = mscc_result.reverse_sum
                    mappable_ccbins = mscc_result.ccbins
            elif isinstance(genome_wide_result, BothGenomeWideResult):
                mscc_result = genome_wide_result.mappable_chroms.get(ref)
                if mscc_result:
                    mappable_len = mscc_result.mappable_len
                    mappable_forward_sum = mscc_result.forward_sum
                    mappable_reverse_sum = mscc_result.reverse_sum
                    mappable_ccbins = mscc_result.ccbins

            # Create UnifiedStats for this chromosome
            chromosome_stats[ref] = UnifiedStats(
                read_len=read_length,
                mv_avr_filter_len=mv_avr_filter_len,
                expected_library_len=expected_library_len,
                genomelen=genomelen,
                forward_sum=ncc_forward_sum or 0,
                reverse_sum=ncc_reverse_sum or 0,
                ccbins=ccbins,
                mappable_len=mappable_len,
                mappable_forward_sum=mappable_forward_sum,
                mappable_reverse_sum=mappable_reverse_sum,
                mappable_ccbins=mappable_ccbins,
                filter_mask_len=filter_mask_len,
                min_calc_width=min_calc_width
            )

        return chromosome_stats

    @staticmethod
    def _calculate_genome_wide_stats_from_containers(
        genome_wide_result: GenomeWideResult,
        references: List[str],
        lengths: List[int],
        read_length: int,
        mv_avr_filter_len: int,
        expected_library_len: Optional[int],
        filter_mask_len: int,
        min_calc_width: int
    ) -> UnifiedStats:
        """Calculate genome-wide statistics from GenomeWideResult containers."""
        # Get per-chromosome stats first
        ref2stats = GenomeWideStatisticsResult._convert_to_chromosome_stats(
            genome_wide_result, references, lengths, read_length,
            mv_avr_filter_len, expected_library_len, filter_mask_len, min_calc_width
        )

        # Calculate totals
        total_genomelen = sum(lengths)
        total_forward = 0
        total_reverse = 0

        if isinstance(genome_wide_result, (NCCGenomeWideResult, BothGenomeWideResult)):
            total_forward = genome_wide_result.forward_sum
            total_reverse = genome_wide_result.reverse_sum

        # Aggregate correlation arrays using Fisher z-transformation
        merged_cc = None
        if isinstance(genome_wide_result, (NCCGenomeWideResult, BothGenomeWideResult)):
            ncc_data = [(lengths[i], ref2stats[ref].cc.cc)
                       for i, ref in enumerate(references)
                       if ref in ref2stats and ref2stats[ref].cc is not None]
            if ncc_data:
                genome_lengths, correlation_arrays = zip(*ncc_data)
                from PyMaSC.core.constants import MERGED_CC_CONFIDENCE_INTERVAL
                from PyMaSC.utils.stats_utils import CrossCorrelationMerger
                merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
                merged_cc, _, _ = merger.merge_correlations(
                    np.array(genome_lengths), correlation_arrays, read_length
                )

        merged_mscc = None
        total_mappable_len = None
        total_mappable_forward = None
        total_mappable_reverse = None

        if isinstance(genome_wide_result, (MSCCGenomeWideResult, BothGenomeWideResult)):
            # Extract MSCC data for aggregation
            mscc_data = []
            mappable_lens = []
            mappable_forwards = []
            mappable_reverses = []

            for i, ref in enumerate(references):
                if ref in ref2stats and ref2stats[ref].masc is not None:
                    if isinstance(genome_wide_result, MSCCGenomeWideResult):
                        mscc_result = genome_wide_result.chroms.get(ref)
                    else:  # BothGenomeWideResult
                        mscc_result = genome_wide_result.mappable_chroms.get(ref)

                    if mscc_result and mscc_result.mappable_len is not None:
                        # Use the final mappable length value
                        mappable_len_value = mscc_result.mappable_len[read_length - 1] if len(mscc_result.mappable_len) > read_length - 1 else mscc_result.mappable_len[-1]
                        mscc_data.append((mappable_len_value, ref2stats[ref].masc.cc))
                        mappable_lens.append(mscc_result.mappable_len)
                        mappable_forwards.append(mscc_result.forward_sum)
                        mappable_reverses.append(mscc_result.reverse_sum)

            if mscc_data:
                mappable_lengths, correlation_arrays = zip(*mscc_data)
                from PyMaSC.core.constants import MERGED_CC_CONFIDENCE_INTERVAL
                from PyMaSC.utils.stats_utils import CrossCorrelationMerger, ArrayAggregator
                merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
                merged_mscc, _, _ = merger.merge_correlations(
                    np.array(mappable_lengths), correlation_arrays, read_length
                )

                # Aggregate mappability arrays
                if mappable_lens:
                    total_mappable_len = ArrayAggregator.safe_array_sum(mappable_lens)
                if mappable_forwards:
                    total_mappable_forward = ArrayAggregator.safe_array_sum(mappable_forwards)
                if mappable_reverses:
                    total_mappable_reverse = ArrayAggregator.safe_array_sum(mappable_reverses)

        return UnifiedStats(
            read_len=read_length,
            mv_avr_filter_len=mv_avr_filter_len,
            expected_library_len=expected_library_len,
            genomelen=total_genomelen,
            forward_sum=total_forward,
            reverse_sum=total_reverse,
            cc=merged_cc,
            mappable_len=total_mappable_len,
            mappable_forward_sum=total_mappable_forward,
            mappable_reverse_sum=total_mappable_reverse,
            masc=merged_mscc,
            filter_mask_len=filter_mask_len,
            warning=True,
            min_calc_width=min_calc_width
        )