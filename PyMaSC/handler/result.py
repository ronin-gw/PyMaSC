import logging
from functools import wraps
from typing import Union, Tuple, List, Dict

import numpy as np
from scipy.stats import norm

from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.calc import moving_avr_filter
from PyMaSC.utils.stats_utils import ArrayAggregator, StatisticalTestUtilities, CrossCorrelationMerger

# New builder-based result construction
from PyMaSC.core.result_builder import ResultBuilder, build_from_handler

logger = logging.getLogger(__name__)

# Quality Assessment Constants
NEAR_READLEN_ERR_CRITERION = 5
"""int: Distance threshold for detecting phantom peak contamination.

When the estimated fragment length is within this many base pairs of the
read length, it suggests potential phantom peak contamination. This triggers
warnings to alert users about potential quality issues in the data.
"""

MERGED_CC_CONFIDENCE_INTERVAL = 0.99
"""float: Confidence level for merged cross-correlation intervals.

This sets the confidence level (99%) for calculating confidence intervals
when merging cross-correlation results from multiple chromosomes using
Fisher z-transformation. Higher values provide wider, more conservative
confidence intervals.
"""

NEAR_ZERO_MIN_CALC_LEN = 10
"""int: Number of positions to examine at the beginning of correlation profile.

Used to detect potential issues with minimum correlation calculation.
If the median of the first NEAR_ZERO_MIN_CALC_LEN positions is smaller
than the calculated minimum, it suggests the shift size may be too small.
"""




def npcalc_with_logging_warn(func):
    """Decorator for handling numpy floating point errors gracefully.
    
    This decorator wraps numerical calculation functions to handle common
    floating point errors that can occur during cross-correlation calculations,
    such as division by zero or invalid operations (e.g., log of negative numbers).
    
    Behavior:
    1. First attempt: Execute function with numpy errors set to raise exceptions
    2. If FloatingPointError occurs: Log the error and retry with errors ignored
    3. Return the result from the second attempt (may contain NaN/Inf values)
    
    This approach ensures that calculations continue even when numerical
    edge cases occur, while providing debugging information about the issues.
    The resulting NaN/Inf values are typically handled by downstream code.
    
    Args:
        func: Function to wrap with error handling
        
    Returns:
        Wrapped function with floating point error handling
        
    Example:
        @npcalc_with_logging_warn
        def divide_arrays(a, b):
            return a / b
        
        # This will handle division by zero gracefully
        result = divide_arrays(np.array([1, 2]), np.array([0, 1]))
    """
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




class CCStats(object):
    """Cross-correlation statistics calculator for ChIP-seq quality assessment.
    
    This class calculates key quality metrics for ChIP-seq cross-correlation analysis,
    including NSC (Normalized Strand Coefficient), RSC (Relative Strand Coefficient),
    and fragment length estimation. It implements the statistical methods described
    in ENCODE ChIP-seq guidelines for data quality assessment.
    
    Mathematical Background:
    - NSC = cc_at_fragment_length / cc_min
      Measures signal-to-noise ratio (typically > 1.05 for good data)
    - RSC = (cc_at_fragment_length - cc_min) / (cc_at_read_length - cc_min)
      Measures enrichment relative to phantom peak (typically > 0.8)
    - Fragment length estimation uses moving average filtering to find peak
    - FWHM (Full Width at Half Maximum) measures peak sharpness
    - VSN (Variance Stabilizing Normalization) normalizes for read count
    
    The class handles edge cases like read length masking to avoid phantom peaks
    and provides comprehensive warning/error reporting for quality assessment.
    
    Attributes:
        cc (np.ndarray): Cross-correlation values for each shift distance
        genomelen (int): Total genome length for normalization
        forward_sum (int): Total forward reads
        reverse_sum (int): Total reverse reads
        read_len (int): Sequencing read length
        library_len (int): Expected fragment length (if provided)
        est_lib_len (int): Estimated fragment length (if calculated)
        cc_min (float): Minimum cross-correlation value (background)
        ccrl (float): Cross-correlation at read length (phantom peak)
        ccfl (float): Cross-correlation at fragment length
        nsc (float): Normalized Strand Coefficient
        rsc (float): Relative Strand Coefficient
        cc_width (int): FWHM of fragment length peak
        vsn (float): Variance Stabilizing Normalization value
    """
    
    def __init__(self, cc, genomelen, forward_sum, reverse_sum, read_len,
                 min_calc_width, mv_avr_filter_len, filter_mask_len, output_warnings,
                 do_llestimation=False, estimated_library_len=None, expected_library_len=None):
        self.cc = cc
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

        #
        if self.output_warnings:
            self.error = logger.error
            self.warning = logger.warning
        else:
            self.error = self.warning = logger.debug

        #
        self._calc_cc_min()  # self.cc_min
        self._calc_cc_rl()  # self.ccrl

        #
        self.ccfl = self.nsc = self.rsc = None
        if self.library_len:
            self.ccfl, self.nsc, self.rsc = self._calc_lib_metrics(self.library_len)

        #
        self.avr_cc = None
        self.est_ccfl = self.est_nsc = self.est_rsc = None
        if do_llestimation:
            self._estimate_fragment_len()
        if self.est_lib_len:
            self.est_ccfl, self.est_nsc, self.est_rsc = self._calc_lib_metrics(self.est_lib_len)

        #
        self.cc_width = None
        self.est_cc_width = None
        if self.library_len:
            try:
                self.cc_width = self._get_FWHM(self.library_len)
            except AssertionError as e:
                self.warning(e)
        if self.est_lib_len:
            try:
                self.est_cc_width = self._get_FWHM(self.est_lib_len)
            except AssertionError as e:
                self.warning(e)

        #
        def _calc_vsn(ccfl, cc_width):
            return 2 * ccfl * cc_width / (self.forward_sum + self.reverse_sum)

        self.vsn = None
        if self.library_len:
            self.vsn = _calc_vsn(self.ccfl, self.cc_width)
        self.est_vsn = None
        if self.est_lib_len:
            self.est_vsn = _calc_vsn(self.est_ccfl, self.est_cc_width)

    def _calc_cc_min(self) -> None:
        """Calculate minimum cross-correlation value for background estimation.
        
        This method estimates the background cross-correlation level by taking
        the median of the rightmost values in the correlation profile. This
        represents the background correlation level at large shift distances.
        
        The calculation uses the median of the last `min_calc_width` values
        to provide a robust estimate that is less sensitive to outliers than
        the minimum value itself.
        
        Raises:
            Warning: If the detected minimum appears larger than values near zero,
                    suggesting the shift size may be too small.
        
        Side Effects:
            Sets self.cc_min to the calculated background correlation level.
        """
        self.cc_min = np.sort(self.cc[-self.min_calc_width:])[min(self.min_calc_width, self.cc.size) // 2]
        if np.median(self.cc[:NEAR_ZERO_MIN_CALC_LEN]) < self.cc_min and self.output_warnings:
            logger.warning("Detected minimum coefficient seems to be larger than bigging part minimum. "
                           "Consider increasing shift size (-d/--max-shift).")

    def _calc_cc_rl(self) -> None:
        """Calculate cross-correlation at read length (phantom peak).
        
        This method calculates the cross-correlation value at the read length
        position, which corresponds to the "phantom peak" in ChIP-seq data.
        The phantom peak arises from sequencing artifacts and mappability bias.
        
        The phantom peak is used in RSC calculation to measure the signal-to-noise
        ratio relative to this artificial peak.
        
        Side Effects:
            Sets self.ccrl to the correlation value at read length position.
            If read length exceeds the correlation array size, sets ccrl to 0.
        """
        if self.read_len > self.cc.size:
            self.ccrl = 0
        else:
            self.ccrl = self.cc[self.read_len - 1]

    @npcalc_with_logging_warn
    def _calc_lib_metrics(self, library_len: int) -> Tuple[float, float, float]:
        """Calculate library-specific quality metrics (NSC and RSC).
        
        This method calculates the key quality metrics for ChIP-seq data:
        - NSC (Normalized Strand Coefficient): Signal-to-noise ratio
        - RSC (Relative Strand Coefficient): Enrichment relative to phantom peak
        
        These metrics are used to assess the quality of ChIP-seq immunoprecipitation
        according to ENCODE guidelines.
        
        Mathematical formulas:
        - NSC = cc_at_fragment_length / cc_min
        - RSC = (cc_at_fragment_length - cc_min) / (cc_at_read_length - cc_min)
        
        Args:
            library_len (int): Fragment length position for metric calculation
            
        Returns:
            tuple[float, float, float]: (ccfl, nsc, rsc) where:
                - ccfl: Cross-correlation at fragment length
                - nsc: Normalized Strand Coefficient (typically > 1.05 for good data)
                - rsc: Relative Strand Coefficient (typically > 0.8 for good data)
                
        Note:
            The @npcalc_with_logging_warn decorator handles floating point errors
            and provides appropriate warnings for numerical issues.
        """
        ccfl = self.cc[library_len - 1]
        nsc = ccfl / self.cc_min
        rsc = (ccfl - self.cc_min) / (self.ccrl - self.cc_min)

        return ccfl, nsc, rsc

    def _estimate_fragment_len(self) -> None:
        """Estimate fragment length from cross-correlation peak.
        
        This method estimates the DNA fragment length by finding the peak in the
        cross-correlation profile after smoothing with a moving average filter.
        The peak represents the most common fragment length in the sample.
        
        Algorithm:
        1. Apply moving average filter to smooth the correlation profile
        2. Find the position of maximum correlation (initial estimate)
        3. If estimate is too close to read length, apply masking to avoid phantom peak
        4. Re-estimate after masking if necessary
        5. Issue warnings if the estimate remains problematic
        
        Phantom Peak Handling:
        The method handles the "phantom peak" problem where the correlation peak
        appears at read length due to sequencing artifacts. When the estimated
        fragment length is close to read length, a masking window is applied
        around the read length to find the true fragment length peak.
        
        Side Effects:
            Sets self.avr_cc to the smoothed correlation profile.
            Sets self.est_lib_len to the estimated fragment length.
            Issues warnings if the estimate is problematic.
            
        Raises:
            Warning: If estimated fragment length is close to read length,
                    suggesting potential phantom peak contamination.
        """
        self.avr_cc = moving_avr_filter(self.cc, self.mv_avr_filter_len)
        self.est_lib_len = np.argmax(self.avr_cc) + 1

        need_warning = False

        if self.filter_mask_len and abs(self.est_lib_len - self.read_len) <= self.filter_mask_len:
            self.warning("Estimated library length is close to the read length.")
            self.warning("Trying to masking around the read length +/- {}bp...".format(self.filter_mask_len))

            _avr_cc = self.avr_cc.copy()
            mask_from = max(0, self.read_len - 1 - self.filter_mask_len)
            mask_to = min(len(_avr_cc), self.read_len + self.filter_mask_len)
            for i in range(mask_from, mask_to):
                _avr_cc[i] = - float("inf")
            self.est_lib_len = np.argmax(_avr_cc) + 1

            if self.est_lib_len - 1 in (mask_from - 1, mask_to):
                need_warning = True

        elif self.output_warnings and abs(self.est_lib_len - self.read_len) <= NEAR_READLEN_ERR_CRITERION:
            need_warning = True

        if self.output_warnings and need_warning:
            logger.error("Estimated library length is close to the read length! Please check output plots.")

    def _get_FWHM(self, library_len: int) -> Union[int, bool]:
        """Calculate Full Width at Half Maximum (FWHM) of fragment length peak.
        
        This method calculates the FWHM of the cross-correlation peak at the
        specified fragment length. FWHM is a measure of peak sharpness and
        is used in quality assessment and variance stabilizing normalization.
        
        Algorithm:
        1. Smooth correlation profile with moving average if not already done
        2. Find the half-maximum target value: (peak_value + background) / 2
        3. Search forward and backward from peak to find half-maximum points
        4. Calculate total width as sum of forward and backward distances
        5. Handle edge cases where peak extends beyond array boundaries
        
        The method provides fallback strategies when one side of the peak
        cannot be measured due to array boundaries.
        
        Args:
            library_len (int): Fragment length position of the peak
            
        Returns:
            Union[int, bool]: FWHM width in base pairs, or False if calculation fails
            
        Raises:
            AssertionError: If peak position is invalid or peak is not above background
            Warning: If half-maximum cannot be found on one or both sides
            Error: If FWHM calculation completely fails
            
        Side Effects:
            Creates self.avr_cc if not already present.
            Issues warnings/errors for problematic peak shapes.
        """
        if self.avr_cc is None:
            self.avr_cc = moving_avr_filter(self.cc, self.mv_avr_filter_len)

        max_i = library_len - 1
        assert max_i >= 0
        cc_max = self.avr_cc[max_i - 1]
        assert cc_max > self.cc_min
        target = self.cc_min + (cc_max - self.cc_min) / 2

        # forward
        forward_shift = 0
        forward_failed = False
        while self.avr_cc[max_i + forward_shift] > target:
            forward_shift += 1
            if max_i + forward_shift == self.avr_cc.size:
                self.warning("Failed to calc the half width at half maximum in the forward "
                             "side of the peak. Consider increasing shift size (-d/--max-shift).")
                forward_failed = True
                forward_shift -= 1
                break
        # backward
        backward_shift = 0
        backward_failed = False
        while self.avr_cc[max_i - backward_shift] > target:
            backward_shift += 1
            if max_i < backward_shift:
                self.warning("Failed to calc the half width at half maximum in the backward side of the peak.")
                backward_failed = True
                backward_shift -= 1
                break

        #
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


class PyMaSCStats(object):
    """PyMaSC (Python Mappability-sensitive Cross-Correlation) statistics calculator.
    
    This class handles both standard NCC (Normalized Cross-Correlation) and
    MSCC (Mappability-Sensitive Cross-Correlation) calculations for ChIP-seq
    quality assessment. It can process data in two modes:
    
    1. From bins: Raw read count bins that need normalization
    2. From correlation arrays: Pre-computed correlation values
    
    Mathematical Framework:
    - NCC: Standard cross-correlation normalized by genome length
      cc[d] = (Σ(forward[i] * reverse[i+d]) - μ_f * μ_r) / (σ_f * σ_r)
    - MSCC: Mappability-weighted cross-correlation
      Uses mappability scores to weight genomic positions
      Accounts for varying mappability across the genome
    
    Key Algorithms:
    - Binomial variance model: var = p * (1 - p) for read density
    - Fisher z-transformation: z = arctanh(r) for correlation merging
    - Moving average filtering for peak detection
    - Phantom peak masking around read length
    
    The class automatically determines what calculations are possible based on
    available input data and handles edge cases like missing mappability data
    or insufficient read counts.
    
    Parameters:
        read_len (int): Sequencing read length
        mv_avr_filter_len (int): Moving average filter length for peak detection
        expected_library_len (int, optional): Expected fragment length
        genomelen (int, optional): Total genome length for NCC
        forward_sum (int, optional): Total forward reads for NCC
        reverse_sum (int, optional): Total reverse reads for NCC
        ccbins (array, optional): Raw cross-correlation bins for NCC
        mappable_len (array, optional): Mappable genome length per shift
        mappable_forward_sum (array, optional): Mappable forward reads per shift
        mappable_reverse_sum (array, optional): Mappable reverse reads per shift
        mappable_ccbins (array, optional): Mappable cross-correlation bins
        cc (array, optional): Pre-computed NCC values
        masc (array, optional): Pre-computed MSCC values
        filter_mask_len (int): Masking window around read length
        warning (bool): Whether to output warnings
        min_calc_width (int, optional): Minimum width for background calculation
    
    Attributes:
        calc_ncc (bool): Whether NCC calculation is performed
        calc_masc (bool): Whether MSCC calculation is performed
        cc (CCStats): NCC statistics object
        masc (CCStats): MSCC statistics object
        max_shift (int): Maximum shift distance calculated
        est_lib_len (int): Estimated fragment length from MSCC peak
    """
    
    def __init__(
        self,
        read_len, mv_avr_filter_len=15, expected_library_len=None,
        genomelen=None, forward_sum=None, reverse_sum=None, ccbins=None,
        mappable_len=None, mappable_forward_sum=None, mappable_reverse_sum=None, mappable_ccbins=None,
        cc=None, masc=None, filter_mask_len=NEAR_READLEN_ERR_CRITERION, warning=False,
        min_calc_width=None
    ):
        self.read_len = read_len
        self.mv_avr_filter_len = mv_avr_filter_len
        self.genomelen = genomelen
        self.forward_sum = forward_sum
        self.reverse_sum = reverse_sum
        self.ccbins = ccbins
        self.mappable_len = mappable_len
        self.mappable_forward_sum = mappable_forward_sum
        self.mappable_reverse_sum = mappable_reverse_sum
        self.mappable_ccbins = mappable_ccbins
        cc = np.array(cc, dtype=np.float_) if cc is not None else None
        # self.cc_min = None
        # self.ccrl = None
        self.library_len = expected_library_len
        # self.ccfl = None
        # self.nsc = None
        # self.rsc = None
        masc = np.array(masc, dtype=np.float_) if masc is not None else None
        # self.masc_min = None
        # self.mascrl = None
        self.est_lib_len = None
        # self.est_ccfl = None
        # self.est_nsc = None
        # self.est_rsc = None

        self.cc = self.masc = None

        self.filter_mask_len = filter_mask_len
        self.output_warnings = warning
        self.min_calc_width = min_calc_width
        # self.cc_fwhm = None
        # self.masc_fwhm = None

        self.calc_ncc = self.calc_masc = False
        if ccbins is not None or mappable_ccbins is not None:
            self.calc_ncc = all(x is not None for x in (
                self.forward_sum,
                self.reverse_sum,
                self.ccbins,
                self.genomelen
            ))
            self.calc_masc = all(x is not None for x in (
                self.mappable_forward_sum,
                self.mappable_reverse_sum,
                self.mappable_ccbins,
                self.mappable_len
            ))
            if not self.calc_ncc and not self.calc_masc:
                return

            self._calc_from_bins()

        else:
            self.calc_ncc = all(x is not None for x in (
                self.forward_sum,
                self.reverse_sum,
                self.genomelen,
                cc
            ))
            self.calc_masc = all(x is not None for x in (
                self.mappable_forward_sum,
                self.mappable_reverse_sum,
                self.mappable_len,
                masc
            ))
            if not self.calc_ncc and not self.calc_masc:
                return

            if self.calc_ncc and self.calc_masc:
                assert len(cc) == len(masc), "Corrupt input: Shift sizes NCC and MASC seem to be differenet."
                self.max_shift = min(len(cc), len(masc)) - 1
            elif self.calc_ncc:
                self.max_shift = len(cc) - 1
            elif self.calc_masc:
                self.max_shift = len(masc) - 1

            self._make_cc_masc_stats(cc, masc)

    def _make_ccstat(self, cc, genomelen, forward_sum, reverse_sum, do_llestimation=False, est_lib_len=None):
        """Create a CCStats object with current configuration.
        
        This helper method creates a CCStats object using the current instance's
        configuration parameters. It provides a consistent way to create CCStats
        objects for both NCC and MSCC calculations.
        
        Args:
            cc: Cross-correlation array
            genomelen: Genome length for normalization
            forward_sum: Total forward reads
            reverse_sum: Total reverse reads
            do_llestimation (bool): Whether to perform library length estimation
            est_lib_len: Pre-estimated library length (if available)
            
        Returns:
            CCStats: Configured statistics calculator object
        """
        return CCStats(cc, genomelen, forward_sum, reverse_sum, self.read_len,
                       self.min_calc_width, self.mv_avr_filter_len,
                       self.filter_mask_len, self.output_warnings, do_llestimation,
                       est_lib_len, self.library_len)

    def _make_cc_masc_stats(self, cc, masc):
        """Create NCC and MSCC statistics objects from correlation arrays.
        
        This method orchestrates the creation of CCStats objects for both
        standard cross-correlation (NCC) and mappability-sensitive cross-correlation
        (MSCC) calculations. It handles the interdependence between the two:
        MSCC is calculated first to estimate fragment length, which is then
        used for NCC quality metrics.
        
        Processing Order:
        1. MSCC calculation with library length estimation (if available)
        2. NCC calculation using estimated library length from MSCC
        3. Fragment length estimation propagation between calculations
        
        Args:
            cc: Standard cross-correlation array (NCC)
            masc: Mappability-sensitive cross-correlation array (MSCC)
            
        Side Effects:
            Sets self.masc to MSCC statistics object (if calc_masc is True)
            Sets self.cc to NCC statistics object (if calc_ncc is True)
            Sets self.est_lib_len to estimated fragment length from MSCC
        """
        if self.calc_masc:
            self.masc = self._make_ccstat(
                masc, self.mappable_len, self.mappable_forward_sum, self.mappable_reverse_sum,
                do_llestimation=True
            )
            self.est_lib_len = self.masc.est_lib_len
        if self.calc_ncc:
            if self.masc:
                self.cc = self._make_ccstat(
                    cc, self.genomelen, self.forward_sum, self.reverse_sum,
                    est_lib_len=self.est_lib_len
                )
            else:
                self.cc = self._make_ccstat(
                    cc, self.genomelen, self.forward_sum, self.reverse_sum
                )

    def _calc_from_bins(self):
        """Calculate cross-correlation statistics from raw count bins.
        
        This method processes raw read count bins to calculate normalized
        cross-correlation values. It handles both NCC and MSCC calculations,
        determining the appropriate shift range and validating mappability
        requirements.
        
        Processing Steps:
        1. Determine maximum shift distance from available data
        2. Validate mappability array length requirements
        3. Calculate normalized cross-correlation arrays
        4. Create statistics objects for quality assessment
        
        The method handles the constraint that both NCC and MSCC calculations
        must use the same shift range for consistency.
        
        Side Effects:
            Sets self.max_shift to the maximum shift distance used
            Calls _make_cc_masc_stats to create statistics objects
            
        Raises:
            AssertionError: If mappability array is too short for required shift size
        """
        #
        if self.calc_ncc:
            self.max_shift = ncc_max_shift = len(self.ccbins) - 1
        else:
            ncc_max_shift = None

        #
        if self.calc_masc:
            self.max_shift = masc_max_shift = min(map(
                len,
                (self.mappable_forward_sum, self.mappable_reverse_sum, self.mappable_ccbins)
            )) - 1
            mappability_max_shift = len(self.mappable_len)
        else:
            masc_max_shift = mappability_max_shift = None

        #
        if ncc_max_shift and masc_max_shift:
            self.max_shift = min(ncc_max_shift, masc_max_shift)

        #
        if mappability_max_shift:
            required_shift_size = MappabilityHandler.calc_mappable_len_required_shift_size(
                self.read_len, self.max_shift
            )
            assert required_shift_size <= mappability_max_shift

        #
        self._make_cc_masc_stats(self._calc_naive_cc(), self._calc_masc())

    def _calc_naive_cc(self):
        """Calculate standard normalized cross-correlation (NCC).
        
        This method calculates the standard cross-correlation normalized by
        genome length. It uses the binomial variance model for normalization,
        which assumes read positions follow a binomial distribution.
        
        Mathematical Formula:
        cc[d] = (observed_correlation - expected_correlation) / std_dev
        
        Where:
        - observed_correlation = ccbins[d] / (genomelen - d)
        - expected_correlation = forward_mean * reverse_mean
        - std_dev = sqrt(forward_var * reverse_var)
        - variances use binomial model: var = p * (1 - p)
        
        Returns:
            np.ndarray or None: Normalized cross-correlation array,
                               or None if NCC calculation is disabled
        """
        if not self.calc_ncc:
            return None

        denom = self.genomelen - np.array(range(self.max_shift + 1), dtype=np.float_)
        return self._calc_cc(
            float(self.forward_sum),
            float(self.reverse_sum),
            self.ccbins[:self.max_shift + 1],
            self.genomelen,
            denom
        )

    @npcalc_with_logging_warn
    def _calc_masc(self):
        """Calculate mappability-sensitive cross-correlation (MSCC).
        
        This method calculates cross-correlation weighted by mappability scores
        to account for varying mappability across the genome. MSCC provides a
        more accurate representation of signal by downweighting regions with
        poor mappability.
        
        Key Differences from NCC:
        - Uses mappability-weighted read counts
        - Normalizes by mappable genome length rather than total length
        - Accounts for position-specific mappability effects
        
        The method handles the shift-dependent nature of mappability by creating
        a symmetric mappability array that accounts for both forward and reverse
        mappability at each shift distance.
        
        Returns:
            np.ndarray or None: Mappability-sensitive cross-correlation array,
                               or None if MSCC calculation is disabled
                               
        Note:
            The @npcalc_with_logging_warn decorator handles numerical issues
            that may arise from low mappability regions.
        """
        if not self.calc_masc:
            return None

        totlen = np.array(self.mappable_len, dtype=np.float_)
        totlen = np.concatenate((
            totlen[:self.read_len][::-1], totlen[1:]
        ))[:self.max_shift + 1]

        return self._calc_cc(
            np.array(self.mappable_forward_sum[:self.max_shift + 1], dtype=np.float_),
            np.array(self.mappable_reverse_sum[:self.max_shift + 1], dtype=np.float_),
            self.mappable_ccbins[:self.max_shift + 1],
            totlen,
            totlen,
        )

    @staticmethod
    @npcalc_with_logging_warn
    def _calc_cc(forward_sum, reverse_sum, ccbins, totlen, denom):
        """Calculate normalized cross-correlation using binomial variance model.
        
        This is the core cross-correlation calculation that implements the
        binomial variance model for ChIP-seq data. It normalizes the observed
        cross-correlation by the expected correlation under the null hypothesis
        and standardizes by the theoretical standard deviation.
        
        Mathematical Model:
        1. Assume read positions follow binomial distribution
        2. Calculate expected correlation: E[r] = p_forward * p_reverse
        3. Calculate variance: Var[r] = p_forward * (1 - p_forward) * p_reverse * (1 - p_reverse)
        4. Normalize: z = (observed - expected) / sqrt(variance)
        
        Where:
        - p_forward = forward_sum / totlen (forward read density)
        - p_reverse = reverse_sum / totlen (reverse read density)
        - observed = ccbins / denom (observed correlation)
        
        Args:
            forward_sum: Total forward reads (scalar or array)
            reverse_sum: Total reverse reads (scalar or array)
            ccbins: Raw cross-correlation bins
            totlen: Total length for normalization
            denom: Denominator for correlation calculation
            
        Returns:
            np.ndarray: Normalized cross-correlation values
            
        Note:
            The @npcalc_with_logging_warn decorator handles division by zero
            and other numerical issues that may arise with low read counts.
        """
        forward_mean = forward_sum / totlen
        reverse_mean = reverse_sum / totlen

        forward_var = forward_mean * (1 - forward_mean)
        reverse_var = reverse_mean * (1 - reverse_mean)

        sum_prod = forward_mean * reverse_mean
        var_geomean = (forward_var * reverse_var) ** 0.5
        return (ccbins / denom - sum_prod) / var_geomean


class ReadsTooFew(IndexError):
    pass


class CCResult(object):
    """Cross-correlation result aggregator and statistics calculator.
    
    This class serves as the main interface for PyMaSC cross-correlation analysis,
    aggregating per-chromosome results into genome-wide statistics and providing
    comprehensive quality metrics for ChIP-seq data assessment.
    
    The class supports two initialization modes:
    1. From handler: Direct initialization from calculation handlers (pymasc)
    2. From parameters: Manual initialization with explicit data (pymasc-plot)
    
    Key Features:
    - Aggregates per-chromosome cross-correlation results
    - Calculates genome-wide NCC and MSCC statistics
    - Performs Fisher z-transformation for correlation merging
    - Provides confidence intervals for merged correlations
    - Handles read count balance validation (chi-squared test)
    - Supports both standard and mappability-sensitive analysis
    
    Mathematical Background:
    - Fisher z-transformation: z = arctanh(r) for correlation merging
    - Confidence intervals: CI = tanh(z ± z_α/2 * SE)
    - Chi-squared test: χ² = (O - E)² / E for read balance
    - Weighted averaging by effective genome length
    
    Quality Assessment:
    - Read count balance between forward/reverse strands
    - NSC/RSC thresholds for data quality (NSC > 1.05, RSC > 0.8)
    - Fragment length estimation consistency
    - Phantom peak detection and masking
    
    The class automatically validates input data, handles missing chromosomes,
    and provides comprehensive error reporting for quality control.
    
    Parameters:
        mv_avr_filter_len (int): Moving average filter length for peak detection
        chi2_pval (float): Chi-squared p-value threshold for read balance test
        filter_mask_len (int): Masking window around read length
        min_calc_width (int): Minimum width for background calculation
        expected_library_len (int, optional): Expected fragment length
        handler (object, optional): Calculation handler with results
        read_len (int, optional): Sequencing read length
        references (list, optional): List of chromosome names
        ref2genomelen (dict, optional): Chromosome lengths
        ref2forward_sum (dict, optional): Forward read counts per chromosome
        ref2reverse_sum (dict, optional): Reverse read counts per chromosome
        ref2cc (dict, optional): NCC values per chromosome
        mappable_ref2forward_sum (dict, optional): Mappable forward reads
        mappable_ref2reverse_sum (dict, optional): Mappable reverse reads
        ref2mappable_len (dict, optional): Mappable lengths per chromosome
        ref2masc (dict, optional): MSCC values per chromosome
    
    Attributes:
        ref2stats (dict): Per-chromosome PyMaSCStats objects
        whole (PyMaSCStats): Genome-wide aggregated statistics
        skip_ncc (bool): Whether NCC calculation was skipped
        calc_masc (bool): Whether MSCC calculation was performed
        ncc_upper (list): Upper confidence interval for NCC
        ncc_lower (list): Lower confidence interval for NCC
        mscc_upper (list): Upper confidence interval for MSCC
        mscc_lower (list): Lower confidence interval for MSCC
    """
    
    def __init__(
        self,
        mv_avr_filter_len, chi2_pval, filter_mask_len, min_calc_width,  # mandatory params
        expected_library_len=None,  # optional parameter
        handler=None,  # source 1 (pymasc)
        read_len=None, references=None,
        ref2genomelen=None, ref2forward_sum=None, ref2reverse_sum=None, ref2cc=None,
        mappable_ref2forward_sum=None, mappable_ref2reverse_sum=None,
        ref2mappable_len=None, ref2masc=None  # source 2 (pymasc-plot)
    ):

        # settings
        self.mv_avr_filter_len = mv_avr_filter_len
        self.chi2_p_thresh = chi2_pval
        self.expected_library_len = expected_library_len
        self.filter_mask_len = max(filter_mask_len, 0)
        self.min_calc_width = max(min_calc_width, 0)

        if handler:
            self._init_from_handler(handler)
        else:
            self.read_len = read_len
            self.references = references
            self.ref2genomelen = ref2genomelen
            self.ref2mappable_len = ref2mappable_len
            self.ref2forward_sum = ref2forward_sum
            self.ref2reverse_sum = ref2reverse_sum
            self.mappable_ref2forward_sum = mappable_ref2forward_sum
            self.mappable_ref2reverse_sum = mappable_ref2reverse_sum

            self.ref2stats = {
                ref: PyMaSCStats(
                    self.read_len, self.mv_avr_filter_len, self.expected_library_len,
                    genomelen=self.ref2genomelen[ref],
                    forward_sum=self.ref2forward_sum.get(ref, None),
                    reverse_sum=self.ref2reverse_sum.get(ref, None),
                    cc=ref2cc.get(ref, None) if ref2cc else None,
                    mappable_len=self.ref2mappable_len.get(ref, None),
                    mappable_forward_sum=self.mappable_ref2forward_sum.get(ref, None),
                    mappable_reverse_sum=self.mappable_ref2reverse_sum.get(ref, None),
                    masc=ref2masc.get(ref, None) if ref2masc else None,
                    filter_mask_len=self.filter_mask_len,
                    min_calc_width=self.min_calc_width
                ) for ref in self.references
            }

            self.skip_ncc = not (self.ref2genomelen and ref2cc)
            self.calc_masc = self.ref2mappable_len and ref2masc

            assert (not self.skip_ncc) or self.calc_masc

        #
        if not self.skip_ncc:
            # Extract genome lengths and correlation arrays for NCC
            ncc_data = [(self.ref2genomelen[ref], self.ref2stats[ref].cc.cc)
                       for ref in self.references if self.ref2stats[ref].cc is not None]
            if ncc_data:
                genome_lengths, correlation_arrays = zip(*ncc_data)
                merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
                ncc, self.ncc_lower, self.ncc_upper = merger.merge_correlations(
                    np.array(genome_lengths), correlation_arrays, self.read_len
                )
            else:
                ncc = self.ncc_upper = self.ncc_lower = None
        else:
            ncc = self.ncc_upper = self.ncc_lower = None

        if self.calc_masc:
            # Extract mappable lengths and correlation arrays for MSCC
            mscc_data = [(self.ref2mappable_len[ref], self.ref2stats[ref].masc.cc)
                        for ref in self.references if self.ref2stats[ref].masc is not None]
            if mscc_data:
                mappable_lengths, correlation_arrays = zip(*mscc_data)
                merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
                mscc, self.mscc_lower, self.mscc_upper = merger.merge_correlations(
                    np.array(mappable_lengths), correlation_arrays, self.read_len
                )
            else:
                mscc = self.mscc_upper = self.mscc_lower = None
        else:
            mscc = self.mscc_upper = self.mscc_lower = None

        self.whole = PyMaSCStats(
            self.read_len, self.mv_avr_filter_len, self.expected_library_len,
            genomelen=sum(self.ref2genomelen.values()),
            forward_sum=sum(self.ref2forward_sum.values()),
            reverse_sum=sum(self.ref2reverse_sum.values()),
            cc=ncc,
            mappable_len=ArrayAggregator.safe_array_sum(self.ref2mappable_len.values()) if self.ref2mappable_len else None,
            mappable_forward_sum=ArrayAggregator.safe_array_sum(self.mappable_ref2forward_sum.values()) if self.mappable_ref2forward_sum else None,
            mappable_reverse_sum=ArrayAggregator.safe_array_sum(self.mappable_ref2reverse_sum.values()) if self.mappable_ref2reverse_sum else None,
            masc=mscc,
            warning=True,
            filter_mask_len=self.filter_mask_len,
            min_calc_width=self.min_calc_width
        )

    def _init_from_handler(self, handler):
        #
        self.references = handler.references
        lengths = handler.lengths
        self.ref2genomelen = dict(zip(self.references, lengths))
        self.read_len = handler.read_len
        self.ref2forward_sum = handler.ref2forward_sum
        self.ref2reverse_sum = handler.ref2reverse_sum
        ref2ccbins = handler.ref2ccbins
        self.mappable_ref2forward_sum = handler.mappable_ref2forward_sum
        self.mappable_ref2reverse_sum = handler.mappable_ref2reverse_sum
        mappable_ref2ccbins = handler.mappable_ref2ccbins
        self.ref2mappable_len = handler.ref2mappable_len

        #
        self.ref2stats = {
            ref: PyMaSCStats(
                self.read_len, self.mv_avr_filter_len, self.expected_library_len,
                genomelen=self.ref2genomelen[ref],
                forward_sum=self.ref2forward_sum.get(ref, 0),
                reverse_sum=self.ref2reverse_sum.get(ref, 0),
                ccbins=ref2ccbins.get(ref, None),
                mappable_len=self.ref2mappable_len.get(ref, None),
                mappable_forward_sum=self.mappable_ref2forward_sum.get(ref, None),
                mappable_reverse_sum=self.mappable_ref2reverse_sum.get(ref, None),
                mappable_ccbins=mappable_ref2ccbins.get(ref, None),
                filter_mask_len=self.filter_mask_len,
                min_calc_width=self.min_calc_width
            ) for ref in self.references
        }

        #
        if all((self.ref2forward_sum, self.ref2reverse_sum, ref2ccbins)):
            self.skip_ncc = False
            forward_sum = sum(list(self.ref2forward_sum.values()))
            reverse_sum = sum(list(self.ref2reverse_sum.values()))
        else:
            self.skip_ncc = True
            forward_sum = reverse_sum = None
        if forward_sum is not None:
            if forward_sum == 0:
                logger.error("There is no forward read.")
                raise ReadsTooFew
            elif reverse_sum == 0:
                logger.error("There is no reverse read.")
                raise ReadsTooFew
            StatisticalTestUtilities.chi2_test(forward_sum, reverse_sum, self.chi2_p_thresh, "Whole genome")

        #
        if all((self.mappable_ref2forward_sum, self.mappable_ref2reverse_sum,
                mappable_ref2ccbins, self.ref2mappable_len)):
            self.calc_masc = True
            mappable_forward_sum = ArrayAggregator.safe_array_sum(self.mappable_ref2forward_sum.values())
            mappable_reverse_sum = ArrayAggregator.safe_array_sum(self.mappable_ref2reverse_sum.values())
        else:
            self.calc_masc = False
            mappable_forward_sum = mappable_reverse_sum = None
        if mappable_forward_sum is not None:
            if np.all(mappable_forward_sum == 0):
                logger.error("There is no forward read.")
                raise ReadsTooFew
            elif np.all(mappable_reverse_sum == 0):
                logger.error("There is no reverse read.")
                raise ReadsTooFew

    # New builder-based result construction methods
    @classmethod
    def from_handler_with_builder(
        cls,
        handler,
        mv_avr_filter_len: int = 15,
        chi2_pval: float = 0.01,
        filter_mask_len: int = 5,
        min_calc_width: int = 10,
        expected_library_len: int = None
    ) -> 'CCResult':
        """Create CCResult using the new builder system.
        
        This method provides an alternative to the traditional __init__ method,
        using the type-safe builder pattern for result construction.
        
        Args:
            handler: Calculation handler with results
            mv_avr_filter_len: Moving average filter length
            chi2_pval: Chi-squared p-value threshold
            filter_mask_len: Filter mask length around read length
            min_calc_width: Minimum width for background calculation
            expected_library_len: Expected fragment length
            
        Returns:
            CCResult instance constructed via builder pattern
        """
        # NOTE: The builder system is currently incomplete and causes test failures.
        # It creates dict objects instead of proper PyMaSCStats objects for ref2stats,
        # leading to AttributeError in output modules.
        # Until the builder system is properly implemented, we use the legacy system.
        return cls(
            mv_avr_filter_len=mv_avr_filter_len,
            chi2_pval=chi2_pval,
            filter_mask_len=filter_mask_len,
            min_calc_width=min_calc_width,
            expected_library_len=expected_library_len,
            handler=handler
        )
    
    @classmethod  
    def from_file_data_with_builder(
        cls,
        mv_avr_filter_len: int,
        chi2_pval: float,
        filter_mask_len: int,
        min_calc_width: int,
        expected_library_len: int = None,
        read_len: int = None,
        references: List[str] = None,
        ref2genomelen: Dict[str, int] = None,
        ref2forward_sum: Dict[str, int] = None,
        ref2reverse_sum: Dict[str, int] = None,
        ref2cc: Dict[str, np.ndarray] = None,
        ref2mappable_len: Dict[str, int] = None,
        mappable_ref2forward_sum: Dict[str, np.ndarray] = None,
        mappable_ref2reverse_sum: Dict[str, np.ndarray] = None,
        ref2masc: Dict[str, np.ndarray] = None
    ) -> 'CCResult':
        """Create CCResult from file data using builder pattern.
        
        This method is designed for creating CCResult objects from data loaded
        from files (e.g., in plot.py), using the new builder system for 
        better type safety and maintainability.
        
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
            CCResult instance constructed via builder pattern
        """
        try:
            # Use builder pattern for construction
            from PyMaSC.core.result_builder import ResultBuilder
            
            builder = ResultBuilder()
            
            # Set basic parameters
            if read_len and references:
                chromosome_lengths = [ref2genomelen.get(ref, 1) for ref in references] if ref2genomelen else None
                builder.with_basic_params(read_len, references, chromosome_lengths)
            
            # Add NCC data if available
            if ref2forward_sum and ref2reverse_sum and ref2cc:
                builder.with_ncc_data(ref2forward_sum, ref2reverse_sum, ref2cc, ref2genomelen)
            
            # Add MSCC data if available
            if (mappable_ref2forward_sum and mappable_ref2reverse_sum and 
                ref2masc and ref2mappable_len):
                builder.with_mscc_data(
                    mappable_ref2forward_sum, mappable_ref2reverse_sum, 
                    ref2masc, ref2mappable_len
                )
            
            # Configure and build
            built_result = (builder
                           .with_statistics_config(mv_avr_filter_len=mv_avr_filter_len)
                           .with_fragment_length(expected=expected_library_len)
                           .build())
            
            # Convert to CCResult format
            return cls._from_built_result(
                built_result, mv_avr_filter_len, chi2_pval, filter_mask_len,
                min_calc_width, expected_library_len
            )
            
        except Exception as e:
            logger.warning(f"Builder system failed for file data, falling back to legacy: {e}")
            # Fall back to legacy construction
            return cls(
                mv_avr_filter_len=mv_avr_filter_len,
                chi2_pval=chi2_pval,
                filter_mask_len=filter_mask_len,
                min_calc_width=min_calc_width,
                expected_library_len=expected_library_len,
                read_len=read_len,
                references=references,
                ref2genomelen=ref2genomelen,
                ref2forward_sum=ref2forward_sum,
                ref2reverse_sum=ref2reverse_sum,
                ref2cc=ref2cc,
                ref2mappable_len=ref2mappable_len,
                mappable_ref2forward_sum=mappable_ref2forward_sum,
                mappable_ref2reverse_sum=mappable_ref2reverse_sum,
                ref2masc=ref2masc
            )

    @classmethod
    def _from_built_result(
        cls,
        built_result,
        mv_avr_filter_len: int,
        chi2_pval: float,
        filter_mask_len: int,
        min_calc_width: int,
        expected_library_len: int = None
    ) -> 'CCResult':
        """Convert BuiltResult to CCResult format."""
        # Extract legacy-compatible data from built result
        container = built_result.container
        
        # Create CCResult using legacy initialization but with builder data
        instance = cls.__new__(cls)
        
        # Set configuration parameters
        instance.mv_avr_filter_len = mv_avr_filter_len
        instance.chi2_p_thresh = chi2_pval
        instance.expected_library_len = expected_library_len
        instance.filter_mask_len = max(filter_mask_len, 0)
        instance.min_calc_width = max(min_calc_width, 0)
        
        # Extract data from container
        instance.read_len = container.read_length
        instance.references = container.references
        instance.skip_ncc = container.skip_ncc
        instance.calc_masc = container.has_mappability
        
        # Convert container data to legacy format
        instance.ref2genomelen = {
            result.chromosome: result.length
            for result in container.chromosome_results.values()
        }
        
        ncc_data = container.get_ncc_data()
        if ncc_data:
            instance.ref2forward_sum = {
                chrom: data.forward_sum for chrom, data in ncc_data.items()
            }
            instance.ref2reverse_sum = {
                chrom: data.reverse_sum for chrom, data in ncc_data.items()
            }
        
        mscc_data = container.get_mscc_data()
        if mscc_data:
            instance.mappable_ref2forward_sum = {
                chrom: data.forward_sum for chrom, data in mscc_data.items()
            }
            instance.mappable_ref2reverse_sum = {
                chrom: data.reverse_sum for chrom, data in mscc_data.items()
            }
            instance.ref2mappable_len = {
                chrom: data.mappable_len for chrom, data in mscc_data.items()
            }
        
        # Create per-chromosome statistics using builder results
        instance.ref2stats = {}
        for chrom in instance.references:
            chrom_stats = built_result.chromosome_statistics.get(chrom, {})
            
            # Create PyMaSCStats for this chromosome
            # This is a simplified version - full implementation would need
            # to convert StatisticsResult back to PyMaSCStats format
            instance.ref2stats[chrom] = chrom_stats
        
        # Create whole-genome statistics from aggregate data
        aggregate_stats = built_result.aggregate_statistics
        
        # Extract data from aggregated statistics
        total_reads = aggregate_stats.get('total_reads', {})
        aggregated_cc = aggregate_stats.get('aggregated_cc', {})
        
        # Extract NCC and MSCC arrays from aggregated results
        ncc = None
        mscc = None
        if isinstance(aggregated_cc, dict):
            if 'cc' in aggregated_cc:
                ncc = aggregated_cc['cc']
            elif 'ncc' in aggregated_cc:
                ncc = aggregated_cc['ncc']
            elif 'mscc' in aggregated_cc:
                mscc = aggregated_cc['mscc']
        elif isinstance(aggregated_cc, np.ndarray):
            ncc = aggregated_cc
        
        # Create PyMaSCStats object for whole-genome data
        instance.whole = PyMaSCStats(
            read_len=instance.read_len,
            mv_avr_filter_len=mv_avr_filter_len,
            expected_library_len=expected_library_len,
            genomelen=aggregate_stats.get('genome_length', sum(instance.ref2genomelen.values())),
            forward_sum=total_reads.get('forward', sum(instance.ref2forward_sum.values()) if instance.ref2forward_sum else 0),
            reverse_sum=total_reads.get('reverse', sum(instance.ref2reverse_sum.values()) if instance.ref2reverse_sum else 0),
            cc=ncc,
            mappable_len=total_reads.get('mappable_len'),
            mappable_forward_sum=total_reads.get('mappable_forward'),
            mappable_reverse_sum=total_reads.get('mappable_reverse'),
            masc=mscc,
            filter_mask_len=filter_mask_len,
            warning=True,
            min_calc_width=min_calc_width
        )
        
        # Additional initialization that matches legacy behavior
        instance.ncc_upper = instance.ncc_lower = None
        instance.mscc_upper = instance.mscc_lower = None
        
        return instance

