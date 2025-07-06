"""Cross-correlation result processing and statistical analysis.

This module contains classes and functions for processing cross-correlation
calculation results, computing statistical metrics, and performing fragment
length estimation. It handles both naive cross-correlation (NCC) and
mappability-sensitive cross-correlation (MSCC) results.

Key functionality includes:
- Cross-correlation statistics calculation (NSC, RSC, etc.)
- Fragment length estimation with peak detection
- Quality metrics computation and validation
- Result aggregation across multiple chromosomes
- Statistical testing for strand bias detection
"""
import logging
from functools import wraps
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple, TypeVar, Union

import numpy as np
from scipy.stats.distributions import chi2
from scipy.stats import norm

from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.calc import moving_avr_filter

logger = logging.getLogger(__name__)

NEAR_READLEN_ERR_CRITERION = 5
MERGED_CC_CONFIDENCE_INTERVAL = 0.99
NEAR_ZERO_MIN_CALC_LEN = 10


def _skip_none(i: Iterable[Optional[Any]]) -> List[Any]:
    """Filter None values from iterable.

    Args:
        i: Iterable containing values and None entries

    Returns:
        List with None values removed
    """
    return [x for x in i if x is not None]


F = TypeVar('F', bound=Callable[..., Any])

def npcalc_with_logging_warn(func: F) -> F:
    """Decorator for numpy calculation error handling.

    Wraps numpy calculations to catch floating point errors and warnings,
    logging them for debugging while allowing calculations to continue.

    Args:
        func: Function to wrap for error handling

    Returns:
        Wrapped function with error handling
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
    return _inner  # type: ignore


def chi2_test(a: Union[int, float], b: Union[int, float], chi2_p_thresh: float, label: str) -> None:
    """Chi-squared test for strand bias detection.

    Performs a chi-squared test to detect significant imbalance between
    forward and reverse read counts, which may indicate strand bias.

    Args:
        a: Forward read count
        b: Reverse read count
        chi2_p_thresh: P-value threshold for significance
        label: Label for logging output

    Note:
        Logs warnings if significant strand bias is detected
    """
    sum_ = a + b
    chi2_val = (((a - sum_ / 2.) ** 2) + ((b - sum_ / 2.) ** 2)) / sum_
    chi2_p = chi2.sf(chi2_val, 1)

    if chi2_p <= chi2_p_thresh:
        logger.warning("{} Forward/Reverse read count imbalance.".format(label))
        logger.warning("+/- = {} / {}, Chi-squared test p-val = {} <= {}".format(
            a, b, chi2_p, chi2_p_thresh
        ))
    else:
        logger.info("{} Forward/Reverse read count +/- = {} / {}".format(label, a, b))
        logger.info("Chi-squared test p-val = {} > {}".format(chi2_p, chi2_p_thresh))


class CCStats(object):
    """Cross-correlation statistics calculator.

    Computes comprehensive statistical metrics from cross-correlation data,
    including NSC (Normalized Strand Coefficient), RSC (Relative Strand Coefficient),
    fragment length estimation, and peak width calculations.

    This class handles the core statistical analysis of cross-correlation results,
    providing quality metrics commonly used in ChIP-seq analysis.

    Attributes:
        cc: Cross-correlation values array
        genomelen: Total mappable genome length
        forward_sum: Total forward reads
        reverse_sum: Total reverse reads
        read_len: Read length
        cc_min: Minimum cross-correlation value
        ccrl: Cross-correlation at read length
        ccfl: Cross-correlation at fragment length
        nsc: Normalized Strand Coefficient
        rsc: Relative Strand Coefficient
        est_lib_len: Estimated library/fragment length
        cc_width: Full Width at Half Maximum of CC peak
        vsn: Variance Stabilizing Normalization factor
    """
    def __init__(self, cc: np.ndarray, genomelen: int, forward_sum: int, reverse_sum: int, read_len: int,
                 min_calc_width: int, mv_avr_filter_len: int, filter_mask_len: Optional[int], output_warnings: bool,
                 do_llestimation: bool = False, estimated_library_len: Optional[int] = None, expected_library_len: Optional[int] = None) -> None:
        """Initialize cross-correlation statistics calculator.

        Args:
            cc: Cross-correlation values array
            genomelen: Total mappable genome length
            forward_sum: Total forward reads count
            reverse_sum: Total reverse reads count
            read_len: Read length in base pairs
            min_calc_width: Width for minimum value calculation
            mv_avr_filter_len: Moving average filter length
            filter_mask_len: Mask length around read length
            output_warnings: Whether to output warning messages
            do_llestimation: Whether to perform library length estimation
            estimated_library_len: Pre-computed estimated library length
            expected_library_len: User-specified expected library length
        """
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
        self.avr_cc: Optional[np.ndarray] = None
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
            if cc_width is None or cc_width is False:
                return None
            return 2 * ccfl * cc_width / (self.forward_sum + self.reverse_sum)

        self.vsn = None
        if self.library_len:
            self.vsn = _calc_vsn(self.ccfl, self.cc_width)
        self.est_vsn = None
        if self.est_lib_len:
            self.est_vsn = _calc_vsn(self.est_ccfl, self.est_cc_width)

    def _calc_cc_min(self) -> None:
        self.cc_min = np.sort(self.cc[-self.min_calc_width:])[min(self.min_calc_width, self.cc.size) // 2]
        if np.median(self.cc[:NEAR_ZERO_MIN_CALC_LEN]) < self.cc_min and self.output_warnings:
            logger.warning("Detected minimum coefficient seems to be larger than bigging part minimum. "
                           "Consider increasing shift size (-d/--max-shift).")

    def _calc_cc_rl(self) -> None:
        if self.read_len > self.cc.size:
            self.ccrl = 0
        else:
            self.ccrl = self.cc[self.read_len - 1]

    @npcalc_with_logging_warn
    def _calc_lib_metrics(self, library_len: int) -> Tuple[float, float, float]:
        ccfl = self.cc[library_len - 1]
        nsc = ccfl / self.cc_min
        rsc = (ccfl - self.cc_min) / (self.ccrl - self.cc_min)

        return ccfl, nsc, rsc

    def _estimate_fragment_len(self) -> None:
        self.avr_cc = moving_avr_filter(self.cc, self.mv_avr_filter_len)
        self.est_lib_len = int(np.argmax(self.avr_cc) + 1)

        need_warning = False

        if self.filter_mask_len and self.est_lib_len is not None and abs(self.est_lib_len - self.read_len) <= self.filter_mask_len:
            self.warning("Estimated library length is close to the read length.")
            self.warning("Trying to masking around the read length +/- {}bp...".format(self.filter_mask_len))

            _avr_cc = self.avr_cc.copy()
            mask_from = max(0, self.read_len - 1 - self.filter_mask_len)
            mask_to = min(len(_avr_cc), self.read_len + self.filter_mask_len)
            for i in range(mask_from, mask_to):
                _avr_cc[i] = - float("inf")
            self.est_lib_len = int(np.argmax(_avr_cc) + 1)

            if self.est_lib_len is not None and self.est_lib_len - 1 in (mask_from - 1, mask_to):
                need_warning = True

        elif self.output_warnings and self.est_lib_len is not None and abs(self.est_lib_len - self.read_len) <= NEAR_READLEN_ERR_CRITERION:
            need_warning = True

        if self.output_warnings and need_warning:
            logger.error("Estimated library length is close to the read length! Please check output plots.")

    def _get_FWHM(self, library_len: int) -> Union[int, bool]:
        if self.avr_cc is None:
            self.avr_cc = moving_avr_filter(self.cc, self.mv_avr_filter_len)

        assert self.avr_cc is not None  # Help mypy understand avr_cc is not None
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
    """Main PyMaSC statistics container and calculator.

    This class serves as the primary container for PyMaSC analysis results,
    handling both naive cross-correlation (NCC) and mappability-sensitive
    cross-correlation (MSCC) statistics. It can initialize from either
    raw calculation data or pre-computed correlation arrays.

    The class automatically determines which analyses are possible based
    on the provided data and creates appropriate CCStats objects for
    detailed statistical calculations.

    Attributes:
        read_len: Read length in base pairs
        calc_ncc: Whether naive cross-correlation was calculated
        calc_masc: Whether MSCC was calculated
        cc: CCStats object for naive cross-correlation
        masc: CCStats object for MSCC
        max_shift: Maximum shift distance used in analysis
    """
    def __init__(
        self,
        read_len: int, mv_avr_filter_len: int = 15, expected_library_len: Optional[int] = None,
        genomelen: Optional[int] = None, forward_sum: Optional[Union[int, np.ndarray]] = None, reverse_sum: Optional[Union[int, np.ndarray]] = None, ccbins: Optional[np.ndarray] = None,
        mappable_len: Optional[int] = None, mappable_forward_sum: Optional[Union[int, np.ndarray]] = None, mappable_reverse_sum: Optional[Union[int, np.ndarray]] = None, mappable_ccbins: Optional[np.ndarray] = None,
        cc: Optional[Union[List[float], np.ndarray]] = None, masc: Optional[Union[List[float], np.ndarray]] = None, filter_mask_len: int = NEAR_READLEN_ERR_CRITERION, warning: bool = False,
        min_calc_width: Optional[int] = None
    ) -> None:
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
        self.est_lib_len: Optional[int] = None
        # self.est_ccfl = None
        # self.est_nsc = None
        # self.est_rsc = None

        self.cc: Optional[CCStats] = None
        self.masc: Optional[CCStats] = None

        self.filter_mask_len = filter_mask_len
        self.output_warnings = warning
        self.min_calc_width = min_calc_width
        # self.cc_fwhm = None
        # self.masc_fwhm = None

        self.calc_ncc = self.calc_masc = False
        if ccbins is not None or mappable_ccbins is not None:
            # Check if we have NCC scalars or MSCC arrays
            has_ncc_data = all(x is not None for x in (
                self.forward_sum,
                self.reverse_sum,
                self.ccbins,
                self.genomelen
            ))

            if has_ncc_data:
                # Check if forward_sum is a scalar (NCC) or array (MSCC)
                try:
                    if hasattr(self.forward_sum, '__len__') and not isinstance(self.forward_sum, str):
                        # forward_sum is an array, this is MSCC data, don't calculate NCC
                        self.calc_ncc = False
                    else:
                        # forward_sum is a scalar, this is NCC data
                        self.calc_ncc = True
                except:
                    # On any error, skip NCC to be safe
                    self.calc_ncc = False
            else:
                self.calc_ncc = False
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
                assert cc is not None and masc is not None, "Corrupt input: NCC and MASC arrays are None"
                assert len(cc) == len(masc), "Corrupt input: Shift sizes NCC and MASC seem to be differenet."
                self.max_shift = min(len(cc), len(masc)) - 1
            elif self.calc_ncc:
                assert cc is not None, "Corrupt input: NCC array is None"
                self.max_shift = len(cc) - 1
            elif self.calc_masc:
                assert masc is not None, "Corrupt input: MASC array is None"
                self.max_shift = len(masc) - 1

            self._make_cc_masc_stats(cc, masc)

    def _make_ccstat(self, cc: np.ndarray, genomelen: int, forward_sum: Union[int, float, np.ndarray], reverse_sum: Union[int, float, np.ndarray], do_llestimation: bool = False, est_lib_len: Optional[int] = None) -> CCStats:
        # Convert arrays to scalars if needed
        if isinstance(forward_sum, np.ndarray):
            forward_sum_val = int(np.sum(forward_sum))
        else:
            forward_sum_val = int(forward_sum)

        if isinstance(reverse_sum, np.ndarray):
            reverse_sum_val = int(np.sum(reverse_sum))
        else:
            reverse_sum_val = int(reverse_sum)

        return CCStats(cc, genomelen, forward_sum_val, reverse_sum_val, self.read_len,
                       self.min_calc_width or 0, self.mv_avr_filter_len,
                       self.filter_mask_len, self.output_warnings, do_llestimation,
                       est_lib_len, self.library_len)

    def _make_cc_masc_stats(self, cc: Optional[np.ndarray], masc: Optional[np.ndarray]) -> None:
        if self.calc_masc:
            assert masc is not None and self.mappable_len is not None
            assert self.mappable_forward_sum is not None and self.mappable_reverse_sum is not None
            self.masc = self._make_ccstat(
                masc, self.mappable_len, self.mappable_forward_sum, self.mappable_reverse_sum,
                do_llestimation=True
            )
            if self.masc and self.masc.est_lib_len is not None:
                self.est_lib_len = self.masc.est_lib_len
        if self.calc_ncc:
            assert cc is not None and self.genomelen is not None
            assert self.forward_sum is not None and self.reverse_sum is not None
            if self.masc:
                self.cc = self._make_ccstat(
                    cc, self.genomelen, self.forward_sum, self.reverse_sum,
                    est_lib_len=self.est_lib_len
                )
            else:
                self.cc = self._make_ccstat(
                    cc, self.genomelen, self.forward_sum, self.reverse_sum
                )

    def _calc_from_bins(self) -> None:
        #
        if self.calc_ncc:
            assert self.ccbins is not None
            self.max_shift = ncc_max_shift = len(self.ccbins) - 1
        else:
            ncc_max_shift = None

        #
        if self.calc_masc:
            # Filter out None values before calculating length
            mappable_arrays = [
                arr for arr in (self.mappable_forward_sum, self.mappable_reverse_sum, self.mappable_ccbins)
                if arr is not None
            ]

            # Debug: check types of mappable arrays
            for i, arr in enumerate(mappable_arrays):
                if not hasattr(arr, '__len__'):
                    # If any array doesn't have len(), treat as no valid mappable data
                    mappable_arrays = []
                    break

            if mappable_arrays:
                # Type-safe length calculation
                lengths = []
                for arr in mappable_arrays:
                    if hasattr(arr, '__len__'):
                        lengths.append(len(arr))
                if lengths:
                    self.max_shift = masc_max_shift = min(lengths) - 1
                else:
                    masc_max_shift = 0
            else:
                masc_max_shift = 0  # Fallback if no valid arrays

            if self.mappable_len is not None:
                if isinstance(self.mappable_len, (list, np.ndarray)):
                    mappability_max_shift = len(self.mappable_len)
                else:
                    mappability_max_shift = None  # Scalar mappable_len
            else:
                mappability_max_shift = None
        else:
            masc_max_shift = None
            mappability_max_shift = None

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

    def _calc_naive_cc(self) -> Optional[np.ndarray]:
        if not self.calc_ncc:
            return None

        assert self.genomelen is not None
        assert self.forward_sum is not None
        assert self.reverse_sum is not None
        assert self.ccbins is not None

        denom = self.genomelen - np.array(range(self.max_shift + 1), dtype=np.float_)
        return self._calc_cc(
            float(self.forward_sum),
            float(self.reverse_sum),
            self.ccbins,
            self.genomelen,
            denom
        )

    @npcalc_with_logging_warn
    def _calc_masc(self) -> Optional[np.ndarray]:
        if not self.calc_masc:
            return None

        assert self.mappable_len is not None
        assert self.mappable_forward_sum is not None
        assert self.mappable_reverse_sum is not None
        assert self.mappable_ccbins is not None

        totlen = np.array(self.mappable_len, dtype=np.float_)
        totlen = np.concatenate((
            totlen[:self.read_len][::-1], totlen[1:]
        ))[:self.max_shift + 1]

        return self._calc_cc(
            np.array(self.mappable_forward_sum, dtype=np.float_),
            np.array(self.mappable_reverse_sum, dtype=np.float_),
            self.mappable_ccbins,
            totlen,
            totlen,
        )

    @staticmethod
    @npcalc_with_logging_warn
    def _calc_cc(forward_sum: Union[float, np.ndarray], reverse_sum: Union[float, np.ndarray], ccbins: np.ndarray, totlen: Union[int, float, np.ndarray], denom: Union[float, np.ndarray]) -> np.ndarray:
        forward_mean = forward_sum / totlen
        reverse_mean = reverse_sum / totlen

        forward_var = forward_mean * (1 - forward_mean)
        reverse_var = reverse_mean * (1 - reverse_mean)

        sum_prod = forward_mean * reverse_mean
        var_geomean = (forward_var * reverse_var) ** 0.5
        result = (ccbins / denom - sum_prod) / var_geomean
        return np.asarray(result, dtype=np.float64)


class ReadsTooFew(IndexError):
    """Exception raised when insufficient reads are available for analysis.

    This exception is raised when the number of reads is too low to perform
    reliable cross-correlation analysis or statistical calculations.
    """
    pass


class CCResult(object):
    """Complete result container and processor for PyMaSC analysis.

    This class serves as the main result container for PyMaSC cross-correlation
    analysis, aggregating results across all chromosomes and providing
    high-level statistical summaries. It can be initialized either from
    a calculation handler or from pre-loaded data tables.

    The class handles:
    - Aggregation of per-chromosome statistics
    - Genome-wide statistical summaries
    - Strand bias testing
    - Quality metric calculation
    - Result validation and error checking

    Attributes:
        read_len: Read length in base pairs
        references: List of analyzed chromosome names
        skip_ncc: Whether naive cross-correlation was skipped
        calc_masc: Whether MSCC was calculated
        ref2stats: Dictionary mapping chromosomes to PyMaSCStats objects
        merged_stats: Genome-wide aggregated statistics
    """
    def __init__(
        self,
        mv_avr_filter_len: int, chi2_pval: float, filter_mask_len: int, min_calc_width: int,  # mandatory params
        expected_library_len: Optional[int] = None,  # optional parameter
        handler: Optional[Any] = None,  # source 1 (pymasc)
        read_len: Optional[int] = None, references: Optional[List[str]] = None,
        ref2genomelen: Optional[Dict[str, int]] = None, ref2forward_sum: Optional[Dict[str, int]] = None, ref2reverse_sum: Optional[Dict[str, int]] = None, ref2cc: Optional[Dict[str, np.ndarray]] = None,
        mappable_ref2forward_sum: Optional[Dict[str, int]] = None, mappable_ref2reverse_sum: Optional[Dict[str, int]] = None,
        ref2mappable_len: Optional[Dict[str, int]] = None, ref2masc: Optional[Dict[str, np.ndarray]] = None  # source 2 (pymasc-plot)
    ) -> None:
        """Initialize CCResult from handler or pre-loaded data.

        Args:
            mv_avr_filter_len: Moving average filter length for smoothing
            chi2_pval: P-value threshold for chi-squared strand bias test
            filter_mask_len: Mask length around read length for peak detection
            min_calc_width: Width for minimum value calculation
            expected_library_len: User-specified expected fragment length
            handler: Calculation handler object (for direct initialization)
            read_len: Read length (for data table initialization)
            references: List of chromosome names (for data table initialization)
            ref2genomelen: Chromosome lengths dictionary
            ref2forward_sum: Forward read counts by chromosome
            ref2reverse_sum: Reverse read counts by chromosome
            ref2cc: Cross-correlation data by chromosome
            mappable_ref2forward_sum: Mappable forward read counts
            mappable_ref2reverse_sum: Mappable reverse read counts
            ref2mappable_len: Mappable lengths by chromosome
            ref2masc: MSCC data by chromosome

        Raises:
            ReadsTooFew: If insufficient reads for reliable analysis
        """

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

            # Type check required attributes
            assert self.read_len is not None
            assert self.references is not None
            assert self.ref2genomelen is not None

            self.ref2stats = {
                ref: PyMaSCStats(
                    self.read_len, self.mv_avr_filter_len, self.expected_library_len,
                    genomelen=self.ref2genomelen[ref],
                    forward_sum=self.ref2forward_sum.get(ref, None) if self.ref2forward_sum else None,
                    reverse_sum=self.ref2reverse_sum.get(ref, None) if self.ref2reverse_sum else None,
                    cc=ref2cc.get(ref, None) if ref2cc else None,
                    mappable_len=self.ref2mappable_len.get(ref, None) if self.ref2mappable_len else None,
                    mappable_forward_sum=self.mappable_ref2forward_sum.get(ref, None) if self.mappable_ref2forward_sum else None,
                    mappable_reverse_sum=self.mappable_ref2reverse_sum.get(ref, None) if self.mappable_ref2reverse_sum else None,
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
            assert self.ref2genomelen is not None
            assert self.references is not None
            ncc_data = []
            for ref in self.references:
                stats = self.ref2stats[ref]
                if stats.cc is not None:
                    ncc_data.append((self.ref2genomelen[ref], stats.cc.cc))

            if ncc_data:
                ncc, self.ncc_upper, self.ncc_lower = self._merge_cc(*zip(*ncc_data))
            else:
                ncc = self.ncc_upper = self.ncc_lower = None
        else:
            ncc = self.ncc_upper = self.ncc_lower = None

        if self.calc_masc:
            assert self.ref2mappable_len is not None
            assert self.references is not None
            mscc_data = []
            for ref in self.references:
                stats = self.ref2stats[ref]
                if stats.masc is not None:
                    mscc_data.append((self.ref2mappable_len[ref], stats.masc.cc))

            if mscc_data:
                mscc, self.mscc_upper, self.mscc_lower = self._merge_cc(*zip(*mscc_data))
            else:
                mscc = self.mscc_upper = self.mscc_lower = None
        else:
            mscc = self.mscc_upper = self.mscc_lower = None

        # Check if forward_sum values are scalars or arrays
        sample_forward = next(iter(self.ref2forward_sum.values()), None) if self.ref2forward_sum else None
        has_scalar_sums = sample_forward is None or not (hasattr(sample_forward, '__len__') and not isinstance(sample_forward, str))

        if has_scalar_sums:
            # Normal case: scalar sums for NCC
            genomelen_sum = sum(self.ref2genomelen.values()) if self.ref2genomelen else None
            forward_sum_total = sum(self.ref2forward_sum.values()) if self.ref2forward_sum else None
            reverse_sum_total = sum(self.ref2reverse_sum.values()) if self.ref2reverse_sum else None
        else:
            # MSCC case: arrays instead of scalars, use None for NCC data
            genomelen_sum = None
            forward_sum_total = None
            reverse_sum_total = None

        assert self.read_len is not None

        self.whole = PyMaSCStats(
            self.read_len, self.mv_avr_filter_len, self.expected_library_len,
            genomelen=genomelen_sum,
            forward_sum=forward_sum_total,
            reverse_sum=reverse_sum_total,
            cc=ncc,
            mappable_len=np.sum(tuple(self.ref2mappable_len.values()), axis=0) if self.ref2mappable_len else None,
            mappable_forward_sum=np.sum(tuple(v for v in self.mappable_ref2forward_sum.values() if v is not None), axis=0) if self.mappable_ref2forward_sum else None,
            mappable_reverse_sum=np.sum(tuple(v for v in self.mappable_ref2reverse_sum.values() if v is not None), axis=0) if self.mappable_ref2reverse_sum else None,
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
                mappable_forward_sum=self.mappable_ref2forward_sum.get(ref, 0),
                mappable_reverse_sum=self.mappable_ref2reverse_sum.get(ref, 0),
                mappable_ccbins=mappable_ref2ccbins.get(ref, None),
                filter_mask_len=self.filter_mask_len,
                min_calc_width=self.min_calc_width
            ) for ref in self.references
        }

        #
        if all((self.ref2forward_sum, self.ref2reverse_sum, ref2ccbins)):
            # Check if we have NCC scalars or MSCC arrays
            try:
                # Try to get a sample value to check type
                sample_forward = next(iter(self.ref2forward_sum.values()), None)
                if sample_forward is not None and hasattr(sample_forward, '__len__') and not isinstance(sample_forward, str):
                    # We have arrays (MSCC data), skip NCC
                    self.skip_ncc = True
                    forward_sum = reverse_sum = None
                else:
                    # We have scalars (NCC data)
                    self.skip_ncc = False
                    forward_sum = sum(list(self.ref2forward_sum.values()))
                    reverse_sum = sum(list(self.ref2reverse_sum.values()))
            except:
                # On any error, skip NCC to be safe
                self.skip_ncc = True
                forward_sum = reverse_sum = None
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
            chi2_test(forward_sum, reverse_sum, self.chi2_p_thresh, "Whole genome")

        #
        if all((self.mappable_ref2forward_sum, self.mappable_ref2reverse_sum,
                mappable_ref2ccbins, self.ref2mappable_len)):
            self.calc_masc = True
            mappable_forward_sum = np.sum(_skip_none(self.mappable_ref2forward_sum.values()), axis=0)
            mappable_reverse_sum = np.sum(_skip_none(self.mappable_ref2reverse_sum.values()), axis=0)
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

    def _merge_cc(self, ns, ccs):
        ns = np.array(ns)

        merged_r = []
        interval_upper = []
        interval_lower = []

        for i, _ccs in enumerate(zip(*ccs)):
            nans = np.isnan(_ccs)
            _ccs = np.array(_ccs)[~nans]

            if len(ns.shape) == 1:
                _ns = ns[~nans] - 3
            else:
                # Ensure index is within bounds
                index = abs(self.read_len - i - 1)  # Convert to 0-based indexing
                if index >= ns.shape[1]:
                    index = ns.shape[1] - 1  # Use last valid index if out of bounds
                _ns = ns[~nans, index] - 3

            zs = np.arctanh(_ccs)
            infs = np.isinf(zs)
            zs = zs[~infs]
            _ns = _ns[~infs]
            avr_z = np.average(zs, weights=_ns)

            z_interval = norm.ppf(1 - (1 - MERGED_CC_CONFIDENCE_INTERVAL) / 2) * np.sqrt(1 / np.sum(_ns))
            z_interval_upper = avr_z + z_interval
            z_interval_lower = avr_z - z_interval

            merged_r.append(np.tanh(avr_z))
            interval_upper.append(np.tanh(z_interval_upper))
            interval_lower.append(np.tanh(z_interval_lower))

        return merged_r, interval_lower, interval_upper
