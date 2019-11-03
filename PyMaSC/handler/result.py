import logging
from functools import wraps

import numpy as np
from scipy.stats.distributions import chi2
from scipy.stats import norm

from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.calc import moving_avr_filter

logger = logging.getLogger(__name__)

NEAR_READLEN_ERR_CRITERION = 5
MERGED_CC_CONFIDENCE_INTERVAL = 0.99
NEAR_ZERO_MIN_CALC_LEN = 10


def _skip_none(i):
    return [x for x in i if x is not None]


def npcalc_with_logging_warn(func):
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


def chi2_test(a, b, chi2_p_thresh, label):
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
                self.error(repr(e))
        if self.est_lib_len:
            try:
                self.est_cc_width = self._get_FWHM(self.est_lib_len)
            except AssertionError as e:
                self.error(repr(e))

        #
        def _calc_vsn(ccfl, cc_width):
            return 2 * ccfl * cc_width / (self.forward_sum + self.reverse_sum)

        self.vsn = None
        if self.library_len:
            self.vsn = _calc_vsn(self.ccfl, self.cc_width)
        self.est_vsn = None
        if self.est_lib_len:
            self.est_vsn = _calc_vsn(self.est_ccfl, self.est_cc_width)

    def _calc_cc_min(self):
        self.cc_min = np.sort(self.cc[-self.min_calc_width:])[min(self.min_calc_width, self.cc.size) // 2]
        if np.median(self.cc[:NEAR_ZERO_MIN_CALC_LEN]) < self.cc_min and self.output_warnings:
            logger.warning("Detected minimum coefficient seems to be larger than bigging part minimum. "
                           "Consider increasing shift size (-d/--max-shift).")

    def _calc_cc_rl(self):
        if self.read_len > self.cc.size:
            self.ccrl = 0
        else:
            self.ccrl = self.cc[self.read_len - 1]

    @npcalc_with_logging_warn
    def _calc_lib_metrics(self, library_len):
        ccfl = self.cc[library_len - 1]
        nsc = ccfl / self.cc_min
        rsc = (ccfl - self.cc_min) / (self.ccrl - self.cc_min)

        return ccfl, nsc, rsc

    def _estimate_fragment_len(self):
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

    def _get_FWHM(self, library_len):
        if self.avr_cc is None:
            self.avr_cc = moving_avr_filter(self.cc, self.mv_avr_filter_len)

        max_i = library_len - 1
        assert max_i >= 0

        cc_max = self.avr_cc[max_i]
        if cc_max <= self.cc_min:
            logger.error("Estimated max is less than minimum. Abort FWHM calculation.")
            return False

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
        return CCStats(cc, genomelen, forward_sum, reverse_sum, self.read_len,
                       self.min_calc_width, self.mv_avr_filter_len,
                       self.filter_mask_len, self.output_warnings, do_llestimation,
                       est_lib_len, self.library_len)

    def _make_cc_masc_stats(self, cc, masc):
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
        if not self.calc_ncc:
            return None

        denom = self.genomelen - np.array(range(self.max_shift + 1), dtype=np.float_)
        return self._calc_cc(
            float(self.forward_sum),
            float(self.reverse_sum),
            self.ccbins,
            self.genomelen,
            denom
        )

    @npcalc_with_logging_warn
    def _calc_masc(self):
        if not self.calc_masc:
            return None

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
    def _calc_cc(forward_sum, reverse_sum, ccbins, totlen, denom):
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
            ncc, self.ncc_upper, self.ncc_lower = self._merge_cc(
                *zip(*((self.ref2genomelen[ref], self.ref2stats[ref].cc.cc)
                       for ref in self.references if self.ref2stats[ref].cc is not None))
            )
        else:
            ncc = self.ncc_upper = self.ncc_lower = None

        if self.calc_masc:
            mscc, self.mscc_upper, self.mscc_lower = self._merge_cc(
                *zip(*((self.ref2mappable_len[ref], self.ref2stats[ref].masc.cc)
                       for ref in self.references if self.ref2stats[ref].masc is not None))
            )
        else:
            mscc = self.mscc_upper = self.mscc_lower = None

        self.whole = PyMaSCStats(
            self.read_len, self.mv_avr_filter_len, self.expected_library_len,
            genomelen=sum(self.ref2genomelen.values()),
            forward_sum=sum(self.ref2forward_sum.values()),
            reverse_sum=sum(self.ref2reverse_sum.values()),
            cc=ncc,
            mappable_len=np.sum(tuple(self.ref2mappable_len.values()), axis=0),
            mappable_forward_sum=np.sum(tuple(self.mappable_ref2forward_sum.values()), axis=0),
            mappable_reverse_sum=np.sum(tuple(self.mappable_ref2reverse_sum.values()), axis=0),
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
                _ns = ns[~nans, abs(self.read_len - i)] - 3

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
