import logging
from functools import wraps

import numpy as np
from scipy.stats.distributions import chi2

from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.calc import moving_avr_filter

logger = logging.getLogger(__name__)


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


class PyMaSCStats(object):
    def __init__(
        self,
        read_len, filter_len=15, expected_library_len=None,
        genomelen=None, forward_sum=None, reverse_sum=None, ccbins=None,
        mappable_len=None, mappable_forward_sum=None, mappable_reverse_sum=None, mappable_ccbins=None,
        cc=None, masc=None
    ):
        self.read_len = read_len
        self.filter_len = filter_len
        self.genomelen = genomelen
        self.forward_sum = forward_sum
        self.reverse_sum = reverse_sum
        self.ccbins = ccbins
        self.mappable_len = mappable_len
        self.mappable_forward_sum = mappable_forward_sum
        self.mappable_reverse_sum = mappable_reverse_sum
        self.mappable_ccbins = mappable_ccbins
        self.cc = np.array(cc, dtype=np.float_) if cc is not None else None
        self.cc_min = None
        self.ccrl = None
        self.library_len = expected_library_len
        self.ccfl = None
        self.nsc = None
        self.rsc = None
        self.masc = np.array(masc, dtype=np.float_) if masc is not None else None
        self.masc_min = None
        self.mascrl = None
        self.est_lib_len = None
        self.est_ccfl = None
        self.est_nsc = None
        self.est_rsc = None

        self.calc_ncc = self.calc_masc = False
        if ccbins is not None or mappable_ccbins is not None:
            self.calc_ncc = all(x is not None for x in (
                self.forward_sum,
                self.reverse_sum,
                self.ccbins
            ))
            self.calc_masc = all(x is not None for x in (
                self.mappable_forward_sum,
                self.mappable_reverse_sum,
                self.mappable_ccbins,
                self.mappable_len
            ))
            assert self.calc_ncc or self.calc_masc
            self._calc_from_bins()
        elif cc or masc:
            self.calc_ncc = cc is not None
            self.calc_masc = masc is not None
            self.max_shift = min(len(cc), len(masc)) - 1
            self.calc_metrics()

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
        if self.calc_ncc:
            self._calc_naive_cc()
        if self.calc_masc:
            self._calc_masc()

        #
        self.calc_metrics()

    def _calc_naive_cc(self):
        denom = self.genomelen - np.array(range(self.max_shift + 1), dtype=np.float_)
        self.cc = self._calc_cc(
            float(self.forward_sum),
            float(self.reverse_sum),
            self.ccbins,
            self.genomelen,
            denom
        )

    @npcalc_with_logging_warn
    def _calc_masc(self):
        totlen = np.array(self.mappable_len, dtype=np.float_)
        totlen = np.concatenate((
            totlen[:self.read_len][::-1], totlen[1:]
        ))[:self.max_shift + 1]

        self.masc = self._calc_cc(
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

    def calc_metrics(self):
        if self.cc is not None:
            self.cc_min, self.ccrl = self._calc_rl_metrics(self.cc)
            if self.library_len:
                self.ccfl, self.nsc, self.rsc = self._calc_lib_metrics(
                    self.cc, self.cc_min, self.ccrl, self.library_len
                )

        if self.masc is not None:
            self.est_lib_len = np.argmax(moving_avr_filter(self.masc, self.filter_len)) + 1
            self.masc_min, self.mascrl = self._calc_rl_metrics(self.masc)
            if self.cc is not None:
                self.est_ccfl, self.est_nsc, self.est_rsc = self._calc_lib_metrics(
                    self.cc, self.cc_min, self.ccrl, self.est_lib_len
                )

    def _calc_rl_metrics(self, cc):
        cc_min = min(cc)

        if self.read_len > self.max_shift:
            ccrl = 0
        else:
            ccrl = cc[self.read_len - 1]

        return cc_min, ccrl

    @staticmethod
    @npcalc_with_logging_warn
    def _calc_lib_metrics(cc, cc_min, ccrl, library_len):
        ccfl = cc[library_len - 1]
        nsc = ccfl / cc_min
        rsc = (ccfl - cc_min) / (ccrl - cc_min)

        return ccfl, nsc, rsc


class ReadsTooFew(IndexError):
    pass


class CCResult(object):
    def __init__(
        self, handler=None,
        filter_len=15, chi2_pval=0.05, expected_library_len=None
    ):
        # settings
        self.filter_len = filter_len
        self.chi2_p_thresh = chi2_pval
        self.expected_library_len = expected_library_len

        #
        self.references = references = handler.references
        lengths = handler.lengths
        ref2genomelen = dict(zip(references, lengths))
        read_len = handler.read_len
        ref2forward_sum = handler.ref2forward_sum
        ref2reverse_sum = handler.ref2reverse_sum
        ref2ccbins = handler.ref2ccbins
        mappable_ref2forward_sum = handler.mappable_ref2forward_sum
        mappable_ref2reverse_sum = handler.mappable_ref2reverse_sum
        mappable_ref2ccbins = handler.mappable_ref2ccbins
        ref2mappable_len = handler.ref2mappable_len

        #
        self.ref2stats = {
            ref: PyMaSCStats(
                read_len, filter_len, expected_library_len,
                genomelen=ref2genomelen[ref],
                forward_sum=ref2forward_sum.get(ref, None),
                reverse_sum=ref2reverse_sum.get(ref, None),
                ccbins=ref2ccbins.get(ref, None),
                mappable_len=ref2mappable_len.get(ref, None),
                mappable_forward_sum=mappable_ref2forward_sum.get(ref, None),
                mappable_reverse_sum=mappable_ref2reverse_sum.get(ref, None),
                mappable_ccbins=mappable_ref2ccbins.get(ref, None)
            ) for ref in references
        }

        #
        if all((ref2forward_sum, ref2reverse_sum, ref2ccbins)):
            self.skip_ncc = False
            forward_sum = sum(list(ref2forward_sum.values()))
            reverse_sum = sum(list(ref2reverse_sum.values()))
            ccbins = np.sum(_skip_none(ref2ccbins.values()), axis=0)
        else:
            self.skip_ncc = True
            forward_sum = reverse_sum = ccbins = None
        if forward_sum is not None:
            if forward_sum == 0:
                logger.error("There is no forward read.")
                raise ReadsTooFew
            elif reverse_sum == 0:
                logger.error("There is no reverse read.")
                raise ReadsTooFew
            chi2_test(forward_sum, reverse_sum, self.chi2_p_thresh, "Whole genome")

        #
        if all((mappable_ref2forward_sum, mappable_ref2reverse_sum, mappable_ref2ccbins, ref2mappable_len)):
            self.calc_masc = True
            mappable_forward_sum = np.sum(_skip_none(mappable_ref2forward_sum.values()), axis=0)
            mappable_reverse_sum = np.sum(_skip_none(mappable_ref2reverse_sum.values()), axis=0)
            mappable_ccbins = np.sum(_skip_none(mappable_ref2ccbins.values()), axis=0)
            mappable_len = np.sum(list(ref2mappable_len.values()), axis=0)
        else:
            self.calc_masc = False
            mappable_forward_sum = mappable_reverse_sum = mappable_ccbins = mappable_len = None
        if mappable_forward_sum is not None:
            if np.all(mappable_forward_sum == 0):
                logger.error("There is no forward read.")
                raise ReadsTooFew
            elif np.all(mappable_reverse_sum == 0):
                logger.error("There is no reverse read.")
                raise ReadsTooFew

        #
        self.whole = PyMaSCStats(
            read_len, filter_len, expected_library_len,
            genomelen=sum(lengths),
            forward_sum=forward_sum,
            reverse_sum=reverse_sum,
            ccbins=ccbins,
            mappable_len=mappable_len,
            mappable_forward_sum=mappable_forward_sum,
            mappable_reverse_sum=mappable_reverse_sum,
            mappable_ccbins=mappable_ccbins
        )
