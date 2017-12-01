import logging

import numpy as np
from scipy.stats.distributions import chi2

from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.calc import moving_avr_filter

logger = logging.getLogger(__name__)


class ReadsTooFew(IndexError):
    pass


class CCResult(object):
    @staticmethod
    def _skip_none(i):
        return [x for x in i if x is not None]

    def __init__(
        self,
        references, lengths, read_len,
        ref2forward_sum=None, ref2reverse_sum=None, ref2ccbins=None,
        mappable_ref2forward_sum=None, mappable_ref2reverse_sum=None,
        mappable_ref2ccbins=None, ref2mappable_len=None,
        filter_len=15, chi2_pval=0.05, expected_library_len=None
    ):
        #
        self.skip_ncc = not all((ref2forward_sum, ref2reverse_sum, ref2ccbins))
        self.calc_masc = all((mappable_ref2forward_sum, mappable_ref2reverse_sum,
                              mappable_ref2ccbins, ref2mappable_len))
        assert (not self.skip_ncc) or self.calc_masc

        # settings
        self.filter_len = filter_len
        self.chi2_p_thresh = chi2_pval
        self.expected_library_len = expected_library_len

        #
        ncc_max_shift = masc_max_shift = None
        mappability_max_shift = None
        if not self.skip_ncc:
            self.max_shift = ncc_max_shift = min(map(len, self._skip_none(ref2ccbins.values()))) - 1
        if self.calc_masc:
            self.max_shift = masc_max_shift = min(*[
                min(map(len, self._skip_none(d.values())))
                for d in (mappable_ref2forward_sum, mappable_ref2reverse_sum, mappable_ref2ccbins)
            ]) - 1
            mappability_max_shift = min(map(len, ref2mappable_len.values()))
        if (not self.skip_ncc) and self.calc_masc:
            self.max_shift = min(ncc_max_shift, masc_max_shift)

        #
        self.ref2genomelen = dict(zip(references, lengths))
        self.read_len = read_len
        self.genomelen = sum(lengths)

        #
        self.ref2forward_sum = ref2forward_sum
        self.ref2reverse_sum = ref2reverse_sum
        self.ref2ccbins = ref2ccbins

        self.cc = self.cc_min = self.ccrl = None
        self.ref2cc = {}
        self.ref2cc_min = {}
        self.ref2ccrl = {}

        self.ccfl = self.nsc = self.rsc = None
        self.ref2ccfl = {}
        self.ref2nsc = {}
        self.ref2rsc = {}

        #
        if not self.skip_ncc:
            self.forward_sum = sum(list(self.ref2forward_sum.values()))
            self.reverse_sum = sum(list(self.ref2reverse_sum.values()))
            self.ccbins = np.sum(self._skip_none(self.ref2ccbins.values()), axis=0)
            #
            if self.forward_sum == 0:
                logger.error("There is no forward read.")
                raise ReadsTooFew
            elif self.reverse_sum == 0:
                logger.error("There is no reverse read.")
                raise ReadsTooFew
            #
            self.test_read_balance()
            self._calc_stats()

        #
        self.mappable_ref2forward_sum = mappable_ref2forward_sum
        self.mappable_ref2reverse_sum = mappable_ref2reverse_sum
        self.mappable_ref2ccbins = mappable_ref2ccbins
        self.ref2mappable_len = ref2mappable_len

        self.masc = self.masc_min = self.mascrl = None
        self.estimated_library_len = None
        self.ref2masc = {}
        self.ref2masc_min = {}
        self.ref2mascrl = {}

        self.estimated_ccfl = self.estimated_nsc = self.estimated_rsc = None
        self.ref2library_len = {}
        self.ref2est_ccfl = {}
        self.ref2est_nsc = {}
        self.ref2est_rsc = {}

        #
        if self.calc_masc:
            #
            required_shift_size = MappabilityHandler.calc_mappable_len_required_shift_size(
                self.read_len, self.max_shift
            )
            assert required_shift_size <= mappability_max_shift
            #
            self.mappable_forward_sum = np.sum(self._skip_none(self.mappable_ref2forward_sum.values()), axis=0)
            self.mappable_reverse_sum = np.sum(self._skip_none(self.mappable_ref2reverse_sum.values()), axis=0)
            self.mappable_ccbins = np.sum(self._skip_none(self.mappable_ref2ccbins.values()), axis=0)
            self.mappable_len = np.sum(list(self.ref2mappable_len.values()), axis=0)
            #
            if np.all(self.mappable_forward_sum == 0):
                logger.error("There is no forward read.")
                raise ReadsTooFew
            elif np.all(self.mappable_reverse_sum == 0):
                logger.error("There is no reverse read.")
                raise ReadsTooFew

            #
            self._calc_masc_stats()

    def _calc_cc(self, forward_sum, reverse_sum, ccbins, totlen, denom):
        forward_mean = forward_sum / totlen
        reverse_mean = reverse_sum / totlen

        forward_var = forward_mean * (1 - forward_mean)
        reverse_var = reverse_mean * (1 - reverse_mean)

        sum_prod = forward_mean * reverse_mean
        var_geomean = (forward_var * reverse_var) ** 0.5

        cc = (ccbins / denom - sum_prod) / var_geomean

        cc_min = min(cc)

        if self.read_len > self.max_shift:
            ccrl = 0
        else:
            ccrl = cc[self.read_len - 1]

        return cc, cc_min, ccrl

    def _calc_naive_cc(self, totlen, forward_sum, reverse_sum, ccbins):
        return self._calc_cc(
            float(forward_sum), float(reverse_sum), ccbins, totlen,
            totlen - np.array(range(self.max_shift + 1), dtype=np.float_)
        )

    def _calc_masc(self, totlen, forward_sum, reverse_sum, ccbins):
        totlen = np.array(totlen, dtype=np.float_)
        totlen = np.concatenate((
            totlen[:self.read_len][::-1], totlen[1:]
        ))[:self.max_shift + 1]

        return self._calc_cc(
            np.array(forward_sum, dtype=np.float_),
            np.array(reverse_sum, dtype=np.float_),
            ccbins, totlen, totlen
        )

    @staticmethod
    def _calc_cc_stats_for_each_ref(cc, cc_min, ccrl, calcfun, genomelen,
                                    forward_sum, reverse_sum, ref2ccbins):
        for ref, ccbins in ref2ccbins.items():
            if ccbins:
                cc[ref], cc_min[ref], ccrl[ref] = calcfun(
                    genomelen[ref], forward_sum[ref], reverse_sum[ref], ccbins
                )

    @staticmethod
    def _calc_cc_metrics(cc, fragment_len, ccrl):
        cc_min = min(cc)
        ccfl = cc[fragment_len - 1]
        nsc = ccfl / cc_min
        rsc = (ccfl - cc_min) / (ccrl - cc_min)
        return ccfl, nsc, rsc

    def _calc_cc_metrics_for_each_ref(self, cc, fragment_len, ccrl, ccfl, nsc, rsc):
        if not isinstance(fragment_len, dict):
            fragment_len = {k: fragment_len for k in ccrl}
        for ref in cc:
            if cc[ref] is not None:
                ccfl[ref], nsc[ref], rsc[ref] = self._calc_cc_metrics(
                    cc[ref], fragment_len[ref], ccrl[ref]
                )

    def _calc_stats(self):
        # calc cc
        self.cc, self.cc_min, self.ccrl = self._calc_naive_cc(
            self.genomelen, self.forward_sum, self.reverse_sum, self.ccbins
        )
        self._calc_cc_stats_for_each_ref(
            self.ref2cc, self.ref2cc_min, self.ref2ccrl, self._calc_naive_cc,
            self.ref2genomelen, self.ref2forward_sum, self.ref2reverse_sum, self.ref2ccbins
        )

        #
        if self.expected_library_len:
            self.ccfl, self.nsc, self.rsc = self._calc_cc_metrics(
                self.cc, self.expected_library_len, self.ccrl
            )
            self._calc_cc_metrics_for_each_ref(
                self.ref2cc, self.expected_library_len, self.ref2ccrl,
                self.ref2ccfl, self.ref2nsc, self.ref2rsc
            )

    def _calc_masc_stats(self):
        # calc masc
        self.masc, self.masc_min, self.mascrl = self._calc_masc(
            self.mappable_len, self.mappable_forward_sum, self.mappable_reverse_sum, self.mappable_ccbins
        )
        self._calc_cc_stats_for_each_ref(
            self.ref2masc, self.ref2masc_min, self.ref2mascrl, self._calc_masc,
            self.ref2mappable_len, self.mappable_ref2forward_sum,
            self.mappable_ref2reverse_sum, self.mappable_ref2ccbins
        )

        #
        self.estimated_library_len = np.argmax(moving_avr_filter(self.masc, self.filter_len)) + 1
        if not self.skip_ncc:
            self.estimated_ccfl, self.estimated_nsc, self.estimated_rsc = self._calc_cc_metrics(
                self.cc, self.estimated_library_len, self.ccrl
            )

        #
        for ref, masc in self.ref2masc.items():
            if masc is not None:
                self.ref2library_len[ref] = np.argmax(moving_avr_filter(masc, self.filter_len)) + 1
        self._calc_cc_metrics_for_each_ref(
            self.ref2masc, self.ref2library_len, self.ref2ccrl,
            self.ref2est_ccfl, self.ref2est_nsc, self.ref2est_rsc
        )

    def _chi2_test(self, a, b, label, info_logging=False):
        sum_ = a + b
        chi2_val = (((a - sum_ / 2.) ** 2) + ((b - sum_ / 2.) ** 2)) / sum_
        chi2_p = chi2.sf(chi2_val, 1)

        if chi2_p <= self.chi2_p_thresh:
            logger.warning("{} Forward/Reverse read count imbalance.".format(label))
            logger.warning("+/- = {} / {}, Chi-squared test p-val = {} <= {}".format(
                a, b, chi2_p, self.chi2_p_thresh
            ))
        elif info_logging:
            logger.info("{} Forward/Reverse read count +/- = {} / {}".format(label, a, b))
            logger.info("Chi-squared test p-val = {} > {}".format(chi2_p, self.chi2_p_thresh))

    def test_read_balance(self):
        self._chi2_test(self.forward_sum, self.reverse_sum, "Whole genome", True)
        # for ref in self.ref2genomelen:
        #     self._chi2_test(self.ref2forward_sum[ref], self.ref2reverse_sum[ref],
        #                     ref, info_logging=False)
