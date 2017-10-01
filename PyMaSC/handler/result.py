import logging

import numpy as np
from scipy.stats.distributions import chi2

from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.calc import moving_avr_filter

logger = logging.getLogger(__name__)


class ReadsTooFew(IndexError):
    pass


class CCResult(object):
    default_chi2_p_thresh = 0.05

    CALCULATOR_RESULT_ATTRS = (
        "ref2forward_sum", "ref2reverse_sum", "ref2ccbins",
    )

    @staticmethod
    def _skip_none(i):
        return [x for x in i if x is not None]

    def __init__(self, calc_handler, filter_len, chi2_pval, expected_library_len=None, skip_ncc=False):
        self.filter_len = filter_len
        self.chi2_p_thresh = chi2_pval if chi2_pval else self.default_chi2_p_thresh
        self.expected_library_len = expected_library_len
        self.skip_ncc = skip_ncc

        #
        self.alignfile = calc_handler.path
        self.read_len = calc_handler.read_len
        self.mapq_criteria = calc_handler.mapq_criteria

        self.max_shift = calc_handler.max_shift
        self.ref2genomelen = dict(zip(calc_handler.references, calc_handler.lengths))
        self.genomelen = sum(calc_handler.lengths)

        for attr in self.CALCULATOR_RESULT_ATTRS:
            setattr(self, attr, getattr(calc_handler, attr))
        self.forward_sum = sum(list(self.ref2forward_sum.values()))
        self.reverse_sum = sum(list(self.ref2reverse_sum.values()))
        self.ccbins = np.sum(self._skip_none(self.ref2ccbins.values()), axis=0)

        ###
        self.estimated_library_len = None

        #
        if calc_handler.mappability_handler:
            self.calc_masc = True
            #
            self.regionfile = calc_handler.mappability_handler.path
            required_size = MappabilityHandler.calc_mappable_len_required_shift_size(self.read_len, self.max_shift)
            assert required_size <= calc_handler.mappability_handler.max_shift

            self.ref2mappable_len = calc_handler.mappability_handler.chrom2mappable_len
            for attr in self.CALCULATOR_RESULT_ATTRS:
                setattr(self, "mappable_" + attr, getattr(calc_handler, "mappable_" + attr))
            self.mappable_len = np.sum(list(self.ref2mappable_len.values()), axis=0)
            self.mappable_forward_sum = np.sum(self._skip_none(self.mappable_ref2forward_sum.values()), axis=0)
            self.mappable_reverse_sum = np.sum(self._skip_none(self.mappable_ref2reverse_sum.values()), axis=0)
            self.mappable_ccbins = np.sum(self._skip_none(self.mappable_ref2ccbins.values()), axis=0)
        else:
            self.calc_masc = False
            self.regionfile = None
            self.mappable_forward_sum = np.zeros(1)
            self.mappable_reverse_sum = np.zeros(1)

        #
        if self.forward_sum == 0 and np.all(self.mappable_forward_sum == 0):
            logger.error("There is no forward read.")
            raise ReadsTooFew
        elif self.reverse_sum == 0 and np.all(self.mappable_reverse_sum == 0):
            logger.error("There is no reverse read.")
            raise ReadsTooFew

        #
        self.cc = self.cc_min = self.ccrl = None
        self.ref2cc = {}
        self.ref2cc_min = {}
        self.ref2ccrl = {}

        self.ccfl = self.nsc = self.rsc = None
        self.ref2ccfl = {}
        self.ref2nsc = {}
        self.ref2rsc = {}

        #
        self.masc = self.masc_min = self.mascrl = None
        self.ref2masc = {}
        self.ref2masc_min = {}
        self.ref2mascrl = {}

        self.estimated_ccfl = self.estimated_nsc = self.estimated_rsc = None
        self.ref2library_len = {}
        self.ref2est_ccfl = {}
        self.ref2est_nsc = {}
        self.ref2est_rsc = {}

        self.test_read_balance()

        if not self.skip_ncc:
            self._calc_stats()
        if self.calc_masc:
            self._calc_masc_stats()

    def _calc_cc(self, totlen, forward_sum, reverse_sum, ccbins):
        forward_mean = float(forward_sum) / totlen
        reverse_mean = float(reverse_sum) / totlen

        forward_var = forward_mean * (1 - forward_mean)
        reverse_var = reverse_mean * (1 - reverse_mean)

        sum_prod = forward_mean * reverse_mean
        var_geomean = (forward_var * reverse_var) ** 0.5

        cc = (ccbins / (totlen - np.array(range(self.max_shift + 1), dtype=np.float_)) - sum_prod) / var_geomean

        cc_min = min(cc)

        if self.read_len > self.max_shift:
            ccrl = 0
        else:
            ccrl = cc[self.read_len - 1]

        return cc, cc_min, ccrl

    def _calc_stats(self):
        # calc cc
        self.cc, self.cc_min, self.ccrl = self._calc_cc(
            self.genomelen, self.forward_sum, self.reverse_sum, self.ccbins
        )

        for ref, ccbins in self.ref2ccbins.items():
            if ccbins:
                self.ref2cc[ref], self.ref2cc_min[ref], self.ref2ccrl[ref] = self._calc_cc(
                    self.ref2genomelen[ref], self.ref2forward_sum[ref], self.ref2reverse_sum[ref], ccbins
                )

        if self.expected_library_len:
            self.ccfl = self.cc[self.expected_library_len - 1]
            self.nsc = self.ccfl / self.cc_min
            self.rsc = (self.ccfl - self.cc_min) / (self.ccrl - self.cc_min)

            for ref, ccbins in self.ref2ccbins.items():
                self.ref2ccfl[ref] = self.ref2cc[self.expected_library_len - 1]
                self.ref2nsc[ref] = self.ref2ccfl[ref] / self.ref2cc_min[ref]
                self.ref2rsc[ref] = ((self.ref2ccfl[ref] - self.ref2cc_min[ref]) /
                                     (self.ref2ccrl[ref] - self.ref2cc_min[ref]))

    def _calc_masc(self, totlen, forward_sum, reverse_sum, ccbins):
        totlen = np.array(totlen, dtype=np.float_)
        totlen = np.concatenate((
            totlen[:self.read_len][::-1], totlen[1:]
        ))[:self.max_shift + 1]

        forward_sum = np.array(forward_sum, dtype=np.float_)
        reverse_sum = np.array(reverse_sum, dtype=np.float_)
        forward_mean = forward_sum / totlen
        reverse_mean = reverse_sum / totlen

        forward_var = forward_mean * (1 - forward_mean)
        reverse_var = reverse_mean * (1 - reverse_mean)

        sum_prod = forward_mean * reverse_mean
        var_geomean = (forward_var * reverse_var) ** 0.5

        cc = (ccbins / totlen - sum_prod) / var_geomean

        cc_min = min(cc)

        if self.read_len > self.max_shift:
            ccrl = 0
        else:
            ccrl = cc[self.read_len]

        return cc, cc_min, ccrl

    def _calc_masc_stats(self):
        # calc masc
        self.masc, self.masc_min, self.mascrl = self._calc_masc(
            self.mappable_len, self.mappable_forward_sum, self.mappable_reverse_sum, self.mappable_ccbins
        )

        for ref, ccbins in self.mappable_ref2ccbins.items():
            if ccbins is not None and ref in self.ref2mappable_len:
                self.ref2masc[ref], self.ref2masc_min[ref], self.ref2mascrl[ref] = self._calc_masc(
                    self.ref2mappable_len[ref], self.mappable_ref2forward_sum[ref],
                    self.mappable_ref2reverse_sum[ref], ccbins
                )

        # calc additional stats
        self.estimated_library_len = np.argmax(moving_avr_filter(self.masc, self.filter_len)) + 1
        if not self.skip_ncc:
            self.estimated_ccfl = self.cc[self.estimated_library_len - 1]
            self.estimated_nsc = self.estimated_ccfl / self.cc_min
            self.estimated_rsc = (self.estimated_ccfl - self.cc_min) / (self.ccrl - self.cc_min)

        for ref, masc in self.ref2masc.items():
            if masc is not None:
                self.ref2library_len[ref] = np.argmax(moving_avr_filter(masc, self.filter_len)) + 1

            if masc is not None and not self.skip_ncc:
                self.ref2est_ccfl[ref] = self.ref2cc[ref][self.ref2library_len[ref] - 1]
                self.ref2est_nsc[ref] = self.ref2est_ccfl[ref] / self.ref2cc_min[ref]
                self.ref2est_rsc[ref] = ((self.ref2est_ccfl[ref] - self.ref2cc_min[ref]) /
                                         (self.ref2ccrl[ref] - self.ref2cc_min[ref]))

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
        if not self.skip_ncc:
            self._chi2_test(self.forward_sum, self.reverse_sum, "Whole genome", True)
            # for ref in self.ref2genomelen:
            #     self._chi2_test(self.ref2forward_sum[ref], self.ref2reverse_sum[ref],
            #                     ref, info_logging=False)
