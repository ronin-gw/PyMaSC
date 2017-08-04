import logging

import numpy as np
from scipy.stats.distributions import chi2

from PyMaSC.utils.calc import moving_avr_filter

logger = logging.getLogger(__name__)


class ReadsTooFew(IndexError):
    pass


class CCResult(object):
    CALCULATOR_RESULT_ATTRS = (
        "max_shift", "calc_masc", "chi2_p_thresh",
        "genomelen", "ref2genomelen",

        "forward_sum", "reference2forward_sum", "forward_read_len_sum",
        "reverse_sum", "reference2reverse_sum", "reverse_read_len_sum",
        "ccbins", "reference2ccbins",

        "mappable_forward_sum", "mappable_reference2forward_sum", "mappable_forward_len_sum",
        "mappable_reverse_sum", "mappable_reference2reverse_sum", "mappable_reverse_len_sum",
        "mappable_ccbins", "mappable_reference2ccbins",
    )

    def __init__(self, ccc, alignpath, mapq_criteria, mappability, filter_len, expected_library_len=None):
        self.alignfile = alignpath
        self.mapq_criteria = mapq_criteria
        self.filter_len = filter_len
        self.expected_library_len = None

        #
        for attr in self.CALCULATOR_RESULT_ATTRS:
            setattr(self, attr, getattr(ccc, attr))

        #
        if mappability is None:
            self.regionfile = self.alignable_len = self.ref2library_len = None
        else:
            self.regionfile = mappability.path
            assert mappability.is_called
            assert self.max_shift <= mappability.max_shift
            self.alignable_len = mappability.alignable_len
            self.ref2alignable_len = mappability.chrom2alignable_len

        #
        if ccc.forward_sum == 0:
            logger.error("There is no forward read.")
            raise ReadsTooFew
        elif ccc.reverse_sum == 0:
            logger.error("There is no reverse read.")
            raise ReadsTooFew

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

        cc = (ccbins / (totlen - np.array(range(self.max_shift)) - 1.) - sum_prod) / var_geomean

        cc_min = min(cc)

        if self.estimated_read_len > self.max_shift:
            ccrl = 0
        else:
            ccrl = cc[self.estimated_read_len - 1]

        return cc, cc_min, ccrl

    def _calc_stats(self):
        # calc read len
        self.forward_read_mean_len = self.forward_read_len_sum / float(self.forward_sum)
        self.reverse_read_mean_len = self.reverse_read_len_sum / float(self.reverse_sum)
        self.read_mean_len = (
            (self.forward_read_len_sum + self.reverse_read_len_sum) /
            (float(self.forward_sum) + self.reverse_sum)
        )
        self.estimated_read_len = int(round(self.read_mean_len))
        if self.estimated_read_len > self.max_shift:
            logger.warning("Observed read lengh ({}) longer than shift size ({}).".format(
                self.read_mean_len, self.max_shift
            ))

        # calc cc
        self.cc, self.cc_min, self.ccrl = self._calc_cc(
            self.genomelen, self.forward_sum, self.reverse_sum, self.ccbins
        )
        self.ref2cc = {}
        self.ref2cc_min = {}
        self.ref2ccrl = {}
        for ref, ccbins in self.ref2ccbins.items():
            if ccbins:
                self.ref2cc[ref], self.ref2cc_min[ref], self.ref2ccrl[ref] = self._calc_cc(
                    self.ref2genomelen[ref], self.ref2forward_sum[ref], self.ref2reverse_sum[ref], ccbins
                )
            else:
                self.ref2cc[ref] = self.ref2cc_min[ref] = self.ref2ccrl[ref] = None

        if self.expected_library_len:
            self.ccfl = self.cc[self.expected_library_len - 1]
            self.nsc = self.ccfl / self.cc_min
            self.rsc = (self.ccfl - self.cc_min) / (self.ccrl - self.cc_min)

            self.ref2ccfl = {}
            self.ref2nsc = {}
            self.ref2rsc = {}
            for ref, ccbins in self.ref2ccbins.items():
                self.ref2ccfl[ref] = self.ref2cc[self.expected_library_len - 1]
                self.ref2nsc[ref] = self.ref2ccfl[ref] / self.ref2cc_min[ref]
                self.ref2rsc[ref] = ((self.ref2ccfl[ref] - self.ref2cc_min[ref]) /
                                     (self.ref2ccrl[ref] - self.ref2cc_min[ref]))
        else:
            self.ccfl = self.nsc = self.rsc = None
            self.ref2ccfl = self.ref2nsc = self.ref2rsc = None

    def _calc_masc(self, totlen, forward_sum, reverse_sum, ccbins):
        totlen = np.array(totlen[1:], dtype=float)
        forward_sum = np.array(forward_sum, dtype=float)
        reverse_sum = np.array(reverse_sum, dtype=float)
        forward_mean = forward_sum / totlen
        reverse_mean = reverse_sum / totlen

        forward_var = forward_mean * (1 - forward_mean)
        reverse_var = reverse_mean * (1 - reverse_mean)

        sum_prod = forward_mean * reverse_mean
        var_geomean = (forward_var * reverse_var) ** 0.5

        cc = (ccbins / totlen - sum_prod) / var_geomean

        cc_min = min(cc)

        if self.estimated_read_len > self.max_shift:
            ccrl = 0
        else:
            ccrl = cc[self.estimated_read_len]

        return cc, cc_min, ccrl

    def _calc_masc_stats(self):
        # calc masc
        self.masc, self.masc_min, self.mascrl_y = self._calc_masc(
            self.alignable_len, self.mappable_forward_sum, self.mappable_reverse_sum, self.mappable_ccbins
        )
        self.ref2masc = {}
        self.ref2masc_min = {}
        self.ref2mascrl = {}
        for ref, ccbins in self.mappable_reference2ccbins.items():
            if ccbins is not None and ref in self.ref2alignable_len:
                self.ref2masc[ref], self.ref2masc_min[ref], self.ref2mascrl[ref] = self._calc_masc(
                    self.ref2alignable_len[ref], self.mappable_reference2forward_sum[ref],
                    self.mappable_reference2reverse_sum[ref], ccbins
                )
            else:
                self.ref2masc[ref] = self.ref2masc_min[ref] = self.ref2mascrl[ref] = None

        # calc additional stats
        self.estimated_library_len = np.argmax(moving_avr_filter(self.masc, self.filter_len)) + 1
        self.estimated_ccfl = self.cc[self.estimated_library_len - 1]
        self.estimated_nsc = self.estimated_ccfl / self.cc_min
        self.estimated_rsc = (self.estimated_ccfl - self.cc_min) / (self.ccrl - self.cc_min)

        self.ref2library_len = {}
        self.ref2est_ccfl = {}
        self.ref2est_nsc = {}
        self.ref2est_rsc = {}
        for ref, masc in self.ref2masc.items():
            if masc is not None:
                self.ref2library_len[ref] = np.argmax(moving_avr_filter(masc, self.filter_len)) + 1
                self.ref2est_ccfl[ref] = self.ref2cc[ref][self.ref2library_len[ref] - 1]
                self.ref2est_nsc[ref] = self.ref2est_ccfl[ref] / self.ref2cc_min[ref]
                self.ref2est_rsc[ref] = ((self.ref2est_ccfl[ref] - self.ref2cc_min[ref]) /
                                         (self.ref2ccrl[ref] - self.ref2cc_min[ref]))
            else:
                self.ref2library_len[ref] = None
                self.ref2est_ccfl[ref] = None
                self.ref2est_nsc[ref] = None
                self.ref2est_rsc[ref] = None

    def _chi2_test(self, a, b, label):
        sum_ = a + b
        chi2_val = (((a - sum_ / 2.) ** 2) + ((b - sum_ / 2.) ** 2)) / sum_
        chi2_p = chi2.sf(chi2_val, 1)

        if chi2_p < self.chi2_p_thresh:
            logger.warning("{} Forward/Reverse read count imbalance.".format(label))
            logger.warning("+/- = {} / {}, Chi-squared test p-val = {} < {}".format(
                a, b, chi2_p, self.chi2_p_thresh
            ))

    def test_read_balance(self):
        self._chi2_test(self.forward_sum, self.reverse_sum, "Whole genome")
        for ref in self.ref2genomelen:
            self._chi2_test(self.ref2forward_sum[ref],
                            self.ref2reverse_sum[ref],
                            ref)
        if self.calc_masc:
            self._chi2_test(self.mappable_forward_sum, self.mappable_reverse_sum, "Mappable region")
            for ref in self.ref2genomelen:
                self._chi2_test(self.mappable_reference2forward_sum[ref],
                                self.mappable_reference2reverse_sum[ref],
                                ref + " mappable region")
