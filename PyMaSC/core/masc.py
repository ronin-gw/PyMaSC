import logging
import types

import pyBigWig
import numpy as np
from scipy.stats.distributions import chi2

from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


class ReadUnsortedError(IndexError):
    pass


class ReadsTooFew(IndexError):
    pass


class CCCalculator(object):
    chi2_p_thresh = 0.05
    region_file_path = None
    region_file = None
    is_region_file_closed = False

    @classmethod
    def set_mappable_region(cls, region_file_path):
        if not region_file_path:
            logger.info("Mappable region file is not specified. Calc naive cross-correlation only.")
            cls.region_file = False
        else:
            cls.region_file_path = region_file_path
            logger.info("Open mappable region file: {}".format(region_file_path))
            cls.region_file = pyBigWig.open(region_file_path)
            cls.is_region_file_closed = False

    @classmethod
    def close_region_file(cls):
        if not cls.region_file:
            raise IOError("There is no region file to close.")
        logger.info("Cose mappable region file.")
        cls.region_file.close()
        cls.is_region_file_closed = True

    def next(self):
        pass

    @staticmethod
    def _get_shift_update_fun(obj, forward_buff_attr_name, reverse_buff_attr_name, ccbins_attr_name):
        def update_func(self, offset, update_forward=False):
            forward_buff = getattr(self, forward_buff_attr_name)
            reverse_buff = getattr(self, reverse_buff_attr_name)
            ccbins = getattr(self, ccbins_attr_name)

            remain_buff_len = forward_buff.size - offset - 1
            if remain_buff_len > 0:
                rotate = offset + (0 if update_forward else 1)
                i = -1
                for i, rev_state in enumerate(reverse_buff[:rotate]):
                    if rev_state:
                        ccbins[i:] += forward_buff[i:][::-1]
                forward_buff = np.concatenate((forward_buff[rotate:], np.repeat(0, rotate)))

                if update_forward:
                    if reverse_buff[offset:offset+1]:
                        ccbins += forward_buff[::-1]
                    forward_buff = np.concatenate((forward_buff[1:], [1]))
            else:
                for i, rev_state in enumerate(reverse_buff[:forward_buff.size]):
                    if rev_state:
                        ccbins[i:] += forward_buff[i:][::-1]
                forward_buff.fill(0)

                if update_forward:
                    forward_buff[-1] = 1

            setattr(self, forward_buff_attr_name, forward_buff)
            setattr(self, reverse_buff_attr_name, reverse_buff[offset+1:])
            setattr(self, ccbins_attr_name, ccbins)

        return types.MethodType(update_func, obj)

    def _get_pass_func(self):
        def pass_func(self, *args, **kwargs):
            pass

        return types.MethodType(pass_func, self)

    def __init__(self, max_shift, references, lengths):
        """
        self._ccbins: summation bins to calc cc
        self._forward_buff: forward strand read mappability buff.
                            fixed length: `self.max_shift`
                            The last base points `self._deq_tail_pos.`
        self._reverse_buff: reverse strand read mappability buff.
                            The first base points `self._deq_tail_pos + 1`
        """
        if self.is_region_file_closed:
            logger.warning("Mappable region file already closed. Try to reopen...")
            self.set_mappable_region(self.region_file_path)

        self.max_shift = max_shift
        self.ref2genomelen = dict(zip(references, lengths))
        self.genomelen = sum(lengths)

        self._chr = None
        # self._ccbins = [0] * max_shift
        self._ccbins = np.repeat(0, max_shift)
        self._forward_buff = np.repeat(0, self.max_shift)
        self._reverse_buff = np.array([], dtype=np.int64)

        self.forward_sum = self.reverse_sum = 0
        self.forward_read_len_sum = self.reverse_read_len_sum = 0

        self.cc = []

        self._shift_with_update_ncc = self._get_shift_update_fun(
            self, "_forward_buff", "_reverse_buff", "_ccbins"
        )

        if self.region_file:
            self._m_ccbins = np.repeat(0, max_shift)
            self._m_forward_buff = np.repeat(0, self.max_shift)
            self._m_reverse_buff = np.array([], dtype=np.int64)

            self.mappable_forward_sum = self.mappable_reverse_sum = 0
            self.mappable_forward_len_sum = self.mappable_reverse_len_sum = 0

            self.mappable_cc = []

            self._shift_with_update_mscc = self._get_shift_update_fun(
                self, "_m_forward_buff", "_m_reverse_buff", "_m_ccbins"
            )
        else:
            self._update_cc_mscc = self._get_pass_func()
            self._shift_with_update_mscc = self._get_pass_func()

        self._init_pos_buff()
        self._progress = ProgressBar()

    def _init_pos_buff(self):
        # self._forward_buff = [0] * self.max_shift
        self._forward_buff.fill(0)
        self._reverse_buff.fill(0)
        if self.region_file:
            self._m_forward_buff.fill(0)
            self._m_reverse_buff.fill(0)

        self._deq_tail_pos = self.max_shift
        self._last_pos = 0
        self._last_forward_pos = self._last_reverse_pos = 0

    def _flush(self):
        if self.region_file:
            self._shift_with_update_mscc(self._m_reverse_buff.size - 1)
        self._shift_with_update_ncc(self._reverse_buff.size - 1)
        self._init_pos_buff()
        self._progress.clean()

    def get_mappability(self, pos):
        try:
            return self.region_file.values(self._chr, pos - 1, pos)[0] == 1
        except (RuntimeError, AttributeError):
            return False

    def _check_pos(self, chrom, pos):
        if chrom != self._chr:
            if self._chr is not None:
                self._flush()
            self._chr = chrom

            logger.info("Sum up {}...".format(chrom))

            self._progress.set(self.ref2genomelen[chrom])

        if pos < self._last_pos:
            raise ReadUnsortedError

        self._progress.update(pos)

        self._last_pos = pos

    def feed_forward_read(self, chrom, pos, readlen):
        logger.debug("Got forward read {}".format(pos))
        self._check_pos(chrom, pos)
        if self._last_forward_pos == pos:
            logger.debug("Read is duplicated")
            return None
        self._last_forward_pos = pos

        logger.debug(self._reverse_buff)
        logger.debug(self._m_reverse_buff)
        logger.debug(self._m_forward_buff)
        logger.debug(self._forward_buff)

        self.forward_sum += 1
        self.forward_read_len_sum += readlen

        is_mappabile = self.get_mappability(pos)
        if is_mappabile:
            self.mappable_forward_sum += 1

        offset = pos - self._deq_tail_pos - 1
        if offset < 0:
            logger.debug("Update buff {}".format(offset))
            self._forward_buff[offset] = 1
            if is_mappabile:
                self._m_forward_buff[offset] = 1
        else:
            logger.debug("Shift buff {}".format(offset))
            self._shift_with_update(offset, update_forward=True, is_mappabile=is_mappabile)

        logger.debug(self._reverse_buff)
        logger.debug(self._m_reverse_buff)
        logger.debug(self._m_forward_buff)
        logger.debug(self._forward_buff)

    def feed_reverse_read(self, chrom, pos, readlen):
        logger.debug("Got reverse read {}".format(pos))
        self._check_pos(chrom, pos)
        if self._last_reverse_pos == pos:
            logger.debug("Read is duplicated")
            return None
        self._last_reverse_pos = pos

        logger.debug(self._reverse_buff)
        logger.debug(self._m_reverse_buff)
        logger.debug(self._m_forward_buff)
        logger.debug(self._forward_buff)

        self.reverse_sum += 1
        self.reverse_read_len_sum += readlen
        offset = pos - self._deq_tail_pos - 1
        rev_pos = pos + readlen - 1

        is_mappabile = self.get_mappability(rev_pos)
        if is_mappabile:
            self.mappable_reverse_sum += 1

        rev_offset = rev_pos - self._deq_tail_pos
        if rev_offset < 1:
            logger.debug("Update immediately {}".format(offset))
            self._ccbins[:rev_offset-1] += self._forward_buff[rev_offset-2::-1]
            if is_mappabile:
                self._m_ccbins[:rev_offset-1] += self._m_forward_buff[rev_offset-2::-1]
        else:
            if offset > 0:
                logger.debug("Shift buff {}".format(offset))
                self._shift_with_update(offset - 1, is_mappabile=is_mappabile)
                rev_offset -= offset
            try:
                self._reverse_buff[rev_offset - 1] = 1
                if is_mappabile:
                    self._m_reverse_buff[rev_offset - 1] = 1
            except IndexError:
                self._reverse_buff = np.append(
                    self._reverse_buff, [0] * (rev_offset - self._reverse_buff.size - 1) + [1]
                )
                if is_mappabile:
                    self._m_reverse_buff = np.append(
                        self._m_reverse_buff, [0] * (rev_offset - self._m_reverse_buff.size - 1) + [1]
                    )

        logger.debug(self._reverse_buff)
        logger.debug(self._m_reverse_buff)
        logger.debug(self._m_forward_buff)
        logger.debug(self._forward_buff)

    def _shift_with_update(self, offset, update_forward=False, is_mappabile=False):
        self._shift_with_update_ncc(offset, update_forward)
        if is_mappabile:
            self._shift_with_update_mscc(offset, update_forward)

        self._deq_tail_pos += offset + 1

    def _update_cc(self, rev_offset, is_mappabile=False):
        self._update_ncc(rev_offset)
        if is_mappabile:
            self._update_mscc(rev_offset)

    def finishup_calculation(self):
        self._flush()

        #
        if self.forward_sum == 0:
            logger.error("There is no forward read.")
            raise ReadsTooFew
        elif self.reverse_sum == 0:
            logger.error("There is no reverse read.")
            raise ReadsTooFew

        # Chi2 test
        self.read_count = self.forward_sum + self.reverse_sum
        self.chi2 = (((self.forward_sum - self.read_count / 2.) ** 2) / self.read_count +
                     ((self.reverse_sum - self.read_count / 2.) ** 2) / self.read_count)

        self.chi2_p = chi2.sf(self.chi2, 1)
        if self.chi2_p < self.chi2_p_thresh:
            logger.warning("Forward/Reverse read count imbalance.")
            logger.warning("+/- = {} / {}, Chi-squared test p-val = {} < {}".format(
                self.forward_sum, self.reverse_sum, self.chi2_p, self.chi2_p_thresh
            ))
        else:
            logger.info("Read counts: +/- = {} / {}, Chi-squared test p-val = {} > {}".format(
                self.forward_sum, self.reverse_sum, self.chi2_p, self.chi2_p_thresh
            ))

        # calc read len
        self.forward_read_mean_len = self.forward_read_len_sum / float(self.forward_sum)
        self.reverse_read_mean_len = self.reverse_read_len_sum / float(self.reverse_sum)
        self.read_mean_len = (
            (self.forward_read_len_sum + self.reverse_read_len_sum) /
            (float(self.forward_sum) + self.reverse_sum)
        )

        # calc cc
        self.forward_mean = float(self.forward_sum) / self.genomelen
        self.reverse_mean = float(self.reverse_sum) / self.genomelen

        self.forward_var = self.forward_mean * (1 - self.forward_mean)
        self.reverse_var = self.reverse_mean * (1 - self.reverse_mean)

        sum_prod = self.forward_mean * self.reverse_mean
        var_geomean = (self.forward_var * self.reverse_var) ** 0.5

        self.cc = (
            self._ccbins / (self.genomelen - np.array(range(self.max_shift)) - 1.) - sum_prod
        ) / var_geomean

        # calc other stats
        self.cc_min = min(self.cc)

        self.ccrl_x = int(round(self.read_mean_len))
        if self.ccrl_x > self.max_shift:
            logger.warning("Observed read lengh ({}) longer than shift size ({}).".format(
                self.read_mean_len, self.max_shift
            ))
            self.ccrl_y = 0
        else:
            self.ccrl_y = self.cc[self.ccrl_x - 1]

        # TODO: implement calculate total mappable region length and calc stats for mscc
        # if self.region_file:
        #     self._finishup_calculation_mappable()

    def _finishup_calculation_mappable(self):
        if self.mappable_forward_sum == 0:
            logger.warning("There is no forward read in mappable region.")
        elif self.mappable_reverse_sum == 0:
            logger.error("There is no reverse read in mappable region.")
        else:
            # calc cc
            self.forward_mean = float(self.forward_sum) / self.genomelen
            self.reverse_mean = float(self.reverse_sum) / self.genomelen

            self.forward_var = self.forward_mean * (1 - self.forward_mean)
            self.reverse_var = self.reverse_mean * (1 - self.reverse_mean)

            sum_prod = self.forward_mean * self.reverse_mean
            var_geomean = (self.forward_var * self.reverse_var) ** 0.5

            self.mscc = (
                self._m_ccbins / (self.genomelen - np.array(range(self.max_shift)) - 1.) - sum_prod
            ) / var_geomean
