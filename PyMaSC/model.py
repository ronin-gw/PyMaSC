import logging
import sys

import numpy
from scipy.stats.distributions import chi2

logger = logging.getLogger(__name__)


class ReadUnsortedError(IndexError):
    pass


class ReadsTooFew(IndexError):
    pass


class CCCalculator(object):
    chi2_p_thresh = 0.05

    need_progress_bar = False
    PROGRESS_FORMAT = "\r>{:<72}<"
    PROGRESS_BAR = "<1II1>" * 12

    def __init__(self, max_shift, references, lengths):
        """
        self._ccbins: summation bins to calc cc
        self._forward_buff: forward strand read mappability buff.
                            fixed length: `self.max_shift`
                            The last base points `self._deq_tail_pos.`
        self._reverse_buff: reverse strand read mappability buff.
                            The first base points `self._deq_tail_pos + 1`
        """
        self.max_shift = max_shift
        self.ref2genomelen = dict(zip(references, lengths))
        self.genomelen = sum(lengths)

        self._chr = None
        # self._ccbins = [0] * max_shift
        self._ccbins = numpy.array([0] * max_shift)
        self._forward_buff = numpy.array([0] * self.max_shift)
        self._reverse_buff = numpy.array([], dtype=numpy.int64)
        self._init_pos_buff()

        self.forward_read_len_sum = self.reverse_read_len_sum = 0
        self.forward_sum = self.reverse_sum = 0
        self.cc = []

    def _init_pos_buff(self):
        # self._forward_buff = [0] * self.max_shift
        self._forward_buff.fill(0)
        self._reverse_buff.fill(0)

        self._deq_tail_pos = self.max_shift
        self._last_pos = 0
        self._last_forward_pos = self._last_reverse_pos = 0

    def _flush(self):
        self._shift_with_update(len(self._reverse_buff) - 1)
        self._init_pos_buff()
        if self.need_progress_bar:
            sys.stderr.write("\r\033[K")
            sys.stderr.flush()

    def _check_pos(self, chrom, pos):
        if chrom != self._chr:
            if self._chr is not None:
                self._flush()
            self._chr = chrom

            logger.info("Sum up {}...".format(chrom))

            if self.need_progress_bar:
                self._unit = self.ref2genomelen[chrom] / len(self.PROGRESS_BAR)
                self._pg_pos = 0
                self._pg_base_pos = self._unit
                sys.stderr.write(self.PROGRESS_FORMAT.format(''))

        if pos < self._last_pos:
            raise ReadUnsortedError

        if self.need_progress_bar and pos > self._pg_base_pos:
            while pos > self._pg_base_pos:
                self._pg_pos += 1
                self._pg_base_pos += self._unit
            sys.stderr.write(self.PROGRESS_FORMAT.format(self.PROGRESS_BAR[:self._pg_pos]))

        self._last_pos = pos

    def feed_forward_read(self, chrom, pos, readlen):
        logger.debug("Got forward read {}".format(pos))
        self._check_pos(chrom, pos)
        if self._last_forward_pos == pos:
            logger.debug("Read is duplicated")
            return None
        self._last_forward_pos = pos

        logger.debug(self._reverse_buff)
        logger.debug(self._forward_buff)

        self.forward_sum += 1
        self.forward_read_len_sum += readlen
        offset = pos - self._deq_tail_pos - 1
        if offset < 0:
            logger.debug("Update buff {}".format(offset))
            self._forward_buff[offset] = 1
        else:
            logger.debug("Shift buff {}".format(offset))
            self._shift_with_update(offset, update_forward=True)

        logger.debug(self._reverse_buff)
        logger.debug(self._forward_buff)

    def feed_reverse_read(self, chrom, pos, readlen):
        logger.debug("Got reverse read {}".format(pos))
        self._check_pos(chrom, pos)
        if self._last_reverse_pos == pos:
            logger.debug("Read is duplicated")
            return None
        self._last_reverse_pos = pos

        logger.debug(self._reverse_buff)
        logger.debug(self._forward_buff)

        self.reverse_sum += 1
        self.reverse_read_len_sum += readlen
        offset = pos - self._deq_tail_pos - 1
        rev_pos = pos + readlen - 1
        rev_offset = rev_pos - self._deq_tail_pos
        if rev_offset < 1:
            logger.debug("Update immediately {}".format(offset))
            self._update_cc(rev_offset)
        else:
            if offset > 0:
                logger.debug("Shift buff {}".format(offset))
                self._shift_with_update(offset - 1)
                rev_offset -= offset
            try:
                self._reverse_buff[rev_offset - 1] = 1
            except IndexError:
                self._reverse_buff = numpy.append(
                    self._reverse_buff, [0] * (rev_offset - len(self._reverse_buff) - 1) + [1]
                )

        logger.debug(self._reverse_buff)
        logger.debug(self._forward_buff)

    def _shift_with_update(self, offset, update_forward=False):
        remain_buff_len = len(self._reverse_buff) - offset - 1

        if remain_buff_len > 0:
            for rev_state in self._reverse_buff[:offset + (1 if update_forward else 0)]:
                if rev_state:
                    self._update_cc()
                # self._forward_buff = self._forward_buff[1:] + [0]
                self._forward_buff = numpy.concatenate((self._forward_buff[1:], [0]))

            if update_forward:
                if self._reverse_buff[offset]:
                    self._update_cc()
                self._forward_buff = numpy.concatenate((self._forward_buff[1:], [1]))

            self._reverse_buff = self._reverse_buff[offset+1:]
        else:
            for rev_state in self._reverse_buff:
                if rev_state:
                    self._update_cc()
                # self._forward_buff = self._forward_buff[1:] + [0]
                self._forward_buff = numpy.concatenate((self._forward_buff[1:], [0]))

            self._forward_buff.fill(0)
            if update_forward:
                self._forward_buff[0] = 1
            self._reverse_buff = numpy.array([], dtype=numpy.int64)

        self._deq_tail_pos += offset + 1

    def _update_cc(self, offset=0):
        if offset:
            for i, state in enumerate(self._forward_buff[offset-1::-1]):
                self._ccbins[i] += state
        else:
            self._ccbins += self._forward_buff[::-1]

    def finishup_calculation(self):
        self._flush()

        if self.forward_sum == 0:
            logger.error("There is no forward read.")
            raise ReadsTooFew
        elif self.reverse_sum == 0:
            logger.error("There is no reverse read.")
            raise ReadsTooFew

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

        self.forward_read_mean_len = self.forward_read_len_sum / float(self.forward_sum)
        self.reverse_read_mean_len = self.reverse_read_len_sum / float(self.reverse_sum)
        self.read_mean_len = (
            (self.forward_read_len_sum + self.reverse_read_len_sum) /
            (float(self.forward_sum) + self.reverse_sum)
        )

        self.forward_mean = float(self.forward_sum) / self.genomelen
        self.reverse_mean = float(self.reverse_sum) / self.genomelen

        self.forward_var = self.forward_mean * (1 - self.forward_mean)
        self.reverse_var = self.reverse_mean * (1 - self.reverse_mean)

        sum_prod = self.forward_mean * self.reverse_mean
        var_geomean = (self.forward_var * self.reverse_var) ** 0.5

        # for i in xrange(self.max_shift):
        #     self.cc.append(
        #         ((self._ccbins[i] / float(self.genomelen - i - 1)) - sum_prod) / var_geomean
        #     )
        self.cc = (
            self._ccbins / (self.genomelen - numpy.array(range(self.max_shift)) - 1.) - sum_prod
        ) / var_geomean

        self.cc_min = min(self.cc)

        self.ccrl_x = int(round(self.read_mean_len))
        self.ccrl_y = self.cc[self.ccrl_x - 1]
