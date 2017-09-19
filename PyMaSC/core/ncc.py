import logging

import numpy as np

logger = logging.getLogger(__name__)


class ReadUnsortedError(IndexError):
    pass


class NaiveCCCalculator(object):
    def __init__(self, max_shift, references, lengths):
        """
        self._ccbins: summation bins to calc cc; shift from 0 to max_shift
        self._forward_buff: forward strand read position buff.
                            fixed length: `self.max_shift + 1`
                            `self._fb_tail_pos` points the last element position as 1-based position.
        self._reverse_buff: reverse strand read position buff.
                            reverse position means `position + read_length - 1`.
                            `self._fb_tail_pos - self.max_shift` points the first element position as 1-based position.
        """
        self.max_shift = max_shift
        self.ref2genomelen = dict(zip(references, lengths))
        self.genomelen = sum(lengths)
        # stats
        self.forward_sum = self.reverse_sum = 0
        self.ref2forward_sum = {ref: 0 for ref in references}
        self.ref2reverse_sum = {ref: 0 for ref in references}
        self.forward_read_len_sum = self.reverse_read_len_sum = 0
        self.ccbins = []
        self.ref2ccbins = {ref: None for ref in references}
        # internal buff
        self._chr = None
        self._forward_sum = self._reverse_sum = 0
        self._ccbins = np.zeros(max_shift + 1, dtype=np.int64)
        self._forward_buff = np.zeros(max_shift + 1, dtype=np.int64)
        self._reverse_buff = np.zeros(max_shift + 1, dtype=np.int64)

        self._init_pos_buff()

    def _init_pos_buff(self):
        self._forward_sum = self._reverse_sum = 0
        self._ccbins.fill(0)

        self._forward_buff.fill(0)
        self._fb_tail_pos = self.max_shift + 1

        self._reverse_buff.fill(0)

        self._last_pos = 0
        self._last_forward_pos = self._last_reverse_pos = 0

    def _flush(self):
        self._shift_with_update(self.max_shift + 1)

        self.forward_sum += self._forward_sum
        self.reverse_sum += self._reverse_sum
        self.ref2forward_sum[self._chr] = self._forward_sum
        self.ref2reverse_sum[self._chr] = self._reverse_sum
        self.ref2ccbins[self._chr] = tuple(self._ccbins)

    def _check_pos(self, chrom, pos):
        if chrom != self._chr:
            if self._chr is not None:
                self._flush()
            self._chr = chrom
            self._init_pos_buff()

            logger.info("Calc {}...".format(chrom))

        if pos < self._last_pos:
            raise ReadUnsortedError

        self._last_pos = pos

    def feed_forward_read(self, chrom, pos, readlen):
        self._check_pos(chrom, pos)

        if self._last_forward_pos == pos:
            return None  # duplicated read
        self._last_forward_pos = pos

        self._forward_sum += 1
        self.forward_read_len_sum += readlen

        offset = pos - self._fb_tail_pos

        if offset < 1:
            self._forward_buff[offset - 1] = 1
        else:
            self._shift_with_update(offset, update_forward=True)

    def feed_reverse_read(self, chrom, pos, readlen):
        self._check_pos(chrom, pos)

        if self._last_reverse_pos == pos:
            return None  # duplicated read
        self._last_reverse_pos = pos

        self._reverse_sum += 1
        self.reverse_read_len_sum += readlen
        offset = pos - self._fb_tail_pos
        revbuff_pos = self.max_shift + readlen - 1

        if offset > 0:
            self._shift_with_update(offset)
        else:
            revbuff_pos += offset
        try:
            self._reverse_buff[revbuff_pos] = 1
        except IndexError:
            self._reverse_buff = np.append(
                self._reverse_buff, [0] * (revbuff_pos - self._reverse_buff.size) + [1]
            )

    def _shift_with_update(self, offset, update_forward=False):
        # logger.debug(tuple(self._ccbins))

        # impl 1
        # b = np.zeros(self.max_shift + 1, dtype=np.int64)
        for i, forward_state in enumerate(self._forward_buff[:offset]):
            if forward_state:
                for j, reverse_state in enumerate(self._reverse_buff[i:i + self.max_shift + 1]):
                    if reverse_state:
                        self._ccbins[j] += 1
                        # b[j] += 1

        # impl 2
        # shift = min(offset, self.max_shift + 1)
        # self._ccbins += np.correlate(
        #     self._forward_buff[:shift], self._reverse_buff[:self.max_shift + shift], "full"
        # )[- shift:- shift - self.max_shift - 1:-1]
        # a = np.correlate(
        #     self._forward_buff[:shift], self._reverse_buff[:self.max_shift + shift], "full"
        # )[- shift:- shift - self.max_shift - 1:-1]
        #
        # self._ccbins += a

        self._forward_buff = np.concatenate((
            self._forward_buff[offset:],
            np.zeros(min(self.max_shift + 1, offset), dtype=np.int64)
        ))
        if update_forward:
            self._forward_buff[-1] = 1
        self._reverse_buff = self._reverse_buff[offset:]
        # self._reverse_buff = np.concatenate((
        #     self._reverse_buff[offset:],
        #     np.zeros(max(self.max_shift + 1 - (self._reverse_buff.size - offset), 0), dtype=np.int64)
        # ))

        self._fb_tail_pos += offset

        # logger.debug(tuple(self._ccbins))

    def finishup_calculation(self):
        self._flush()

        for bins in zip(*[v for v in self.ref2ccbins.values() if v]):
            self.ccbins.append(sum(bins))
