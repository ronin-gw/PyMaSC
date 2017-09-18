import logging

import numpy as np

from PyMaSC.core.ncc import NaiveCCCalculator, ReadUnsortedError
from PyMaSC.core.alignability import ContinueCalculation
from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


class MSCCCalculator(NaiveCCCalculator):
    def __init__(self, max_shift, read_len, references, lengths, bwfeeder):
        #    _last_pos      pos + read_len - 1
        #        |-----------------| read_len
        #        /=================>
        #    <---------------------| max_shift
        #    |<--| _max_shift_from_f

        self.max_shift = max_shift
        self.ref2genomelen = dict(zip(references, lengths))
        self.genomelen = sum(lengths)
        self.read_len = read_len
        # _max_shift_from_f = relative pos from forward read pos to max_shift
        self._max_shift_from_f = -(read_len - 1 - max_shift)

        # stats
        # {forward,reverse}_sum[d[0...i]] = Number of reads in Md
        self.forward_sum = np.repeat(0, max_shift + 1)
        self.reverse_sum = np.repeat(0, max_shift + 1)
        self.ref2forward_sum = {ref: None for ref in references}
        self.ref2reverse_sum = {ref: None for ref in references}
        #
        self.forward_read_len_sum = self.reverse_read_len_sum = 0
        # ccbins[d[1...i]] = Number of doubly mapped region in Md
        self.ccbins = []
        self.ref2ccbins = {ref: None for ref in references}
        # internal buff
        self._chr = None
        self._forward_sum = np.zeros(max_shift + 1, dtype=np.int64)
        self._reverse_sum = np.zeros(max_shift + 1, dtype=np.int64)
        self._ccbins = np.zeros(max_shift + 1, dtype=np.int64)
        # None for read absence
        self._forward_buff = [None for _ in range(self.max_shift + 1)]
        self._reverse_buff = [None for _ in range(self.max_shift + 1)]

        # mappability
        self.bwfeeder = bwfeeder
        self.bwfeeder.disable_progress_bar()
        self._bwiter_stopped = False
        # _bwbuff[i] = mappability of pos i
        self._bwbuff = np.array([], dtype=np.int64)
        # _bwbuff_head = 1-based position for first element in _bwbuff
        self._bwbuff_head = 1
        # _extra_size = extra buff size from _last_pos save to _bwbuff
        self._extra_size = max(0, self._max_shift_from_f)
        # print self._max_shift_from_f, self._extra_size

        self._init_pos_buff()
        self._progress = ProgressBar()
        self._progress.disable_bar()
        self.bwfeeder.disable_progress_bar()

    def _init_pos_buff(self):
        self._forward_sum.fill(0)
        self._reverse_sum.fill(0)
        self._ccbins.fill(0)

        for i in range(len(self._forward_buff)):
            self._forward_buff[i] = None
        self._fb_tail_pos = self.max_shift + 1

        for i in range(len(self._reverse_buff)):
            self._reverse_buff[i] = None

        self._last_pos = 0
        self._last_forward_pos = self._last_reverse_pos = 0

    def _init_bw(self):
        self._bwbuff = np.array([], dtype=np.int64)
        self._bwbuff_head = 1
        self._bwiter_stopped = False
        try:
            self._feeder = self.bwfeeder.feed(self._chr)
        except KeyError as e:
            logger.info("Mappability for '{}' not found. "
                        "Skip calc mappability sensitive CC.".format(e.args[0]))
            self._bwiter_stopped = True

    def _flush(self):
        if self._reverse_buff:
            self._shift_with_update(self.max_shift + 1)

        self.forward_sum += self._forward_sum
        self.reverse_sum += self._reverse_sum
        self.ref2forward_sum[self._chr] = tuple(self._forward_sum)
        self.ref2reverse_sum[self._chr] = tuple(self._reverse_sum)
        self.ref2ccbins[self._chr] = tuple(self._ccbins)

        # update bwfeeder
        if not self._bwiter_stopped:
            try:
                self._feeder.throw(ContinueCalculation)
            except StopIteration:
                pass

        self._init_pos_buff()

    def _check_pos(self, chrom, pos):
        if chrom != self._chr:
            if self._chr is None:
                self._init_pos_buff()
            else:
                self._flush()
            self._chr = chrom
            self._init_bw()

            logger.info("Calc {}...".format(chrom))

            self._progress.set(self.ref2genomelen[chrom])

        if pos < self._last_pos:
            raise ReadUnsortedError

        self._progress.update(pos)

        self._last_pos = pos

    def feed_forward_read(self, chrom, pos, readlen):
        self._check_pos(chrom, pos)

        if self._last_forward_pos == pos:
            return None
        self._last_forward_pos = pos

        offset = pos - self._fb_tail_pos

        try:
            read_mappability, mappability = self._get_forward_mappabiility()
        except StopIteration:
            return None
        if not read_mappability:
            return None

        # logger.debug("Got forward read {}".format(pos))
        self._forward_sum += mappability
        self.forward_read_len_sum += readlen

        if offset < 1:
            self._forward_buff[offset - 1] = mappability
        else:
            self._shift_with_update(offset, forward_mappability=mappability)

    def feed_reverse_read(self, chrom, pos, readlen):
        self._check_pos(chrom, pos)

        if self._last_reverse_pos == pos:
            return None
        self._last_reverse_pos = pos

        try:
            mappability = self._get_reverse_read_mappability()
        except StopIteration:
            return None

        # logger.debug("Got reverse read {} ({})".format(pos, offset))
        self._reverse_sum += mappability
        self.reverse_read_len_sum += readlen
        offset = pos - self._fb_tail_pos
        revbuff_pos = self.max_shift + self.read_len - 1

        if offset > 0:
            self._shift_with_update(offset)
        else:
            revbuff_pos += offset
        try:
            # self._reverse_buff[revbuff_pos] = mappability
            self._reverse_buff[revbuff_pos] = 1
        except IndexError:
            # self._reverse_buff += [None] * (revbuff_pos - len(self._reverse_buff)) + [mappability]
            self._reverse_buff += [None] * (revbuff_pos - len(self._reverse_buff)) + [1]

    def _get_forward_mappabiility(self):
        to = self._last_pos + self.read_len - 1
        from_ = min(to - self.max_shift - 1, self._last_pos - 1)

        try:
            if from_ < 0:
                mappability = np.concatenate((
                    np.zeros(-from_, dtype=np.int64), self._get_bwbuff(0, to)
                ))
            else:
                mappability = self._get_bwbuff(from_, to)
        except StopIteration:
            return 0, np.zeros(self.max_shift + 1, dtype=np.int64)

        return mappability[self._extra_size], mappability[self.max_shift::-1]

    def _get_reverse_read_mappability(self):
        forward_to = self._last_pos + self.read_len - 1
        forward_from = forward_to - self.max_shift - 1

        if forward_from < 0:
            forward_mappability = np.concatenate((
                np.zeros(-forward_from, dtype=np.int64), self._get_bwbuff(0, forward_to)
            ))[::-1]
        else:
            forward_mappability = self._get_bwbuff(forward_from, forward_to)[::-1]

        reverse_to = forward_to + self.read_len - 1
        reverse_from = reverse_to - self.max_shift * 2 - 1

        if reverse_from < 0:
            reverse_mappability = np.concatenate((
                np.zeros(-reverse_from, dtype=np.int64),
                self._get_bwbuff(0, reverse_to)
            ))[::-2]
        else:
            reverse_mappability = self._get_bwbuff(reverse_from, reverse_to)
            if reverse_mappability.size != self.max_shift * 2 + 1:
                reverse_mappability = np.concatenate((
                    reverse_mappability,
                    np.zeros(self.max_shift*2 + 1 - reverse_mappability.size, dtype=np.int64)
                ))[::-2]
            else:
                reverse_mappability = reverse_mappability[::-2]

        return forward_mappability * reverse_mappability

    def _get_bwbuff(self, from_, to):
        if self._bwiter_stopped and not self._bwbuff.size:
            raise StopIteration

        bwbuff_head = max(0, self._last_pos - self.max_shift * 2 - 1)
        self._bwbuff = self._bwbuff[bwbuff_head - self._bwbuff_head:]
        self._bwbuff_head = bwbuff_head

        from_ -= bwbuff_head
        to -= bwbuff_head

        while self._bwbuff.size < to:
            try:
                _chrom, start, end, _val = next(self._feeder)
            except StopIteration:
                self._bwiter_stopped = True
                m = self._bwbuff[from_:to]
                return np.concatenate((
                    m, np.zeros(max(0, self.max_shift - m.size + 1), dtype=np.int64)
                ))

            bwbuff_tail = bwbuff_head + self._bwbuff.size
            if bwbuff_tail <= start:
                gap = start - bwbuff_tail
                self._bwbuff = np.concatenate((
                    self._bwbuff, np.zeros(gap, dtype=np.int64),
                    np.ones(end - start, dtype=np.int64)
                ))
            elif bwbuff_head < end:
                self._bwbuff[max(0, start - bwbuff_head):(end - bwbuff_head)] = 1
                self._bwbuff = np.concatenate(
                    (self._bwbuff, np.ones(max(0, end - bwbuff_tail), dtype=np.int64))
                )

        return self._bwbuff[from_:to]

    def _shift_with_update(self, offset, forward_mappability=None):
        for i, forward_state in enumerate(self._forward_buff[:offset]):
            if forward_state is None:
                continue

            for j, reverse_state in enumerate(self._reverse_buff[i:i + self.max_shift + 1]):
                if reverse_state is None:
                    continue

                # assert forward_state[j] == reverse_state[j]
                self._ccbins[j] += forward_state[j]

        self._forward_buff = self._forward_buff[offset:] + [None] * min(offset, self.max_shift + 1)
        if forward_mappability is not None:
            self._forward_buff[-1] = forward_mappability
        self._reverse_buff = self._reverse_buff[offset:]
        self._fb_tail_pos += offset

    def finishup_calculation(self):
        super(MSCCCalculator, self).finishup_calculation()

        self.bwfeeder.enable_progress_bar()
        if not self._bwiter_stopped:
            try:
                self.bwfeeder.throw(ContinueCalculation)
            except StopIteration:
                pass
        self.bwfeeder._calc_alignability()
