import logging, sys, time

from memory_profiler import profile
import numpy as np

from reader.bigwig import BigWigFile

logger = logging.getLogger(__name__)


class BWFeederWithAlignableRegoinSum(BigWigFile):
    MAPPABILITY_THRESHOLD = 1

    def __init__(self, path, max_shift=0, chrom_size=None):
        super(BWFeederWithAlignableRegoinSum, self).__init__(path, chrom_size)
        self.max_shift = max_shift
        self.chrom2is_called = {c: False for c in self.chromsizes}
        self.chrom2alignable_len = {}
        self.alignable_len = [0] * (max_shift + 1)
        self.is_called = False

        self._sum_bin_size = self.max_shift + 1
        self._chr = None

        self._sumbins = np.repeat(0, self.max_shift + 1)
        self._buff_tail_pos = 0
        self._buff = np.repeat(0, self.max_shift)

    def _init_buff(self, chrom):
        self._sumbins.fill(0)
        self._buff_tail_pos = 0
        self._buff.fill(0)

        logger.info("Process {}...".format(chrom))
        self._chr = chrom

    def _flush(self):
        if self._chr:
            self.chrom2alignable_len[self._chr] = tuple(self._sumbins)
            for i in xrange(self.max_shift + 1):
                self.alignable_len[i] += self._sumbins[i]
            self.chrom2is_called[self._chr] = True

            if all(self.chrom2is_called.values()):
                self.is_called = True

    def feed(self, chrom):
        if self.chrom2is_called[chrom]:
            return self._yield(chrom)
        else:
            return self._yield_with_calc(chrom)

    def _yield(self, chrom):
        for wig in super(BWFeederWithAlignableRegoinSum, self).fetch(chrom):
            if wig[3] < self.MAPPABILITY_THRESHOLD:
                continue
            yield wig

    def _yield_with_calc_alignability(self, chrom):
        stop_yield = False

        self._init_buff(chrom)
        for wig in super(BWFeederWithAlignableRegoinSum, self).fetch(chrom):
            if wig[3] < self.MAPPABILITY_THRESHOLD:
                continue

            self._feed_track(wig[1], wig[2])

            if not stop_yield:
                try:
                    yield wig
                except GeneratorExit:
                    stop_yield = True

        self._flush()

    def _calc_alignability(self, chrom=None):
        if chrom is None:
            chroms = [c for c, b in self.chrom2is_called.items() if b is False]
        else:
            chroms = [chrom]

        for c in chroms:
            self._init_buff(c)
            for wig in super(BWFeederWithAlignableRegoinSum, self).fetch(c):
                if wig[3] < self.MAPPABILITY_THRESHOLD:
                    continue
                self._feed_track(wig[1], wig[2])

            self._flush()

    def _feed_track(self, begin, end):
        gap_len = begin - self._buff_tail_pos - 1
        remain_buff_len = self.max_shift - gap_len
        track_len = end - begin

        if remain_buff_len > track_len:
            track = np.repeat(1, track_len)

            self._sumbins[gap_len + 1:] += np.correlate(
                self._buff[- remain_buff_len:], track, "full"
            )[:- 1 - remain_buff_len:-1]

            for i in range(track_len):
                self._sumbins[i] += track_len - i

            self._buff = np.concatenate((
                self._buff[track_len - remain_buff_len:], np.repeat(0, gap_len), track
            ))
        else:
            if remain_buff_len > 0:
                self._sumbins[gap_len + 1:] += np.correlate(
                    self._buff[- remain_buff_len:], np.repeat(1, remain_buff_len), "full"
                )[:- 1 - remain_buff_len:-1]

            if track_len < self.max_shift:
                for i in range(track_len):
                    self._sumbins[i] += track_len - i

                self._buff = np.concatenate((
                    np.repeat(0, self.max_shift - track_len), np.repeat(1, track_len)
                ))
            else:
                for i in range(self._sum_bin_size):
                    self._sumbins[i] += track_len - i

                self._buff = np.repeat(1, self.max_shift)

        self._buff_tail_pos = end - 1

    def get_alignable_len(self, chrom=None, shift_from=None, shift_to=None, force=False):
        if chrom is not None:
            if chrom not in self.chrom2is_called:
                return None
            if self.chrom2is_called[chrom]:
                return self.chrom2alignable_len[chrom][shift_from:shift_to]
            elif force:
                self._calc_alignability(chrom)
                return self.chrom2alignable_len[chrom][shift_from:shift_to]
            else:
                raise KeyError("Alignable length for '{}' is not calculated yet.".format(chrom))

        if self.is_called:
            return self.alignable_len[shift_from:shift_to]

        self._calc_alignability()
        return self.alignable_len[shift_from:shift_to]

    @profile
    def _test_del(self):
        del self._sumbins
        time.sleep(5)
        del self._buff
        time.sleep(5)
        del self.chrom2alignable_len
        time.sleep(5)
        del self.alignable_len
        time.sleep(5)
        return None
