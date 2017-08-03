import logging

import numpy as np

from PyMaSC.reader.bigwig.bigwig import BigWigFile
from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


class AlignableLengthCalculator(BigWigFile):
    MAPPABILITY_THRESHOLD = 1

    def __init__(self, path, max_shift=0, chrom_size=None):
        super(AlignableLengthCalculator, self).__init__(path, chrom_size)
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

        self._progress = ProgressBar()

    def enable_progress_bar(self):
        self._progress.enable_bar()

    def disable_progress_bar(self):
        self._progress.disable_bar()

    def _init_buff(self, chrom):
        self._sumbins.fill(0)
        self._buff_tail_pos = 0
        self._buff.fill(0)

        logger.info("Calculate mappable length for {}...".format(chrom))
        self._chr = chrom
        self._progress.set(self.chromsizes[chrom])

    def _flush(self):
        if self._chr:
            self.chrom2alignable_len[self._chr] = tuple(self._sumbins)
            for i in xrange(self.max_shift + 1):
                self.alignable_len[i] += self._sumbins[i]
            self.chrom2is_called[self._chr] = True

            if all(self.chrom2is_called.values()):
                self.is_called = True

            self._progress.clean()

    def _calc_alignability(self, chrom=None):
        if chrom is None:
            chroms = [c for c, b in self.chrom2is_called.items() if b is False]
        elif self.chrom2is_called[chrom]:
            return None
        else:
            chroms = [chrom]

        for c in chroms:
            self._init_buff(c)
            for wig in self.fetch(self.MAPPABILITY_THRESHOLD, c):
                self._progress.update(wig[2])
                self._feed_track(wig[1], wig[2])

            self._flush()

    def _feed_track(self, begin, end):
        gap_len = begin - self._buff_tail_pos
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


class BWFeederWithAlignableRegionSum(AlignableLengthCalculator):
    def feed(self, chrom):
        if self.chrom2is_called[chrom]:
            return self._yield(chrom)
        else:
            return self._yield_with_calc_alignability(chrom)

    def _yield(self, chrom):
        for wig in self.fetch(self.MAPPABILITY_THRESHOLD, chrom):
            try:
                yield wig
            except GeneratorExit:
                break

    def _yield_with_calc_alignability(self, chrom):
        stop_yield = False

        self._init_buff(chrom)
        for wig in self.fetch(self.MAPPABILITY_THRESHOLD, chrom):
            self._progress.update(wig[2])
            self._feed_track(wig[1], wig[2])

            if not stop_yield:
                try:
                    yield wig
                except GeneratorExit:
                    logger.info("Continue calc alignability for {}...".format(chrom))
                    stop_yield = True

        self._flush()
