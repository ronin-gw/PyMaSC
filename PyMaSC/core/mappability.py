import logging

import numpy as np

from PyMaSC.reader.bigwig.bigwig import BigWigFile
from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


class ContinueCalculation(Exception):
    pass


class MappableLengthCalculator(BigWigFile):
    MAPPABILITY_THRESHOLD = 1

    def __init__(self, path, max_shift=0, chrom_size=None):
        super(MappableLengthCalculator, self).__init__(path, chrom_size)
        self.max_shift = max_shift
        self.chrom2is_called = {c: False for c in self.chromsizes}
        self.chrom2mappable_len = {}
        self.mappable_len = [0] * (max_shift + 1)
        self.is_called = False

        self._sum_bin_size = self.max_shift + 1
        self._chr = None

        self._sumbins = np.zeros(self.max_shift + 1, dtype=np.int64)
        self._buff_tail_pos = -1
        self._buff = np.zeros(self.max_shift, dtype=np.int64)

        self._progress = ProgressBar()

    def enable_progress_bar(self):
        self._progress.enable_bar()

    def disable_progress_bar(self):
        self._progress.disable_bar()

    def _init_buff(self, chrom):
        self._sumbins.fill(0)
        self._buff_tail_pos = -1
        self._buff = np.zeros(self.max_shift, dtype=np.int64)

        logger.info("Calc {} mappable length...".format(chrom))
        self._chr = chrom
        self._progress.set(self.chromsizes[chrom])

    def _flush(self):
        if self._chr:
            self.chrom2mappable_len[self._chr] = tuple(self._sumbins)
            for i in xrange(self.max_shift + 1):
                self.mappable_len[i] += self._sumbins[i]
            self.chrom2is_called[self._chr] = True

            if all(self.chrom2is_called.values()):
                self.is_called = True

            self._progress.clean()

    def calc_mappability(self, chrom=None):
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
        gap_len = begin - self._buff_tail_pos - 1
        track_len = end - begin

        for i in xrange(min(self.max_shift + 1, track_len)):
            self._sumbins[i] += track_len - i

        if gap_len < self.max_shift:
            overlap_len = self.max_shift - gap_len

            self._sumbins[gap_len+1:] += np.correlate(
                np.ones(min(overlap_len, track_len), dtype=np.int64),
                self._buff[-(overlap_len):], "full"
            )[:overlap_len]

            self._buff = np.concatenate((
                self._buff, np.zeros(gap_len, dtype=np.int64),
                np.ones(track_len, dtype=np.int64)
            ))[-self.max_shift:]
        elif track_len < self.max_shift:
            self._buff = np.concatenate((
                np.zeros(self.max_shift - track_len, dtype=np.int64),
                np.ones(track_len, dtype=np.int64)
            ))
        else:
            self._buff = np.ones(self.max_shift, dtype=np.int64)

        self._buff_tail_pos = end - 1

    def get_mappable_len(self, chrom=None, shift_from=None, shift_to=None, force=False):
        if chrom is not None:
            if chrom not in self.chrom2is_called:
                return None
            if self.chrom2is_called[chrom]:
                return self.chrom2mappable_len[chrom][shift_from:shift_to]
            elif force:
                self.calc_mappability(chrom)
                return self.chrom2mappable_len[chrom][shift_from:shift_to]
            else:
                raise KeyError("Mappable length for '{}' is not calculated yet.".format(chrom))

        if self.is_called:
            return self.mappable_len[shift_from:shift_to]

        self.calc_mappability()
        return self.mappable_len[shift_from:shift_to]


class BWFeederWithMappableRegionSum(MappableLengthCalculator):
    def feed(self, chrom):
        if self.chrom2is_called[chrom]:
            return self._yield(chrom)
        else:
            return self._yield_with_calc_mappability(chrom)

    def _yield(self, chrom):
        for wig in self.fetch(self.MAPPABILITY_THRESHOLD, chrom):
            try:
                yield wig
            except ContinueCalculation:
                break

    def _yield_with_calc_mappability(self, chrom):
        stop_yield = False

        self._init_buff(chrom)
        for wig in self.fetch(self.MAPPABILITY_THRESHOLD, chrom):
            self._progress.update(wig[2])
            self._feed_track(wig[1], wig[2])

            if not stop_yield:
                try:
                    yield wig
                except ContinueCalculation:
                    logger.info("Continue calc mappability for {}...".format(chrom))
                    stop_yield = True
                    self.enable_progress_bar()

        self._flush()
