# cython: profile=True
import logging

cimport numpy as np
from PyMaSC.reader.bigwig.bigwig cimport BigWigFile
from PyMaSC.reader.bigwig.interval cimport Interval
from PyMaSC.reader.bigwig.kentlib.types cimport bits32
from cython cimport boundscheck, wraparound

import numpy as np
from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


class ContinueCalculation(Exception):
    pass


cdef class MappableLengthCalculator(BigWigFile):
    cdef:
        public double MAPPABILITY_THRESHOLD
        readonly bits32 max_shift
        public dict chrom2is_called
        public dict chrom2mappable_len
        public list mappable_len
        public bint is_called

        unsigned int _sum_bin_size
        char* _chr

        # np.ndarray[np.long] _sumbins
        bits32 _buff_tail_pos
        # np.ndarray[np.long] _buff

    def __init__(self, char* path, unsigned int max_shift=0, chrom_size=None):
        self.MAPPABILITY_THRESHOLD = 1.

        super(MappableLengthCalculator, self).__init__(path, chrom_size)
        self.max_shift = max_shift
        self.chrom2is_called = {c: False for c in self.chromsizes}
        self.chrom2mappable_len = {}
        self.mappable_len = [0] * (max_shift + 1)
        self.is_called = False

        self._sum_bin_size = self.max_shift + 1
        self._chr = NULL

        self._sumbins = np.zeros(self.max_shift + 1, dtype=np.long)
        self._buff_tail_pos = -1
        self._buff = np.zeros(self.max_shift, dtype=np.long)

        self._progress = ProgressBar()

    def enable_progress_bar(self):
        self._progress.enable_bar()

    def disable_progress_bar(self):
        self._progress.disable_bar()

    cdef _init_buff(self, chrom):
        self._sumbins.fill(0)
        self._buff_tail_pos = -1
        self._buff = np.zeros(self.max_shift, dtype=np.long)

        logger.info("Calc {} mappable length...".format(chrom))
        self._chr = chrom
        self._progress.set(self.chromsizes[chrom])

    cdef _flush(self):
        if self._chr:
            self.chrom2mappable_len[self._chr] = tuple(self._sumbins)
            for i in xrange(self.max_shift + 1):
                self.mappable_len[i] += self._sumbins[i]
            print self.mappable_len
            self.chrom2is_called[self._chr] = True

            if all(self.chrom2is_called.values()):
                self.is_called = True

            self._progress.clean()

    def calc_mappability(self, chrom=None):
        if not chrom:
            chroms = [c for c, b in self.chrom2is_called.items() if b is False]
        elif self.chrom2is_called[chrom]:
            return None
        else:
            chroms = [chrom]

        cdef char* _chrom
        cdef bits32 begin, end
        cdef double _val

        for c in chroms:
            self._init_buff(c)
            for _chrom, begin, end, _val in self.fetch(self.MAPPABILITY_THRESHOLD, c):
                self._progress.update(end)
                self._feed_track(begin, end)

            self._flush()

    @wraparound(False)
    @boundscheck(False)
    cdef _feed_track(self, bits32 begin, bits32 end):
        cdef bits32 i
        cdef bits32 track_len = end - begin
        cdef bits32 gap_len = begin - self._buff_tail_pos - 1
        cdef overlap_len, remain_buff_len

        for i in range(min(self.max_shift + 1, track_len)):
            self._sumbins[i] += track_len - i

        if gap_len < self.max_shift:
            overlap_len = self.max_shift - gap_len

            self._sumbins[gap_len+1:] += np.correlate(
                np.ones(min(overlap_len, track_len), dtype=np.long),
                self._buff[self.max_shift - overlap_len:], "full"
            )[:overlap_len]

        if track_len + gap_len < self.max_shift:
            # `overlap_len` must be defined @line 115 under this condition.
            remain_buff_len = overlap_len - track_len
            self._buff[:remain_buff_len] = self._buff[self.max_shift - remain_buff_len:]
            self._buff[remain_buff_len:remain_buff_len + gap_len] = 0
            self._buff[remain_buff_len + gap_len:] = 1
        elif track_len < self.max_shift:
            self._buff[:self.max_shift - track_len] = 0
            self._buff[self.max_shift - track_len:] = 1
        else:
            self._buff.fill(1)

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
            return self._yield_withcalc_mappability(chrom)

    def _yield(self, chrom):
        for wig in self.fetch(self.MAPPABILITY_THRESHOLD, chrom):
            try:
                yield wig
            except ContinueCalculation:
                break

    def _yield_withcalc_mappability(self, chrom):
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
