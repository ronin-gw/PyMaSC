import logging

cimport numpy as np
from bx.bbi.types cimport bits32

from PyMaSC.reader.bigwig cimport BigWigReader
from PyMaSC.reader.bx.bbi_file cimport interval
from PyMaSC.reader.bx.bigwig_file cimport BigWigFile
from cython cimport boundscheck, wraparound

from .utils cimport bits32_min

import numpy as np
from PyMaSC.utils.progress import ProgressBar
from PyMaSC.utils.compatible import tostr

logger = logging.getLogger(__name__)


class ContinueCalculation(Exception):
    pass


cdef class MappableLengthCalculator(BigWigReader):
    cdef:
        public double MAPPABILITY_THRESHOLD
        readonly bits32 max_shift
        public dict chrom2is_called
        public dict chrom2mappable_len
        public list mappable_len
        public bint is_called

        unsigned int _sum_bin_size
        str _chr

        np.ndarray _sumbins
        bits32 _buff_tail_pos
        np.ndarray _buff

        public object _progress

        object logger_lock

    def __init__(self, str path, unsigned int max_shift=0, logger_lock=None):
        self.MAPPABILITY_THRESHOLD = 1.

        super(MappableLengthCalculator, self).__init__(path)
        self.max_shift = max_shift
        self.chrom2is_called = {c: False for c in self.chromsizes}
        self.chrom2mappable_len = {}
        self.mappable_len = [0] * (max_shift + 1)
        self.is_called = False

        self._sum_bin_size = self.max_shift + 1
        self._chr = ''

        self._sumbins = np.zeros(self.max_shift + 1, dtype=np.long)
        self._buff_tail_pos = -1
        self._buff = np.zeros(self.max_shift, dtype=np.long)

        self._progress = ProgressBar()
        self.logger_lock = logger_lock

    def enable_progress_bar(self):
        self._progress.enable_bar()

    def disable_progress_bar(self):
        self._progress.disable_bar()

    cdef _init_buff(self, chrom):
        self._sumbins.fill(0)
        self._buff_tail_pos = -1
        self._buff = np.zeros(self.max_shift, dtype=np.long)

        if self.logger_lock:
            self.logger_lock.acquire()
        logger.info("Calc {} mappable length...".format(chrom))
        if self.logger_lock:
            self.logger_lock.release()
        self._chr = chrom
        self._progress.set(chrom, self.chromsizes[chrom])

    cpdef flush(self):
        if self._chr:
            self.chrom2mappable_len[self._chr] = tuple(self._sumbins)
            for i in xrange(self.max_shift + 1):
                self.mappable_len[i] += self._sumbins[i]
            self.chrom2is_called[self._chr] = True

            if all(self.chrom2is_called.values()):
                self.is_called = True

            self._progress.clean()

    def calc_mappability(self, chrom=None):
        if not chrom:
            chroms = [tostr(c) for c, b in self.chrom2is_called.items() if b is False]
        elif self.chrom2is_called[chrom]:
            return None
        else:
            chroms = [chrom]

        cdef char* _chrom
        cdef BigWigFile gen
        cdef interval i

        for c in chroms:
            self._init_buff(c)
            gen = self.fetch(self.MAPPABILITY_THRESHOLD, c)
            while True:
                i = gen.next()
                if i.end == 0:
                    break
                self._progress.update(i.end)
                self._feed_track(i.begin, i.end)

            self.flush()

    @wraparound(False)
    @boundscheck(False)
    cdef inline _feed_track(self, bits32 begin, bits32 end):
        cdef bits32 i
        cdef bits32 track_len = end - begin
        cdef bits32 gap_len = begin - self._buff_tail_pos - 1
        cdef overlap_len, remain_buff_len

        for i in range(bits32_min(self.max_shift + 1, track_len)):
            self._sumbins[i] += track_len - i

        if gap_len < self.max_shift:
            overlap_len = self.max_shift - gap_len

            self._sumbins[gap_len+1:] += np.correlate(
                np.ones(bits32_min(overlap_len, track_len), dtype=np.long),
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


cdef class BWFeederWithMappableRegionSum(MappableLengthCalculator):
    def __init__(self, path, max_shift, chrom2mappable_len=None):
        super(BWFeederWithMappableRegionSum, self).__init__(path, max_shift)

        if chrom2mappable_len:
            for chrom, l in chrom2mappable_len.items():
                if chrom in self.chrom2is_called:
                    self.chrom2is_called[chrom] = True
                    self.chrom2mappable_len[chrom] = l
                    for i in range(self.max_shift + 1):
                        self.mappable_len[i] += l[i]

            if all(self.chrom2is_called.values()):
                self.is_called = True

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
            self._progress.update(wig[1])
            self._feed_track(wig[0], wig[1])

            if not stop_yield:
                try:
                    yield wig
                except ContinueCalculation:
                    if self.logger_lock:
                        self.logger_lock.acquire()
                    logger.info("Continue calc mappability for {}...".format(chrom))
                    if self.logger_lock:
                        self.logger_lock.release()
                    stop_yield = True
                    self.enable_progress_bar()

        self.flush()
