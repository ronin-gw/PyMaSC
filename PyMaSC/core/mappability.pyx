# cython: language_level=3
"""Mappability calculation and BigWig data processing.

This module provides Cython-optimized implementations for calculating
mappability statistics from BigWig files. Mappability data is essential
for MSCC calculations to correct for genomic regions with variable
mappability that can confound cross-correlation analysis.

Key components:
- MappableLengthCalculator: Core mappability calculation engine
- BWFeederWithMappableRegionSum: Specialized feeder for on-demand calculation
- ContinueCalculation: Exception for flow control in generators

The module efficiently processes BigWig mappability data and calculates
shift-dependent mappability statistics required for MSCC analysis.
"""
import logging

cimport numpy as np
from PyMaSC.reader.bigwig cimport BigWigReader, BigWigFileIterator, interval
from cython cimport boundscheck, wraparound


import numpy as np
from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


class ContinueCalculation(Exception):
    """Exception for flow control in mappability calculation generators.

    Used to signal that a generator should stop yielding data points
    but continue internal calculation processes. This allows for
    flexible control flow in mappability data processing.
    """
    pass


cdef class MappableLengthCalculator(BigWigReader):
    """Core mappability calculation engine for BigWig data.

    Calculates mappability statistics from BigWig files for use in
    MSCC cross-correlation analysis. Processes mappability data
    efficiently using Cython optimization and maintains statistics
    for different shift distances.

    The calculator processes BigWig intervals and computes cumulative
    mappability statistics required for correcting cross-correlation
    analysis in regions with variable mappability.

    Attributes:
        MAPPABILITY_THRESHOLD: Minimum mappability value to consider
        max_shift: Maximum shift distance for calculations
        chrom2is_called: Per-chromosome calculation status
        chrom2mappable_len: Per-chromosome mappability arrays
        mappable_len: Genome-wide mappability statistics
        is_called: Whether all chromosomes have been processed
    """
    cdef:
        public double MAPPABILITY_THRESHOLD
        readonly unsigned int max_shift
        public dict chrom2is_called
        public dict chrom2mappable_len
        public list mappable_len
        public bint is_called

        unsigned int _sum_bin_size
        str _chr

        np.ndarray _sumbins
        unsigned int _buff_tail_pos
        np.ndarray _buff

        public object _progress

        object logger_lock

    def __init__(self, str path, unsigned int max_shift=0, logger_lock=None):
        """Initialize mappability calculator.

        Args:
            path: Path to BigWig mappability file
            max_shift: Maximum shift distance for calculations
            logger_lock: Optional lock for thread-safe logging
        """
        self.MAPPABILITY_THRESHOLD = 1.

        super(MappableLengthCalculator, self).__init__(path)
        self.max_shift = max_shift
        self.chrom2is_called = {c: False for c in self.chromsizes}
        self.chrom2mappable_len = {}
        self.mappable_len = [0] * (max_shift + 1)
        self.is_called = False

        self._sum_bin_size = self.max_shift + 1
        self._chr = ''

        self._sumbins = np.zeros(self.max_shift + 1, dtype=np.int64)
        self._buff_tail_pos = -1
        self._buff = np.zeros(self.max_shift, dtype=np.int64)

        self._progress = ProgressBar()
        self.logger_lock = logger_lock

    def enable_progress_bar(self):
        self._progress.enable_bar()

    def disable_progress_bar(self):
        self._progress.disable_bar()

    cdef _init_buff(self, chrom):
        self._sumbins.fill(0)
        self._buff_tail_pos = -1
        self._buff = np.zeros(self.max_shift, dtype=np.int64)

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
        """Calculate mappability statistics for specified chromosome(s).

        Processes BigWig mappability data and computes shift-dependent
        mappability statistics. Can process a single chromosome or
        all unprocessed chromosomes.

        Args:
            chrom: Chromosome name to process (None for all unprocessed)

        Note:
            Updates chrom2mappable_len and mappable_len attributes
            with calculated statistics
        """
        if not chrom:
            chroms = [c for c, b in self.chrom2is_called.items() if b is False]
        elif self.chrom2is_called[chrom]:
            return None
        else:
            chroms = [chrom]

        cdef char* _chrom
        cdef BigWigFileIterator gen
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
    cdef inline _feed_track(self, unsigned int begin, unsigned int end):
        cdef unsigned int i
        cdef unsigned int track_len = end - begin
        cdef unsigned int gap_len = begin - self._buff_tail_pos - 1
        cdef overlap_len, remain_buff_len

        for i in range(min(self.max_shift + 1, track_len)):
            self._sumbins[i] += track_len - i

        if gap_len < self.max_shift:
            overlap_len = self.max_shift - gap_len

            self._sumbins[gap_len+1:] += np.correlate(
                np.ones(min(overlap_len, track_len), dtype=np.int64),
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
        """Retrieve mappability statistics for specified parameters.

        Returns mappability statistics for a chromosome or genome-wide,
        optionally filtered by shift distance range. Can trigger
        calculation if not already computed.

        Args:
            chrom: Chromosome name (None for genome-wide)
            shift_from: Start of shift range (None for beginning)
            shift_to: End of shift range (None for end)
            force: Whether to force calculation if not computed

        Returns:
            List of mappability values for specified range

        Raises:
            KeyError: If chromosome data not calculated and force=False
        """
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
    """Specialized BigWig feeder with on-demand mappability calculation.

    Extends MappableLengthCalculator to provide streaming access to
    BigWig data while calculating mappability statistics on-demand.
    This class is used in MSCC calculations where mappability data
    is needed incrementally during cross-correlation computation.

    The feeder can operate in two modes:
    1. Using pre-calculated mappability statistics (if available)
    2. Calculating mappability on-demand while streaming data

    Attributes:
        Inherits all attributes from MappableLengthCalculator
    """

    def __init__(self, path, max_shift, chrom2mappable_len=None):
        """Initialize BigWig feeder with mappability calculation.

        Args:
            path: Path to BigWig mappability file
            max_shift: Maximum shift distance for calculations
            chrom2mappable_len: Pre-calculated mappability statistics (optional)
        """
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
        """Feed BigWig data for specified chromosome.

        Returns a generator that yields BigWig intervals for the
        specified chromosome. If mappability statistics haven't
        been calculated for this chromosome, performs on-demand
        calculation while yielding data.

        Args:
            chrom: Chromosome name to process

        Returns:
            Generator yielding (start, end, value) tuples for BigWig intervals

        Note:
            Automatically chooses between pre-calculated and on-demand
            calculation modes based on available data
        """
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
