# cython: language_level=3
"""Mappability-Sensitive Cross-Correlation (MSCC) calculation implementation.

This module implements the MSCC algorithm, which extends naive cross-correlation
by incorporating mappability information from BigWig files. The algorithm
corrects for non-uniform mappability across the genome, addressing the
'phantom peak' problem in ChIP-seq analysis.

Key components:
- MSCCCalculator: Cython-optimized MSCC calculator with mappability correction
- Integration with BigWig mappability tracks
- Sliding window algorithm with mappability weighting
- Memory-efficient processing of large genomic datasets

The MSCC algorithm provides more accurate fragment length estimation by
accounting for regions of variable mappability that can confound traditional
cross-correlation analysis.
"""
import logging

cimport numpy as np
from libc.stdint cimport int_fast64_t as int64
from libc.string cimport strcmp, strcpy
from cython cimport boundscheck, wraparound

from .utils cimport int64_min, int64_max

import numpy as np
from ..exceptions import ReadUnsortedError
from ..mappability import ContinueCalculation
from PyMaSC.result import MSCCResult, MSCCGenomeWideResult, EmptyMSCCResult

logger = logging.getLogger(__name__)


cdef class MSCCCalculator(object):
    """Cython implementation of Mappability-Sensitive Cross-Correlation.

    This class implements the MSCC algorithm using Cython for performance.
    It extends naive cross-correlation by incorporating mappability weights
    from BigWig files to correct for genomic regions with variable mappability.

    The algorithm processes reads using a sliding window approach while
    consulting mappability tracks to weight cross-correlation contributions
    appropriately. This helps eliminate phantom peaks caused by mappability
    biases in repetitive or low-mappability regions.

    Attributes:
        max_shift: Maximum shift distance for correlation calculation
        read_len: Read length for mappability window calculation
        ref2genomelen: Dictionary mapping chromosome names to lengths
        genomelen: Total genome length
        logger_lock: Thread lock for logging (optional)

    Note:
        Result data (forward/reverse sums, mappable lengths, MSCCResult objects) are stored internally
        and accessed via get_genome_wide_result() method.
    """
    cdef readonly:
        int64 max_shift
        int64 read_len
        object logger_lock
    cdef:
        # Result storage (internal)
        dict ref2genomelen
        int64 genomelen
        int64 forward_read_len_sum, reverse_read_len_sum
        dict ref2mscc_result
        # Calculation state
        str _chr
        int64 _max_shift_from_f, _forward_buff_size, _forward_read_len_sum, _reverse_read_len_sum
        np.ndarray _forward_sum, _reverse_sum, _ccbins
        list _forward_buff, _reverse_buff, _solved_chr
        int64 _fb_tail_pos, _last_pos, _last_forward_pos
        # Mappability handling
        object _bwfeeder, _feeder
        list references
        np.ndarray _bwbuff
        bint _bwiter_stopped
        int64 _bwbuff_head, _extra_size, _double_shift_len

    def __init__(self, max_shift, read_len, references, lengths, bwfeeder, logger_lock=None):
        """Initialize MSCC calculator with mappability support.

        Args:
            max_shift: Maximum shift distance for correlation calculation
            read_len: Read length for mappability window calculation
            references: List of chromosome names
            lengths: List of chromosome lengths
            bwfeeder: BigWig reader for mappability data
            logger_lock: Optional lock for thread-safe logging
        """
        #    _last_pos      pos + read_len - 1
        #        |-----------------| read_len
        #        /=================>
        #    <---------------------| max_shift
        #    |<--| _max_shift_from_f

        self.max_shift = max_shift
        self.ref2genomelen = dict(zip(references, lengths))
        self.genomelen = sum(lengths)
        self.read_len = read_len

        self._forward_buff_size = self.max_shift + 1
        # _max_shift_from_f = relative pos from forward read pos to max_shift
        self._max_shift_from_f = -(read_len - self._forward_buff_size)

        # stats
        self.ref2mscc_result = {}
        self.references = references
        self.forward_read_len_sum = self.reverse_read_len_sum = 0
        # internal buff
        self._chr = ''
        self._forward_sum = np.zeros(self._forward_buff_size, dtype=np.int64)
        self._reverse_sum = np.zeros(self._forward_buff_size, dtype=np.int64)
        self._ccbins = np.zeros(self._forward_buff_size, dtype=np.int64)
        # None for read absence
        self._forward_buff = [None for _ in range(self._forward_buff_size)]
        self._reverse_buff = [None for _ in range(self._forward_buff_size)]
        #
        self._solved_chr = []

        # mappability
        self._bwfeeder = bwfeeder
        self._bwfeeder.disable_progress_bar()
        self._bwiter_stopped = False
        # _bwbuff[i] = mappability of pos i
        self._bwbuff = np.array([], dtype=np.int64)
        # _bwbuff_head = 1-based position for first element in _bwbuff
        self._bwbuff_head = 1
        # _extra_size = extra buff size from _last_pos save to _bwbuff
        self._extra_size = int64_max(0, self._max_shift_from_f)
        self._double_shift_len = self.max_shift * 2 + 1

        self._init_pos_buff()
        self._bwfeeder.disable_progress_bar()

        self.logger_lock = logger_lock

    def _init_pos_buff(self):
        """Initialize position buffers for a new chromosome.

        Resets all internal buffers, counters, and mappability state
        for processing a new chromosome in the MSCC calculation.
        """
        self._forward_sum.fill(0)
        self._reverse_sum.fill(0)
        self._ccbins.fill(0)
        self._forward_read_len_sum = self._reverse_read_len_sum = 0

        for i in range(len(self._forward_buff)):
            self._forward_buff[i] = None
        self._fb_tail_pos = self._forward_buff_size

        for i in range(len(self._reverse_buff)):
            self._reverse_buff[i] = None

        self._last_pos = 0
        self._last_forward_pos = 0

    def _init_bw(self):
        """Initialize BigWig mappability reader for current chromosome.

        Sets up the BigWig iterator and buffers for reading mappability
        values for the current chromosome. Handles cases where mappability
        data is not available for certain chromosomes.

        Raises:
            ContinueCalculation: If mappability data is unavailable
        """
        self._bwbuff = np.array([], dtype=np.int64)
        self._bwbuff_head = 1
        self._bwiter_stopped = False
        try:
            self._feeder = self._bwfeeder.feed(self._chr)
        except KeyError as e:
            if self.logger_lock:
                self.logger_lock.acquire()
            logger.info("Mappability for '{}' not found. "
                        "Skip calc mappability sensitive CC.".format(e.args[0]))
            if self.logger_lock:
                self.logger_lock.release()
            self._bwiter_stopped = True

    def flush(self, chrom: Optional[str] = None):
        """Finalize MSCC calculation for current chromosome.

        Completes the mappability-weighted cross-correlation calculation
        for the current chromosome and stores results as an MSCCResult object.
        Handles final mappability corrections before storing.
        """
        if self._chr == '':
            return

        if chrom is not None and chrom != self._chr:
            self._chr = chrom

        if self._reverse_buff:
            self._shift_with_update(self._forward_buff_size)

        self.forward_read_len_sum += self._forward_read_len_sum
        self.reverse_read_len_sum += self._reverse_read_len_sum

        # update bwfeeder
        if not self._bwiter_stopped and self._feeder:
            try:
                self._feeder.throw(ContinueCalculation)
            except StopIteration:
                pass
            self._bwiter_stopped = True

        # Get mappable_len
        mappable_len = self._bwfeeder.get_mappable_len(self._chr)
        if mappable_len is None:
            self._init_pos_buff()
            return

        # Create MSCCResult for this chromosome
        result = self.ref2mscc_result[self._chr] = MSCCResult(
            max_shift=self.max_shift,
            read_len=self.read_len,
            genomelen=self.ref2genomelen[self._chr],
            forward_sum=self._forward_sum.copy(),
            reverse_sum=self._reverse_sum.copy(),
            forward_read_len_sum=self._forward_read_len_sum,
            reverse_read_len_sum=self._reverse_read_len_sum,
            ccbins=self._ccbins.copy(),
            mappable_len=mappable_len
        )
        result.calc_cc()

        self._init_pos_buff()

    cdef inline _check_pos(self, str chrom, int64 pos):
        """Validate read position and handle chromosome transitions.

        Ensures reads are sorted by position and handles switching between
        chromosomes. Initializes mappability data for new chromosomes.

        Args:
            chrom: Chromosome name
            pos: Read position

        Raises:
            ReadUnsortedError: If reads are not sorted by position
        """
        if chrom != self._chr:
            if self._chr != '':
                if chrom in self._solved_chr:
                    raise ReadUnsortedError
                self._solved_chr.append(self._chr)
                self.flush()
            else:
                self._init_pos_buff()
            self._chr = chrom
            self._init_bw()
            if not self._bwiter_stopped:
                if self.logger_lock:
                    self.logger_lock.acquire()
                logger.info("Calculate cross-correlation for {}...".format(chrom))
                if self.logger_lock:
                    self.logger_lock.release()

        if pos < self._last_pos:
            raise ReadUnsortedError

        self._last_pos = pos

    @boundscheck(False)
    def feed_forward_read(self, str chrom, int64 pos, int64 readlen):
        """Process a forward strand read with mappability weighting.

        Incorporates a forward strand read into the MSCC calculation,
        applying mappability weights from BigWig data to correct for
        regional mappability biases.

        Args:
            chrom: Chromosome name
            pos: Read start position (1-based)
            readlen: Read length in base pairs

        Note:
            Duplicate reads at the same position are ignored
        """
        cdef int64 offset, read_mappability
        cdef np.ndarray[np.int64_t] mappability

        self._check_pos(chrom, pos)

        if self._last_forward_pos == pos:
            return None
        self._last_forward_pos = pos

        offset = pos - self._fb_tail_pos

        read_mappability, mappability = self._get_forward_mappabiility()
        if not read_mappability:
            return None

        self._forward_sum += mappability
        self.forward_read_len_sum += readlen

        if offset < 1:
            self._forward_buff[offset - 1] = mappability
        else:
            self._shift_with_update(offset)
            self._forward_buff[self.max_shift] = mappability

    @wraparound(False)
    @boundscheck(False)
    def feed_reverse_read(self, str chrom, int64 pos, int64 readlen):
        """Process a reverse strand read with mappability weighting.

        Incorporates a reverse strand read into the MSCC calculation,
        applying appropriate mappability weights and position adjustments
        for reverse strand processing.

        Args:
            chrom: Chromosome name
            pos: Read start position (1-based)
            readlen: Read length in base pairs

        Note:
            Reverse read positions are adjusted for 3' end calculation
        """
        cdef int64 offset, revbuff_pos
        cdef np.ndarray[np.int64_t] mappability

        self._check_pos(chrom, pos)

        offset = pos - self._fb_tail_pos
        revbuff_pos = self.max_shift + readlen - 1  # revbuff_tail_pos

        if offset > 0:
            self._shift_with_update(offset)
        else:
            revbuff_pos += offset

        if revbuff_pos < len(self._reverse_buff) and self._reverse_buff[revbuff_pos] == 1:
            return None

        mappability = self._get_reverse_read_mappability(readlen)
        if not np.any(mappability):
            return None

        self._reverse_sum += mappability
        self.reverse_read_len_sum += readlen

        if revbuff_pos < len(self._reverse_buff):
            self._reverse_buff[revbuff_pos] = 1
        else:
            self._reverse_buff += [None] * (revbuff_pos - len(self._reverse_buff)) + [1]

    @wraparound(False)
    @boundscheck(False)
    cdef inline tuple _get_forward_mappabiility(self):
        if self._bwiter_stopped and not self._bwbuff.size:
            return 0, np.zeros(self._forward_buff_size, dtype=np.int64)

        cdef int64 from_, to
        cdef np.ndarray[np.int64_t] mappability

        to = self._last_pos + self.read_len - 1
        from_ = int64_min(to - self._forward_buff_size, self._last_pos - 1)

        if from_ < 0:
            mappability = np.concatenate((
                np.zeros(-from_, dtype=np.int64), self._get_bwbuff(0, to)
            ))
        else:
            mappability = self._get_bwbuff(from_, to)

        return mappability[self._extra_size], mappability[self.max_shift::-1]

    @wraparound(False)
    @boundscheck(False)
    cdef inline np.ndarray[np.int64_t] _get_reverse_read_mappability(self, int64 readlen):
        if self._bwiter_stopped and not self._bwbuff.size:
            return np.zeros(self._forward_buff_size, dtype=np.int64)

        cdef int64 forward_from, forward_to, reverse_from, reverse_to
        cdef np.ndarray[np.int64_t] forward_mappability, reverse_mappability

        forward_to = self._last_pos + readlen - 1
        forward_from = forward_to - self._forward_buff_size

        if forward_from < 0:
            forward_mappability = np.concatenate((
                self._get_bwbuff(0, forward_to)[::-1], np.zeros(-forward_from, dtype=np.int64)
            ))
        else:
            forward_mappability = self._get_bwbuff(forward_from, forward_to)[::-1]

        reverse_to = forward_to + self.read_len - 1
        reverse_from = reverse_to - self._double_shift_len

        if reverse_from < 0:
            reverse_mappability = np.concatenate((
                self._get_bwbuff(0, reverse_to)[::-2],
                np.zeros(self._forward_buff_size - (reverse_to + 1) / 2, dtype=np.int64)
            ))
        else:
            reverse_mappability = self._get_bwbuff(reverse_from, reverse_to)[::-2]

        return forward_mappability * reverse_mappability

    @wraparound(False)
    @boundscheck(False)
    cdef inline np.ndarray[np.int64_t] _get_bwbuff(self, int64 from_, int64 to):
        cdef int64 bwbuff_head = int64_max(0, self._last_pos - self._double_shift_len)
        cdef int64 bwbuff_tail, gap

        cdef int64 start, end
        cdef float _val

        cdef np.ndarray[np.int64_t] m

        self._bwbuff = self._bwbuff[bwbuff_head - self._bwbuff_head:]
        self._bwbuff_head = bwbuff_head

        from_ -= bwbuff_head
        to -= bwbuff_head

        while self._bwbuff.size < to:
            try:
                start, end, _val = next(self._feeder)
            except StopIteration:
                self._bwiter_stopped = True
                m = self._bwbuff[from_:to]
                return np.concatenate((
                    m, np.zeros(int64_max(0, to - from_ - m.size), dtype=np.int64)
                ))

            bwbuff_tail = bwbuff_head + self._bwbuff.size
            if bwbuff_tail <= start:
                gap = start - bwbuff_tail
                self._bwbuff = np.concatenate((
                    self._bwbuff, np.zeros(gap, dtype=np.int64),
                    np.ones(end - start, dtype=np.int64)
                ))
            elif bwbuff_head < end:
                self._bwbuff[int64_max(0, start - bwbuff_head):(end - bwbuff_head)] = 1
                self._bwbuff = np.concatenate(
                    (self._bwbuff, np.ones(int64_max(0, end - bwbuff_tail), dtype=np.int64))
                )

        return self._bwbuff[from_:to]

    @wraparound(False)
    @boundscheck(False)
    cdef inline _shift_with_update(self, int64 offset):
        cdef int64 i, j
        cdef np.ndarray[np.int64_t] forward_state

        for i in range(int64_min(int64_min(offset, self._forward_buff_size), len(self._reverse_buff))):
            forward_state = self._forward_buff[i]
            if forward_state is None:
                continue

            for j in range(int64_min(self._forward_buff_size, len(self._reverse_buff) - i)):
                if self._reverse_buff[i + j] is None:
                    continue

                self._ccbins[j] += forward_state[j]

        self._forward_buff = self._forward_buff[offset:] + [None] * int64_min(offset, self._forward_buff_size)
        self._reverse_buff = self._reverse_buff[offset:]
        self._fb_tail_pos += offset

    def finishup_calculation(self):
        """Complete MSCC calculation for all chromosomes.

        Finalizes the mappability-sensitive cross-correlation calculation
        by flushing any remaining data and ensuring all results are
        properly stored with mappability corrections applied.
        """
        self.flush()

        if not self._bwiter_stopped and self._feeder is not None:
            try:
                self._feeder.throw(ContinueCalculation)
            except StopIteration:
                pass
        self._bwfeeder.calc_mappability()

        # Initialize any unprocessed chromosomes with zero-filled results
        cdef int64 zero_bins_size = self._forward_buff_size
        cdef np.ndarray zero_bins = np.zeros(zero_bins_size, dtype=np.int64)
        for ref in self.references:
            if ref in self._bwfeeder.chrom2mappable_len and ref not in self.ref2mscc_result:
                result = self.ref2mscc_result[ref] = EmptyMSCCResult.create_empty(
                    self.ref2genomelen[ref], self.max_shift, self.read_len
                )
                # Update mappable_len with actual mappability data if available
                result.mappable_len = self._bwfeeder.chrom2mappable_len[ref]

    def get_result(self, chrom: str) -> MSCCResult:
        return self.ref2mscc_result[chrom]

    def get_whole_result(self):
        """Get the genome-wide MSCC calculation results.

        Returns:
            MSCCGenomeWideResult: Complete results for all chromosomes including
                per-chromosome MSCCResult objects and genome-wide statistics.
        """
        return MSCCGenomeWideResult(
            genomelen=self.genomelen,
            chroms=self.ref2mscc_result.copy(),
            forward_read_len_sum=self.forward_read_len_sum,
            reverse_read_len_sum=self.reverse_read_len_sum
        )
