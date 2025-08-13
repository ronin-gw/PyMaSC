# cython: language_level=3
"""Naive cross-correlation calculation implementation.

This module implements the core naive cross-correlation (NCC) algorithm for
ChIP-seq fragment length estimation. The algorithm calculates cross-correlation
between forward and reverse strand reads at different shift distances to
identify the characteristic fragment length peak.

Key components:
- NaiveCCCalculator: Cython-optimized cross-correlation calculator
- ReadUnsortedError: Exception for unsorted input data
- Sliding window algorithm for memory-efficient processing

The implementation uses Cython for performance optimization and handles
large genomic datasets efficiently through buffered processing.
"""
import logging

cimport numpy as np
from libc.stdint cimport int_fast64_t as int64
from libc.string cimport strcmp, strcpy
from cython cimport boundscheck, wraparound

from .utils cimport int64_min

import numpy as np

from ..exceptions import ReadUnsortedError
from PyMaSC.result import NCCResult, NCCGenomeWideResult, EmptyNCCResult

logger = logging.getLogger(__name__)


cdef class NaiveCCCalculator(object):
    """Cython implementation of naive cross-correlation calculation.

    This class implements the core NCC algorithm using Cython for performance.
    It processes reads from forward and reverse strands using a sliding window
    approach to calculate cross-correlation at different shift distances.

    The algorithm maintains buffers for forward and reverse reads and updates
    cross-correlation bins as reads are processed, enabling memory-efficient
    computation for large genomic datasets.

    Attributes:
        max_shift: Maximum shift distance for correlation calculation
        ref2genomelen: Dictionary mapping chromosome names to lengths
        genomelen: Total genome length
        logger_lock: Thread lock for logging (optional)

    Note:
        Result data (forward/reverse sums, NCCResult objects) are stored internally
        and accessed via get_genome_wide_result() method.
    """
    cdef readonly:
        int64 max_shift
        object logger_lock
    cdef:
        # Result storage (internal)
        int64 genomelen
        dict ref2genomelen
        int64 forward_sum, reverse_sum
        int64 forward_read_len_sum, reverse_read_len_sum
        dict ref2ncc_result
        # Calculation state
        str _chr
        int64 _forward_buff_size, _forward_sum, _reverse_sum, _forward_read_len_sum, _reverse_read_len_sum
        np.ndarray _ccbins
        list _forward_buff, _reverse_buff, _solved_chr
        int64 _fb_tail_pos, _last_pos, _last_forward_pos

    def __init__(self, int64 max_shift, references, lengths, logger_lock=None):
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
        self.forward_read_len_sum = self.reverse_read_len_sum = 0
        self.ref2ncc_result = {}
        # internal buff
        self._forward_buff_size = max_shift + 1
        self._chr = ''
        self._forward_sum = self._reverse_sum = 0
        self._ccbins = np.zeros(self._forward_buff_size, dtype=np.int64)
        self._forward_buff = [0] * self._forward_buff_size
        self._reverse_buff = [0] * self._forward_buff_size
        self._solved_chr = []

        self._init_pos_buff()

        self.logger_lock = logger_lock

    def _init_pos_buff(self):
        """Initialize position buffers for a new chromosome.

        Resets all internal buffers and counters for processing a new
        chromosome. This includes clearing cross-correlation bins,
        forward/reverse read buffers, and position tracking variables.
        """
        self._forward_sum = self._reverse_sum = 0
        self._forward_read_len_sum = self._reverse_read_len_sum = 0
        self._ccbins.fill(0)

        self._forward_buff = [0] * self._forward_buff_size
        self._fb_tail_pos = self._forward_buff_size

        self._reverse_buff = [0] * self._forward_buff_size

        self._last_pos = 0
        self._last_forward_pos = 0

    def flush(self, chrom: Optional[str] = None):
        """Finalize cross-correlation calculation for current chromosome.

        Completes the sliding window calculation for the current chromosome
        and stores the results as an NCCResult object.
        Called when switching to a new chromosome or finishing calculation.
        """
        if self._chr == '':
            return

        if chrom is not None and chrom != self._chr:
            self._chr = chrom

        self._shift_with_update(self._forward_buff_size)

        self.forward_sum += self._forward_sum
        self.reverse_sum += self._reverse_sum
        self.forward_read_len_sum += self._forward_read_len_sum
        self.reverse_read_len_sum += self._reverse_read_len_sum

        # Create NCCResult for this chromosome
        result = self.ref2ncc_result[self._chr] = NCCResult(
            max_shift=self.max_shift,
            read_len=self._forward_buff_size,
            genomelen=self.ref2genomelen[self._chr],
            forward_sum=self._forward_sum,
            reverse_sum=self._reverse_sum,
            forward_read_len_sum=self._forward_read_len_sum,
            reverse_read_len_sum=self._reverse_read_len_sum,
            ccbins=np.array(self._ccbins, dtype=np.int64)
        )
        result.calc_cc()
        self._init_pos_buff()

    cdef inline _check_pos(self, str chrom, int64 pos):
        """Validate read position and handle chromosome transitions.

        Ensures reads are sorted by position and handles switching between
        chromosomes. Raises ReadUnsortedError if reads are out of order.

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
            self._chr = chrom
            self._init_pos_buff()

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
        """Process a forward strand read for cross-correlation calculation.

        Incorporates a forward strand read into the sliding window algorithm.
        Updates forward read buffers and triggers cross-correlation updates
        as the window slides along the chromosome.

        Args:
            chrom: Chromosome name
            pos: Read start position (1-based)
            readlen: Read length in base pairs

        Note:
            Duplicate reads at the same position are ignored
        """
        self._check_pos(chrom, pos)

        if self._last_forward_pos == pos:
            return None  # duplicated read
        self._last_forward_pos = pos

        self._forward_sum += 1
        self._forward_read_len_sum += readlen

        offset = pos - self._fb_tail_pos

        if offset < 1:
            self._forward_buff[offset - 1] = 1
        else:
            self._shift_with_update(offset, update_forward=True)

    @wraparound(False)
    @boundscheck(False)
    def feed_reverse_read(self, str chrom, int64 pos, int64 readlen):
        """Process a reverse strand read for cross-correlation calculation.

        Incorporates a reverse strand read into the sliding window algorithm.
        The reverse read position is adjusted by read length to represent
        the 3' end position for proper cross-correlation calculation.

        Args:
            chrom: Chromosome name
            pos: Read start position (1-based)
            readlen: Read length in base pairs

        Note:
            Reverse read positions are converted to 3' end positions
            before processing in the sliding window
        """
        cdef int64 offset, revbuff_pos
        # cdef np.ndarray[int64] appendarray
        cdef int64 appendtail
        cdef list appendlist

        self._check_pos(chrom, pos)

        offset = pos - self._fb_tail_pos
        revbuff_pos = self.max_shift + readlen - 1

        if offset > 0:
            self._shift_with_update(offset)
        else:
            revbuff_pos += offset

        if len(self._reverse_buff) > revbuff_pos and self._reverse_buff[revbuff_pos] == 1:
            return None  # duplicated read

        self._reverse_sum += 1
        self._reverse_read_len_sum += readlen

        if len(self._reverse_buff) > revbuff_pos:
            self._reverse_buff[revbuff_pos] = 1
        else:
            appendtail = revbuff_pos - len(self._reverse_buff)
            appendlist = [0] * (appendtail + 1)
            appendlist[appendtail] = 1
            self._reverse_buff += appendlist

    @wraparound(False)
    @boundscheck(False)
    cdef inline _shift_with_update(self, int64 offset, bint update_forward=False):
        """Slide the window and update cross-correlation calculations.

        Advances the sliding window by the specified offset and updates
        cross-correlation bins based on overlapping read positions.
        This is the core operation of the sliding window algorithm.

        Args:
            offset: Distance to slide the window forward
            update_forward: Whether to update forward buffer during slide

        Note:
            This method performs the actual cross-correlation computation
            as the window slides across the chromosome
        """
        cdef int64 buff_max_shift = int64_min(self._forward_buff_size, offset)
        cdef int64 i, j

        for i in range(buff_max_shift):
            if self._forward_buff[i]:
                for j in range(int64_min(self._forward_buff_size, len(self._reverse_buff) - i)):
                    if self._reverse_buff[i + j]:
                        self._ccbins[j] += 1

        self._forward_buff = self._forward_buff[offset:] + [0] * buff_max_shift

        if update_forward:
            self._forward_buff[self.max_shift] = 1
        self._reverse_buff = self._reverse_buff[offset:]

        self._fb_tail_pos += offset

    def finishup_calculation(self):
        """Complete cross-correlation calculation for all chromosomes.

        Finalizes the calculation by flushing any remaining data in buffers
        and ensuring all cross-correlation results are properly stored.
        Creates empty results for chromosomes with no data to maintain
        genome length consistency.
        Called at the end of the analysis workflow.
        """
        self.flush()

        # Create empty results for chromosomes with no data
        for ref in self.ref2genomelen.keys():
            if ref not in self.ref2ncc_result:
                empty_result = EmptyNCCResult.create_empty(
                    self.ref2genomelen[ref], self.max_shift, self._forward_buff_size
                )
                self.ref2ncc_result[ref] = empty_result

    def get_result(self, chrom: str) -> NCCResult:
        return self.ref2ncc_result[chrom]

    def get_whole_result(self):
        """Get the genome-wide NCC calculation results.

        Returns:
            NCCGenomeWideResult: Complete results for all chromosomes including
                per-chromosome NCCResult objects and genome-wide statistics.
        """
        return NCCGenomeWideResult(
            genomelen=self.genomelen,
            forward_sum=self.forward_sum,
            reverse_sum=self.reverse_sum,
            chroms=self.ref2ncc_result.copy(),
            forward_read_len_sum=self.forward_read_len_sum,
            reverse_read_len_sum=self.reverse_read_len_sum
        )
