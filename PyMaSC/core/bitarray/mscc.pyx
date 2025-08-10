# cython: language_level=3
"""Mappability-Sensitive Cross-Correlation using BitArray operations.

This module provides Cython-optimized implementation of MSCC calculation
using bit array operations for maximum performance. The CCBitArrayCalculator
class efficiently processes read data and mappability information to compute
both naive and mappability-sensitive cross-correlation statistics.

Key functionality:
- High-performance bit array operations for genomic position tracking
- Simultaneous calculation of NCC and MSCC statistics
- Memory-efficient streaming processing of large datasets
- Integration with BigWig mappability data

This implementation is critical for processing large ChIP-seq datasets
where performance and memory efficiency are essential.
"""
import logging

cimport numpy as np
from libc.stdint cimport int_fast64_t as int64
from libc.string cimport strcmp, strcpy
from cython cimport boundscheck, wraparound

from .bitarray cimport bit_index_t, bitarray
from PyMaSC.reader.bigwig cimport BigWigFileIterator

import numpy as np

from ..exceptions import ReadUnsortedError
from PyMaSC.utils.progress import ProgressBar
from PyMaSC.result import (
    NCCResult, MSCCResult, BothChromResult,
    NCCGenomeWideResult, MSCCGenomeWideResult, BothGenomeWideResult,
    EmptyNCCResult, EmptyMSCCResult, EmptyBothChromResult
)

logger = logging.getLogger(__name__)


cdef class CCBitArrayCalculator(object):
    """High-performance MSCC calculator using bit array operations.

    Implements mappability-sensitive cross-correlation calculation using
    optimized bit array operations for maximum performance. Capable of
    computing both naive cross-correlation (NCC) and MSCC statistics
    simultaneously while efficiently managing memory usage.

    The calculator processes sequencing reads and mappability data to
    generate cross-correlation profiles that account for genomic regions
    with variable mappability. This addresses the 'phantom peak' problem
    in ChIP-seq quality assessment.

    Key features:
    - Bit array operations for O(1) position marking and correlation
    - Streaming processing to handle large datasets efficiently
    - Simultaneous NCC and MSCC calculation when mappability data available
    - Memory-efficient chromosome-by-chromosome processing
    - Thread-safe operation for multiprocessing environments

    Attributes:
        MAPPABILITY_THRESHOLD: Minimum mappability value to consider
        EXTRA_ALLOCATE_SIZE: Extra buffer space for variable read lengths
        max_shift: Maximum shift distance for cross-correlation
        read_len: Expected read length for calculations
        ref2genomelen: Chromosome name to length mapping
        genomelen: Total genome length
        forward_sum, reverse_sum: Total read counts by strand
        skip_ncc: Whether to skip naive cross-correlation calculation
    """
    cdef public:
        double MAPPABILITY_THRESHOLD
        int64 EXTRA_ALLOCATE_SIZE

    cdef readonly:
        int64 max_shift
        int64 read_len
        bint skip_ncc
        object logger_lock
    cdef:
        list references
        dict ref2genomelen
        int64 genomelen
        int64 forward_sum, reverse_sum
        int64 forward_read_len_sum, reverse_read_len_sum
        dict ref2ncc_result, ref2mscc_result

        str _chr
        list _solved_chr
        bitarray _forward_array, _reverse_array
        object _bwfeeder, _progress

        int64 _last_pos, _last_forward_pos, _last_reverse_pos, _array_extend_size
        int64 _forward_read_len_sum, _reverse_read_len_sum
        bint _buff_flashed, _allow_bitarray_progress

    def __init__(self, int64 max_shift, int64 read_len, references, lengths,
                 bwfeeder=None, skip_ncc=False, logger_lock=None, progress_bar=None):
        """Initialize MSCC calculator with specified parameters.

        Sets up the calculator for processing sequencing data with optional
        mappability correction. Initializes bit arrays and statistics
        containers for efficient cross-correlation computation.

        Args:
            max_shift: Maximum shift distance for cross-correlation
            read_len: Expected read length for calculations
            references: List of chromosome names
            lengths: List of chromosome lengths
            bwfeeder: BigWig feeder for mappability data (optional)
            skip_ncc: Whether to skip naive cross-correlation calculation
            logger_lock: Lock for thread-safe logging (optional)
            progress_bar: Progress reporting object (optional)
        """
        self.MAPPABILITY_THRESHOLD = 1.
        # bitarary extra allocate size for reads longer than self.read_len
        self.EXTRA_ALLOCATE_SIZE = 100

        self.max_shift = max_shift
        self.read_len = read_len
        self.references = references
        self.ref2genomelen = dict(zip(references, lengths))
        self.genomelen = sum(lengths)
        # stats
        self.ref2ncc_result = {}
        self.ref2mscc_result = {}
        self.forward_sum = self.reverse_sum = 0
        self.forward_read_len_sum = self.reverse_read_len_sum = 0

        # internal buff
        self._chr = ''
        self._solved_chr = []
        self._buff_flashed = False
        self._array_extend_size = self.read_len + self.max_shift + self.EXTRA_ALLOCATE_SIZE

        #
        self._bwfeeder = bwfeeder
        self.skip_ncc = skip_ncc
        self.logger_lock = logger_lock
        self._progress = progress_bar if progress_bar else ProgressBar()
        # Issue: bitarray progress stacks in multiprocessing
        self._allow_bitarray_progress = progress_bar is None

    def _logging_info(self, msg):
        if self.logger_lock:
            self.logger_lock.acquire()
        logger.info(msg)
        if self.logger_lock:
            self.logger_lock.release()

    def _init_mappability_progress(self):
        self._progress.enable_bar()
        self._progress.reset_progress("@@@@@^" * 12, '^', ' ', self._chr, self.ref2genomelen[self._chr])

    def _init_bitarray_progress(self):
        if self._allow_bitarray_progress:
            self._progress.reset_progress("xxxxxX" * 12, 'X', ' ', self._chr, self.max_shift)
        else:
            self._progress.disable_bar()

    def _init_buff(self):
        # To access by 1-based index, allocate arrays with 1 more size.
        # Additionally, extra spaces same or more than (read_length - 1)
        # is needed to contain reverse strand reads.
        chromsize = self.ref2genomelen[self._chr] + self._array_extend_size
        self._forward_array = bitarray(chromsize)
        self._reverse_array = bitarray(chromsize)

        self._last_pos = 0
        self._last_forward_pos = self._last_reverse_pos = 0
        self._forward_read_len_sum = self._reverse_read_len_sum = 0

    def flush(self, chrom: Optional[str] = None):
        if self._chr != '' and not self._buff_flashed:
            self._calc_correlation()
        if chrom is not None:
            self._fill_result(chrom)

        self._buff_flashed = True

    cdef inline int _fill_result(self, chrom: str) except -1:
        cdef np.ndarray zero_bins = np.zeros(self.max_shift + 1, dtype=np.int64)

        self._chr = chrom

        #
        if self._chr not in self.ref2ncc_result:
            result = self.ref2ncc_result[self._chr] = EmptyNCCResult.create_empty(
                self.ref2genomelen[self._chr], self.max_shift, self.read_len
            )

        #
        if not self._bwfeeder or self._chr in self.ref2mscc_result:
            return 0

        result = self.ref2mscc_result[self._chr] = EmptyMSCCResult.create_empty(
            self.ref2genomelen[self._chr], self.max_shift, self.read_len
        )

        try:
            mappability = self._load_mappability()
            if mappability is None:
                raise KeyError(self._chr)
        except KeyError:
            return 0

        mappable_len = []
        forward_mappability = mappability
        reverse_mappability = mappability.clone()

        self._logging_info("Calc {} mappable length...".format(self._chr))
        for i in xrange(self.max_shift + 1):
            mappable_len.append(forward_mappability.acount(reverse_mappability))
            reverse_mappability.rshift(1, 0)
        result.mappable_len = tuple(mappable_len)

    cdef inline int _calc_correlation(self) except -1:
        chromsize = self.ref2genomelen[self._chr] + 1

        cdef bit_index_t forward_sum, reverse_sum
        cdef list ccbin
        cdef list mappable_len, mappable_forward_sum, mappable_reverse_sum, mascbin
        cdef bitarray mappability, forward_mappability, reverse_mappability
        cdef bitarray mappable_forward = bitarray(chromsize)
        cdef bitarray mappable_reverse = bitarray(chromsize)
        cdef bitarray doubly_mappable_region = bitarray(chromsize)
        cdef list reverse_mappability_buff
        cdef int i

        #
        self.forward_read_len_sum += self._forward_read_len_sum
        self.reverse_read_len_sum += self._reverse_read_len_sum

        # naive cross-correlation
        if not self.skip_ncc:
            forward_sum = self._forward_array.count()
            reverse_sum = self._reverse_array.count()
            self.forward_sum += forward_sum
            self.reverse_sum += reverse_sum

            self.ref2ncc_result[self._chr] = NCCResult(
                max_shift=self.max_shift,
                read_len=self.read_len,
                genomelen=self.ref2genomelen[self._chr],
                forward_sum=forward_sum,
                reverse_sum=reverse_sum,
                forward_read_len_sum=self._forward_read_len_sum,
                reverse_read_len_sum=self._reverse_read_len_sum,
                ccbins=[]  # Assign later
            )
            ccbin = self.ref2ncc_result[self._chr].ccbins

        # mappability-sensitive cross-correlation
        try:
            mappability = self._load_mappability()
        except KeyError as e:
            self._logging_info("Mappability for '{}' not found. "
                               "Skip calc mappability sensitive CC.".format(e.args[0]))
            mappability = None
        if mappability:
            # Attributes will be assigned later
            result = self.ref2mscc_result[self._chr] = MSCCResult(
                max_shift=self.max_shift,
                read_len=self.read_len,
                genomelen=self.ref2genomelen[self._chr],
                forward_sum=[],
                reverse_sum=[],
                forward_read_len_sum=self._forward_read_len_sum,
                reverse_read_len_sum=self._reverse_read_len_sum,
                ccbins=[],
                mappable_len=[None] * self.read_len
            )

            mappable_len = result.mappable_len
            mappable_forward_sum = result.forward_sum
            mappable_reverse_sum = result.reverse_sum
            mascbin = result.ccbins

            forward_mappability = mappability
            reverse_mappability = mappability.clone()
            reverse_mappability_buff = reverse_mappability[:self.read_len - 1]
            reverse_mappability.rshift(self.read_len - 1, 0)

        # calc cross-correlation
        self._logging_info("Calculate cross-correlation for {}...".format(self._chr))
        self._init_bitarray_progress()

        for i in xrange(self.max_shift + 1):
            # calc mappability-sensitive cross-correlation
            if mappability:
                doubly_mappable_region.alloc_and(forward_mappability, reverse_mappability)
                if i < self.read_len:
                    mappable_len[self.read_len - i - 1] = doubly_mappable_region.count()
                elif i < self.read_len * 2 - 1:
                    # assert mappable_len[i - self.read_len + 1] == doubly_mappable_region.count()
                    pass
                else:
                    mappable_len.append(doubly_mappable_region.count())

                mappable_forward.alloc_and(self._forward_array, doubly_mappable_region)
                mappable_reverse.alloc_and(self._reverse_array, doubly_mappable_region)

                mappable_forward_sum.append(mappable_forward.count())
                mappable_reverse_sum.append(mappable_reverse.count())
                mascbin.append(mappable_forward.acount(mappable_reverse))

                reverse_mappability.lshift(
                    1,
                    reverse_mappability_buff.pop() if reverse_mappability_buff else 0
                )

            # calc naive cross-correlation
            if not self.skip_ncc:
                ccbin.append(self._forward_array.acount(self._reverse_array))

            self._reverse_array.rshift(1, 0)
            self._progress.update(i)

        #
        if mappability:
            self.ref2mscc_result[self._chr].calc_cc()
        if not self.skip_ncc:
            self.ref2ncc_result[self._chr].calc_cc()

        self._progress.clean()

    cdef inline bitarray _load_mappability(self):
        cdef bitarray mappability
        cdef BigWigFileIterator feeder
        cdef float _val

        if not self._bwfeeder:
            return None

        # KeyError raises if there is no mappability tracks for _chr
        feeder = self._bwfeeder.fetch(self.MAPPABILITY_THRESHOLD, self._chr)
        self._logging_info("Loading {} mappability to bit array...".format(self._chr))

        # Allocate same size as forward/reverse read array
        chromsize = self.ref2genomelen[self._chr] + self._array_extend_size
        mappability = bitarray(chromsize)
        self._init_mappability_progress()
        for begin, end, _val in feeder:
            mappability.set(begin + 1, end)
            self._progress.update(end)
        self._progress.update(self.ref2genomelen[self._chr])
        self._progress.clean()

        return mappability

    cdef inline _check_pos(self, str chrom, int64 pos):
        if chrom != self._chr:
            if self._chr != '':
                if chrom in self._solved_chr:
                    raise ReadUnsortedError
                self._solved_chr.append(self._chr)
                self.flush()
                self._buff_flashed = False
            self._chr = chrom
            self._init_buff()
            self._logging_info("Loading {} reads to bit array...".format(chrom))

        if pos < self._last_pos:
            raise ReadUnsortedError

        self._last_pos = pos

    @wraparound(False)
    @boundscheck(False)
    def feed_forward_read(self, str chrom, int64 pos, int64 readlen):
        """Process a forward strand read at specified position.

        Marks the read's 5' end position in the forward strand bit array
        for cross-correlation calculation. Handles chromosome switching
        and position validation automatically.

        Args:
            chrom: Chromosome name
            pos: 1-based genomic position of read's 5' end
            readlen: Length of the read in base pairs

        Note:
            Automatically triggers correlation calculation when switching
            between chromosomes. Ignores duplicate reads at same position.
        """
        self._check_pos(chrom, pos)

        if self._last_forward_pos == pos:
            return None  # duplicated read
        self._last_forward_pos = pos

        self._forward_read_len_sum += readlen
        self._forward_array[pos] = 1

    @wraparound(False)
    @boundscheck(False)
    def feed_reverse_read(self, str chrom, int64 pos, int64 readlen):
        """Process a reverse strand read at specified position.

        Marks the read's 3' end position (corresponding to 5' end of
        the reverse complement) in the reverse strand bit array for
        cross-correlation calculation.

        Args:
            chrom: Chromosome name
            pos: 1-based genomic position of read's 5' end
            readlen: Length of the read in base pairs

        Note:
            For reverse reads, marks position at pos + readlen - 1 to
            represent the 5' end of the reverse complement strand.
            Prevents duplicate counting at the same effective position.
        """
        self._check_pos(chrom, pos)

        if self._reverse_array[pos + readlen - 1] == 0:
            self._reverse_array[pos + readlen - 1] = 1
            self._reverse_read_len_sum += readlen

    def finishup_calculation(self):
        """Complete cross-correlation calculation for all chromosomes.

        Finalizes the calculation by processing any remaining buffered
        data and calculating mappability statistics for chromosomes
        with no aligned reads. This ensures complete genome-wide
        statistics are available.

        Note:
            Must be called after all reads have been processed to
            ensure complete and accurate cross-correlation results.
        """
        cdef bitarray mappability, forward_mappability, reverse_mappability
        cdef int i

        self.flush(self._chr)

        #
        for chrom in self.references:
            self._fill_result(chrom)

    def get_result(self, chrom: str) -> BothChromResult:
        if chrom not in self.ref2ncc_result and chrom not in self.ref2mscc_result:
            raise KeyError(chrom)
        return BothChromResult(
            chrom=self.ref2ncc_result.get(chrom),
            mappable_chrom=self.ref2mscc_result.get(chrom)
        )

    def get_whole_result(self):
        """Get the genome-wide NCC and MSCC calculation results.

        Returns:
            MSCCGenomeWideResult: Complete results for all chromosomes including
                per-chromosome MSCCResult objects and genome-wide statistics.
        """
        if not self.ref2mscc_result:
            assert self.ref2ncc_result, \
                "No results available for either NCC or MSCC."
            return NCCGenomeWideResult(
                genomelen=self.genomelen,
                forward_sum=self.forward_sum,
                reverse_sum=self.reverse_sum,
                chroms=self.ref2ncc_result.copy(),
                forward_read_len_sum=self.forward_read_len_sum,
                reverse_read_len_sum=self.reverse_read_len_sum
            )
        elif not self.ref2ncc_result:
            return MSCCGenomeWideResult(
                genomelen=self.genomelen,
                chroms=self.ref2mscc_result.copy(),
                forward_read_len_sum=self.forward_read_len_sum,
                reverse_read_len_sum=self.reverse_read_len_sum
            )
        else:
            return BothGenomeWideResult(
                genomelen=self.genomelen,
                forward_sum=self.forward_sum,
                reverse_sum=self.reverse_sum,
                chroms=self.ref2ncc_result.copy(),
                mappable_chroms=self.ref2mscc_result.copy(),
                forward_read_len_sum=self.forward_read_len_sum,
                reverse_read_len_sum=self.reverse_read_len_sum
            )
