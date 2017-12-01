import logging

from libc.stdint cimport int_fast64_t as int64
from libc.string cimport strcmp, strcpy
from cython cimport boundscheck, wraparound

from bitarray import bitarray

from PyMaSC.core.ncc import ReadUnsortedError
from PyMaSC.utils.compatible import tostr

logger = logging.getLogger(__name__)


cdef class NCCBitArrayCalculator(object):
    cdef readonly:
        int64 max_shift
        dict ref2genomelen
        int64 genomelen
        int64 forward_sum, reverse_sum
        dict ref2forward_sum, ref2reverse_sum
        int64 forward_read_len_sum, reverse_read_len_sum
        dict ref2ccbins
    cdef:
        str _chr
        list _solved_chr
        object _forward_array, _reverse_array, logger_lock

        int64 _last_pos, _last_forward_pos, _last_reverse_pos
        bint _buff_flashed

    def __init__(self, int64 max_shift, references, lengths, logger_lock=None):
        """
        self._ccbins: summation bins to calc cc; shift from 0 to max_shift
        """
        self.max_shift = max_shift
        self.ref2genomelen = dict(zip(references, lengths))
        self.genomelen = sum(lengths)
        # stats
        self.forward_sum = self.reverse_sum = 0
        self.ref2forward_sum = {ref: 0 for ref in references}
        self.ref2reverse_sum = {ref: 0 for ref in references}
        self.forward_read_len_sum = self.reverse_read_len_sum = 0
        self.ref2ccbins = {ref: None for ref in references}
        # internal buff
        self._chr = ''
        self._solved_chr = []
        self._buff_flashed = False

        self.logger_lock = logger_lock

    def _init_buff(self):
        # To access by 1-based index, allocate arrays with 1 more size.
        chromsize = self.ref2genomelen[self._chr] + 1
        self._forward_array = bitarray(chromsize)
        self._reverse_array = bitarray(chromsize)
        self._forward_array.setall(0)
        self._reverse_array.setall(0)

        self._last_pos = 0
        self._last_forward_pos = self._last_reverse_pos = 0

    def flush(self):
        if self._buff_flashed:
            return None

        if self.logger_lock:
            self.logger_lock.acquire()
        logger.info("Calculate cross-correlation for {}...".format(self._chr))
        if self.logger_lock:
            self.logger_lock.release()

        forward_sum = self._forward_array.count()
        reverse_sum = self._reverse_array.count()
        ccbins = [(self._forward_array & self._reverse_array).count()]
        ccbins += [
            (self._forward_array[:-i] & self._reverse_array[i:]).count()
            for i in xrange(1, self.max_shift + 1)
        ]
        print ccbins

        self.forward_sum += forward_sum
        self.reverse_sum += reverse_sum
        self.ref2forward_sum[self._chr] = forward_sum
        self.ref2reverse_sum[self._chr] = reverse_sum
        self.ref2ccbins[self._chr] = tuple(ccbins)

        self._buff_flashed = True

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

            if self.logger_lock:
                self.logger_lock.acquire()
            logger.info("Loading {} reads to bit array...".format(chrom))
            if self.logger_lock:
                self.logger_lock.release()

        if pos < self._last_pos:
            raise ReadUnsortedError

        self._last_pos = pos

    @wraparound(False)
    @boundscheck(False)
    def feed_forward_read(self, str chrom, int64 pos, int64 readlen):
        self._check_pos(chrom, pos)

        if self._last_forward_pos == pos:
            return None  # duplicated read
        self._last_forward_pos = pos

        self.forward_read_len_sum += readlen
        self._forward_array[pos] = 1

    @wraparound(False)
    @boundscheck(False)
    def feed_reverse_read(self, str chrom, int64 pos, int64 readlen):
        self._check_pos(chrom, pos)

        if self._last_reverse_pos == pos:
            return None  # duplicated read
        self._last_reverse_pos = pos

        self.reverse_read_len_sum += readlen
        self._reverse_array[pos + readlen - 1] = 1

    def finishup_calculation(self):
        self.flush()
