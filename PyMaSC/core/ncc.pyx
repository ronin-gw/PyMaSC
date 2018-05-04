import logging

cimport numpy as np
from libc.stdint cimport int_fast64_t as int64
from libc.string cimport strcmp, strcpy
from cython cimport boundscheck, wraparound

from .utils cimport int64_min

import numpy as np

logger = logging.getLogger(__name__)


class ReadUnsortedError(IndexError):
    pass


cdef class NaiveCCCalculator(object):
    # cdef readonly:
    #     int64 max_shift
    #     dict ref2genomelen
    #     int64 genomelen
    #     int64 forward_sum, reverse_sum
    #     dict ref2forward_sum, ref2reverse_sum
    #     int64 forward_read_len_sum, reverse_read_len_sum
    #     dict ref2ccbins
    # cdef:
    #     char _chr[1024]
    #     int64 _forward_buff_size, _forward_sum, _reverse_sum
    #     np.ndarray _ccbins
    #     list _forward_buff, _reverse_buff
    #
    #     int64 _fb_tail_pos, _last_pos, _last_forward_pos

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
        self.ref2forward_sum = {ref: 0 for ref in references}
        self.ref2reverse_sum = {ref: 0 for ref in references}
        self.forward_read_len_sum = self.reverse_read_len_sum = 0
        self.ref2ccbins = {ref: None for ref in references}
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
        self._forward_sum = self._reverse_sum = 0
        self._ccbins.fill(0)

        self._forward_buff = [0] * self._forward_buff_size
        self._fb_tail_pos = self._forward_buff_size

        self._reverse_buff = [0] * self._forward_buff_size

        self._last_pos = 0
        self._last_forward_pos = 0

    def flush(self):
        self._shift_with_update(self._forward_buff_size)

        self.forward_sum += self._forward_sum
        self.reverse_sum += self._reverse_sum
        self.ref2forward_sum[self._chr] = self._forward_sum
        self.ref2reverse_sum[self._chr] = self._reverse_sum
        self.ref2ccbins[self._chr] = tuple(self._ccbins)

    cdef inline _check_pos(self, str chrom, int64 pos):
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
        self._check_pos(chrom, pos)

        if self._last_forward_pos == pos:
            return None  # duplicated read
        self._last_forward_pos = pos

        self._forward_sum += 1
        self.forward_read_len_sum += readlen

        offset = pos - self._fb_tail_pos

        if offset < 1:
            self._forward_buff[offset - 1] = 1
        else:
            self._shift_with_update(offset, update_forward=True)

    @wraparound(False)
    @boundscheck(False)
    def feed_reverse_read(self, str chrom, int64 pos, int64 readlen):
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
        self.reverse_read_len_sum += readlen

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
        self.flush()
