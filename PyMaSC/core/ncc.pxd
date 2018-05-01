from libc.stdint cimport int_fast64_t as int64
cimport numpy as np


cdef class NaiveCCCalculator(object):
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
        int64 _forward_buff_size, _forward_sum, _reverse_sum
        np.ndarray _ccbins
        list _forward_buff, _reverse_buff, _solved_chr
        int64 _fb_tail_pos, _last_pos, _last_forward_pos
        object logger_lock

    # def __init__(self, int64 max_shift, references, lengths)
    # def _init_pos_buff(self)
    # def _flush(self)
    cdef inline _check_pos(self, str chrom, int64 pos)
    # def feed_forward_read(self, char* chrom, int64 pos, int64 readlen):
    # def feed_reverse_read(self, char* chrom, int64 pos, int64 readlen)
    cdef inline _shift_with_update(self, int64 offset, bint update_forward=False)
    # def finishup_calculation(self)
