from libc.stdint cimport int_fast64_t as int64


cdef class NaiveCCCalculator(object):
    cdef:
        readonly int64 max_shift
        readonly dict ref2genomelen
        readonly int64 genomelen
        readonly int64 forward_sum, reverse_sum
        readonly dict ref2forward_sum, ref2reverse_sum
        readonly int64 forward_read_len_sum, reverse_read_len_sum
        readonly list ccbins
        readonly dict ref2ccbins

        char _chr[1024]
        int64 _forward_sum, _reverse_sum
        object _ccbins, _forward_buff, _reverse_buff
        int64 _fb_tail_pos, _last_pos, _last_forward_pos, _last_reverse_pos

    # def __init__(self, int64 max_shift, references, lengths)
    # def _init_pos_buff(self)
    # def _flush(self)
    cdef inline _check_pos(self, char* chrom, int64 pos)
    # def feed_forward_read(self, char* chrom, int64 pos, int64 readlen):
    # def feed_reverse_read(self, char* chrom, int64 pos, int64 readlen)
    cdef inline _shift_with_update(self, int64 offset, bint update_forward=False)
    # def finishup_calculation(self)
