from libc.stdint cimport uint64_t


cdef extern from "bitcnt.h":
    cdef uint64_t popcnt(uint64_t n)


cdef extern from "bit_array.h":
    ctypedef uint64_t word_t
    ctypedef uint64_t word_addr_t
    ctypedef uint64_t bit_index_t

    struct BIT_ARRAY:
        word_t* words
        bit_index_t num_of_bits
        word_addr_t num_of_words
        word_addr_t capacity_in_words

    # Constructor and Deconstructor
    BIT_ARRAY* bit_array_create(bit_index_t nbits)
    void bit_array_free(BIT_ARRAY* bitarray);
    #
    void bit_array_set_bit(BIT_ARRAY* bitarr, bit_index_t b)
    void bit_array_set_region(BIT_ARRAY* bitarr, bit_index_t start, bit_index_t len)
    void bit_array_clear_all(BIT_ARRAY* bitarr)
    #
    bit_index_t bit_array_num_bits_set(const BIT_ARRAY* bitarr)
    #
    BIT_ARRAY* bit_array_clone(const BIT_ARRAY* bitarr)
    #
    void bit_array_and(BIT_ARRAY* dest, const BIT_ARRAY* src1, const BIT_ARRAY* src2)
    #
    void bit_array_shift_right(BIT_ARRAY* bitarr, bit_index_t shift_dist, char fill)
    void bit_array_shift_left(BIT_ARRAY* bitarr, bit_index_t shift_dist, char fill)
    #
    char bit_array_get_bit(BIT_ARRAY* bitarr, bit_index_t b);

cdef class bitarray(object):
    cdef:
        BIT_ARRAY *array

    cdef void set(self, bit_index_t from_, bit_index_t to)
    cdef void clear(self)
    cdef bit_index_t count(self)
    cdef bit_index_t acount(self, bitarray other)
    cdef bitarray clone(self)
    cdef void rshift(self, bit_index_t shift_dist, char fill)
    cdef void lshift(self, bit_index_t shift_dist, char fill)
    cdef void alloc_and(self, bitarray src1, bitarray src2)
