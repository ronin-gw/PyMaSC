from libc.stdint cimport uint64_t


cdef inline word_addr_t min(word_addr_t a, word_addr_t b):
    return a if a < b else b


cdef class bitarray(object):
    def __init__(self, bit_index_t nbits):
        if nbits > 0:
            self.array = bit_array_create(nbits)

    def __getitem__(self, obj):
        if isinstance(obj, slice):
            if obj.start:
                if obj.step:
                    r = range(obj.start, obj.stop, obj.step)
                else:
                    r = range(obj.start, obj.stop)
            else:
                r = range(obj.stop)
            return [bit_array_get_bit(self.array, i) for i in r]
        return bit_array_get_bit(self.array, obj)

    def __setitem__(self, bit_index_t key, value):
        bit_array_set_bit(self.array, key)

    def __dealloc__(self):
        bit_array_free(self.array)

    cdef void set(self, bit_index_t from_, bit_index_t to):
        bit_array_set_region(self.array, from_, to - from_ + 1)

    cdef void clear(self):
        bit_array_clear_all(self.array)

    cdef bit_index_t count(self):
        return bit_array_num_bits_set(self.array)

    cdef bit_index_t acount(self, bitarray other):
        cdef word_t b1, b2
        cdef word_addr_t i
        cdef bit_index_t num_of_bits_set = 0

        cdef word_addr_t min_words = min(self.array.num_of_words, other.array.num_of_words)

        for i in range(min_words):
            b1 = self.array.words[i]
            b2 = other.array.words[i]
            if b1 > 0 and b2 > 0:
                num_of_bits_set += popcnt(b1 & b2)

        return num_of_bits_set

    cdef bitarray clone(self):
        cdef bitarray a = bitarray(0)
        cdef BIT_ARRAY *clone = bit_array_clone(self.array)
        a.array = clone
        return a

    cdef void rshift(self, bit_index_t shift_dist, char fill):
        bit_array_shift_right(self.array, shift_dist, fill)

    cdef void lshift(self, bit_index_t shift_dist, char fill):
        bit_array_shift_left(self.array, shift_dist, fill)

    cdef void alloc_and(bitarray self, bitarray src1, bitarray src2):
        bit_array_and(self.array, src1.array, src2.array)
