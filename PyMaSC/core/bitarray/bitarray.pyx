# cython: language_level=3
"""BitArray operations for efficient genomic position tracking.

This module provides Cython-optimized BitArray data structures for
efficient tracking of genomic positions in PyMaSC cross-correlation
calculations. The BitArray approach offers significant performance
benefits for large-scale genomic analysis.

Key components:
- bitarray: Cython wrapper for C BitArray library
- Bit manipulation operations (set, clear, count, shift)
- Efficient memory usage for genomic coordinate tracking
- High-performance bitwise operations for correlation calculations

The BitArray implementation interfaces with external C library routines
for maximum performance in genomic position tracking and manipulation.
"""
from libc.stdint cimport uint64_t


cdef inline word_addr_t min(word_addr_t a, word_addr_t b):
    return a if a < b else b


cdef class bitarray(object):
    """Cython wrapper for C BitArray library.

    Provides efficient bit manipulation operations for tracking genomic
    positions in PyMaSC cross-correlation calculations. The class wraps
    the external C BitArray library to provide high-performance bitwise
    operations on large arrays.

    The BitArray data structure is particularly effective for genomic
    applications where memory efficiency and fast bit operations are
    critical for processing large datasets.

    Attributes:
        array: Pointer to underlying C BitArray structure
    """
    def __init__(self, bit_index_t nbits):
        """Initialize BitArray with specified number of bits.

        Args:
            nbits: Number of bits to allocate in the array

        Note:
            Creates underlying C BitArray structure if nbits > 0
        """
        if nbits > 0:
            self.array = bit_array_create(nbits)

    def __getitem__(self, obj):
        """Get bit value(s) at specified index or slice.

        Args:
            obj: Integer index or slice object

        Returns:
            Boolean value for single index, list of values for slice
        """
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
        """Set bit at specified index.

        Args:
            key: Bit index to set
            value: Value to set (currently unused - always sets bit)
        """
        bit_array_set_bit(self.array, key)

    def __dealloc__(self):
        """Free underlying C BitArray memory.

        Called automatically when BitArray object is garbage collected.
        """
        bit_array_free(self.array)

    cdef void set(self, bit_index_t from_, bit_index_t to):
        """Set bits in specified range.

        Args:
            from_: Starting bit index (inclusive)
            to: Ending bit index (exclusive)
        """
        bit_array_set_region(self.array, from_, to - from_ + 1)

    cdef void clear(self):
        """Clear all bits in the array."""
        bit_array_clear_all(self.array)

    cdef bit_index_t count(self):
        """Count number of set bits in the array.

        Returns:
            Number of bits set to 1
        """
        return bit_array_num_bits_set(self.array)

    cdef bit_index_t acount(self, bitarray other):
        """Count bits set in both arrays (AND operation).

        Performs bitwise AND between this array and another array,
        then counts the number of set bits in the result.

        Args:
            other: BitArray to AND with

        Returns:
            Number of bits set in both arrays
        """
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
        """Create a copy of this BitArray.

        Returns:
            New BitArray instance with identical bit pattern
        """
        cdef bitarray a = bitarray(0)
        cdef BIT_ARRAY *clone = bit_array_clone(self.array)
        a.array = clone
        return a

    cdef void rshift(self, bit_index_t shift_dist, char fill):
        """Right shift bits by specified distance.

        Args:
            shift_dist: Number of positions to shift right
            fill: Fill value for shifted-in bits
        """
        bit_array_shift_right(self.array, shift_dist, fill)

    cdef void lshift(self, bit_index_t shift_dist, char fill):
        """Left shift bits by specified distance.

        Args:
            shift_dist: Number of positions to shift left
            fill: Fill value for shifted-in bits
        """
        bit_array_shift_left(self.array, shift_dist, fill)

    cdef void alloc_and(bitarray self, bitarray src1, bitarray src2):
        """Allocate and perform AND operation between two BitArrays.

        Allocates memory for this array and stores the result of
        bitwise AND operation between src1 and src2.

        Args:
            src1: First source BitArray
            src2: Second source BitArray
        """
        bit_array_and(self.array, src1.array, src2.array)
