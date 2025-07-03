from libc.stdint cimport uint32_t, int_fast64_t as int64


cdef inline uint32_t bits32_min(uint32_t a, uint32_t b):
    return a if a < b else b


cdef inline int64 int64_min(int64 a, int64 b):
    return a if a < b else b


cdef inline int64 int64_max(int64 a, int64 b):
    return a if a > b else b
