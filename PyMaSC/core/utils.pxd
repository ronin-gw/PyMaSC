from bx.bbi.types cimport bits32
from libc.stdint cimport int_fast64_t as int64


cdef inline bits32 bits32_min(bits32 a, bits32 b):
    return a if a < b else b


cdef inline int64 int64_min(int64 a, int64 b):
    return a if a < b else b


cdef inline int64 int64_max(int64 a, int64 b):
    return a if a > b else b
