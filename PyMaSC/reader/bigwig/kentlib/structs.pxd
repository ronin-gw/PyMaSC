from PyMaSC.reader.bigwig.kentlib.types cimport bits64


cdef extern from "common.h":
    struct fileOffsetSize:
        fileOffsetSize *next
        bits64	offset
        bits64	size
