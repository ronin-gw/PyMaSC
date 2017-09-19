from .kentlib.files cimport bbiFile
from .kentlib.types cimport bits32


cdef extern from "bbiFile.h":
    struct bbiChromInfo:
        bbiChromInfo *next
        char *name
        bits32 id
        bits32 size


cdef class BigWigFile:
    cdef:
        readonly char path[1024]
        readonly bint closed
        readonly dict chromsizes

        bbiFile *file
        bbiChromInfo *chrominfo
