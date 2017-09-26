from PyMaSC.reader.bx.bigwig_file cimport BigWigFile


cdef class BigWigReader(object):
    cdef readonly:
        str path
        object file
        bint closed
        BigWigFile bigwig
        dict chromsizes

    cpdef BigWigFile fetch(self, float valfilter, chrom)
