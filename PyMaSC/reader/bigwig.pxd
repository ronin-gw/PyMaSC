from PyMaSC.reader.bx.bigwig_file cimport BigWigFile


cdef class BigWigReader(object):
    cdef:
        object file
        bint closed
        BigWigFile bigwig
        readonly dict chromsizes

    cpdef BigWigFile fetch(self, float valfilter, chrom)
