# Define interval structure to match original API
cdef struct interval:
    unsigned int begin, end
    float value


cdef class BigWigFileIterator(object):
    cdef:
        list _intervals
        int _index

    cdef interval next(self)


cdef class BigWigReader(object):
    cdef:
        public str path
        object file
        bint closed
        public dict chromsizes

    cpdef BigWigFileIterator fetch(self, float valfilter, chrom)