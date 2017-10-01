from bx.bbi.types cimport bits32
from .bbi_file cimport BBIFile


cdef class BigWigFile( BBIFile ):
    cpdef get_chrom_sizes( self )
    cpdef get( self, str chrom, bits32 start, bits32 end, float valfilter )
