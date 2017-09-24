from bx.bbi.bbi_file cimport BlockHandler
from bx.bbi.types cimport bits32, bits16, UBYTE
from .bbi_file cimport BBIFile

from .bbi_file cimport interval


cdef class BigWigFile( BBIFile ):
    cpdef get_chrom_sizes( self )
    cpdef get( self, char * chrom, bits32 start, bits32 end, float valfilter )
