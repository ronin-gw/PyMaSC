from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t


cdef class BinaryFileReader( object ):
    cdef:
        bint is_little_endian, byteswap_needed
        object file, endian_code

    cpdef read_c_string( self )
    cpdef read_raw_array( self, dtype, size )
    cpdef tell( self )
    cpdef skip( self, count )
    cpdef uint8_t read_uint8( self )
    cpdef uint16_t read_uint16( self )
    cpdef uint32_t read_uint32( self )
    cpdef uint64_t read_uint64( self )
    cpdef float read_float( self )
