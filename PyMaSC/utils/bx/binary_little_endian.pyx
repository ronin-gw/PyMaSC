from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t

import struct
import sys

import numpy


cdef union c2uint16_t:
    uint8_t cdata[2]
    uint16_t unpack

cdef union c4uint32_t:
    uint8_t cdata[4]
    uint32_t unpack

cdef union c8uint64_t:
    uint8_t cdata[8]
    uint64_t unpack

cdef union c4float:
    uint8_t cdata[4]
    float unpack


class BadMagicNumber( IOError ):
    pass


cdef class BinaryFileReader( object ):
    """
    Wrapper for doing binary reads on any file like object.
    """
    def __init__( self, file, magic = None, is_little_endian = False ):
        self.is_little_endian = is_little_endian
        self.file = file
        if magic is not None:
            # Attempt to read magic number and chuck endianess
            bytes = file.read( 4 )
            if struct.unpack( ">I", bytes )[0] == magic:
                pass
            elif struct.unpack( "<I", bytes )[0] == magic:
                self.is_little_endian = True
            else:
                raise BadMagicNumber( "File does not have expected magic number: %x != %x or %x" \
                        % ( magic, struct.unpack( ">I", bytes )[0], struct.unpack( "<I", bytes )[0] ) )
        # Set endian code
        if self.is_little_endian:
            self.endian_code = "<"
            self.byteswap_needed = ( sys.byteorder != "little" )
        else:
            self.endian_code = ">"
            self.byteswap_needed = ( sys.byteorder != "big" )

    def unpack( self, format, buffer, byte_count=None ):
        pattern = "%s%s" % ( self.endian_code, format )
        if byte_count is None:
            byte_count = struct.calcsize( pattern )
        return struct.unpack( pattern, buffer )

    def read_and_unpack( self, str format, byte_count=None ):
        """
        Read enough bytes to unpack according to `format` and return the
        tuple of unpacked values.
        """
        cdef str pattern

        pattern = "%s%s" % ( self.endian_code, format )
        if byte_count is None:
            byte_count = struct.calcsize( pattern )
        return struct.unpack( pattern, self.file.read( byte_count ) )

    cpdef read_c_string( self ):
        """
        Read a zero terminated (C style) string
        """
        rval = []
        while 1:
            ch = self.file.read(1)
            assert len( ch ) == 1, "Unexpected end of file"
            if ch == '\0':
                break
            rval.append( ch )
        return ''.join( rval )

    cpdef read_raw_array( self, dtype, size ):
        a = numpy.fromfile( self.file, dtype=dtype, count=size )
        if self.byteswap_needed:
            a.byteswap()
        return a

    def read( self, byte_count=1 ):
        return self.file.read( byte_count )

    cpdef tell( self ):
        return self.file.tell()

    cpdef skip( self, count ):
        self.file.seek( count, 1 )

    def seek( self, pos, whence=0 ):
        return self.file.seek( pos, whence )

    cpdef uint8_t read_uint8( self ):
        cdef bytes cdata = self.file.read(1)
        return (<uint8_t*>(<char*>cdata))[0]

    cpdef uint16_t read_uint16( self ):
        cdef c2uint16_t data

        cdef bytes cdata = self.file.read(2)
        data.cdata = <char*>cdata
        return data.unpack

    cpdef uint32_t read_uint32( self ):
        cdef c4uint32_t data

        cdef bytes cdata = self.file.read(4)
        data.cdata = <char*>cdata
        return data.unpack

    cpdef uint64_t read_uint64( self ):
        cdef c8uint64_t data

        cdef bytes cdata = self.file.read(8)
        data.cdata = <char*>cdata
        return data.unpack

    cpdef float read_float( self ):
        cdef c4float data

        cdef bytes cdata = self.file.read(4)
        data.cdata = <char*>cdata
        return data.unpack
