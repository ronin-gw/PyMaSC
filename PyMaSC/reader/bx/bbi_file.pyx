from libc.stdint cimport uint64_t

from bx.bbi.cirtree_file cimport CIRTreeFile
from bx.bbi.types cimport bits32

from .bpt_file cimport BPTFile
from PyMaSC.utils.bx.binary_file cimport BinaryFileReader

import zlib


cdef class BlockHandler(object):
    cdef handle_block( self, bytes block_data, BBIFile bbi_file ):
        pass

    cdef set_handle_block( self, bytes block_data, BBIFile bbi_file ):
        pass

    cdef interval next(self):
        pass


cdef class BBIFile(object):
    def __init__( self, file=None, expected_sig=None, type_name=None ):
        if file is not None:
            self.open( file, expected_sig, type_name )

    def open( self, file, expected_sig, type_name ):
        """
        Initialize from an existing bbi file, signature (magic) must be passed
        in since this is generic.
        """
        assert expected_sig is not None
        self.file = file
        # Open the file in a BinaryFileReader, handles magic and byteswapping
        self.reader = reader = BinaryFileReader( file, expected_sig )
        # self.magic = expected_sig
        # self.is_byteswapped = self.reader.byteswap_needed
        # Read header stuff
        # self.version = reader.read_uint16()
        reader.read_uint16()
        self.zoom_levels = reader.read_uint16()
        self.chrom_tree_offset = reader.read_uint64()
        # self.unzoomed_data_offset = reader.read_uint64()
        reader.read_uint64()
        self.unzoomed_index_offset = reader.read_uint64()
        # self.field_count = reader.read_uint16()
        reader.read_uint16()
        # self.defined_field_count = reader.read_uint16()
        reader.read_uint16()
        # self.as_offset = reader.read_uint64()
        reader.read_uint64()
        # self.total_summary_offset = reader.read_uint64()
        reader.read_uint64()
        self.uncompress_buf_size = reader.read_uint32()
        # Skip reserved
        reader.seek( 64 )
        # Read zoom headers
        # self.level_list = []
        for i from 0 <= i < self.zoom_levels:
            # level = ZoomLevel()
            # level.bbi_file = self
            # level.reduction_level = reader.read_uint32()
            reader.read_uint32()
            # level.reserved = reader.read_uint32()
            reader.read_uint32()
            # level.data_offset = reader.read_uint64()
            reader.read_uint64()
            # level.index_offset = reader.read_uint64()
            reader.read_uint64()
            # self.level_list.append( level )
        # Initialize and attach embedded BPTFile containing chromosome names and ids
        reader.seek( self.chrom_tree_offset )
        self.chrom_bpt = BPTFile( file=self.file )

    cdef _get_chrom_id_and_size( self, char * chrom ):
        """
        Lookup id and size from the chromosome named `chrom`
        """
        bytes = self.chrom_bpt.find( chrom )
        if bytes is not None:
            # The value is two 32 bit uints, use the BPT's reader for checking byteswapping
            assert len( bytes ) == 8
            chrom_id, chrom_size = self.chrom_bpt.reader.unpack( "II", bytes )
            return chrom_id, chrom_size
        else:
            return None, None

    cdef set_visitor_blocks_in_region( self, bits32 chrom_id, bits32 start, bits32 end, BlockHandler handler ):
        """
        Visit each block from the full data that overlaps a specific region
        """
        cdef CIRTreeFile ctf
        reader = self.reader
        reader.seek( self.unzoomed_index_offset )
        ctf = CIRTreeFile( reader.file )
        self._block_list = ctf.find_overlapping_blocks( chrom_id, start, end )
        self._handler = handler
        self._set_next_block()

    cdef inline int _set_next_block(self) except *:
        cdef uint64_t offset, size
        cdef bytes block_data

        try:
            offset, size = self._block_list.pop(0)
        except IndexError:
            return 1

        # Seek to and read all data for the block
        self.reader.seek( offset )
        block_data = self.reader.read( size )
        # Might need to uncompress
        if self.uncompress_buf_size > 0:
            block_data = zlib.decompress( block_data )
        self._handler.set_handle_block( block_data, self )

        return 0

    def __iter__(self):
        return self

    cdef interval next(self) except *:
        cdef interval result
        cdef int status

        while True:
            result = self._handler.next()
            if result.end == 0:
                status = self._set_next_block()
                if status == 1:
                    return interval(0, 0, 0)
            else:
                return result

    def __next__(self):
        cdef interval result = self.next()
        if result.end == 0:
            raise StopIteration
        return result.begin, result.end, result.value

    def _yield_all_chrom_key_id_and_size( self ):
        """
        """
        cdef bytes key_bytes
        for key_bytes, val_bytes in self.chrom_bpt.traverse():
            chrom = key_bytes.decode("utf-8").rstrip('\0')
            # The value is two 32 bit uints, use the BPT's reader for checking byteswapping
            # assert len( val_bytes ) == 8
            chrom_id, chrom_size = self.chrom_bpt.reader.unpack( "II", val_bytes )
            yield chrom, chrom_id, chrom_size
