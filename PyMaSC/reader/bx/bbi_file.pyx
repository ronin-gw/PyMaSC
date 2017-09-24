#cython: profile=True
from libc.stdint cimport uint64_t

from bx.bbi.bbi_file cimport BlockHandler as bxBlockHandler, BBIFile as bxBBIFile
from bx.bbi.cirtree_file cimport CIRTreeFile
from bx.bbi.types cimport bits32

from .bpt_file cimport BPTFile

import zlib


cdef class BlockHandler(bxBlockHandler):
    cdef set_handle_block( self, str block_data, BBIFile bbi_file ):
        pass

    cdef interval next(self):
        pass


cdef class BBIFile(bxBBIFile):
    def open( self, file, expected_sig, type_name ):
        super(BBIFile, self).open( file, expected_sig, type_name )
        # Initialize and attach embedded BPTFile containing chromosome names and ids
        self.reader.seek( self.chrom_tree_offset )
        self.chrom_bpt = BPTFile( file=self.file )

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
        cdef str block_data

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
        for key_bytes, val_bytes in self.chrom_bpt.traverse():
            chrom = key_bytes.rstrip('\0')
            # The value is two 32 bit uints, use the BPT's reader for checking byteswapping
            # assert len( val_bytes ) == 8
            chrom_id, chrom_size = self.chrom_bpt.reader.unpack( "II", val_bytes )
            yield chrom, chrom_id, chrom_size
