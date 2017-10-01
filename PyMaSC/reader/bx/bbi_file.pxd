from bx.bbi.types cimport bits16, bits32, bits64

from .bpt_file cimport BPTFile
from PyMaSC.utils.bx.binary_file cimport BinaryFileReader


cdef struct interval:
    bits32 begin, end
    float value

cdef class BlockHandler(object):
    cdef handle_block( self, bytes block_data, BBIFile bbi_file )
    cdef set_handle_block( self, bytes block_data, BBIFile bbi_file )
    cdef interval next(self)


cdef class BBIFile(object):
    cdef:
        object file
        BinaryFileReader reader
        public bits16 zoom_levels
        chrom_tree_offset
        bits64 unzoomed_index_offset
        bits32 uncompress_buf_size

        BPTFile chrom_bpt

        list _block_list
        BlockHandler _handler

    cdef _get_chrom_id_and_size( self, char * chrom )
    cdef set_visitor_blocks_in_region( self, bits32 chrom_id, bits32 start, bits32 end, BlockHandler handler )
    cdef int _set_next_block(self) except *
    cdef interval next(self) except *
