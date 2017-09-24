from bx.bbi.bbi_file cimport BlockHandler as bxBlockHandler, BBIFile as bxBPTFile
from bx.bbi.types cimport bits32


cdef struct interval:
    bits32 begin, end
    float value

cdef class BlockHandler(bxBlockHandler):
    cdef set_handle_block( self, str block_data, BBIFile bbi_file )
    cdef interval next(self)


cdef class BBIFile(bxBPTFile):
    cdef:
        list _block_list
        BlockHandler _handler

    cdef set_visitor_blocks_in_region( self, bits32 chrom_id, bits32 start, bits32 end, BlockHandler handler )
    cdef int _set_next_block(self) except *
    cdef interval next(self) except *
