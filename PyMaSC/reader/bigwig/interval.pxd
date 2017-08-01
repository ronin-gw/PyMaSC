from .kentlib.types cimport bits32, bits16, boolean, UBYTE
from .kentlib.structs cimport fileOffsetSize
from .kentlib.files cimport bbiFile, udcFile


cdef struct Interval:
    bits32 start, end
    double val


cdef extern from "bwgInternal.h":
    struct bwgSectionHead:
        bits32 chromId
        bits32 start, end
        bits32 itemStep
        bits32 itemSpan
        UBYTE type
        UBYTE reserved
        bits16 itemCount


cdef class BWIntervalGenerator:
    cdef:
        bbiFile *file
        char chrom[256]
        bits32 start
        bits32 end

        fileOffsetSize *block
        fileOffsetSize *blockList

        fileOffsetSize *beforeGap
        fileOffsetSize *afterGap
        udcFile *udc
        boolean isSwapped

        char *uncomp_buff

        # contig blocks
        char *mergedBuf
        char *blockBuf

        # block
        char *block_point
        char *blockEnd
        bwgSectionHead head

        # iter
        int yield_count
        bits32 s, e, clippedS, clippedE

    cdef init(self, bbiFile *bigwig, char *chrom, bits32 start, bits32 end)
    cdef _init_contigblocks(self)
    cdef int _init_block(self)

    cdef Interval _iter_block(self)

    cdef _update_block(self)
    cdef _clean_contig(self)
    cdef _clean_up(self)
