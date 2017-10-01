cimport numpy
from bx.bbi.types cimport bits32, bits16, UBYTE

from PyMaSC.utils.bx.binary_file cimport BinaryFileReader

from .bbi_file cimport BlockHandler, BBIFile, interval

DEF big_wig_sig = 0x888FFC26
DEF bwg_bed_graph = 1
DEF bwg_variable_step = 2
DEF bwg_fixed_step = 3

import numpy
from PyMaSC.utils.compatible import StringIO


cdef class BigWigBlockGenerationHandler( BlockHandler ):
    """
    BlockHandler that parses the block into a series of wiggle records, and calls `handle_interval_value` for each.
    """
    cdef:
        bits32 start
        bits32 end

        BinaryFileReader block_reader
        bits32 _b_start
        bits32 _b_item_span
        bits16 _b_item_count
        UBYTE _b_type
        bits32 s, e

        float valfilter
        bits16 _current_item_count

    def __init__( self, bits32 start, bits32 end, float valfilter ):
        BlockHandler.__init__( self )
        self.start = start
        self.end = end
        self.valfilter = valfilter

    cdef set_handle_block( self, bytes block_data, BBIFile bbi_file ):
        cdef bits32 b_chrom_id, b_end, b_valid_count
        cdef bits32 b_item_step
        # Now we parse the block, first the header
        self.block_reader = block_reader = BinaryFileReader(
            StringIO( block_data ), is_little_endian=bbi_file.reader.is_little_endian
        )
        b_chrom_id = block_reader.read_uint32()
        self._b_start = block_reader.read_uint32()
        b_end = block_reader.read_uint32()
        b_item_step = block_reader.read_uint32()
        self._b_item_span = block_reader.read_uint32()
        self._b_type = block_reader.read_uint8()
        block_reader.skip(1)
        self._b_item_count = block_reader.read_uint16()
        self._current_item_count = 0

    cdef interval next(self):
        cdef float val
        cdef bits16 i

        while True:
            if self._current_item_count == self._b_item_count:
                return interval(0, 0, 0)

            # Depending on the type, s and e are either read or
            # generate using header, val is always read
            if self._b_type == bwg_bed_graph:
                self.s = self.block_reader.read_uint32()
                self.e = self.block_reader.read_uint32()
                val = self.block_reader.read_float()
            elif self._b_type == bwg_variable_step:
                self.s = self.block_reader.read_uint32()
                self.e = self.s + self._b_item_span
                val = self.block_reader.read_float()
            elif self._b_type == bwg_fixed_step:
                self.s = self._b_start + ( self._current_item_count * self._b_item_span )
                self.e = self.s + self._b_item_span
                val = self.block_reader.read_float()

            self._current_item_count += 1

            if self.s < self.start:
                self.s = self.start
            if self.e > self.end:
                self.e = self.end
            if self.s < self.e and self.valfilter <= val:
                break

        return interval(self.s, self.e, val)


cdef class BigWigFile( BBIFile ):
    """
    A "big binary indexed" file whose raw data is in wiggle format.
    """
    def __init__( self, file=None ):
        BBIFile.__init__( self, file, big_wig_sig, "bigwig" )

    cpdef get_chrom_sizes( self ):
        return {chrom: size for chrom, _, size in self._yield_all_chrom_key_id_and_size()}

    cpdef get( self, str chrom, bits32 start, bits32 end, float valfilter ):
        """
        Gets all data points over the regions `chrom`:`start`-`end`.
        """
        if start >= end:
            return None
        chrom_id, chrom_size = self._get_chrom_id_and_size( chrom.encode("utf-8") )
        if chrom_id is None:
            return None

        gen = BigWigBlockGenerationHandler( start, end, valfilter )
        self.set_visitor_blocks_in_region( chrom_id, start, end, gen )
        return self
