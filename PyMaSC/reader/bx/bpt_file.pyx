from bx.bbi.bpt_file cimport BPTFile as bxBPTFile
from bx.bbi.types cimport UBYTE, bits16, bits64



cdef class BPTFile(bxBPTFile):
    def r_traverse(self, bits64 block_start):
        """
        """
        cdef UBYTE is_leaf
        cdef bits16 child_count
        cdef bits64 offset
        self.reader.seek( block_start )
        # Block header
        is_leaf = self.reader.read_uint8()
        self.reader.read_uint8()
        child_count = self.reader.read_uint16()
        if is_leaf:
            for i from 0 <= i < child_count:
                node_key = self.reader.read( self.key_size )
                node_value = self.reader.read( self.value_size )
                yield node_key, node_value
        else:
            # Read and discard first key, store offset
            self.reader.read( self.key_size )
            offset = self.reader.read_uint64()
            # Loop until correct subtree is found
            for i from 0 <= i < child_count:
                self.reader.skip( self.key_size )
                offset = self.reader.read_uint64()
            yield self.traverse( offset )

    def traverse(self):
        """
        """
        return self.r_traverse(self.root_offset)