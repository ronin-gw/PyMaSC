from bx.bbi.bpt_file cimport BPTFile as bxBPTFile
from bx.bbi.types cimport UBYTE, bits16, bits64



cdef class BPTFile(bxBPTFile):
    def r_find( self, bits64 block_start, key ):
        """
        Recursively seek the value matching key under the subtree starting
        at file offset `block_start`
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
                if node_key == key:
                    return node_value
            return None
        else:
            # Read and discard first key, store offset
            self.reader.read( self.key_size )
            offset = self.reader.read_uint64()
            # Loop until correct subtree is found
            for i from 0 <= i < child_count - 1:
                node_key = self.reader.read( self.key_size )
                if node_key > key:
                    break
                offset = self.reader.read_uint64()
            return self.r_find( offset, key )

    def find( self, key ):
        """
        Find the value matching `key` (a bytestring). Returns the matching
        value as a bytestring if found, or None
        """
        # Key is greater than key_size, must not be a match
        if len(key) > self.key_size:
            return None
        # Key is less than key_size, right pad with 0 bytes
        if len(key) < self.key_size:
            key += ( b'\0' * ( self.key_size - len(key) ) )
        # Call the recursive finder
        return self.r_find( self.root_offset, key )

    def r_traverse(self, bits64 block_start):
        """
        """
        cdef UBYTE is_leaf
        cdef bits16 child_count
        cdef bits64 offset, backpoint
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
            # Traverse all subtree
            for i from 0 <= i < child_count:
                self.reader.skip( self.key_size )
                offset = self.reader.read_uint64()
                if offset == 0:
                    break
                backpoint = self.reader.tell()
                yield from self.r_traverse( offset )
                self.reader.seek( backpoint )

    def traverse(self):
        """
        """
        return self.r_traverse(self.root_offset)
