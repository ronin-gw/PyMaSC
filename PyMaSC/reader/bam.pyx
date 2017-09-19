from pysam.libcalignmentfile cimport AlignmentFile
from pysam.libcalignedsegment cimport AlignedSegment


class BamReadGenerator(AlignmentFile):
    def next(self):
        cdef AlignedSegment read

        while True:
            read = super(BamReadGenerator, self).next()
            # https://github.com/pysam-developers/pysam/issues/268
            try:
                return (read.is_reverse, read.is_duplicate, read.reference_name,
                        read.reference_start + 1, read.mapping_quality, read.query_length)
            except ValueError:
                pass
