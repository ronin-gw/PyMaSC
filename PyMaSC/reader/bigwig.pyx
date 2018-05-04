from bx.bbi.types cimport bits32
from PyMaSC.reader.bx.bbi_file cimport interval
from PyMaSC.reader.bx.bigwig_file cimport BigWigFile

import os
import sys


cdef class BigWigReader(object):
    def __init__(self, str path):
        if not os.path.exists(path) and path != '-':
            raise IOError("input file '{0}' dose not exist.".format(path))

        self.path = path
        if path == '-':
            self.file = sys.stdin
        else:
            self.file = open(path, "rb")
        self.closed = False

        self.bigwig = BigWigFile(self.file)
        self.chromsizes = self.bigwig.get_chrom_sizes()

    cpdef BigWigFile fetch(self, float valfilter, chrom):
        cdef bits32 begin, end

        if chrom not in self.chromsizes:
            raise KeyError(chrom)

        begin = 0
        end = self.chromsizes[chrom]

        return self.bigwig.get(chrom, begin, end, valfilter)

    def close(self):
        if not self.closed:
            self.file.close()
            self.closed = True
