import sys
import os

from PyMaSC.reader.bx.bigwig_file import BigWigFile


class BigWigReader(object):
    def __init__(self, path, bw_readsize=None):
        if not os.path.exists(path) and path != '-':
            raise IOError("input file '{0}' dose not exist.".format(path))

        if path == '-':
            self.file = sys.stdin
        else:
            self.file = open(path)
        self.closed = False

        self.bigwig = BigWigFile(self.file)
        self.chromsizes = self.bigwig.get_chrom_sizes()
        self.bw_readsize = bw_readsize

    def fetch(self, valfilter, chrom=None, begin=None, end=None):
        if chrom is None:
            chroms = self.chromsizes.keys()
        else:
            if chrom not in self.chromsizes:
                raise AttributeError("Reference name '{}' not found.".format(chrom))
            chroms = [chrom]
        if begin is None or begin < 0:
            begin = 0

        for chrom in chroms:
            if end is None:
                end = self.chromsizes[chrom]

            for s, e, val in self.bigwig.get(chrom, begin, end, valfilter):
                yield s, e, val

    def close(self):
        if not self.closed:
            self.file.close()
            self.closed = True

    def __del__(self):
        self.close()
