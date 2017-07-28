import os

import wWigIO


class BigWigFile(object):
    @staticmethod
    def wigToBigWig(wigfile, sizefile, bwfile):
        wWigIO.wigToBigWig(wigfile, sizefile, bwfile)

    @staticmethod
    def bigWigToWig(bwfile, wigfile):
        wWigIO.bigWigToWig(bwfile, wigfile)

    def __init__(self, path, chrom_size=None):
        if not os.path.exists(path) and path != '-':
            raise IOError("input file '{0}' dose not exist.".format(path))
        elif path == '-':
            path = "stdin"

        prefix, ext = os.path.splitext(path)[0]
        if ext == "wig":
            bwfile = prefix + ".bw"
            if os.path.exists(bwfile):
                self.path = bwfile
            else:
                if chrom_size is None:
                    raise IOError("Failed to convet wig to bigwig. 'chrom_size' file required.")
                BigWigFile.wigToBigWig(path, chrom_size, bwfile)
            self.path = bwfile
        else:
            self.path = path

        wWigIO.open(self.path)
        self.set_chromsizes()

        self.closed = False

    def set_chromsizes(self):
        self.chromsizes = wWigIO.getChromSize(self.path)

    def _getIntervals(self, chrom, begin, end):
        wigs = wWigIO.getIntervals(self.path, chrom, begin, end)

        if wigs == 1:
            raise ValueError("wWigIO.getIntervals doesn't have correct parameters.")
        if wigs == 2:
            raise ValueError("Fail to open BigWig file.")

        return wigs

    def fetch(self, chrom=None, begin=None, end=None):
        if chrom is None:
            chroms = self.chromsizes.keys()
        else:
            chroms = [chrom]
        if begin is None or begin < 0:
            begin = 0
        if end is None:
            end = 0

        for chrom in chroms:
            for wig in self._getIntervals(chrom, begin, end):
                yield chrom, wig[0], wig[1], wig[2]

    def close(self):
        if not self.closed:
            wWigIO.close(self.infile)
            self.closed = True

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()

    def __del__(self):
        self.close()
