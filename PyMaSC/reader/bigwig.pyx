# cython: language_level=3
import os
import sys
import pyBigWig


# Define interval structure to match original API
cdef struct interval:
    unsigned int begin, end
    float value


cdef class BigWigFileIterator(object):
    """Iterator wrapper to match the original BigWigFile API."""

    def __init__(self, intervals):
        self._intervals = intervals if intervals else []
        self._index = 0

    cdef interval next(self):
        cdef interval result

        if self._index >= len(self._intervals):
            # Signal end of iteration with end=0
            result.begin = 0
            result.end = 0
            result.value = 0.0
            return result

        entry = self._intervals[self._index]
        result.begin = entry[0]  # start
        result.end = entry[1]    # end
        result.value = entry[2]  # value
        self._index += 1

        return result

    def __iter__(self):
        """Python iterator protocol - return self."""
        return self

    def __next__(self):
        """Python iterator protocol - return next interval as tuple."""
        if self._index >= len(self._intervals):
            raise StopIteration

        entry = self._intervals[self._index]
        self._index += 1
        return (entry[0], entry[1], entry[2])  # (begin, end, value)


cdef class BigWigReader(object):
    """pyBigWig-based BigWig reader with API compatibility."""

    def __init__(self, str path):
        if not os.path.exists(path) and path != '-':
            raise IOError("input file '{0}' dose not exist.".format(path))

        self.path = path
        self.closed = False

        if path == '-':
            raise ValueError("pyBigWig does not support stdin input")

        # Open BigWig file with pyBigWig
        self.file = pyBigWig.open(path)
        if self.file is None:
            raise IOError("Failed to open BigWig file: {}".format(path))

        # Get chromosome sizes
        self.chromsizes = self.file.chroms()

    cpdef BigWigFileIterator fetch(self, float valfilter, chrom):
        """Fetch intervals for chromosome with value filter."""
        if chrom not in self.chromsizes:
            raise KeyError(chrom)

        chrom_size = self.chromsizes[chrom]

        # Get intervals from pyBigWig - returns list of (start, end, value) tuples
        intervals = self.file.intervals(chrom, 0, chrom_size)

        # Filter by value threshold if needed
        if intervals and valfilter > 0:
            intervals = [interval for interval in intervals if interval[2] >= valfilter]

        return BigWigFileIterator(intervals)

    def close(self):
        """Close the BigWig file."""
        if not self.closed:
            if self.file:
                self.file.close()
            self.closed = True