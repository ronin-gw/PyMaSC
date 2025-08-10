# cython: language_level=3
"""BigWig file reader interface for mappability data.

This module provides a Cython-optimized interface for reading BigWig files
containing mappability data used in MSCC calculations. It wraps the pyBigWig
library to provide efficient access to mappability tracks.

Key components:
- BigWigFileIterator: Iterator for BigWig intervals
- BigWigReader: Main BigWig file reader interface
- Efficient interval data structures for genomic regions
- Compatible API for seamless integration with MSCC algorithms

The module replaces the original bx-python implementation with pyBigWig
for improved performance and maintainability.
"""
import os
import sys
import pyBigWig


# Define interval structure to match original API
cdef struct interval:
    unsigned int begin, end
    float value


cdef class BigWigFileIterator(object):
    """Iterator for BigWig intervals with original API compatibility.

    Provides iteration over BigWig intervals while maintaining compatibility
    with the original BigWigFile API. The iterator yields interval structures
    containing genomic coordinates and associated values.

    Attributes:
        _intervals: List of intervals to iterate over
        _index: Current iteration position
    """

    def __init__(self, intervals):
        """Initialize iterator with interval data.

        Args:
            intervals: List of interval tuples (start, end, value) or None
        """
        self._intervals = intervals if intervals else []
        self._index = 0

    cdef interval next(self):
        """Get next interval in Cython-optimized format.

        Returns:
            interval structure with begin, end, and value fields
            Returns end=0 when iteration is complete
        """
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
        """Return iterator object for Python iteration protocol.

        Returns:
            Self as the iterator object
        """
        """Python iterator protocol - return self."""
        return self

    def __next__(self):
        """Get next interval for Python iteration protocol.

        Returns:
            Tuple of (start, end, value) for the next interval

        Raises:
            StopIteration: When no more intervals are available
        """
        """Python iterator protocol - return next interval as tuple."""
        if self._index >= len(self._intervals):
            raise StopIteration

        entry = self._intervals[self._index]
        self._index += 1
        return (entry[0], entry[1], entry[2])  # (begin, end, value)


cdef class BigWigReader(object):
    """pyBigWig-based BigWig file reader with original API compatibility.

    Provides efficient access to BigWig files containing mappability data
    for MSCC calculations. The class wraps pyBigWig functionality while
    maintaining compatibility with the original bx-python-based API.

    Attributes:
        path: Path to the BigWig file
        closed: Whether the file has been closed
        _bw: pyBigWig file object for data access
        chromsizes: Dictionary mapping chromosome names to lengths
    """

    def __init__(self, path):
        """Initialize BigWig reader.

        Opens the BigWig file and extracts chromosome size information
        for efficient data access during MSCC calculations.

        Args:
            path: Path to BigWig file to read (Union[str, os.PathLike[str]])

        Raises:
            IOError: If the specified file does not exist
            RuntimeError: If the BigWig file cannot be opened
        """
        # Convert Path objects to string for Cython compatibility
        path_str = str(path)

        if not os.path.exists(path_str) and path_str != '-':
            raise IOError("input file '{0}' dose not exist.".format(path_str))

        self.path = path_str
        self.closed = False

        if path_str == '-':
            raise ValueError("pyBigWig does not support stdin input")

        # Open BigWig file with pyBigWig
        self.file = pyBigWig.open(path_str)
        if self.file is None:
            raise IOError("Failed to open BigWig file: {}".format(path_str))

        # Get chromosome sizes
        self.chromsizes = self.file.chroms()

    cpdef BigWigFileIterator fetch(self, float valfilter, chrom):
        """Fetch intervals from BigWig file for specified chromosome.

        Retrieves all intervals from the specified chromosome that meet
        the value filter criteria. Returns an iterator for efficient
        processing of large datasets.

        Args:
            valfilter: Minimum value threshold for intervals
            chrom: Chromosome name to fetch data from

        Returns:
            BigWigFileIterator containing filtered intervals

        Note:
            Only intervals with values >= valfilter are included
        """
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

    def disable_progress_bar(self):
        """Disable progress bar for compatibility with MSCC calculations.

        This method provides compatibility with the original bx-python API
        but has no actual effect since pyBigWig doesn't have progress bars.
        """
        # No-op for pyBigWig compatibility
        pass

    def close(self):
        """Close the BigWig file and release resources.

        Closes the underlying pyBigWig file object and marks the
        reader as closed to prevent further operations.
        """
        """Close the BigWig file."""
        if not self.closed:
            if self.file:
                self.file.close()
            self.closed = True