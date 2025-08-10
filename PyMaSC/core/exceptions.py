"""Exceptions for PyMaSC cross-correlation analysis.

This module centralizes PyMaSC-specific exceptions that were originally
defined in result.py.
"""


class ReadUnsortedError(IndexError):
    """Exception raised when input reads are not sorted by position.

    This exception is raised when the algorithm detects that input reads
    are not properly sorted by genomic position, which is required for
    the sliding window cross-correlation calculation.
    """
    pass


class NothingToCalc(Exception):
    """Exception raised when no chromosomes are available for calculation."""
    pass


class InputUnseekable(Exception):
    """Exception raised when input stream cannot be rewound."""
    pass


class ReadsTooFew(IndexError):
    """Exception raised when there are insufficient reads for analysis.

    This exception indicates that the number of reads in the input data
    is too low to perform reliable cross-correlation analysis. It typically
    occurs when the BAM file contains very few mapped reads or when
    strict filtering criteria remove most reads.

    Inherits from IndexError to maintain compatibility with existing
    exception handling code.
    """
    pass
