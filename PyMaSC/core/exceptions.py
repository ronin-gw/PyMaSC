"""Exceptions for PyMaSC cross-correlation analysis.

This module centralizes PyMaSC-specific exceptions that were originally
defined in result.py.
"""


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