"""Custom exceptions for PyMaSC analysis."""


class ReadUnsortedError(IndexError):
    """Raised when input reads are not sorted by genomic position."""
    pass


class NothingToCalc(Exception):
    """Exception raised when no chromosomes are available for calculation."""
    pass


class InputUnseekable(Exception):
    """Exception raised when input stream cannot be rewound."""
    pass


class ReadsTooFew(IndexError):
    """Raised when input data contains insufficient reads for analysis."""
    pass
