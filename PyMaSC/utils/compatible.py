"""Python compatibility utilities.

This module provides compatibility functions and imports to ensure
PyMaSC works consistently across different Python environments.
Currently focused on Python 3 support with legacy compatibility
stubs for previously supported Python 2 functionality.

Key functionality:
- String conversion utilities (tostr)
- Iterator compatibility imports (zip_longest, xrange)
- I/O compatibility (StringIO)

These utilities abstract away version-specific differences to maintain
clean code throughout the PyMaSC codebase.
"""


def tostr(s):
    """Convert input to string (Python 3 compatibility stub).
    
    In Python 3, this function simply returns the input unchanged.
    Previously provided Python 2/3 compatibility for string handling.
    
    Args:
        s: Input object to convert to string
        
    Returns:
        The input object unchanged (Python 3 strings are already Unicode)
        
    Note:
        This function is maintained for backward compatibility but
        is essentially a pass-through in Python 3 environments
    """
    return s

from itertools import zip_longest
from io import BytesIO as StringIO
xrange = range
