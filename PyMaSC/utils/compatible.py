"""Python compatibility utilities.

This module provides compatibility functions and imports to ensure
PyMaSC works consistently across different Python environments.
Currently focused on Python 3 support with legacy compatibility
stubs for previously supported Python 2 functionality.

Key functionality:
- Iterator compatibility imports (zip_longest)

These utilities abstract away version-specific differences to maintain
clean code throughout the PyMaSC codebase.
"""

from itertools import zip_longest
