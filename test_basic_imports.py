#!/usr/bin/env python
"""
Basic import test for PyMaSC across Python versions
"""

import sys
print(f"Python {sys.version}")

try:
    import PyMaSC
    print(f"✓ PyMaSC version: {PyMaSC.VERSION}")
except ImportError as e:
    print(f"✗ Failed to import PyMaSC: {e}")
    sys.exit(1)

try:
    from PyMaSC.core.ncc import NaiveCCCalculator
    print("✓ NaiveCCCalculator imported")
except ImportError as e:
    print(f"✗ Failed to import NaiveCCCalculator: {e}")
    sys.exit(1)

try:
    from PyMaSC.bacore.bitarray import bitarray
    print("✓ bitarray imported")
except ImportError as e:
    print(f"✗ Failed to import bitarray: {e}")
    sys.exit(1)

try:
    import PyMaSC.plot
    print("✓ Plot module imported")
except ImportError as e:
    print(f"✗ Failed to import plot module: {e}")

print("\nAll basic imports successful!")