#!/usr/bin/env python
"""Simple test runner script for PyMaSC."""

import sys
import subprocess


def run_tests():
    """Run PyMaSC test suite."""
    print("Running PyMaSC Test Suite")
    print("=" * 50)
    
    # Run pytest with coverage
    cmd = [
        sys.executable, "-m", "pytest", 
        "tests/", 
        "--verbose",
        "--tb=short",
        "--cov=PyMaSC",
        "--cov-report=term-missing",
        "--cov-report=html"
    ]
    
    try:
        result = subprocess.run(cmd, check=True)
        print("\n✅ All tests passed!")
        return 0
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Tests failed with exit code {e.returncode}")
        return e.returncode
    except KeyboardInterrupt:
        print("\n⚠️  Tests interrupted by user")
        return 1


if __name__ == "__main__":
    sys.exit(run_tests())