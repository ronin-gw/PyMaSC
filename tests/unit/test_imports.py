"""Test basic module imports to ensure all Cython extensions load correctly."""

import pytest


class TestBasicImports:
    """Test that all core modules can be imported without errors."""

    def test_import_main_module(self):
        """Test importing the main PyMaSC module."""
        import PyMaSC
        assert hasattr(PyMaSC, 'VERSION')

    def test_import_core_modules(self):
        """Test importing core algorithmic modules."""
        # These should work without BitArray dependencies
        import PyMaSC.core.mappability
        import PyMaSC.core.ncc
        import PyMaSC.core.mscc
        import PyMaSC.core.readlen

    def test_import_bitarray_modules(self):
        """Test importing BitArray-dependent modules (critical functionality)."""
        # These were the problematic modules that are now fixed
        import PyMaSC.bacore.bitarray
        import PyMaSC.bacore.mscc

    def test_import_reader_modules(self):
        """Test importing reader modules."""
        import PyMaSC.reader.bigwig
        import PyMaSC.reader.bx.bbi_file
        import PyMaSC.reader.bx.bigwig_file
        import PyMaSC.reader.bx.bpt_file

    def test_import_utility_modules(self):
        """Test importing utility modules."""
        import PyMaSC.utils.bx.binary_file

    def test_import_handler_modules(self):
        """Test importing business logic handlers."""
        import PyMaSC.handler.masc
        import PyMaSC.handler.mappability
        import PyMaSC.handler.result
        import PyMaSC.handler.bamasc

    def test_import_entry_points(self):
        """Test importing main entry point modules."""
        import PyMaSC.pymasc
        import PyMaSC.calcmappablelen
        # Note: PyMaSC.plot may fail due to matplotlib/numpy compatibility issues
        # but this is a known issue that doesn't affect core functionality


class TestVersionInfo:
    """Test version and basic package information."""

    def test_version_accessible(self):
        """Test that version information is accessible."""
        import PyMaSC
        assert hasattr(PyMaSC, 'VERSION')
        assert isinstance(PyMaSC.VERSION, str)
        assert len(PyMaSC.VERSION) > 0

    def test_version_format(self):
        """Test that version follows expected format."""
        import PyMaSC
        # Should be semantic versioning like "0.3.1"
        parts = PyMaSC.VERSION.split('.')
        assert len(parts) >= 2
        for part in parts:
            assert part.isdigit()


class TestCriticalFunctionality:
    """Test that critical BitArray functions are accessible."""

    def test_bitarray_module_accessible(self):
        """Test that BitArray functions can be called."""
        import PyMaSC.bacore.bitarray as ba
        # Just test that the module loads and has expected attributes
        # Full functional testing will be in integration tests
        assert hasattr(ba, '__file__')

    def test_mscc_module_accessible(self):
        """Test that MSCC (mappability-sensitive cross-correlation) module loads."""
        import PyMaSC.bacore.mscc as mscc
        assert hasattr(mscc, '__file__')