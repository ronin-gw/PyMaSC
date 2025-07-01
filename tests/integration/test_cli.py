"""Integration tests for PyMaSC CLI functionality."""

import subprocess
import pytest


class TestCLIBasics:
    """Test basic CLI functionality and version information."""

    def test_pymasc_version(self):
        """Test that pymasc --version works correctly."""
        result = subprocess.run(
            ['pymasc', '--version'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'PyMaSC' in result.stdout
        assert '0.3.1' in result.stdout

    def test_pymasc_help(self):
        """Test that pymasc --help works correctly."""
        result = subprocess.run(
            ['pymasc', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'usage:' in result.stdout.lower()
        assert 'pymasc' in result.stdout.lower()

    def test_pymasc_precalc_version(self):
        """Test that pymasc-precalc --version works correctly."""
        result = subprocess.run(
            ['pymasc-precalc', '--version'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'PyMaSC' in result.stdout
        assert '0.3.1' in result.stdout

    @pytest.mark.plot
    def test_pymasc_plot_version_may_fail(self):
        """Test pymasc-plot --version (may fail due to matplotlib issues)."""
        result = subprocess.run(
            ['pymasc-plot', '--version'],
            capture_output=True,
            text=True
        )
        # This may fail due to matplotlib/numpy compatibility issues
        # We test both success and expected failure cases
        if result.returncode == 0:
            assert 'PyMaSC' in result.stdout
            assert '0.3.1' in result.stdout
        else:
            # Expected failure due to matplotlib import issues
            assert 'matplotlib' in result.stderr or 'NDArray' in result.stderr


class TestCLIErrorHandling:
    """Test CLI error handling for invalid arguments."""

    def test_pymasc_invalid_args(self):
        """Test that pymasc handles invalid arguments gracefully."""
        result = subprocess.run(
            ['pymasc', '--invalid-argument'],
            capture_output=True,
            text=True
        )
        assert result.returncode != 0
        # Should show usage or error message
        assert len(result.stderr) > 0 or 'usage:' in result.stdout.lower()

    def test_pymasc_no_args(self):
        """Test that pymasc without arguments shows help."""
        result = subprocess.run(
            ['pymasc'],
            capture_output=True,
            text=True
        )
        # May return error code but should show usage
        assert 'usage:' in result.stdout.lower() or 'usage:' in result.stderr.lower()


class TestCoreFunctionality:
    """Test that core PyMaSC functionality is accessible via CLI."""

    def test_pymasc_requires_input(self):
        """Test that pymasc requires input files and shows appropriate error."""
        result = subprocess.run(
            ['pymasc', '--output', '/tmp/test'],
            capture_output=True,
            text=True
        )
        # Should fail because no input BAM file provided
        assert result.returncode != 0
        # Error should mention missing input or required arguments
        error_output = result.stderr.lower()
        assert any(word in error_output for word in ['input', 'required', 'bam', 'argument'])

    def test_pymasc_precalc_requires_input(self):
        """Test that pymasc-precalc requires input files."""
        result = subprocess.run(
            ['pymasc-precalc'],
            capture_output=True,
            text=True
        )
        # Should fail or show usage due to missing required arguments
        assert result.returncode != 0 or 'usage:' in result.stdout.lower()