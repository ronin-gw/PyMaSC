"""Integration tests for pymasc-precalc CLI functionality."""

import subprocess
import tempfile
import pytest
import json
from pathlib import Path
from PyMaSC import VERSION


class TestPymascPrecalcBasics:
    """Test basic pymasc-precalc CLI functionality."""

    def test_pymasc_precalc_version(self):
        """Test that pymasc-precalc --version works correctly."""
        result = subprocess.run(
            ['pymasc-precalc', '--version'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'PyMaSC' in result.stdout
        assert VERSION in result.stdout

    def test_pymasc_precalc_help(self):
        """Test that pymasc-precalc --help works correctly."""
        result = subprocess.run(
            ['pymasc-precalc', '--help'],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert 'usage:' in result.stdout.lower()
        assert 'pymasc-precalc' in result.stdout.lower()


class TestPymascPrecalcFunctionality:
    """Test pymasc-precalc mappability calculation functionality."""

    @pytest.fixture
    def test_data_paths(self):
        """Provide paths to test data."""
        base_dir = Path(__file__).parent.parent.parent
        return {
            'bigwig': base_dir / 'tests' / 'data' / 'hg19_36mer-test.bigwig',
            'golden_json': base_dir / 'tests' / 'data' / 'hg19_36mer-test_mappability.json'
        }

    @pytest.fixture
    def expected_json_data(self, test_data_paths):
        """Load expected JSON data from golden file."""
        with open(test_data_paths['golden_json'], 'r') as f:
            return json.load(f)

    def test_precalc_single_process(self, test_data_paths, expected_json_data):
        """Test pymasc-precalc with single process (default)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_json = Path(temp_dir) / 'test_mappability.json'
            
            cmd = [
                'pymasc-precalc',
                '-m', str(test_data_paths['bigwig']),
                '--mappability-stats', str(output_json),
                '-d', '1000',  # Use default max_shift
                '-r', '36'     # Use read length 36
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            assert result.returncode == 0, f"pymasc-precalc failed: {result.stderr}"
            
            # Check that JSON file was created
            assert output_json.exists(), "Output JSON file not generated"
            
            # Load and compare JSON data
            with open(output_json, 'r') as f:
                actual_data = json.load(f)
            
            # Verify JSON structure and compare overlapping content
            self._compare_json_structure_and_values(expected_json_data, actual_data)

    def test_precalc_multi_process(self, test_data_paths, expected_json_data):
        """Test pymasc-precalc with multi-process mode (-p 2)."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_json = Path(temp_dir) / 'test_mappability_multi.json'
            
            cmd = [
                'pymasc-precalc',
                '-m', str(test_data_paths['bigwig']),
                '--mappability-stats', str(output_json),
                '-p', '2',     # Use 2 processes
                '-d', '1000',  # Use default max_shift
                '-r', '36'     # Use read length 36
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            assert result.returncode == 0, f"pymasc-precalc failed: {result.stderr}"
            
            # Check that JSON file was created
            assert output_json.exists(), "Output JSON file not generated"
            
            # Load and compare JSON data
            with open(output_json, 'r') as f:
                actual_data = json.load(f)
            
            # Verify JSON structure and compare overlapping content
            self._compare_json_structure_and_values(expected_json_data, actual_data)

    def test_precalc_default_output_path(self, test_data_paths):
        """Test pymasc-precalc with default output path."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Copy bigwig to temp dir to test default output naming
            temp_bigwig = Path(temp_dir) / 'test.bigwig'
            import shutil
            shutil.copy(test_data_paths['bigwig'], temp_bigwig)
            
            cmd = [
                'pymasc-precalc',
                '-m', str(temp_bigwig),
                # No --mappability-stats specified, should use default
                '-d', '100',
                '-r', '36'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            assert result.returncode == 0, f"pymasc-precalc failed: {result.stderr}"
            
            # Check that default JSON file was created
            # The default naming is basename + _mappability.json
            default_json = Path(temp_dir) / 'test_mappability.json'
            assert default_json.exists(), f"Default output JSON not generated at {default_json}"
            
            # Verify it's valid JSON
            with open(default_json, 'r') as f:
                data = json.load(f)
                assert '__whole__' in data
                assert 'max_shift' in data
                assert 'references' in data

    def test_precalc_different_parameters(self, test_data_paths):
        """Test pymasc-precalc with different max_shift and max_readlen parameters."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_json = Path(temp_dir) / 'test_custom_params.json'
            
            cmd = [
                'pymasc-precalc',
                '-m', str(test_data_paths['bigwig']),
                '--mappability-stats', str(output_json),
                '-d', '100',  # Different max_shift
                '-r', '50'    # Different read length
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            assert result.returncode == 0, f"pymasc-precalc failed: {result.stderr}"
            
            # Check that JSON file was created
            assert output_json.exists(), "Output JSON file not generated"
            
            # Load and verify structure
            with open(output_json, 'r') as f:
                data = json.load(f)
            
            # Verify basic structure
            assert '__whole__' in data, "Missing __whole__ key in JSON"
            assert 'max_shift' in data, "Missing max_shift key in JSON"
            assert 'references' in data, "Missing references key in JSON"
            
            # The actual max_shift may be limited by the data (e.g., 51 instead of 100)
            # when using different read length (50 instead of 36)
            actual_max_shift = data['max_shift']
            assert actual_max_shift > 0, f"Invalid max_shift: {actual_max_shift}"
            assert actual_max_shift <= 100, f"max_shift should not exceed requested value: {actual_max_shift}"
            
            # Verify array lengths match actual max_shift
            assert len(data['__whole__']) == actual_max_shift + 1, \
                f"Expected {actual_max_shift + 1} values in __whole__, got {len(data['__whole__'])}"

    def test_precalc_log_levels(self, test_data_paths):
        """Test pymasc-precalc with different log levels."""
        with tempfile.TemporaryDirectory() as temp_dir:
            for log_level in ['DEBUG', 'INFO', 'WARNING', 'ERROR']:
                output_json = Path(temp_dir) / f'test_{log_level}.json'
                
                cmd = [
                    'pymasc-precalc',
                    '-m', str(test_data_paths['bigwig']),
                    '--mappability-stats', str(output_json),
                    '-v', log_level,
                    '-d', '50',  # Small value for speed
                    '-r', '36'
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                assert result.returncode == 0, f"pymasc-precalc failed with log level {log_level}: {result.stderr}"
                assert output_json.exists(), f"Output not generated for log level {log_level}"

    def test_precalc_disable_progress(self, test_data_paths):
        """Test pymasc-precalc with progress bar disabled."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_json = Path(temp_dir) / 'test_no_progress.json'
            
            cmd = [
                'pymasc-precalc',
                '-m', str(test_data_paths['bigwig']),
                '--mappability-stats', str(output_json),
                '--disable-progress',
                '-d', '50',  # Small value for speed
                '-r', '36'
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            assert result.returncode == 0, f"pymasc-precalc failed: {result.stderr}"
            assert output_json.exists(), "Output JSON file not generated"

    def _compare_json_structure_and_values(self, expected, actual):
        """Helper to compare JSON data structures, handling different max_shift values."""
        # Check required keys exist
        assert '__whole__' in actual, "Missing __whole__ key in actual JSON"
        assert 'max_shift' in actual, "Missing max_shift key in actual JSON"
        assert 'references' in actual, "Missing references key in actual JSON"
        
        # The actual max_shift may be less than requested due to data limitations
        # For the test data, it caps at 230 even if we request more
        actual_max_shift = actual['max_shift']
        expected_max_shift = expected['max_shift']
        
        # Check that we got a reasonable max_shift value
        assert actual_max_shift > 0, f"Invalid max_shift: {actual_max_shift}"
        assert len(actual['__whole__']) == actual_max_shift + 1, \
            f"__whole__ length ({len(actual['__whole__'])}) doesn't match max_shift + 1 ({actual_max_shift + 1})"
        
        # Compare overlapping values (up to the minimum of the two max_shifts)
        min_shift = min(actual_max_shift, expected_max_shift)
        for i in range(min(min_shift + 1, len(expected['__whole__']), len(actual['__whole__']))):
            assert actual['__whole__'][i] == expected['__whole__'][i], \
                f"__whole__[{i}] mismatch: expected {expected['__whole__'][i]}, got {actual['__whole__'][i]}"
        
        # Check references
        assert set(actual['references'].keys()) == set(expected['references'].keys()), \
            f"References keys mismatch: expected {set(expected['references'].keys())}, got {set(actual['references'].keys())}"
        
        # Compare reference values (up to overlapping range)
        for chrom in expected['references']:
            assert chrom in actual['references'], f"Missing chromosome {chrom} in actual references"
            
            exp_ref = expected['references'][chrom]
            act_ref = actual['references'][chrom]
            
            # Compare overlapping values
            for i in range(min(min_shift + 1, len(exp_ref), len(act_ref))):
                assert act_ref[i] == exp_ref[i], \
                    f"Reference {chrom}[{i}] mismatch: expected {exp_ref[i]}, got {act_ref[i]}"


class TestPymascPrecalcErrorHandling:
    """Test error handling in pymasc-precalc."""

    def test_precalc_missing_input(self):
        """Test pymasc-precalc with missing required input."""
        result = subprocess.run(
            ['pymasc-precalc'],
            capture_output=True,
            text=True
        )
        # Should fail or show usage due to missing required arguments
        assert result.returncode != 0 or 'usage:' in result.stdout.lower()

    def test_precalc_invalid_bigwig(self):
        """Test pymasc-precalc with non-existent BigWig file."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_json = Path(temp_dir) / 'test.json'
            
            result = subprocess.run(
                ['pymasc-precalc',
                 '-m', '/nonexistent/path/file.bigwig',
                 '--mappability-stats', str(output_json)],
                capture_output=True,
                text=True
            )
            assert result.returncode != 0
            # Should show error about missing file
            assert 'error' in result.stderr.lower() or 'not found' in result.stderr.lower()

    def test_precalc_invalid_option(self):
        """Test pymasc-precalc with invalid option."""
        result = subprocess.run(
            ['pymasc-precalc', '--invalid-option'],
            capture_output=True,
            text=True
        )
        assert result.returncode != 0
        # Should show usage or error
        assert len(result.stderr) > 0 or 'usage:' in result.stdout.lower()

    def test_precalc_invalid_process_count(self):
        """Test pymasc-precalc with invalid process count."""
        with tempfile.TemporaryDirectory() as temp_dir:
            result = subprocess.run(
                ['pymasc-precalc',
                 '-m', 'dummy.bigwig',
                 '-p', '0'],  # Invalid process count
                capture_output=True,
                text=True
            )
            # Should fail with invalid process count
            assert result.returncode != 0 or 'error' in result.stderr.lower()