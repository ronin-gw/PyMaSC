"""Golden test for exact PyMaSC output validation.

This module tests that PyMaSC produces identical numerical results
compared to baseline 'golden' outputs generated from known test data.
"""

import pytest
import os
import tempfile
import shutil
import json
import subprocess
from pathlib import Path
import numpy as np
import csv


class TestGoldenOutputs:
    """Test PyMaSC outputs against golden baselines."""

    @pytest.fixture
    def test_data_paths(self):
        """Provide paths to test data."""
        base_dir = Path(__file__).parent.parent.parent
        return {
            'bam': base_dir / 'tests' / 'data' / 'ENCFF000RMB-test.bam',
            'bigwig': base_dir / 'tests' / 'data' / 'hg19_36mer-test.bigwig',
            'golden_dir': base_dir / 'tests' / 'golden'
        }

    @pytest.fixture
    def golden_files(self, test_data_paths):
        """Provide paths to golden output files."""
        golden_dir = test_data_paths['golden_dir']
        return {
            'stats': golden_dir / 'ENCFF000RMB-test_stats.tab',
            'cc': golden_dir / 'ENCFF000RMB-test_cc.tab',
            'mscc': golden_dir / 'ENCFF000RMB-test_mscc.tab',
            'nreads': golden_dir / 'ENCFF000RMB-test_nreads.tab'
        }

    def test_golden_files_exist(self, golden_files):
        """Test that all golden baseline files exist."""
        for name, path in golden_files.items():
            assert path.exists(), f"Golden file missing: {name} at {path}"

    def test_exact_stats_output(self, test_data_paths, golden_files):
        """Test that stats output matches golden baseline exactly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Run PyMaSC with identical parameters
            cmd = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-o', temp_dir,
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots'  # Skip PDF generation for speed
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            assert result.returncode == 0, f"PyMaSC failed: {result.stderr}"

            # Compare stats output exactly
            output_stats = Path(temp_dir) / 'ENCFF000RMB-test_stats.tab'
            assert output_stats.exists(), "Output stats file not generated"

            # Read both files and compare line by line
            with open(golden_files['stats'], 'r') as golden_f, \
                 open(output_stats, 'r') as output_f:
                golden_lines = golden_f.readlines()
                output_lines = output_f.readlines()

                assert len(golden_lines) == len(output_lines), \
                    f"Different number of lines: golden={len(golden_lines)}, output={len(output_lines)}"

                for i, (golden_line, output_line) in enumerate(zip(golden_lines, output_lines)):
                    assert golden_line == output_line, \
                        f"Line {i+1} differs:\\nGolden: {golden_line.strip()}\\nOutput: {output_line.strip()}"

    def test_exact_cc_output(self, test_data_paths, golden_files):
        """Test that cross-correlation output matches exactly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Run PyMaSC
            cmd = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-o', temp_dir,
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots'
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            assert result.returncode == 0, f"PyMaSC failed: {result.stderr}"

            # Compare CC output
            output_cc = Path(temp_dir) / 'ENCFF000RMB-test_cc.tab'
            assert output_cc.exists(), "Output CC file not generated"

            # Read as CSV data for precise comparison
            def read_tsv_file(filepath):
                data = []
                with open(filepath, 'r') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    headers = reader.fieldnames
                    for row in reader:
                        data.append(row)
                return headers, data

            golden_headers, golden_data = read_tsv_file(golden_files['cc'])
            output_headers, output_data = read_tsv_file(output_cc)

            # Check structure
            assert golden_headers == output_headers, \
                f"Column names differ: {golden_headers} vs {output_headers}"
            assert len(golden_data) == len(output_data), \
                f"Row count differs: {len(golden_data)} vs {len(output_data)}"

            # Check numerical values with exact precision
            for i, (golden_row, output_row) in enumerate(zip(golden_data, output_data)):
                for col in golden_headers:
                    golden_val = golden_row[col]
                    output_val = output_row[col]

                    # Try to convert to float for numerical comparison
                    try:
                        golden_float = float(golden_val)
                        output_float = float(output_val)

                        # Use high precision comparison for floats
                        np.testing.assert_almost_equal(
                            golden_float, output_float,
                            decimal=15,
                            err_msg=f"Row {i+1}, Column {col}: {golden_val} vs {output_val}"
                        )
                    except ValueError:
                        # Non-numeric values should match exactly
                        assert golden_val == output_val, \
                            f"Row {i+1}, Column {col}: {golden_val} vs {output_val}"

    def test_exact_mscc_output(self, test_data_paths, golden_files):
        """Test that MSCC output matches exactly."""
        with tempfile.TemporaryDirectory() as temp_dir:
            # Run PyMaSC
            cmd = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-o', temp_dir,
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots'
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            assert result.returncode == 0, f"PyMaSC failed: {result.stderr}"

            # Compare MSCC output
            output_mscc = Path(temp_dir) / 'ENCFF000RMB-test_mscc.tab'
            assert output_mscc.exists(), "Output MSCC file not generated"

            # Read and compare numerically using CSV reader
            def read_tsv_file(filepath):
                data = []
                with open(filepath, 'r') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    headers = reader.fieldnames
                    for row in reader:
                        data.append(row)
                return headers, data

            golden_headers, golden_data = read_tsv_file(golden_files['mscc'])
            output_headers, output_data = read_tsv_file(output_mscc)

            # Structural checks
            assert golden_headers == output_headers
            assert len(golden_data) == len(output_data)

            # Numerical comparison
            for i, (golden_row, output_row) in enumerate(zip(golden_data, output_data)):
                for col in golden_headers:
                    golden_val = golden_row[col]
                    output_val = output_row[col]

                    try:
                        golden_float = float(golden_val)
                        output_float = float(output_val)

                        np.testing.assert_almost_equal(
                            golden_float, output_float,
                            decimal=15,
                            err_msg=f"MSCC Row {i+1}, Column {col}: {golden_val} vs {output_val}"
                        )
                    except ValueError:
                        assert golden_val == output_val, \
                            f"MSCC Row {i+1}, Column {col}: {golden_val} vs {output_val}"

    def test_command_reproducibility(self, test_data_paths):
        """Test that running the same command twice produces identical results."""
        with tempfile.TemporaryDirectory() as temp_dir1, \
             tempfile.TemporaryDirectory() as temp_dir2:

            cmd_base = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots'
            ]

            # Run twice with different output directories
            cmd1 = cmd_base + ['-o', temp_dir1]
            cmd2 = cmd_base + ['-o', temp_dir2]

            result1 = subprocess.run(cmd1, capture_output=True, text=True, timeout=60)
            result2 = subprocess.run(cmd2, capture_output=True, text=True, timeout=60)

            assert result1.returncode == 0, f"First run failed: {result1.stderr}"
            assert result2.returncode == 0, f"Second run failed: {result2.stderr}"

            # Compare outputs
            files_to_compare = ['ENCFF000RMB-test_stats.tab', 'ENCFF000RMB-test_cc.tab', 'ENCFF000RMB-test_mscc.tab']

            for filename in files_to_compare:
                file1 = Path(temp_dir1) / filename
                file2 = Path(temp_dir2) / filename

                assert file1.exists() and file2.exists(), f"Missing file: {filename}"

                # Read and compare content
                with open(file1, 'r') as f1, open(file2, 'r') as f2:
                    content1 = f1.read()
                    content2 = f2.read()
                    assert content1 == content2, f"File {filename} differs between runs"

    def test_key_statistics_values(self, golden_files):
        """Test that key statistical values are within expected ranges."""
        # Read the golden stats file
        stats_data = {}
        with open(golden_files['stats'], 'r') as f:
            for line in f:
                if '\t' in line:
                    key, value = line.strip().split('\t')
                    try:
                        stats_data[key] = float(value)
                    except ValueError:
                        stats_data[key] = value

        # Validate key statistics
        assert stats_data['Read length'] == 36, "Read length should be 36"
        assert stats_data['Forward reads'] == 622, "Forward reads count unexpected"
        assert stats_data['Reverse reads'] == 670, "Reverse reads count unexpected"
        assert stats_data['Genome length'] == 3137454505, "Genome length unexpected"

        # Check that estimated library length is reasonable for ChIP-seq
        assert 50 <= stats_data['Estimated library length'] <= 300, \
            f"Estimated library length {stats_data['Estimated library length']} outside expected range"

        # Check NSC and RSC values are reasonable (use estimated values)
        assert stats_data['Estimated NSC'] > 1.0, "Estimated NSC should be > 1"
        assert stats_data['Estimated RSC'] > 1.0, "Estimated RSC should be > 1"
        assert stats_data['Estimated MSCC NSC'] > 1.0, "Estimated MSCC NSC should be > 1"
        assert stats_data['Estimated MSCC RSC'] > 1.0, "Estimated MSCC RSC should be > 1"

    def test_cross_correlation_peak_detection(self, golden_files):
        """Test that cross-correlation has reasonable peak structure."""
        # Read CC data using CSV reader
        cc_data = []
        with open(golden_files['cc'], 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                cc_data.append({
                    'shift': int(row['shift']),
                    'whole': float(row['whole']),
                    'chr1': float(row['chr1'])
                })

        # Find peak in whole genome CC
        peak_row = max(cc_data, key=lambda x: x['whole'])
        peak_shift = peak_row['shift']
        peak_value = peak_row['whole']

        # Validate peak properties
        assert 50 <= peak_shift <= 200, f"CC peak at shift {peak_shift} outside expected range"
        assert peak_value > 0.1, f"CC peak value {peak_value} too low"

        # Check that peak is higher than shift=0 (should show fragment length signal)
        cc_at_zero = next(row['whole'] for row in cc_data if row['shift'] == 0)
        assert peak_value > cc_at_zero, "Peak should be higher than zero-shift correlation"

        print(f"Cross-correlation peak detected at shift {peak_shift} with value {peak_value:.6f}")

    def test_parallel_vs_single_process_golden_consistency(self, test_data_paths):
        """Test that parallel processing produces identical results to single-process mode.
        
        This critical test ensures that the MSCC parallel processing bug doesn't regress.
        It validates that parallel mode (-p option) produces exactly the same golden
        results as single-process mode.
        """
        with tempfile.TemporaryDirectory() as single_dir, \
             tempfile.TemporaryDirectory() as parallel_dir:

            # Base command
            cmd_base = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots'
            ]

            # Single-process command
            cmd_single = cmd_base + ['-o', single_dir]
            
            # Parallel-process command
            cmd_parallel = cmd_base + ['-o', parallel_dir, '-p', '2']

            # Run both commands
            result_single = subprocess.run(cmd_single, capture_output=True, text=True, timeout=120)
            result_parallel = subprocess.run(cmd_parallel, capture_output=True, text=True, timeout=120)

            # Both should succeed
            assert result_single.returncode == 0, f"Single-process golden test failed: {result_single.stderr}"
            assert result_parallel.returncode == 0, f"Parallel-process golden test failed: {result_parallel.stderr}"

            # Files to compare
            files_to_compare = [
                'ENCFF000RMB-test_stats.tab',
                'ENCFF000RMB-test_cc.tab', 
                'ENCFF000RMB-test_mscc.tab',
                'ENCFF000RMB-test_nreads.tab'
            ]

            for filename in files_to_compare:
                single_file = Path(single_dir) / filename
                parallel_file = Path(parallel_dir) / filename

                assert single_file.exists(), f"Single-process {filename} not generated"
                assert parallel_file.exists(), f"Parallel-process {filename} not generated"

                # For tabular data, compare with high precision
                if filename.endswith('.tab'):
                    self._compare_tabular_files_exactly(single_file, parallel_file, filename)

        print("✅ Parallel processing produces identical golden results to single-process mode")

    def _compare_tabular_files_exactly(self, file1, file2, filename):
        """Compare two tabular files with exact numerical precision."""
        # Read as CSV for proper parsing
        def read_tab_file(filepath):
            data = []
            with open(filepath, 'r') as f:
                if filename.endswith('_stats.tab'):
                    # Stats file format: key\tvalue
                    for line in f:
                        line = line.strip()
                        if line and '\t' in line:
                            key, value = line.split('\t', 1)
                            data.append({'key': key, 'value': value})
                else:
                    # Table format with headers
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        data.append(row)
            return data

        data1 = read_tab_file(file1)
        data2 = read_tab_file(file2)

        assert len(data1) == len(data2), \
            f"{filename}: Different number of rows: {len(data1)} vs {len(data2)}"

        # Compare each row
        for i, (row1, row2) in enumerate(zip(data1, data2)):
            for key in row1.keys():
                val1 = row1[key]
                val2 = row2[key]

                # Try numerical comparison first
                try:
                    float1 = float(val1)
                    float2 = float(val2)
                    
                    # Use high precision comparison
                    np.testing.assert_almost_equal(
                        float1, float2, 
                        decimal=15,
                        err_msg=f"{filename} Row {i+1}, Column {key}: {val1} vs {val2}"
                    )
                except ValueError:
                    # Non-numeric values must match exactly
                    assert val1 == val2, \
                        f"{filename} Row {i+1}, Column {key}: '{val1}' vs '{val2}'"


class TestGoldenDataIntegrity:
    """Test integrity and properties of golden baseline data."""

    @pytest.fixture
    def golden_files(self):
        """Provide paths to golden files."""
        base_dir = Path(__file__).parent.parent.parent
        golden_dir = base_dir / 'tests' / 'golden'
        return {
            'stats': golden_dir / 'ENCFF000RMB-test_stats.tab',
            'cc': golden_dir / 'ENCFF000RMB-test_cc.tab',
            'mscc': golden_dir / 'ENCFF000RMB-test_mscc.tab'
        }

    def test_golden_stats_format(self, golden_files):
        """Test that golden stats file has correct format."""
        with open(golden_files['stats'], 'r') as f:
            lines = f.readlines()

        assert len(lines) > 30, "Stats file should have multiple entries"

        # Check that each line (except empty ones) has key-value format
        for i, line in enumerate(lines):
            line = line.strip()
            if line:  # Non-empty lines
                assert '\t' in line, f"Line {i+1} missing tab separator: {line}"
                parts = line.split('\t')
                assert len(parts) == 2, f"Line {i+1} should have exactly 2 parts: {line}"

    def test_golden_cc_format(self, golden_files):
        """Test that golden CC file has correct format."""
        # Read CC data using CSV reader
        cc_data = []
        with open(golden_files['cc'], 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            headers = reader.fieldnames
            for row in reader:
                cc_data.append({
                    'shift': int(row['shift']),
                    'whole': float(row['whole']),
                    'chr1': float(row['chr1'])
                })

        # Check required columns
        assert 'shift' in headers, "CC file missing 'shift' column"
        assert 'whole' in headers, "CC file missing 'whole' column"
        assert 'chr1' in headers, "CC file missing 'chr1' column"

        # Check shift values go from 0 to 300
        shifts = [row['shift'] for row in cc_data]
        assert min(shifts) == 0, "CC shifts should start at 0"
        assert max(shifts) == 300, "CC shifts should go to 300"
        assert len(cc_data) == 301, "Should have 301 shift values (0-300)"

        # Check that correlation values are reasonable
        whole_values = [row['whole'] for row in cc_data]
        assert min(whole_values) >= 0, "Correlation values should be >= 0"
        assert max(whole_values) <= 1, "Correlation values should be <= 1"

    def test_golden_mscc_format(self, golden_files):
        """Test that golden MSCC file has correct format."""
        # Read MSCC data using CSV reader
        mscc_data = []
        with open(golden_files['mscc'], 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            headers = reader.fieldnames
            for row in reader:
                mscc_data.append({
                    'shift': int(row['shift']),
                    'whole': float(row['whole']),
                    'chr1': float(row['chr1'])
                })

        # Check structure similar to CC
        assert 'shift' in headers
        assert 'whole' in headers
        assert 'chr1' in headers

        # Check shift range
        shifts = [row['shift'] for row in mscc_data]
        assert min(shifts) == 0
        assert max(shifts) == 300
        assert len(mscc_data) == 301

        # Check correlation values (MSCC can be negative)
        whole_values = [row['whole'] for row in mscc_data]
        assert min(whole_values) >= -1, "MSCC values should be >= -1"
        assert max(whole_values) <= 1, "MSCC values should be <= 1"

    def test_golden_data_consistency(self, golden_files):
        """Test consistency between different golden files."""
        # Read data
        with open(golden_files['stats'], 'r') as f:
            stats_lines = f.readlines()

        # Read CC and MSCC data using CSV reader
        def read_correlation_data(filepath):
            data = []
            with open(filepath, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    data.append({
                        'shift': int(row['shift']),
                        'whole': float(row['whole']),
                        'chr1': float(row['chr1'])
                    })
            return data

        cc_data = read_correlation_data(golden_files['cc'])
        mscc_data = read_correlation_data(golden_files['mscc'])

        # Extract read counts from stats
        forward_reads = None
        reverse_reads = None
        for line in stats_lines:
            if line.startswith('Forward reads\t'):
                forward_reads = int(line.split('\t')[1])
            elif line.startswith('Reverse reads\t'):
                reverse_reads = int(line.split('\t')[1])

        assert forward_reads is not None, "Forward reads count not found in stats"
        assert reverse_reads is not None, "Reverse reads count not found in stats"

        # Verify these counts are reasonable for our test data
        assert forward_reads > 0, "Should have forward reads"
        assert reverse_reads > 0, "Should have reverse reads"
        assert forward_reads + reverse_reads > 1000, "Total reads should be > 1000"

        print(f"Golden data validated: {forward_reads} forward, {reverse_reads} reverse reads")


# Test markers
pytestmark = [
    pytest.mark.integration,
    pytest.mark.golden_test,
    pytest.mark.slow  # These tests run actual PyMaSC commands
]

# Add parallel processing tests to the same markers
pytest.mark.parallel = pytest.mark.integration