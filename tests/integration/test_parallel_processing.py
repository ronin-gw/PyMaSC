"""Parallel processing tests to prevent MSCC bugs.

This module ensures that parallel processing (-p option) produces identical
results to single-process mode, particularly for MSCC calculations which
were previously failing silently in parallel mode.
"""

import pytest
import tempfile
import subprocess
from pathlib import Path
import csv
import numpy as np


@pytest.mark.xdist_group(name="process_consistency_tests")
class TestParallelProcessingConsistency:
    """Test that parallel processing produces identical results to single-process mode."""

    @pytest.fixture
    def test_data_paths(self):
        """Provide paths to test data."""
        base_dir = Path(__file__).parent.parent.parent
        return {
            'bam': base_dir / 'tests' / 'data' / 'ENCFF000RMB-test.bam',
            'bigwig': base_dir / 'tests' / 'data' / 'hg19_36mer-test.bigwig'
        }

    def test_single_vs_parallel_stats_identical(self, test_data_paths):
        """Test that single-process and parallel mode produce identical stats."""
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

            # Parallel-process command (2 workers)
            cmd_parallel = cmd_base + ['-o', parallel_dir, '-p', '2']

            # Run both commands
            result_single = subprocess.run(cmd_single, capture_output=True, text=True, timeout=120)
            result_parallel = subprocess.run(cmd_parallel, capture_output=True, text=True, timeout=120)

            # Both should succeed
            assert result_single.returncode == 0, f"Single-process failed: {result_single.stderr}"
            assert result_parallel.returncode == 0, f"Parallel-process failed: {result_parallel.stderr}"

            # Compare stats files
            single_stats = Path(single_dir) / 'ENCFF000RMB-test_stats.tab'
            parallel_stats = Path(parallel_dir) / 'ENCFF000RMB-test_stats.tab'

            assert single_stats.exists(), "Single-process stats file missing"
            assert parallel_stats.exists(), "Parallel-process stats file missing"

            # Read and compare stats line by line
            with open(single_stats, 'r') as f_single, \
                 open(parallel_stats, 'r') as f_parallel:

                single_lines = f_single.readlines()
                parallel_lines = f_parallel.readlines()

                assert len(single_lines) == len(parallel_lines), \
                    f"Different number of lines: single={len(single_lines)}, parallel={len(parallel_lines)}"

                for i, (single_line, parallel_line) in enumerate(zip(single_lines, parallel_lines)):
                    assert single_line == parallel_line, \
                        f"Stats line {i+1} differs:\\nSingle: {single_line.strip()}\\nParallel: {parallel_line.strip()}"

    def test_mscc_parallel_processing_critical(self, test_data_paths):
        """Critical test: Ensure MSCC results are generated in parallel mode.

        This test specifically addresses the bug where MSCC metrics were 'nan'
        in parallel mode, preventing proper mappability-sensitive analysis.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            cmd = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-o', temp_dir,
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots',
                '-p', '2'  # Parallel mode
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            assert result.returncode == 0, f"Parallel MSCC test failed: {result.stderr}"

            # Read stats and verify MSCC metrics are NOT nan
            stats_file = Path(temp_dir) / 'ENCFF000RMB-test_stats.tab'
            assert stats_file.exists(), "Stats file not generated in parallel mode"

            stats_data = {}
            with open(stats_file, 'r') as f:
                for line in f:
                    if '\t' in line:
                        key, value = line.strip().split('\t')
                        stats_data[key] = value

            # Critical MSCC metrics that were failing before the fix
            critical_mscc_metrics = [
                'Minimum MSCC',
                'MSCC at read length',
                'MSCC at estimated library length',
                'Estimated MSCC NSC',
                'Estimated MSCC RSC',
                'Estimated MSCC FWHM',
                'Estimated MSCC VSN'
            ]

            for metric in critical_mscc_metrics:
                assert metric in stats_data, f"Missing MSCC metric: {metric}"
                value = stats_data[metric]

                # The critical test: MSCC values must NOT be 'nan'
                assert value != 'nan', f"CRITICAL BUG DETECTED: {metric} is 'nan' in parallel mode"

                # For numeric metrics, verify they are valid numbers
                if metric in ['Minimum MSCC', 'MSCC at read length', 'Estimated MSCC NSC', 'Estimated MSCC RSC']:
                    try:
                        float_value = float(value)
                        assert not np.isnan(float_value), f"{metric} is NaN: {value}"
                        assert float_value > 0, f"{metric} should be positive: {value}"
                    except ValueError:
                        pytest.fail(f"{metric} is not a valid number: {value}")

            print("✅ MSCC parallel processing working correctly")

    def test_single_vs_parallel_cc_tables_identical(self, test_data_paths):
        """Test that CC tables are identical between single and parallel modes."""
        with tempfile.TemporaryDirectory() as single_dir, \
             tempfile.TemporaryDirectory() as parallel_dir:

            # Run both modes
            cmd_base = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots'
            ]

            # Single and parallel commands
            subprocess.run(cmd_base + ['-o', single_dir], check=True, timeout=120)
            subprocess.run(cmd_base + ['-o', parallel_dir, '-p', '2'], check=True, timeout=120)

            # Compare CC tables
            single_cc = Path(single_dir) / 'ENCFF000RMB-test_cc.tab'
            parallel_cc = Path(parallel_dir) / 'ENCFF000RMB-test_cc.tab'

            # Read and compare using CSV reader for numerical precision
            def read_cc_data(filepath):
                data = []
                with open(filepath, 'r') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        data.append(row)
                return data

            single_data = read_cc_data(single_cc)
            parallel_data = read_cc_data(parallel_cc)

            assert len(single_data) == len(parallel_data), \
                f"CC table length differs: single={len(single_data)}, parallel={len(parallel_data)}"

            # Compare each row numerically
            for i, (single_row, parallel_row) in enumerate(zip(single_data, parallel_data)):
                for col in single_row.keys():
                    single_val = single_row[col]
                    parallel_val = parallel_row[col]

                    # Compare numerically if possible
                    try:
                        single_float = float(single_val)
                        parallel_float = float(parallel_val)

                        np.testing.assert_almost_equal(
                            single_float, parallel_float,
                            decimal=15,
                            err_msg=f"CC Row {i+1}, Column {col}: single={single_val} vs parallel={parallel_val}"
                        )
                    except ValueError:
                        # Non-numeric values should match exactly
                        assert single_val == parallel_val, \
                            f"CC Row {i+1}, Column {col}: single={single_val} vs parallel={parallel_val}"

    def test_single_vs_parallel_mscc_tables_identical(self, test_data_paths):
        """Test that MSCC tables are identical between single and parallel modes.

        This is the most critical test for the bug that was fixed.
        """
        with tempfile.TemporaryDirectory() as single_dir, \
             tempfile.TemporaryDirectory() as parallel_dir:

            # Run both modes
            cmd_base = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots'
            ]

            # Single and parallel commands
            subprocess.run(cmd_base + ['-o', single_dir], check=True, timeout=120)
            subprocess.run(cmd_base + ['-o', parallel_dir, '-p', '2'], check=True, timeout=120)

            # Compare MSCC tables
            single_mscc = Path(single_dir) / 'ENCFF000RMB-test_mscc.tab'
            parallel_mscc = Path(parallel_dir) / 'ENCFF000RMB-test_mscc.tab'

            assert single_mscc.exists(), "Single-process MSCC file missing"
            assert parallel_mscc.exists(), "Parallel-process MSCC file missing"

            # Read and compare MSCC data
            def read_mscc_data(filepath):
                data = []
                with open(filepath, 'r') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        data.append(row)
                return data

            single_data = read_mscc_data(single_mscc)
            parallel_data = read_mscc_data(parallel_mscc)

            assert len(single_data) == len(parallel_data), \
                f"MSCC table length differs: single={len(single_data)}, parallel={len(parallel_data)}"

            # Compare each row with high numerical precision
            for i, (single_row, parallel_row) in enumerate(zip(single_data, parallel_data)):
                for col in single_row.keys():
                    single_val = single_row[col]
                    parallel_val = parallel_row[col]

                    # Compare numerically
                    try:
                        single_float = float(single_val)
                        parallel_float = float(parallel_val)

                        np.testing.assert_almost_equal(
                            single_float, parallel_float,
                            decimal=15,
                            err_msg=f"MSCC Row {i+1}, Column {col}: single={single_val} vs parallel={parallel_val}"
                        )
                    except ValueError:
                        # Non-numeric values should match exactly
                        assert single_val == parallel_val, \
                            f"MSCC Row {i+1}, Column {col}: single={single_val} vs parallel={parallel_val}"

    def test_parallel_worker_counts(self, test_data_paths):
        """Test different parallel worker counts produce identical results."""
        results = {}

        for nworkers in [1, 2, 4]:
            with tempfile.TemporaryDirectory() as temp_dir:
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

                if nworkers > 1:
                    cmd.extend(['-p', str(nworkers)])

                result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                assert result.returncode == 0, f"Failed with {nworkers} workers: {result.stderr}"

                # Read stats
                stats_file = Path(temp_dir) / 'ENCFF000RMB-test_stats.tab'
                with open(stats_file, 'r') as f:
                    content = f.read()
                    results[nworkers] = content

        # All results should be identical
        base_result = results[1]  # Single-process as baseline
        for nworkers, result_content in results.items():
            assert result_content == base_result, \
                f"Results with {nworkers} workers differ from single-process baseline"

        print(f"✅ Parallel processing consistent across 1, 2, and 4 workers")

    def test_mscc_regression_test(self, test_data_paths):
        """Regression test for the specific MSCC parallel processing bug.

        This test validates the exact fix that was applied and ensures
        the bug doesn't reoccur.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            # Run with parallel processing
            cmd = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-o', temp_dir,
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots',
                '-p', '2'
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            assert result.returncode == 0, f"Regression test failed: {result.stderr}"

            # Validate specific metrics that were failing
            stats_file = Path(temp_dir) / 'ENCFF000RMB-test_stats.tab'

            with open(stats_file, 'r') as f:
                content = f.read()

                # Check that 'nan' doesn't appear in MSCC-related lines that should have values
                # Note: Some MSCC values are legitimately 'nan' when expected library length is unknown
                lines = content.split('\n')

                # Metrics that should NOT be nan in parallel mode (these were the actual bug)
                critical_non_nan_metrics = [
                    'Minimum MSCC',
                    'MSCC at read length',
                    'Estimated MSCC NSC',
                    'Estimated MSCC RSC',
                    'Estimated MSCC FWHM',
                    'Estimated MSCC VSN'
                ]

                for line in lines:
                    if '\t' in line:
                        key, value = line.split('\t', 1)
                        if key in critical_non_nan_metrics:
                            assert value != 'nan', f"REGRESSION: {key} is 'nan' in parallel mode"

            # Check MSCC table exists and has data
            mscc_file = Path(temp_dir) / 'ENCFF000RMB-test_mscc.tab'
            assert mscc_file.exists(), "REGRESSION: MSCC table not generated in parallel mode"

            with open(mscc_file, 'r') as f:
                lines = f.readlines()
                assert len(lines) > 300, "REGRESSION: MSCC table has insufficient data"

            print("✅ MSCC regression test passed - bug fix is working")

    @pytest.mark.successive
    @pytest.mark.parallel
    @pytest.mark.critical
    def test_successive_algorithm_single_vs_parallel_consistency(self, test_data_paths):
        """Test successive algorithm consistency between single and parallel modes.

        This test specifically addresses the successive algorithm parallel processing
        bug where MSCCCalculator receives wrong BigWig reader type in parallel mode.
        """
        with tempfile.TemporaryDirectory() as single_dir, \
             tempfile.TemporaryDirectory() as parallel_dir:

            # Base command for successive algorithm
            cmd_base = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots',
                '--successive'  # Critical: use successive algorithm
            ]

            # Single-process command
            cmd_single = cmd_base + ['-o', single_dir]

            # Parallel-process command
            cmd_parallel = cmd_base + ['-o', parallel_dir, '-p', '2']

            # Run both commands
            result_single = subprocess.run(cmd_single, capture_output=True, text=True, timeout=120)
            result_parallel = subprocess.run(cmd_parallel, capture_output=True, text=True, timeout=120)

            # Both should succeed
            assert result_single.returncode == 0, f"Successive single-process failed: {result_single.stderr}"
            assert result_parallel.returncode == 0, f"Successive parallel-process failed: {result_parallel.stderr}"

            # Compare stats files for exact consistency
            single_stats = Path(single_dir) / 'ENCFF000RMB-test_stats.tab'
            parallel_stats = Path(parallel_dir) / 'ENCFF000RMB-test_stats.tab'

            assert single_stats.exists(), "Successive single-process stats file missing"
            assert parallel_stats.exists(), "Successive parallel-process stats file missing"

            # Read and compare stats line by line
            with open(single_stats, 'r') as f_single, \
                 open(parallel_stats, 'r') as f_parallel:

                single_lines = f_single.readlines()
                parallel_lines = f_parallel.readlines()

                assert len(single_lines) == len(parallel_lines), \
                    f"Successive stats lines differ: single={len(single_lines)}, parallel={len(parallel_lines)}"

                for i, (single_line, parallel_line) in enumerate(zip(single_lines, parallel_lines)):
                    assert single_line == parallel_line, \
                        f"Successive stats line {i+1} differs:\nSingle: {single_line.strip()}\nParallel: {parallel_line.strip()}"

            print("✅ Successive algorithm parallel processing working correctly")

    @pytest.mark.successive
    @pytest.mark.parallel
    @pytest.mark.critical
    def test_successive_algorithm_parallel_feed_method_fix(self, test_data_paths):
        """Critical regression test for successive algorithm parallel processing.

        This test ensures that the BWFeederWithMappableRegionSum fix prevents
        the AttributeError: 'BigWigReader' object has no attribute 'feed' error
        that occurred when using successive algorithm in parallel mode.
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            cmd = [
                'pymasc',
                str(test_data_paths['bam']),
                '-m', str(test_data_paths['bigwig']),
                '-o', temp_dir,
                '-d', '300',
                '-q', '10',
                '-r', '36',
                '--skip-plots',
                '--successive',  # Critical: use successive algorithm
                '-p', '2'  # Critical: use parallel processing
            ]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            # The critical test: this should NOT fail with 'feed' method error
            assert result.returncode == 0, f"Successive parallel processing failed: {result.stderr}"

            # Should not contain the specific error message
            error_output = result.stderr
            assert "'BigWigReader' object has no attribute 'feed'" not in error_output, \
                "REGRESSION: BWFeederWithMappableRegionSum fix not working"

            # Verify output files were generated
            stats_file = Path(temp_dir) / 'ENCFF000RMB-test_stats.tab'
            cc_file = Path(temp_dir) / 'ENCFF000RMB-test_cc.tab'
            mscc_file = Path(temp_dir) / 'ENCFF000RMB-test_mscc.tab'

            assert stats_file.exists(), "Successive parallel processing failed to generate stats"
            # Note: SUCCESSIVE with mappability produces MSCC only (no NCC), which is correct behavior
            assert mscc_file.exists(), "Successive parallel processing failed to generate MSCC table"

            # Verify MSCC data contains actual calculations (not just headers)
            with open(mscc_file, 'r') as f:
                lines = f.readlines()
                assert len(lines) > 300, "Successive parallel MSCC table has insufficient data"

            print("✅ Successive algorithm BWFeederWithMappableRegionSum fix verified")


# Test markers
pytestmark = [
    pytest.mark.integration,
    pytest.mark.parallel,
    pytest.mark.critical  # Mark as critical since this prevents major bugs
]