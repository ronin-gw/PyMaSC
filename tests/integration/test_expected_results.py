"""Test expected results storage and comparison for numerical validation.

This module provides infrastructure for storing and comparing expected
cross-correlation results for regression testing.
"""

import pytest
import json
import numpy as np
from pathlib import Path

# Import PyMaSC components
import pysam
# BigWigFile import removed - no longer used in this test file


class TestResultsStorage:
    """Test infrastructure for storing and retrieving expected results."""

    @pytest.fixture
    def results_dir(self):
        """Create directory for storing expected results."""
        test_dir = Path(__file__).parent
        results_dir = test_dir / "expected_results"
        results_dir.mkdir(exist_ok=True)
        return results_dir

    def test_create_results_directory(self, results_dir):
        """Test that results directory can be created."""
        assert results_dir.exists()
        assert results_dir.is_dir()

    def test_store_encode_data_summary(self, results_dir):
        """Store summary of ENCODE test data for future comparison."""
        base_dir = Path(__file__).parent.parent.parent
        test_data_dir = base_dir / 'tests' / 'data'
        bam_path = test_data_dir / 'ENCFF000RMB-test.bam'
        bigwig_path = test_data_dir / 'hg19_36mer-test.bigwig'

        if not bam_path.exists() or not bigwig_path.exists():
            pytest.skip("Test data files not found")

        # Extract BAM summary
        bam_summary = {}
        with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
            bam_summary['total_reads'] = bam.count()
            bam_summary['references'] = list(bam.references)[:5]  # First 5 only
            bam_summary['chr1_reads'] = bam.count('chr1')

            # Get read position summary
            positions = []
            for read in bam.fetch('chr1'):
                if not read.is_unmapped and read.mapping_quality >= 10:
                    positions.append(read.reference_start)

            if positions:
                bam_summary['chr1_position_range'] = [min(positions), max(positions)]
                bam_summary['chr1_valid_reads'] = len(positions)

        # Store summary
        summary_file = results_dir / 'encode_data_summary.json'
        with open(summary_file, 'w') as f:
            json.dump(bam_summary, f, indent=2)

        print(f"Stored BAM summary: {bam_summary}")
        assert summary_file.exists()

    def test_load_and_validate_summary(self, results_dir):
        """Test loading and validating stored summary."""
        summary_file = results_dir / 'encode_data_summary.json'

        if not summary_file.exists():
            pytest.skip("Summary file not created yet")

        with open(summary_file, 'r') as f:
            summary = json.load(f)

        # Validate summary structure
        required_keys = ['total_reads', 'references', 'chr1_reads', 'chr1_position_range', 'chr1_valid_reads']
        for key in required_keys:
            assert key in summary, f"Missing key in summary: {key}"

        # Validate reasonable values
        assert summary['total_reads'] > 1000, "Too few total reads"
        assert summary['chr1_reads'] > 1000, "Too few chr1 reads"
        assert 'chr1' in summary['references'], "chr1 not in references"

        print(f"Validated summary: {summary}")


class TestCrossCorrelationBaseline:
    """Establish baseline cross-correlation results for regression testing."""

    @pytest.fixture
    def encode_baseline_data(self):
        """Extract consistent baseline data from ENCODE files."""
        base_dir = Path(__file__).parent.parent.parent
        test_data_dir = base_dir / 'tests' / 'data'
        bam_path = test_data_dir / 'ENCFF000RMB-test.bam'

        if not bam_path.exists():
            pytest.skip("ENCODE test BAM file not found")

        # Extract specific region for consistent results
        forward_positions = []
        reverse_positions = []

        with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
            # Use specific genomic region for reproducible results
            region_start = 750000
            region_end = 800000

            for read in bam.fetch('chr1', region_start, region_end):
                if read.is_unmapped or read.mapping_quality < 20:  # Higher threshold
                    continue

                if read.is_reverse:
                    reverse_positions.append(read.reference_end)
                else:
                    forward_positions.append(read.reference_start)

        # Ensure sufficient data
        if len(forward_positions) < 20 or len(reverse_positions) < 20:
            pytest.skip("Insufficient high-quality reads in baseline region")

        return {
            'forward': sorted(forward_positions),
            'reverse': sorted(reverse_positions),
            'region': f"chr1:{region_start}-{region_end}",
            'min_mapq': 20
        }

    def test_calculate_baseline_correlation(self, encode_baseline_data):
        """Calculate baseline cross-correlation for future comparison."""
        forward_pos = np.array(encode_baseline_data['forward'])
        reverse_pos = np.array(encode_baseline_data['reverse'])

        # Calculate correlation at key shifts
        key_shifts = [0, 50, 100, 150, 200, 250]
        correlations = {}

        for shift in key_shifts:
            shifted_reverse = reverse_pos + shift

            # Count overlaps with fixed window
            overlaps = 0
            window_size = 25  # Smaller window for precision

            for fpos in forward_pos:
                overlaps += np.sum((shifted_reverse >= fpos - window_size) &
                                 (shifted_reverse <= fpos + window_size))

            correlations[shift] = int(overlaps)  # Store as int for exact comparison

        # Store baseline results
        baseline = {
            'correlations': correlations,
            'forward_reads': len(forward_pos),
            'reverse_reads': len(reverse_pos),
            'region': encode_baseline_data['region'],
            'window_size': 25,
            'min_mapq': 20
        }

        print(f"Baseline correlation: {correlations}")
        print(f"Baseline data: {baseline}")

        # Basic validation
        assert all(val >= 0 for val in correlations.values()), "Negative correlations found"
        assert max(correlations.values()) > 0, "No positive correlations found"

    def test_correlation_reproducibility(self, encode_baseline_data):
        """Test that correlation calculation is reproducible."""
        # Run calculation twice and compare
        forward_pos = np.array(encode_baseline_data['forward'])
        reverse_pos = np.array(encode_baseline_data['reverse'])

        def calculate_correlation(shift):
            shifted_reverse = reverse_pos + shift
            overlaps = 0
            for fpos in forward_pos:
                overlaps += np.sum((shifted_reverse >= fpos - 25) &
                                 (shifted_reverse <= fpos + 25))
            return overlaps

        # Calculate same shift multiple times
        shift = 100
        result1 = calculate_correlation(shift)
        result2 = calculate_correlation(shift)
        result3 = calculate_correlation(shift)

        # Results should be identical
        assert result1 == result2 == result3, f"Inconsistent results: {result1}, {result2}, {result3}"
        print(f"Reproducible correlation at shift {shift}: {result1}")


class TestRegressionDetection:
    """Test detection of regressions in cross-correlation calculations."""

    def test_detect_significant_changes(self):
        """Test ability to detect significant changes in results."""
        # Mock baseline and current results
        baseline = {
            'correlations': {0: 100, 50: 95, 100: 120, 150: 110, 200: 90},
            'forward_reads': 50,
            'reverse_reads': 55
        }

        # Simulate current results with small variation (acceptable)
        current_small_change = {
            'correlations': {0: 102, 50: 93, 100: 122, 150: 108, 200: 92},
            'forward_reads': 50,
            'reverse_reads': 55
        }

        # Simulate current results with large variation (regression)
        current_large_change = {
            'correlations': {0: 150, 50: 80, 100: 200, 150: 70, 200: 60},
            'forward_reads': 50,
            'reverse_reads': 55
        }

        def calculate_relative_difference(baseline_corr, current_corr):
            """Calculate relative difference between correlation values."""
            max_diff = 0
            for shift in baseline_corr:
                if shift in current_corr:
                    base_val = baseline_corr[shift]
                    curr_val = current_corr[shift]
                    if base_val > 0:
                        diff = abs(curr_val - base_val) / base_val
                        max_diff = max(max_diff, diff)
            return max_diff

        # Test small change detection
        small_diff = calculate_relative_difference(
            baseline['correlations'],
            current_small_change['correlations']
        )
        assert small_diff < 0.15, f"Small change detected as significant: {small_diff}"

        # Test large change detection
        large_diff = calculate_relative_difference(
            baseline['correlations'],
            current_large_change['correlations']
        )
        assert large_diff > 0.3, f"Large change not detected: {large_diff}"

        print(f"Small change: {small_diff:.3f}, Large change: {large_diff:.3f}")

    def test_peak_shift_detection(self):
        """Test detection of peak shift changes."""
        # Mock correlations with peak at different positions
        correlations1 = {0: 80, 50: 90, 100: 120, 150: 110, 200: 85}  # Peak at 100
        correlations2 = {0: 85, 50: 110, 100: 90, 150: 120, 200: 80}  # Peak at 150

        def find_peak(correlations):
            return max(correlations.keys(), key=lambda k: correlations[k])

        peak1 = find_peak(correlations1)
        peak2 = find_peak(correlations2)

        assert peak1 == 100, f"Wrong peak detection: {peak1}"
        assert peak2 == 150, f"Wrong peak detection: {peak2}"

        # Detect peak shift
        peak_shift = abs(peak1 - peak2)
        assert peak_shift == 50, f"Peak shift calculation wrong: {peak_shift}"

        print(f"Peak shift detected: {peak_shift}bp")


# Test markers
pytestmark = [
    pytest.mark.integration,
    pytest.mark.baseline_testing
]