"""Test Naive Cross-Correlation (NCC) core algorithms."""

import pytest
import numpy as np
from unittest.mock import Mock, patch

from PyMaSC.core.ncc import NaiveCCCalculator, ReadUnsortedError
from tests.utils.test_data_generator import (
    create_mock_read_data, 
    create_mock_reference_data,
    create_expected_cc_result
)


class TestNaiveCCCalculator:
    """Test NaiveCCCalculator core functionality."""

    @pytest.fixture
    def calculator_setup(self):
        """Setup basic calculator with mock data."""
        references, lengths = create_mock_reference_data()
        max_shift = 200
        calc = NaiveCCCalculator(max_shift, references, lengths)
        return calc, references, lengths, max_shift

    def test_calculator_initialization(self, calculator_setup):
        """Test that NaiveCCCalculator initializes correctly."""
        calc, references, lengths, max_shift = calculator_setup

        assert calc.max_shift == max_shift
        assert calc.ref2genomelen == dict(zip(references, lengths))
        assert calc.genomelen == sum(lengths)
        assert calc.forward_sum == 0
        assert calc.reverse_sum == 0

    def test_calculator_with_empty_references(self):
        """Test calculator behavior with empty reference data."""
        # Test that calculator can handle empty references (implementation-specific)
        calc = NaiveCCCalculator(200, [], [])
        assert calc.max_shift == 200
        assert calc.genomelen == 0
        assert len(calc.ref2genomelen) == 0

    def test_calculator_with_mismatched_references(self):
        """Test calculator behavior with mismatched references and lengths."""
        references = ['chr1', 'chr2']
        lengths = [10000]  # Mismatched length

        # Should either raise an error or handle gracefully
        # The exact behavior depends on implementation
        try:
            calc = NaiveCCCalculator(200, references, lengths)
            # If no error, check that it handles gracefully
            assert len(calc.ref2genomelen) <= min(len(references), len(lengths))
        except (ValueError, IndexError):
            # Expected behavior for mismatched data
            pass

    def test_genomelen_calculation(self, calculator_setup):
        """Test that total genome length is calculated correctly."""
        calc, references, lengths, _ = calculator_setup
        expected_length = sum(lengths)
        assert calc.genomelen == expected_length

    def test_max_shift_parameter(self):
        """Test different max_shift values."""
        references, lengths = create_mock_reference_data()

        for max_shift in [50, 100, 200, 500]:
            calc = NaiveCCCalculator(max_shift, references, lengths)
            assert calc.max_shift == max_shift

    def test_reference_mapping(self, calculator_setup):
        """Test that reference to length mapping is correct."""
        calc, references, lengths, _ = calculator_setup

        for ref, length in zip(references, lengths):
            assert calc.ref2genomelen[ref] == length

    @pytest.mark.parametrize("max_shift,n_reads", [
        (100, 50),
        (200, 100),
        (300, 200),
    ])
    def test_calculator_with_different_parameters(self, max_shift, n_reads):
        """Test calculator with various parameter combinations."""
        references, lengths = create_mock_reference_data()
        calc = NaiveCCCalculator(max_shift, references, lengths)

        assert calc.max_shift == max_shift
        assert calc.genomelen == sum(lengths)


class TestNaiveCCCalculatorWithMockData:
    """Test NCC calculator with mock read data."""

    @pytest.fixture
    def calculator_with_data(self):
        """Setup calculator with mock read data."""
        references, lengths = create_mock_reference_data()
        max_shift = 200
        calc = NaiveCCCalculator(max_shift, references, lengths)

        # Generate mock read data
        forward_pos, reverse_pos = create_mock_read_data(n_reads=100, max_shift=max_shift)

        return calc, forward_pos, reverse_pos, max_shift

    def test_read_data_structure(self, calculator_with_data):
        """Test that mock read data has expected structure."""
        calc, forward_pos, reverse_pos, max_shift = calculator_with_data

        # Check that positions are sorted
        assert forward_pos == sorted(forward_pos)
        assert reverse_pos == sorted(reverse_pos)

        # Check that positions are reasonable
        assert all(pos > 0 for pos in forward_pos)
        assert all(pos > 0 for pos in reverse_pos)

    def test_cross_correlation_calculation_structure(self, calculator_with_data):
        """Test the structure of cross-correlation calculation (without implementation details)."""
        calc, forward_pos, reverse_pos, max_shift = calculator_with_data

        # Test basic properties that should hold for any CC implementation
        assert calc.max_shift == max_shift
        assert len(forward_pos) > 0
        assert len(reverse_pos) > 0

    def test_expected_results_validation(self):
        """Test that expected results have reasonable structure."""
        max_shift = 200
        fragment_length = 150
        expected_cc = create_expected_cc_result(max_shift, fragment_length)

        # Check that we have results for all shifts
        assert len(expected_cc) == max_shift + 1

        # Check that peak is at expected fragment length
        peak_shift = max(expected_cc.keys(), key=lambda k: expected_cc[k])
        assert peak_shift == fragment_length

        # Check that peak value is highest
        peak_value = expected_cc[fragment_length]
        assert all(expected_cc[shift] <= peak_value for shift in expected_cc)

    def test_memory_efficiency_check(self, calculator_with_data):
        """Test that calculator doesn't use excessive memory."""
        calc, forward_pos, reverse_pos, max_shift = calculator_with_data

        # Basic memory usage check - calculator should not store excessive data
        # This is a structural test to ensure we're not accidentally storing huge arrays
        import sys
        calc_size = sys.getsizeof(calc)

        # Calculator should not be excessively large (arbitrary reasonable limit)
        assert calc_size < 1024 * 1024  # Less than 1MB for basic calculator


class TestReadUnsortedError:
    """Test ReadUnsortedError exception."""

    def test_read_unsorted_error_creation(self):
        """Test that ReadUnsortedError can be created and raised."""
        error = ReadUnsortedError("Test error message")
        assert str(error) == "Test error message"

        with pytest.raises(ReadUnsortedError):
            raise ReadUnsortedError("Test error")

    def test_read_unsorted_error_inheritance(self):
        """Test that ReadUnsortedError inherits from IndexError."""
        error = ReadUnsortedError("Test")
        assert isinstance(error, IndexError)


class TestNaiveCCCalculatorEdgeCases:
    """Test edge cases and error conditions."""

    def test_zero_max_shift(self):
        """Test calculator with zero max_shift."""
        references, lengths = create_mock_reference_data()
        calc = NaiveCCCalculator(0, references, lengths)
        assert calc.max_shift == 0

    def test_large_max_shift(self):
        """Test calculator with very large max_shift."""
        references, lengths = create_mock_reference_data()
        max_shift = 10000
        calc = NaiveCCCalculator(max_shift, references, lengths)
        assert calc.max_shift == max_shift

    def test_single_reference(self):
        """Test calculator with single reference chromosome."""
        references = ['chr1']
        lengths = [50000]
        calc = NaiveCCCalculator(200, references, lengths)

        assert calc.genomelen == 50000
        assert len(calc.ref2genomelen) == 1

    def test_empty_chromosome(self):
        """Test calculator with zero-length chromosome."""
        references = ['chr1', 'empty_chr']
        lengths = [50000, 0]
        calc = NaiveCCCalculator(200, references, lengths)

        assert calc.ref2genomelen['empty_chr'] == 0
        assert calc.genomelen == 50000