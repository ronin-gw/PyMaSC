"""Tests for result aggregation modules in PyMaSC.core.result.

This module tests the result classes and aggregation functions that handle
cross-correlation results from individual chromosomes and their aggregation
into genome-wide statistics.
"""
import unittest
import numpy as np
import numpy.testing as npt
from unittest.mock import Mock

from PyMaSC.result import (
    CorrelationResult, NCCResult, MSCCResult, BothChromResult,
    GenomeWideResult, NCCGenomeWideResult, MSCCGenomeWideResult, BothGenomeWideResult,
    AggregateItems, NCCAggregateItems, MSCCAggregateItems,
    aggregate_results, _aggregate_ncc_results, _aggregate_mscc_results, _aggregate_both_results
)


class TestNCCResult(unittest.TestCase):
    """Test NCCResult class functionality."""

    def test_ncc_result_initialization(self):
        """Test NCCResult initialization."""
        ccbins = [10, 20, 80, 40, 10]
        
        result = NCCResult(
            max_shift=4,
            read_len=36,
            genomelen=3000000000,
            forward_sum=50000,
            reverse_sum=48000,
            ccbins=ccbins,
            forward_read_len_sum=1800000,
            reverse_read_len_sum=1728000
        )
        
        self.assertEqual(result.genomelen, 3000000000)
        self.assertEqual(result.forward_sum, 50000)
        self.assertEqual(result.reverse_sum, 48000)
        self.assertEqual(result.max_shift, 4)
        self.assertEqual(result.read_len, 36)
        self.assertEqual(result.forward_read_len_sum, 1800000)
        self.assertEqual(result.reverse_read_len_sum, 1728000)
        self.assertEqual(result.ccbins, ccbins)

    def test_ncc_result_calc_cc(self):
        """Test NCC cross-correlation calculation."""
        ccbins = [10, 20, 80, 40, 10]
        
        result = NCCResult(
            max_shift=4,
            read_len=36,
            genomelen=1000,  # Small genome for testing
            forward_sum=100,
            reverse_sum=120,
            ccbins=ccbins,
            forward_read_len_sum=3600,
            reverse_read_len_sum=4320
        )
        
        # Calculate cross-correlation
        result.calc_cc()
        
        # Should have calculated cc array
        self.assertTrue(hasattr(result, 'cc'))
        self.assertIsNotNone(result.cc)
        self.assertEqual(len(result.cc), result.max_shift + 1)
        
        # Check that cc values are finite (not NaN or inf)
        self.assertTrue(np.all(np.isfinite(result.cc)))

    def test_ncc_result_calc_cc_zero_ccbins(self):
        """Test NCC calculation with zero ccbins."""
        ccbins = [0, 0, 0, 0, 0]  # All zeros
        
        result = NCCResult(
            max_shift=4,
            read_len=36,
            genomelen=1000,
            forward_sum=100,
            reverse_sum=120,
            ccbins=ccbins,
            forward_read_len_sum=3600,
            reverse_read_len_sum=4320
        )
        
        result.calc_cc()
        
        # Should return NaN for all zero ccbins
        self.assertTrue(np.all(np.isnan(result.cc)))


class TestMSCCResult(unittest.TestCase):
    """Test MSCCResult class functionality."""

    def test_mscc_result_initialization(self):
        """Test MSCCResult initialization."""
        ccbins = [10, 20, 80, 40, 10]
        mappable_len = (100, 200, 300, 400, 500)  # Tuple as specified in interface
        forward_sum = np.array([5, 10, 15, 20, 25], dtype=np.int64)
        reverse_sum = np.array([4, 8, 12, 16, 20], dtype=np.int64)
        
        result = MSCCResult(
            max_shift=4,
            read_len=5,  # Use small read_len to match array sizes
            genomelen=3000000000,
            forward_sum=forward_sum,
            reverse_sum=reverse_sum,
            mappable_len=mappable_len,
            forward_read_len_sum=None,
            reverse_read_len_sum=None,
            ccbins=ccbins
        )
        
        self.assertEqual(result.genomelen, 3000000000)
        self.assertEqual(result.max_shift, 4)
        self.assertEqual(result.read_len, 5)
        self.assertEqual(result.mappable_len, mappable_len)
        npt.assert_array_equal(result.forward_sum, forward_sum)
        npt.assert_array_equal(result.reverse_sum, reverse_sum)
        self.assertEqual(result.ccbins, ccbins)

    def test_mscc_result_calc_cc(self):
        """Test MSCC cross-correlation calculation."""
        ccbins = [10, 20, 80, 40, 10]
        mappable_len = (100, 200, 300, 400, 500)
        forward_sum = np.array([5, 10, 15, 20, 25], dtype=np.int64)
        reverse_sum = np.array([4, 8, 12, 16, 20], dtype=np.int64)
        
        result = MSCCResult(
            max_shift=4,
            read_len=5,  # Small read length for testing
            genomelen=1000,
            forward_sum=forward_sum,
            reverse_sum=reverse_sum,
            mappable_len=mappable_len,
            forward_read_len_sum=None,
            reverse_read_len_sum=None,
            ccbins=ccbins
        )
        
        # Calculate cross-correlation
        result.calc_cc()
        
        # Should have calculated cc array
        self.assertTrue(hasattr(result, 'cc'))
        self.assertIsNotNone(result.cc)
        self.assertEqual(len(result.cc), result.max_shift + 1)

    def test_mscc_result_calc_cc_no_mappable_len(self):
        """Test MSCC calculation fails without mappable_len."""
        ccbins = [10, 20, 80, 40, 10]
        forward_sum = np.array([5, 10, 15, 20, 25], dtype=np.int64)
        reverse_sum = np.array([4, 8, 12, 16, 20], dtype=np.int64)
        
        result = MSCCResult(
            max_shift=4,
            read_len=5,
            genomelen=1000,
            forward_sum=forward_sum,
            reverse_sum=reverse_sum,
            mappable_len=None,  # Missing mappable_len
            forward_read_len_sum=None,
            reverse_read_len_sum=None,
            ccbins=ccbins
        )
        
        # Should raise assertion error
        with self.assertRaises(AssertionError):
            result.calc_cc()


class TestBothChromResult(unittest.TestCase):
    """Test BothChromResult class functionality."""

    def test_both_chrom_result_initialization(self):
        """Test BothChromResult initialization."""
        # Create mock NCC and MSCC results
        mock_ncc = Mock(spec=NCCResult)
        mock_mscc = Mock(spec=MSCCResult)
        
        result = BothChromResult(
            chrom=mock_ncc,
            mappable_chrom=mock_mscc
        )
        
        self.assertEqual(result.chrom, mock_ncc)
        self.assertEqual(result.mappable_chrom, mock_mscc)


class TestGenomeWideResults(unittest.TestCase):
    """Test genome-wide result classes."""

    def test_ncc_genome_wide_result(self):
        """Test NCCGenomeWideResult initialization."""
        mock_chroms = {
            'chr1': Mock(spec=NCCResult),
            'chr2': Mock(spec=NCCResult)
        }
        
        result = NCCGenomeWideResult(
            genomelen=3000000000,
            forward_read_len_sum=1800000,
            reverse_read_len_sum=1728000,
            forward_sum=50000,
            reverse_sum=48000,
            chroms=mock_chroms
        )
        
        self.assertEqual(result.forward_read_len_sum, 1800000)
        self.assertEqual(result.reverse_read_len_sum, 1728000)
        self.assertEqual(result.forward_sum, 50000)
        self.assertEqual(result.reverse_sum, 48000)
        self.assertEqual(result.chroms, mock_chroms)

    def test_mscc_genome_wide_result(self):
        """Test MSCCGenomeWideResult initialization."""
        mock_chroms = {
            'chr1': Mock(spec=MSCCResult),
            'chr2': Mock(spec=MSCCResult)
        }
        
        result = MSCCGenomeWideResult(
            genomelen=3000000000,
            forward_read_len_sum=1800000,
            reverse_read_len_sum=1728000,
            chroms=mock_chroms
        )
        
        self.assertEqual(result.forward_read_len_sum, 1800000)
        self.assertEqual(result.reverse_read_len_sum, 1728000)
        self.assertEqual(result.chroms, mock_chroms)

    def test_both_genome_wide_result(self):
        """Test BothGenomeWideResult initialization."""
        mock_ncc_chroms = {
            'chr1': Mock(spec=NCCResult),
            'chr2': Mock(spec=NCCResult)
        }
        mock_mscc_chroms = {
            'chr1': Mock(spec=MSCCResult),
            'chr2': Mock(spec=MSCCResult)
        }
        
        result = BothGenomeWideResult(
            forward_read_len_sum=1800000,
            reverse_read_len_sum=1728000,
            genomelen=3000000000,
            forward_sum=50000,
            reverse_sum=48000,
            chroms=mock_ncc_chroms,
            mappable_chroms=mock_mscc_chroms
        )
        
        self.assertEqual(result.genomelen, 3000000000)
        self.assertEqual(result.forward_sum, 50000)
        self.assertEqual(result.reverse_sum, 48000)
        self.assertEqual(result.chroms, mock_ncc_chroms)
        self.assertEqual(result.mappable_chroms, mock_mscc_chroms)


class TestAggregateItems(unittest.TestCase):
    """Test aggregate items classes."""

    def test_ncc_aggregate_items_initialization(self):
        """Test NCCAggregateItems initialization."""
        items = NCCAggregateItems()
        
        self.assertEqual(items.genomelen, 0)
        self.assertEqual(items.forward_read_len_sum, 0)
        self.assertEqual(items.reverse_read_len_sum, 0)
        self.assertEqual(items.forward_sum, 0)
        self.assertEqual(items.reverse_sum, 0)

    def test_ncc_aggregate_items_addition(self):
        """Test NCCAggregateItems addition."""
        items = NCCAggregateItems()
        
        # Create mock NCCResult
        mock_result = Mock(spec=NCCResult)
        mock_result.genomelen = 1000000
        mock_result.forward_sum = 500
        mock_result.reverse_sum = 480
        mock_result.forward_read_len_sum = 18000
        mock_result.reverse_read_len_sum = 17280
        
        new_items = items + mock_result
        
        self.assertEqual(new_items.genomelen, 1000000)
        self.assertEqual(new_items.forward_sum, 500)
        self.assertEqual(new_items.reverse_sum, 480)
        self.assertEqual(new_items.forward_read_len_sum, 18000)
        self.assertEqual(new_items.reverse_read_len_sum, 17280)

    def test_mscc_aggregate_items_initialization(self):
        """Test MSCCAggregateItems initialization."""
        items = MSCCAggregateItems()
        
        self.assertEqual(items.genomelen, 0)
        self.assertEqual(items.forward_read_len_sum, 0)
        self.assertEqual(items.reverse_read_len_sum, 0)

    def test_mscc_aggregate_items_addition(self):
        """Test MSCCAggregateItems addition."""
        items = MSCCAggregateItems()
        
        # Create mock MSCCResult
        mock_result = Mock(spec=MSCCResult)
        mock_result.genomelen = 1000000
        mock_result.forward_read_len_sum = 18000
        mock_result.reverse_read_len_sum = 17280
        
        new_items = items + mock_result
        
        self.assertEqual(new_items.genomelen, 1000000)
        self.assertEqual(new_items.forward_read_len_sum, 18000)
        self.assertEqual(new_items.reverse_read_len_sum, 17280)


class TestAggregationFunctions(unittest.TestCase):
    """Test result aggregation functions."""

    def test_aggregate_ncc_results(self):
        """Test NCC results aggregation."""
        # Create mock NCC results
        mock_result1 = Mock(spec=NCCResult)
        mock_result1.genomelen = 1000000
        mock_result1.forward_sum = 500
        mock_result1.reverse_sum = 480
        mock_result1.forward_read_len_sum = 18000
        mock_result1.reverse_read_len_sum = 17280
        
        mock_result2 = Mock(spec=NCCResult)
        mock_result2.genomelen = 2000000
        mock_result2.forward_sum = 600
        mock_result2.reverse_sum = 580
        mock_result2.forward_read_len_sum = 21600
        mock_result2.reverse_read_len_sum = 20880
        
        results = {
            'chr1': mock_result1,
            'chr2': mock_result2
        }
        
        aggregated = _aggregate_ncc_results(results)
        
        self.assertIsInstance(aggregated, NCCGenomeWideResult)
        self.assertEqual(aggregated.genomelen, 3000000)  # Sum of both
        self.assertEqual(aggregated.forward_sum, 1100)
        self.assertEqual(aggregated.reverse_sum, 1060)
        self.assertEqual(aggregated.forward_read_len_sum, 39600)
        self.assertEqual(aggregated.reverse_read_len_sum, 38160)
        self.assertEqual(aggregated.chroms, results)

    def test_aggregate_mscc_results(self):
        """Test MSCC results aggregation."""
        # Create mock MSCC results
        mock_result1 = Mock(spec=MSCCResult)
        mock_result1.genomelen = 1000000
        mock_result1.forward_read_len_sum = 18000
        mock_result1.reverse_read_len_sum = 17280
        
        mock_result2 = Mock(spec=MSCCResult)
        mock_result2.genomelen = 2000000
        mock_result2.forward_read_len_sum = 21600
        mock_result2.reverse_read_len_sum = 20880
        
        results = {
            'chr1': mock_result1,
            'chr2': mock_result2
        }
        
        aggregated = _aggregate_mscc_results(results)
        
        self.assertIsInstance(aggregated, MSCCGenomeWideResult)
        self.assertEqual(aggregated.genomelen, 3000000)
        self.assertEqual(aggregated.forward_read_len_sum, 39600)
        self.assertEqual(aggregated.reverse_read_len_sum, 38160)
        self.assertEqual(aggregated.chroms, results)

    def test_aggregate_both_results(self):
        """Test Both (NCC+MSCC) results aggregation."""
        # Create mock NCC and MSCC results
        mock_ncc1 = Mock(spec=NCCResult)
        mock_ncc1.genomelen = 1000000
        mock_ncc1.forward_sum = 500
        mock_ncc1.reverse_sum = 480
        mock_ncc1.forward_read_len_sum = 18000
        mock_ncc1.reverse_read_len_sum = 17280
        
        mock_mscc1 = Mock(spec=MSCCResult)
        mock_mscc1.genomelen = 1000000
        mock_mscc1.forward_read_len_sum = 18000
        mock_mscc1.reverse_read_len_sum = 17280
        
        both_result1 = BothChromResult(
            chrom=mock_ncc1,
            mappable_chrom=mock_mscc1
        )
        
        results = {'chr1': both_result1}
        
        aggregated = _aggregate_both_results(results)
        
        self.assertIsInstance(aggregated, BothGenomeWideResult)
        self.assertEqual(aggregated.genomelen, 1000000)
        self.assertEqual(aggregated.forward_sum, 500)
        self.assertEqual(aggregated.reverse_sum, 480)

    def test_aggregate_results_dispatch(self):
        """Test aggregate_results function dispatches correctly."""
        # Test NCC results
        mock_ncc = Mock(spec=NCCResult)
        mock_ncc.genomelen = 1000000
        mock_ncc.forward_sum = 500
        mock_ncc.reverse_sum = 480
        mock_ncc.forward_read_len_sum = 18000
        mock_ncc.reverse_read_len_sum = 17280
        
        ncc_results = {'chr1': mock_ncc}
        aggregated = aggregate_results(ncc_results)
        self.assertIsInstance(aggregated, NCCGenomeWideResult)
        
        # Test MSCC results
        mock_mscc = Mock(spec=MSCCResult)
        mock_mscc.genomelen = 1000000
        mock_mscc.forward_read_len_sum = 18000
        mock_mscc.reverse_read_len_sum = 17280
        
        mscc_results = {'chr1': mock_mscc}
        aggregated = aggregate_results(mscc_results)
        self.assertIsInstance(aggregated, MSCCGenomeWideResult)

    def test_aggregate_results_unknown_type(self):
        """Test aggregate_results with unknown result type."""
        # Create mock result of unknown type
        mock_unknown = Mock()
        unknown_results = {'chr1': mock_unknown}
        
        with self.assertRaises(TypeError):
            aggregate_results(unknown_results)


if __name__ == '__main__':
    unittest.main()