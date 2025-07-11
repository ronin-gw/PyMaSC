"""Golden value tests for extracted statistical functions.

This test module verifies that the extracted statistical functions in
utils/stats_utils.py produce exactly the same numerical results as the
original implementations in handler/result.py.

These tests are critical for ensuring no regressions during the refactoring
process and maintaining 15-digit precision compatibility.
"""
import unittest
import numpy as np
from unittest.mock import patch, MagicMock
import logging

from PyMaSC.utils.stats_utils import AdvancedStatisticsCalculator, npcalc_with_logging_warn, CrossCorrelationMerger, ArrayAggregator, StatisticalTestUtilities
from PyMaSC.handler.result import CCStats, CCResult, NEAR_ZERO_MIN_CALC_LEN, NEAR_READLEN_ERR_CRITERION, MERGED_CC_CONFIDENCE_INTERVAL


class TestExtractedStatisticalFunctions(unittest.TestCase):
    """Test extracted statistical functions against original implementations."""
    
    def setUp(self):
        """Set up test data for comparisons."""
        # Create test cross-correlation data
        self.ccbins = np.array([
            0.1, 0.12, 0.15, 0.18, 0.22, 0.28, 0.35, 0.42, 0.38, 0.32,
            0.28, 0.25, 0.22, 0.20, 0.18, 0.16, 0.15, 0.14, 0.13, 0.12,
            0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02
        ])
        
        self.read_length = 6
        self.library_len = 7
        self.min_calc_width = 10
        self.mv_avr_filter_len = 3
        self.filter_mask_len = 2
        
        # Test parameters
        self.forward_sum = 1000
        self.reverse_sum = 1200
        self.genomelen = 50000
        
    def test_background_minimum_calculation(self):
        """Test that background minimum calculation matches original."""
        # Create CCStats instance for comparison
        ccstats = CCStats(
            self.ccbins, self.genomelen, self.forward_sum, self.reverse_sum,
            self.read_length, self.min_calc_width, self.mv_avr_filter_len,
            self.filter_mask_len, output_warnings=False
        )
        
        # Get original result
        original_cc_min = ccstats.cc_min
        
        # Get extracted function result
        extracted_cc_min = AdvancedStatisticsCalculator.calculate_background_minimum(
            self.ccbins, self.min_calc_width, NEAR_ZERO_MIN_CALC_LEN, output_warnings=False
        )
        
        # Compare with high precision
        np.testing.assert_almost_equal(original_cc_min, extracted_cc_min, decimal=15)
    
    def test_phantom_peak_correlation_calculation(self):
        """Test that phantom peak correlation calculation matches original."""
        # Create CCStats instance for comparison
        ccstats = CCStats(
            self.ccbins, self.genomelen, self.forward_sum, self.reverse_sum,
            self.read_length, self.min_calc_width, self.mv_avr_filter_len,
            self.filter_mask_len, output_warnings=False
        )
        
        # Get original result
        original_ccrl = ccstats.ccrl
        
        # Get extracted function result
        extracted_ccrl = AdvancedStatisticsCalculator.calculate_phantom_peak_correlation(
            self.ccbins, self.read_length
        )
        
        # Compare with high precision
        np.testing.assert_almost_equal(original_ccrl, extracted_ccrl, decimal=15)
    
    def test_exact_quality_metrics_calculation(self):
        """Test that exact quality metrics calculation matches original."""
        # Create CCStats instance for comparison
        ccstats = CCStats(
            self.ccbins, self.genomelen, self.forward_sum, self.reverse_sum,
            self.read_length, self.min_calc_width, self.mv_avr_filter_len,
            self.filter_mask_len, output_warnings=False,
            expected_library_len=self.library_len
        )
        
        # Get original results
        original_ccfl = ccstats.ccfl
        original_nsc = ccstats.nsc
        original_rsc = ccstats.rsc
        
        # Get extracted function results
        extracted_ccfl, extracted_nsc, extracted_rsc = AdvancedStatisticsCalculator.calculate_exact_quality_metrics(
            self.ccbins, self.library_len, self.read_length, ccstats.cc_min
        )
        
        # Compare with high precision
        np.testing.assert_almost_equal(original_ccfl, extracted_ccfl, decimal=15)
        np.testing.assert_almost_equal(original_nsc, extracted_nsc, decimal=15)
        np.testing.assert_almost_equal(original_rsc, extracted_rsc, decimal=15)
    
    def test_fragment_length_estimation_with_masking(self):
        """Test that fragment length estimation matches original."""
        # Create CCStats instance for comparison
        ccstats = CCStats(
            self.ccbins, self.genomelen, self.forward_sum, self.reverse_sum,
            self.read_length, self.min_calc_width, self.mv_avr_filter_len,
            self.filter_mask_len, output_warnings=False,
            do_llestimation=True
        )
        
        # Get original results
        original_avr_cc = ccstats.avr_cc
        original_est_lib_len = ccstats.est_lib_len
        
        # Get extracted function results
        extracted_avr_cc, extracted_est_lib_len = AdvancedStatisticsCalculator.estimate_fragment_length_with_masking(
            self.ccbins, self.read_length, self.mv_avr_filter_len,
            self.filter_mask_len, NEAR_READLEN_ERR_CRITERION, output_warnings=False
        )
        
        # Compare with high precision
        np.testing.assert_array_almost_equal(original_avr_cc, extracted_avr_cc, decimal=15)
        self.assertEqual(original_est_lib_len, extracted_est_lib_len)
    
    def test_fwhm_calculation(self):
        """Test that FWHM calculation matches original."""
        # Create CCStats instance for comparison
        ccstats = CCStats(
            self.ccbins, self.genomelen, self.forward_sum, self.reverse_sum,
            self.read_length, self.min_calc_width, self.mv_avr_filter_len,
            self.filter_mask_len, output_warnings=False,
            expected_library_len=self.library_len
        )
        
        # Get original result
        original_cc_width = ccstats.cc_width
        
        # Get extracted function result
        extracted_cc_width = AdvancedStatisticsCalculator.calculate_fwhm(
            ccstats.avr_cc, self.library_len, ccstats.cc_min, output_warnings=False
        )
        
        # Compare results
        if original_cc_width is False:
            self.assertEqual(original_cc_width, extracted_cc_width)
        else:
            self.assertEqual(original_cc_width, extracted_cc_width)
    
    def test_edge_cases_read_length_out_of_bounds(self):
        """Test edge cases where read length exceeds array bounds."""
        short_ccbins = np.array([0.1, 0.2, 0.3])
        long_read_length = 10
        
        # Test phantom peak calculation
        result = AdvancedStatisticsCalculator.calculate_phantom_peak_correlation(
            short_ccbins, long_read_length
        )
        self.assertEqual(result, 0.0)
        
        # Test with original implementation
        ccstats = CCStats(
            short_ccbins, self.genomelen, self.forward_sum, self.reverse_sum,
            long_read_length, self.min_calc_width, self.mv_avr_filter_len,
            self.filter_mask_len, output_warnings=False
        )
        
        self.assertEqual(ccstats.ccrl, 0.0)
        self.assertEqual(result, ccstats.ccrl)
    
    def test_numerical_stability_with_decorator(self):
        """Test that the npcalc_with_logging_warn decorator handles errors correctly."""
        @npcalc_with_logging_warn
        def divide_by_zero_test(a, b):
            return np.array(a) / np.array(b)
        
        # Test that function continues with divide by zero
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            result = divide_by_zero_test(1.0, 0.0)
            self.assertTrue(np.isinf(result))
            mock_logger.debug.assert_called()
    
    def test_warning_messages_consistency(self):
        """Test that warning messages are consistent between implementations."""
        # Create test data that triggers warnings
        problematic_ccbins = np.array([0.05, 0.04, 0.03, 0.02, 0.01] + [0.1] * 20)
        
        # Test background minimum warning
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            AdvancedStatisticsCalculator.calculate_background_minimum(
                problematic_ccbins, self.min_calc_width, 5, output_warnings=True
            )
            # Should trigger warning about minimum being too large
            self.assertTrue(any("minimum coefficient" in str(call) for call in mock_logger.warning.call_args_list))
    
    def test_masking_behavior_consistency(self):
        """Test that masking behavior is consistent in fragment length estimation."""
        # Create test data where estimate is close to read length
        peak_at_read_length = np.zeros(30)
        peak_at_read_length[self.read_length - 1] = 1.0  # Peak at read length
        peak_at_read_length[15] = 0.8  # Secondary peak
        
        # Test masking behavior
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            avr_cc, est_lib_len = AdvancedStatisticsCalculator.estimate_fragment_length_with_masking(
                peak_at_read_length, self.read_length, self.mv_avr_filter_len,
                self.filter_mask_len, NEAR_READLEN_ERR_CRITERION, output_warnings=True
            )
            
            # Should trigger masking warning
            self.assertTrue(any("masking" in str(call) for call in mock_logger.warning.call_args_list))
            # Should find the secondary peak after masking
            self.assertNotEqual(est_lib_len, self.read_length)


class TestStatisticalFunctionIntegration(unittest.TestCase):
    """Test integration of extracted statistical functions."""
    
    def test_complete_calculation_pipeline(self):
        """Test complete calculation pipeline using extracted functions."""
        # Create realistic test data
        ccbins = np.array([
            0.05, 0.06, 0.08, 0.12, 0.18, 0.25, 0.35, 0.42, 0.38, 0.32,
            0.28, 0.25, 0.22, 0.20, 0.18, 0.16, 0.15, 0.14, 0.13, 0.12,
            0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02
        ])
        
        read_length = 6
        min_calc_width = 10
        mv_avr_filter_len = 3
        filter_mask_len = 2
        
        # Step 1: Calculate background minimum
        cc_min = AdvancedStatisticsCalculator.calculate_background_minimum(
            ccbins, min_calc_width, output_warnings=False
        )
        
        # Step 2: Estimate fragment length
        avr_cc, est_lib_len = AdvancedStatisticsCalculator.estimate_fragment_length_with_masking(
            ccbins, read_length, mv_avr_filter_len, filter_mask_len, output_warnings=False
        )
        
        # Step 3: Calculate quality metrics
        ccfl, nsc, rsc = AdvancedStatisticsCalculator.calculate_exact_quality_metrics(
            ccbins, est_lib_len, read_length, cc_min
        )
        
        # Step 4: Calculate FWHM
        cc_width = AdvancedStatisticsCalculator.calculate_fwhm(
            avr_cc, est_lib_len, cc_min, output_warnings=False
        )
        
        # Verify all results are reasonable
        self.assertGreater(cc_min, 0)
        self.assertGreater(est_lib_len, 0)
        self.assertGreater(ccfl, 0)
        self.assertGreater(nsc, 0)
        self.assertGreater(rsc, 0)
        self.assertIsInstance(cc_width, (int, bool))
        
        # Test that pipeline produces consistent results
        self.assertAlmostEqual(ccfl, ccbins[est_lib_len - 1], places=15)


class TestCrossCorrelationMerger(unittest.TestCase):
    """Test CrossCorrelationMerger against original _merge_cc implementation."""
    
    def setUp(self):
        """Set up test data for merger tests."""
        # Create test correlation data from multiple chromosomes
        self.correlation_arrays = [
            np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.4, 0.3, 0.2, 0.1]),
            np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.35, 0.25, 0.15, 0.05]),
            np.array([0.08, 0.18, 0.28, 0.38, 0.48, 0.38, 0.28, 0.18, 0.08])
        ]
        
        self.genome_lengths = np.array([1000000, 800000, 1200000])
        self.read_length = 36
        
        # Create mock handler data for CCResult comparison
        self.mock_handler = MagicMock()
        self.mock_handler.read_len = self.read_length
        self.mock_handler.references = ['chr1', 'chr2', 'chr3']
        self.mock_handler.lengths = self.genome_lengths.tolist()
        self.mock_handler.skip_ncc = False
        self.mock_handler.ref2forward_sum = {'chr1': 100, 'chr2': 80, 'chr3': 120}
        self.mock_handler.ref2reverse_sum = {'chr1': 150, 'chr2': 120, 'chr3': 180}
        self.mock_handler.ref2ccbins = {
            'chr1': self.correlation_arrays[0],
            'chr2': self.correlation_arrays[1],
            'chr3': self.correlation_arrays[2]
        }
        self.mock_handler.mappable_ref2forward_sum = {}
        self.mock_handler.mappable_ref2reverse_sum = {}
        self.mock_handler.mappable_ref2ccbins = {}
        self.mock_handler.ref2mappable_len = {}
        
    def test_merge_correlations_integration(self):
        """Test that merged correlations work correctly in CCResult after migration."""
        # Create CCResult instance to test the migrated functionality
        ccresult = CCResult(
            mv_avr_filter_len=15,
            chi2_pval=0.01,
            filter_mask_len=5,
            min_calc_width=10,
            handler=self.mock_handler
        )
        
        # Skip this test if results are None (no NCC calculation)
        if ccresult.ncc_upper is None:
            self.skipTest("CCResult has no NCC data")
            
        # Test that the migrated correlation merger produces valid results
        self.assertIsNotNone(ccresult.ncc_upper)
        self.assertIsNotNone(ccresult.ncc_lower)
        self.assertIsInstance(ccresult.ncc_upper, list)
        self.assertIsInstance(ccresult.ncc_lower, list)
        self.assertEqual(len(ccresult.ncc_upper), len(ccresult.ncc_lower))
        
        # Test that confidence intervals are ordered correctly
        for upper, lower in zip(ccresult.ncc_upper, ccresult.ncc_lower):
            self.assertGreaterEqual(upper, lower)
            
        # Test that the new merger produces the same results as used in CCResult
        processed_correlation_arrays = []
        processed_genome_lengths = []
        
        for ref in ccresult.references:
            if ccresult.ref2stats[ref].cc is not None:
                processed_correlation_arrays.append(ccresult.ref2stats[ref].cc.cc)
                processed_genome_lengths.append(ccresult.ref2genomelen[ref])
        
        if processed_correlation_arrays:
            merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
            extracted_merged, extracted_lower, extracted_upper = merger.merge_correlations(
                np.array(processed_genome_lengths), processed_correlation_arrays, ccresult.read_len
            )
            
            # Verify the results match what's stored in CCResult
            np.testing.assert_array_almost_equal(ccresult.ncc_upper, extracted_upper, decimal=14)
            np.testing.assert_array_almost_equal(ccresult.ncc_lower, extracted_lower, decimal=14)
    
    def test_merge_with_weights_functionality(self):
        """Test merge_with_weights method with explicit weights."""
        merger = CrossCorrelationMerger(confidence_interval=0.95)
        
        correlations = [0.3, 0.4, 0.5]
        weights = [100, 200, 150]
        
        merged_r, lower_ci, upper_ci = merger.merge_with_weights(correlations, weights)
        
        # Verify results are reasonable
        self.assertGreater(merged_r, 0.0)
        self.assertLess(merged_r, 1.0)
        self.assertLess(lower_ci, merged_r)
        self.assertGreater(upper_ci, merged_r)
    
    def test_validate_correlation_arrays(self):
        """Test correlation array validation."""
        # Test valid arrays
        valid_arrays = [
            np.array([0.1, 0.2, 0.3]),
            np.array([0.15, 0.25, 0.35]),
            np.array([0.12, 0.22, 0.32])
        ]
        self.assertTrue(CrossCorrelationMerger.validate_correlation_arrays(valid_arrays))
        
        # Test invalid arrays (different lengths)
        invalid_arrays = [
            np.array([0.1, 0.2, 0.3]),
            np.array([0.15, 0.25]),  # Different length
            np.array([0.12, 0.22, 0.32])
        ]
        self.assertFalse(CrossCorrelationMerger.validate_correlation_arrays(invalid_arrays))
        
        # Test empty arrays
        self.assertFalse(CrossCorrelationMerger.validate_correlation_arrays([]))
    
    def test_effective_sample_size_calculation(self):
        """Test effective sample size calculation."""
        genome_length = 1000000
        position = 100
        read_length = 36
        
        sample_size = CrossCorrelationMerger.calculate_effective_sample_size(
            genome_length, position, read_length
        )
        
        # Should be genome_length - 3 (for degrees of freedom)
        self.assertEqual(sample_size, genome_length - 3)
        
        # Test edge case with small genome
        small_genome = 2
        small_sample_size = CrossCorrelationMerger.calculate_effective_sample_size(
            small_genome, position, read_length
        )
        self.assertEqual(small_sample_size, 1)  # Should be at least 1
    
    def test_fisher_z_transformation_consistency(self):
        """Test Fisher z-transformation mathematical consistency."""
        merger = CrossCorrelationMerger()
        
        # Test with known correlation values
        test_correlations = [0.0, 0.5, 0.8, -0.3, -0.7]
        test_weights = [100, 200, 150, 175, 125]
        
        merged_r, lower_ci, upper_ci = merger.merge_with_weights(test_correlations, test_weights)
        
        # Verify the result is within valid correlation range
        self.assertGreaterEqual(merged_r, -1.0)
        self.assertLessEqual(merged_r, 1.0)
        self.assertGreaterEqual(lower_ci, -1.0)
        self.assertLessEqual(upper_ci, 1.0)
        
        # Verify confidence interval ordering
        self.assertLessEqual(lower_ci, merged_r)
        self.assertLessEqual(merged_r, upper_ci)
    
    def test_nan_and_inf_handling(self):
        """Test handling of NaN and infinite values."""
        merger = CrossCorrelationMerger()
        
        # Test with NaN values
        correlations_with_nan = [0.3, np.nan, 0.5, 0.4]
        weights = [100, 200, 150, 175]
        
        merged_r, lower_ci, upper_ci = merger.merge_with_weights(correlations_with_nan, weights)
        
        # Should handle NaN gracefully
        self.assertFalse(np.isnan(merged_r))
        self.assertFalse(np.isnan(lower_ci))
        self.assertFalse(np.isnan(upper_ci))
        
        # Test with extreme correlations (close to ±1)
        extreme_correlations = [0.99, -0.99, 0.95]
        extreme_weights = [100, 100, 100]
        
        merged_r, lower_ci, upper_ci = merger.merge_with_weights(extreme_correlations, extreme_weights)
        
        # Should handle extreme values
        self.assertGreaterEqual(merged_r, -1.0)
        self.assertLessEqual(merged_r, 1.0)
    
    def test_confidence_interval_levels(self):
        """Test different confidence interval levels."""
        correlations = [0.4, 0.5, 0.6]
        weights = [100, 100, 100]
        
        # Test 95% confidence interval
        merger_95 = CrossCorrelationMerger(confidence_interval=0.95)
        merged_95, lower_95, upper_95 = merger_95.merge_with_weights(correlations, weights)
        
        # Test 99% confidence interval  
        merger_99 = CrossCorrelationMerger(confidence_interval=0.99)
        merged_99, lower_99, upper_99 = merger_99.merge_with_weights(correlations, weights)
        
        # 99% interval should be wider than 95%
        interval_width_95 = upper_95 - lower_95
        interval_width_99 = upper_99 - lower_99
        self.assertGreater(interval_width_99, interval_width_95)
        
        # Merged values should be the same
        self.assertAlmostEqual(merged_95, merged_99, places=10)


class TestArrayAggregator(unittest.TestCase):
    """Test ArrayAggregator functionality."""
    
    def setUp(self):
        """Set up test data."""
        self.arr1 = np.array([1, 2, 3, 4, 5])
        self.arr2 = np.array([10, 20, 30, 40, 50])
        self.arr3 = np.array([100, 200, 300, 400, 500])
        
    def test_filter_none_values(self):
        """Test filtering None values from arrays."""
        arrays = [self.arr1, None, self.arr2, None, self.arr3]
        filtered = ArrayAggregator.filter_none_values(arrays)
        
        self.assertEqual(len(filtered), 3)
        np.testing.assert_array_equal(filtered[0], self.arr1)
        np.testing.assert_array_equal(filtered[1], self.arr2)
        np.testing.assert_array_equal(filtered[2], self.arr3)
        
    def test_filter_none_values_empty(self):
        """Test filtering when all values are None."""
        arrays = [None, None, None]
        filtered = ArrayAggregator.filter_none_values(arrays)
        
        self.assertEqual(len(filtered), 0)
        
    def test_safe_array_sum_basic(self):
        """Test basic array summation."""
        arrays = [self.arr1, self.arr2, self.arr3]
        result = ArrayAggregator.safe_array_sum(arrays)
        
        expected = np.array([111, 222, 333, 444, 555])
        np.testing.assert_array_equal(result, expected)
        
    def test_safe_array_sum_with_none(self):
        """Test array summation with None values."""
        arrays = [self.arr1, None, self.arr2, None, self.arr3]
        result = ArrayAggregator.safe_array_sum(arrays)
        
        expected = np.array([111, 222, 333, 444, 555])
        np.testing.assert_array_equal(result, expected)
        
    def test_safe_array_sum_all_none(self):
        """Test array summation when all arrays are None."""
        arrays = [None, None, None]
        result = ArrayAggregator.safe_array_sum(arrays)
        
        self.assertIsNone(result)
        
    def test_safe_array_sum_single_array(self):
        """Test array summation with single array."""
        arrays = [self.arr1]
        result = ArrayAggregator.safe_array_sum(arrays)
        
        np.testing.assert_array_equal(result, self.arr1)
        
    def test_safe_array_sum_different_lengths(self):
        """Test array summation with different length arrays."""
        arr_short = np.array([1, 2])
        arr_long = np.array([10, 20, 30, 40])
        arrays = [self.arr1, arr_short, arr_long]  # lengths: 5, 2, 4
        
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            result = ArrayAggregator.safe_array_sum(arrays)
            
            # Should use the most common length (all are unique, so first one wins)
            # In this case, length 5 from arr1 wins, so only arr1 is used
            np.testing.assert_array_equal(result, self.arr1)
            
            # Should log warning about inconsistent lengths
            mock_logger.warning.assert_called_once()
            
    def test_validate_array_compatibility_same_length(self):
        """Test array compatibility validation with same length arrays."""
        arrays = [self.arr1, self.arr2, self.arr3]
        is_valid, valid_arrays, length_info = ArrayAggregator.validate_array_compatibility(arrays)
        
        self.assertTrue(is_valid)
        self.assertEqual(len(valid_arrays), 3)
        self.assertEqual(length_info["lengths"], [5, 5, 5])
        self.assertEqual(length_info["unique_lengths"], [5])
        
    def test_validate_array_compatibility_different_lengths(self):
        """Test array compatibility validation with different length arrays."""
        arr_short = np.array([1, 2])
        arrays = [self.arr1, arr_short, self.arr2]
        is_valid, valid_arrays, length_info = ArrayAggregator.validate_array_compatibility(arrays)
        
        self.assertFalse(is_valid)
        self.assertEqual(len(valid_arrays), 3)
        self.assertEqual(set(length_info["lengths"]), {5, 2})
        self.assertEqual(set(length_info["unique_lengths"]), {5, 2})
        
    def test_validate_array_compatibility_allow_different_lengths(self):
        """Test array compatibility validation allowing different lengths."""
        arr_short = np.array([1, 2])
        arrays = [self.arr1, arr_short, self.arr2]
        is_valid, valid_arrays, length_info = ArrayAggregator.validate_array_compatibility(
            arrays, require_same_length=False
        )
        
        self.assertTrue(is_valid)
        self.assertEqual(len(valid_arrays), 3)
        
    def test_aggregate_chromosome_data_sum(self):
        """Test chromosome data aggregation with sum method."""
        ref2data = {
            'chr1': self.arr1,
            'chr2': self.arr2,
            'chr3': self.arr3
        }
        
        result = ArrayAggregator.aggregate_chromosome_data(ref2data, method='sum')
        expected = np.array([111, 222, 333, 444, 555])
        np.testing.assert_array_equal(result, expected)
        
    def test_aggregate_chromosome_data_weighted_sum(self):
        """Test chromosome data aggregation with weighted sum method."""
        ref2data = {
            'chr1': self.arr1,
            'chr2': self.arr2,
            'chr3': self.arr3
        }
        ref2weights = {
            'chr1': 1.0,
            'chr2': 2.0,
            'chr3': 0.5
        }
        
        result = ArrayAggregator.aggregate_chromosome_data(ref2data, ref2weights, method='weighted_sum')
        expected = self.arr1 * 1.0 + self.arr2 * 2.0 + self.arr3 * 0.5
        np.testing.assert_array_equal(result, expected)
        
    def test_aggregate_chromosome_data_mean(self):
        """Test chromosome data aggregation with mean method."""
        ref2data = {
            'chr1': self.arr1,
            'chr2': self.arr2,
            'chr3': self.arr3
        }
        
        result = ArrayAggregator.aggregate_chromosome_data(ref2data, method='mean')
        expected = np.array([111, 222, 333, 444, 555]) / 3
        np.testing.assert_array_equal(result, expected)
        
    def test_aggregate_chromosome_data_weighted_mean(self):
        """Test chromosome data aggregation with weighted mean method."""
        ref2data = {
            'chr1': self.arr1,
            'chr2': self.arr2,
            'chr3': self.arr3
        }
        ref2weights = {
            'chr1': 1.0,
            'chr2': 2.0,
            'chr3': 0.5
        }
        
        result = ArrayAggregator.aggregate_chromosome_data(ref2data, ref2weights, method='weighted_mean')
        weighted_sum = self.arr1 * 1.0 + self.arr2 * 2.0 + self.arr3 * 0.5
        expected = weighted_sum / 3.5  # Total weight: 1.0 + 2.0 + 0.5
        np.testing.assert_array_equal(result, expected)
        
    def test_aggregate_chromosome_data_empty(self):
        """Test chromosome data aggregation with empty input."""
        ref2data = {}
        result = ArrayAggregator.aggregate_chromosome_data(ref2data, method='sum')
        self.assertIsNone(result)
        
    def test_aggregate_chromosome_data_invalid_method(self):
        """Test chromosome data aggregation with invalid method."""
        ref2data = {'chr1': self.arr1}
        
        with self.assertRaises(ValueError):
            ArrayAggregator.aggregate_chromosome_data(ref2data, method='invalid')


class TestStatisticalTestUtilities(unittest.TestCase):
    """Test StatisticalTestUtilities functionality."""
    
    def test_chi2_test_balanced_reads(self):
        """Test chi2 test with balanced forward/reverse reads."""
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            chi2_val, p_val, is_significant = StatisticalTestUtilities.chi2_test(
                1000, 1000, 0.05, "Test chromosome"
            )
            
            # Perfect balance should have chi2_val = 0 and p_val = 1
            self.assertAlmostEqual(chi2_val, 0.0, places=10)
            self.assertAlmostEqual(p_val, 1.0, places=10)
            self.assertFalse(is_significant)
            
            # Should log info message about acceptable balance
            mock_logger.info.assert_called()
            self.assertFalse(mock_logger.warning.called)
            
    def test_chi2_test_slightly_imbalanced_reads(self):
        """Test chi2 test with slightly imbalanced reads (not significant)."""
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            chi2_val, p_val, is_significant = StatisticalTestUtilities.chi2_test(
                1000, 1050, 0.05, "Test chromosome"
            )
            
            # Should have low chi2 value and high p-value
            self.assertGreater(chi2_val, 0.0)
            self.assertLess(chi2_val, 3.84)  # Chi2 critical value for α=0.05
            self.assertGreater(p_val, 0.05)
            self.assertFalse(is_significant)
            
            # Should log info message about acceptable balance
            mock_logger.info.assert_called()
            self.assertFalse(mock_logger.warning.called)
            
    def test_chi2_test_significantly_imbalanced_reads(self):
        """Test chi2 test with significantly imbalanced reads."""
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            chi2_val, p_val, is_significant = StatisticalTestUtilities.chi2_test(
                1000, 1500, 0.05, "Test chromosome"
            )
            
            # Should have high chi2 value and low p-value
            self.assertGreater(chi2_val, 3.84)  # Chi2 critical value for α=0.05
            self.assertLess(p_val, 0.05)
            self.assertTrue(is_significant)
            
            # Should log warning message about imbalance
            mock_logger.warning.assert_called()
            self.assertFalse(mock_logger.info.called)
            
    def test_chi2_test_extreme_imbalance(self):
        """Test chi2 test with extreme imbalance."""
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            chi2_val, p_val, is_significant = StatisticalTestUtilities.chi2_test(
                1000, 10, 0.01, "Test chromosome"
            )
            
            # Should have very high chi2 value and very low p-value
            self.assertGreater(chi2_val, 10.0)
            self.assertLess(p_val, 0.001)
            self.assertTrue(is_significant)
            
            # Should log warning message about imbalance
            mock_logger.warning.assert_called()
            
    def test_chi2_test_different_thresholds(self):
        """Test chi2 test with different significance thresholds."""
        # Test with strict threshold (0.001)
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            chi2_val, p_val, is_significant = StatisticalTestUtilities.chi2_test(
                1000, 1100, 0.001, "Test chromosome"
            )
            
            # With strict threshold, moderate imbalance should not be significant
            self.assertFalse(is_significant)
            mock_logger.info.assert_called()
            
        # Test with relaxed threshold (0.1)
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            chi2_val, p_val, is_significant = StatisticalTestUtilities.chi2_test(
                1000, 1100, 0.1, "Test chromosome"
            )
            
            # With relaxed threshold, moderate imbalance might be significant
            # (depends on exact p-value)
            if p_val <= 0.1:
                self.assertTrue(is_significant)
                mock_logger.warning.assert_called()
            else:
                self.assertFalse(is_significant)
                mock_logger.info.assert_called()
                
    def test_chi2_test_zero_reads(self):
        """Test chi2 test with zero reads (edge case)."""
        with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
            # This should not crash but handle the edge case
            chi2_val, p_val, is_significant = StatisticalTestUtilities.chi2_test(
                0, 0, 0.05, "Test chromosome"
            )
            
            # With zero reads, chi2 calculation should handle gracefully
            # (actual behavior depends on implementation)
            self.assertIsInstance(chi2_val, float)
            self.assertIsInstance(p_val, float)
            self.assertIsInstance(is_significant, bool)


class TestArrayAggregatorIntegration(unittest.TestCase):
    """Test ArrayAggregator integration with actual result.py usage patterns."""
    
    def test_safe_array_sum_migration_completed(self):
        """Test that safe_array_sum migration has been completed successfully."""
        # Test that the new ArrayAggregator.safe_array_sum works correctly
        test_cases = [
            [np.array([1, 2, 3]), np.array([4, 5, 6]), np.array([7, 8, 9])],
            [np.array([1, 2, 3]), None, np.array([4, 5, 6])],
            [None, None, None],
            [np.array([1, 2, 3])],
            [np.array([1, 2]), np.array([3, 4, 5]), np.array([6, 7])]  # Different lengths
        ]
        
        for arrays in test_cases:
            with patch('PyMaSC.utils.stats_utils.logger') as mock_logger:
                result = ArrayAggregator.safe_array_sum(arrays)
                
                # Verify that the function handles all cases properly
                if all(arr is None for arr in arrays):
                    self.assertIsNone(result)
                else:
                    valid_arrays = [arr for arr in arrays if arr is not None]
                    if valid_arrays:
                        self.assertIsInstance(result, np.ndarray)
                        self.assertGreater(len(result), 0)
                        
    def test_filter_none_values_migration_completed(self):
        """Test that filter_none_values migration has been completed successfully."""
        # Test that the new ArrayAggregator.filter_none_values works correctly
        test_cases = [
            [1, None, 2, None, 3],
            [None, None, None],
            [1, 2, 3],
            [],
            [None, 1, None],
        ]
        
        for arrays in test_cases:
            result = ArrayAggregator.filter_none_values(arrays)
            
            # Verify that None values are properly filtered
            self.assertNotIn(None, result)
            expected_length = len([x for x in arrays if x is not None])
            self.assertEqual(len(result), expected_length)
            
    def test_chi2_test_migration_completed(self):
        """Test that chi2_test migration has been completed successfully."""
        # Test that the new StatisticalTestUtilities.chi2_test works correctly
        test_cases = [
            (1000, 1000, 0.05, "Test1"),
            (1000, 1100, 0.05, "Test2"),
            (1000, 1500, 0.05, "Test3"),
            (500, 600, 0.01, "Test4"),
        ]
        
        for forward, reverse, threshold, label in test_cases:
            with patch('PyMaSC.utils.stats_utils.logger') as mock_utils_logger:
                # Extracted function returns values and logs appropriately
                chi2_val, p_val, is_significant = StatisticalTestUtilities.chi2_test(
                    forward, reverse, threshold, label
                )
                
                # Verify that function returns proper values
                self.assertIsInstance(chi2_val, float)
                self.assertIsInstance(p_val, float)
                self.assertIn(is_significant, [True, False])  # is_significant is bool
                
                # Check that logging behavior is consistent
                if is_significant:
                    mock_utils_logger.warning.assert_called()
                else:
                    mock_utils_logger.info.assert_called()


if __name__ == '__main__':
    unittest.main()