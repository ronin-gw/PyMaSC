"""Tests for refactored UnifiedCalcHandler components.

This module tests the new components created for the refactored
UnifiedCalcHandler to ensure they work correctly and maintain
compatibility with the existing system.
"""
import unittest
import tempfile
import os
from unittest.mock import Mock, patch, MagicMock
import numpy as np

from PyMaSC.core.bam_processor import BAMFileProcessor, BAMValidationError, BAMMetadata
from PyMaSC.core.progress_coordinator import (
    ProgressCoordinator, ProgressMode, ProgressConfig,
    SingleLineProgressReporter, MultiLineProgressReporter,
    ObserverProgressReporter, NullProgressReporter
)
from PyMaSC.utils.stats_utils import ResultAggregator
from PyMaSC.handler.unified import UnifiedCalcHandler
from PyMaSC.core.models import (
    CalculationConfig, ExecutionConfig, MappabilityConfig,
    CalculationTarget, ImplementationAlgorithm, ExecutionMode
)


class TestBAMFileProcessor(unittest.TestCase):
    """Test BAM file processing functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_bam = os.path.join(
            os.path.dirname(__file__), 
            '../data/ENCFF000RMB-test.bam'
        )
    
    def test_bam_processor_creation(self):
        """Test BAM processor can be created."""
        processor = BAMFileProcessor(self.test_bam)
        self.assertEqual(processor.path, self.test_bam)
        self.assertFalse(processor.is_stdin)
    
    def test_bam_validation_success(self):
        """Test successful BAM file validation."""
        processor = BAMFileProcessor(self.test_bam)
        metadata = processor.validate_and_open()
        
        self.assertIsInstance(metadata, BAMMetadata)
        self.assertEqual(metadata.path, self.test_bam)
        self.assertTrue(metadata.has_index)
        self.assertGreater(len(metadata.references), 0)
        self.assertGreater(len(metadata.filtered_references), 0)
    
    def test_bam_validation_with_filter(self):
        """Test BAM validation with chromosome filter."""
        # Test with known chromosome
        processor = BAMFileProcessor(self.test_bam)
        # Use the correct format: List[Tuple[bool, List[str]]]
        metadata = processor.validate_and_open(chromfilter=[(True, ['chr1', 'chr2'])])
        
        self.assertEqual(len(metadata.filtered_references), 2)
        self.assertIn('chr1', metadata.filtered_references)
        self.assertIn('chr2', metadata.filtered_references)
        processor.close()
        
        # Test with non-existent chromosome
        processor2 = BAMFileProcessor(self.test_bam)
        with self.assertRaises(BAMValidationError):
            processor2.validate_and_open(chromfilter=[(True, ['chrNonExistent'])])
    
    def test_bam_validation_error(self):
        """Test BAM validation error handling."""
        with tempfile.NamedTemporaryFile(suffix='.bam') as tmp:
            processor = BAMFileProcessor(tmp.name)
            with self.assertRaises(BAMValidationError):
                processor.validate_and_open()
    
    def test_chromosome_size_validation(self):
        """Test chromosome size validation."""
        processor = BAMFileProcessor(self.test_bam)
        metadata = processor.validate_and_open()
        
        # Test with matching sizes
        external_sizes = {
            'chr1': metadata.filtered_lengths[0]
        }
        updated = processor.validate_chromosome_sizes(external_sizes)
        self.assertEqual(updated['chr1'], metadata.filtered_lengths[0])
        
        # Test with mismatched sizes
        external_sizes = {
            'chr1': metadata.filtered_lengths[0] + 1000
        }
        updated = processor.validate_chromosome_sizes(external_sizes)
        self.assertEqual(updated['chr1'], metadata.filtered_lengths[0] + 1000)
    
    def test_multiprocess_compatibility_check(self):
        """Test multiprocess compatibility checking."""
        processor = BAMFileProcessor(self.test_bam)
        processor.validate_and_open()
        
        # Should pass with index
        self.assertTrue(processor.check_multiprocess_compatibility())
    
    def test_context_manager(self):
        """Test BAM processor as context manager."""
        with BAMFileProcessor(self.test_bam) as processor:
            metadata = processor.validate_and_open()
            self.assertIsNotNone(metadata)


class TestProgressCoordinator(unittest.TestCase):
    """Test progress coordination functionality."""
    
    def test_single_line_reporter(self):
        """Test single-line progress reporter."""
        reporter = SingleLineProgressReporter()
        
        # Should not raise any exceptions
        reporter.start_chromosome('chr1', 1000)
        reporter.update('chr1', 500)
        reporter.finish_chromosome('chr1')
        reporter.cleanup()
    
    def test_multi_line_reporter(self):
        """Test multi-line progress reporter."""
        reporter = MultiLineProgressReporter()
        
        # Should not raise any exceptions
        reporter.start_chromosome('chr1', 1000)
        reporter.update('chr1', 500)
        reporter.finish_chromosome('chr1')
        reporter.cleanup()
    
    def test_null_reporter(self):
        """Test null progress reporter."""
        reporter = NullProgressReporter()
        
        # All methods should be no-ops
        reporter.start_chromosome('chr1', 1000)
        reporter.update('chr1', 500) 
        reporter.finish_chromosome('chr1')
        reporter.cleanup()
    
    def test_progress_coordinator_single_process(self):
        """Test progress coordinator for single process."""
        coordinator = ProgressCoordinator.create_for_single_process()
        reporter = coordinator.setup_single_process()
        
        self.assertIsNotNone(reporter)
        self.assertIsInstance(reporter, SingleLineProgressReporter)
    
    def test_progress_coordinator_with_observer(self):
        """Test progress coordinator with observer mode."""
        mock_observer = Mock()
        coordinator = ProgressCoordinator.create_for_single_process(
            use_observer=True,
            observers=[mock_observer]
        )
        reporter = coordinator.setup_single_process()
        
        self.assertIsInstance(reporter, ObserverProgressReporter)
    
    def test_progress_coordinator_multiprocess(self):
        """Test progress coordinator for multiprocess."""
        coordinator = ProgressCoordinator.create_for_multiprocess()
        
        mock_queue = Mock()
        mock_lock = Mock()
        coordinator.setup_multiprocess(mock_queue, mock_lock)
        
        reporter = coordinator.get_reporter()
        self.assertIsInstance(reporter, MultiLineProgressReporter)
    
    def test_multiprocess_report_handling(self):
        """Test handling of multiprocess progress reports."""
        from multiprocessing import Lock
        
        coordinator = ProgressCoordinator.create_for_multiprocess()
        
        mock_queue = Mock()
        real_lock = Lock()  # Use real lock instead of mock
        coordinator.setup_multiprocess(mock_queue, real_lock)
        
        # Test progress update
        is_progress = coordinator.handle_multiprocess_report(None, ('chr1', 500))
        self.assertTrue(is_progress)
        
        # Test result (not progress)
        is_progress = coordinator.handle_multiprocess_report('chr1', (1, 2, 3))
        self.assertFalse(is_progress)


class TestResultAggregatorExtensions(unittest.TestCase):
    """Test extended ResultAggregator functionality."""
    
    @unittest.skip("aggregate_mappable_results method not implemented yet")
    def test_aggregate_mappable_results(self):
        """Test aggregation of mappable results."""
        # Test with valid data
        mappable_forward = {
            'chr1': np.array([10, 20, 30]),
            'chr2': np.array([5, 10, 15])
        }
        mappable_reverse = {
            'chr1': np.array([15, 25, 35]),
            'chr2': np.array([8, 12, 18])
        }
        mappable_ccbins = {
            'chr1': np.array([100, 200, 300]),
            'chr2': np.array([50, 100, 150])
        }
        
        result = ResultAggregator.aggregate_mappable_results(
            mappable_forward, mappable_reverse, mappable_ccbins
        )
        
        np.testing.assert_array_equal(
            result['mappable_forward_sum'],
            np.array([15, 30, 45])
        )
        np.testing.assert_array_equal(
            result['mappable_reverse_sum'],
            np.array([23, 37, 53])
        )
        np.testing.assert_array_equal(
            result['mappable_ccbins'],
            np.array([150, 300, 450])
        )
    
    @unittest.skip("aggregate_mappable_results method not implemented yet")
    def test_aggregate_mappable_results_with_none(self):
        """Test aggregation with None values."""
        mappable_forward = {
            'chr1': np.array([10, 20, 30]),
            'chr2': None
        }
        
        result = ResultAggregator.aggregate_mappable_results(
            mappable_forward, {}, {}
        )
        
        np.testing.assert_array_equal(
            result['mappable_forward_sum'],
            np.array([10, 20, 30])
        )
        self.assertIsNone(result['mappable_reverse_sum'])
        self.assertIsNone(result['mappable_ccbins'])
    
    @unittest.skip("collect_calculator_results method not implemented yet")
    def test_collect_calculator_results(self):
        """Test collecting results from calculator."""
        # Mock calculator
        mock_calc = Mock()
        mock_calc.ref2forward_sum = {'chr1': 100}
        mock_calc.ref2reverse_sum = {'chr1': 150}
        mock_calc.ref2ccbins = {'chr1': np.array([1, 2, 3])}
        # Ensure no mappability attributes
        mock_calc.ref2mappable_forward_sum = None
        # Ensure no _calculator attribute
        if hasattr(mock_calc, '_calculator'):
            del mock_calc._calculator
        
        results = ResultAggregator.collect_calculator_results(mock_calc)
        
        self.assertEqual(results['ref2forward_sum'], {'chr1': 100})
        self.assertEqual(results['ref2reverse_sum'], {'chr1': 150})
        self.assertIn('chr1', results['ref2ccbins'])
    
    @unittest.skip("collect_calculator_results method not implemented yet")
    def test_collect_calculator_results_with_mappability(self):
        """Test collecting results with mappability data."""
        # Mock calculator with mappability
        mock_calc = Mock()
        mock_calc.ref2forward_sum = {'chr1': 100}
        mock_calc.ref2reverse_sum = {'chr1': 150}
        mock_calc.ref2ccbins = {'chr1': np.array([1, 2, 3])}
        mock_calc.ref2mappable_forward_sum = {'chr1': 90}
        mock_calc.ref2mappable_reverse_sum = {'chr1': 140}
        mock_calc.ref2mascbins = {'chr1': np.array([2, 3, 4])}
        # Mock the inner calculator
        mock_calc._calculator = Mock()
        mock_calc._calculator.ref2mappable_len = Mock()
        mock_calc._calculator.ref2mappable_len.items.return_value = []
        
        results = ResultAggregator.collect_calculator_results(mock_calc)
        
        self.assertEqual(results['mappable_ref2forward_sum'], {'chr1': 90})
        self.assertEqual(results['mappable_ref2reverse_sum'], {'chr1': 140})
        self.assertIn('chr1', results['mappable_ref2ccbins'])


class TestUnifiedCalcHandler(unittest.TestCase):
    """Test refactored UnifiedCalcHandler."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_bam = os.path.join(
            os.path.dirname(__file__), 
            '../data/ENCFF000RMB-test.bam'
        )
        
        self.calc_config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=100,
            mapq_criteria=0,
            skip_ncc=False
        )
        
        self.exec_config = ExecutionConfig(
            mode=ExecutionMode.SINGLE_PROCESS,
            worker_count=1
        )
    
    def test_handler_initialization(self):
        """Test handler can be initialized."""
        handler = UnifiedCalcHandler(
            self.test_bam,
            self.calc_config,
            self.exec_config
        )
        
        self.assertEqual(handler.path, self.test_bam)
        self.assertIsNotNone(handler.references)
        self.assertIsNotNone(handler.lengths)
        self.assertGreater(len(handler.references), 0)
    
    def test_handler_backward_compatibility(self):
        """Test backward compatibility properties."""
        handler = UnifiedCalcHandler(
            self.test_bam,
            self.calc_config,
            self.exec_config
        )
        
        # Test properties
        self.assertEqual(handler.skip_ncc, False)
        self.assertEqual(handler.max_shift, 100)
        self.assertEqual(handler.mapq_criteria, 0)
        self.assertEqual(handler.nworker, 1)
    
    def test_handler_read_length_estimation(self):
        """Test read length estimation."""
        handler = UnifiedCalcHandler(
            self.test_bam,
            self.calc_config,
            self.exec_config
        )
        
        # Set explicit read length
        handler.set_or_estimate_readlen(50)
        self.assertEqual(handler.read_len, 50)
        self.assertEqual(handler.config.read_length, 50)
    
    def test_handler_with_invalid_bam(self):
        """Test handler with invalid BAM file."""
        with tempfile.NamedTemporaryFile(suffix='.bam') as tmp:
            with self.assertRaises(ValueError) as cm:
                UnifiedCalcHandler(
                    tmp.name,
                    self.calc_config,
                    self.exec_config
                )
            self.assertIn("File has no sequences defined", str(cm.exception))
    
    def test_handler_progress_observer_methods(self):
        """Test progress observer methods."""
        handler = UnifiedCalcHandler(
            self.test_bam,
            self.calc_config,
            self.exec_config
        )
        
        mock_observer = Mock()
        
        # Test attach
        handler.attach_progress_observer(mock_observer)
        self.assertIn(mock_observer, handler._progress_observers)
        
        # Test detach
        handler.detach_progress_observer(mock_observer)
        self.assertNotIn(mock_observer, handler._progress_observers)
        
        # Test enable/disable
        handler.enable_observer_progress(False)
        self.assertFalse(handler.use_observer_progress)
        
        handler.enable_observer_progress(True)
        self.assertTrue(handler.use_observer_progress)
    
    @patch('pysam.AlignmentFile')
    def test_handler_multiprocess_fallback(self, mock_alignfile):
        """Test fallback to single process when no index."""
        # Mock BAM file without index
        mock_file = MagicMock()
        mock_file.has_index.return_value = False
        mock_file.references = ['chr1']
        mock_file.lengths = [1000]
        mock_alignfile.return_value = mock_file
        
        exec_config = ExecutionConfig(
            mode=ExecutionMode.MULTI_PROCESS,
            worker_count=4
        )
        
        handler = UnifiedCalcHandler(
            self.test_bam,
            self.calc_config,
            exec_config
        )
        
        # Should fall back to single process
        self.assertEqual(handler.execution_config.mode, ExecutionMode.SINGLE_PROCESS)
        self.assertEqual(handler.execution_config.worker_count, 1)


if __name__ == '__main__':
    unittest.main()