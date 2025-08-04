"""Tests for refactored CalcHandler components.

This module tests the new components created for the refactored
CalcHandler to ensure they work correctly and maintain
compatibility with the existing system.
"""
import unittest
import tempfile
import os
from unittest.mock import Mock, patch, MagicMock
import numpy as np

from PyMaSC.core.bam_processor import BAMFileProcessor, BAMValidationError, BAMMetadata
# ResultAggregator functionality has been integrated into the main codebase
from PyMaSC.handler.calc import CalcHandler
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


# Deprecated test class removed - functionality has been integrated into ResultAggregator


class TestCalcHandler(unittest.TestCase):
    """Test refactored CalcHandler."""

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
        handler = CalcHandler(
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
        handler = CalcHandler(
            self.test_bam,
            self.calc_config,
            self.exec_config
        )

        # Test properties through configuration objects
        self.assertEqual(handler.config.skip_ncc, False)
        self.assertEqual(handler.config.max_shift, 100)
        self.assertEqual(handler.config.mapq_criteria, 0)
        self.assertEqual(handler.execution_config.worker_count, 1)

    def test_handler_read_length_estimation(self):
        """Test read length estimation."""
        handler = CalcHandler(
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
                CalcHandler(
                    tmp.name,
                    self.calc_config,
                    self.exec_config
                )
            self.assertIn("File has no sequences defined", str(cm.exception))


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

        handler = CalcHandler(
            self.test_bam,
            self.calc_config,
            exec_config
        )

        # Should fall back to single process
        self.assertEqual(handler.execution_config.mode, ExecutionMode.SINGLE_PROCESS)
        self.assertEqual(handler.execution_config.worker_count, 1)


if __name__ == '__main__':
    unittest.main()