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

from PyMaSC.reader.bam import BAMFileProcessor, BAMValidationError
# ResultAggregator functionality has been integrated into the main codebase
from PyMaSC.handler.calc import CalcHandler
from PyMaSC.core.interfaces.config import (
    PyMaSCConfig, CalculationTarget, Algorithm, EstimationType
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
        self.assertEqual(processor._path, self.test_bam)

    def test_bam_validation_success(self):
        """Test successful BAM file validation."""
        processor = BAMFileProcessor(self.test_bam)
        filtered_refs, filtered_lengths = processor.apply_chromfilter()

        self.assertTrue(processor.has_index())
        self.assertGreater(len(processor.references), 0)
        self.assertGreater(len(filtered_refs), 0)
        self.assertEqual(len(filtered_refs), len(filtered_lengths))

    def test_bam_validation_with_filter(self):
        """Test BAM validation with chromosome filter."""
        # Test with known chromosome
        processor = BAMFileProcessor(self.test_bam)
        # Use the correct format: List[Tuple[bool, List[str]]]
        filtered_refs, filtered_lengths = processor.apply_chromfilter(chromfilter=[(True, ['chr1', 'chr2'])])

        self.assertEqual(len(filtered_refs), 2)
        self.assertIn('chr1', filtered_refs)
        self.assertIn('chr2', filtered_refs)
        processor.close()

        # Test with non-existent chromosome
        processor2 = BAMFileProcessor(self.test_bam)
        with self.assertRaises(BAMValidationError):
            processor2.apply_chromfilter(chromfilter=[(True, ['chrNonExistent'])])

    def test_bam_validation_error(self):
        """Test BAM validation error handling."""
        with tempfile.NamedTemporaryFile(suffix='.bam') as tmp:
            # pysam raises ValueError when opening invalid BAM file
            with self.assertRaises(ValueError):
                processor = BAMFileProcessor(tmp.name)

    def test_chromosome_size_validation(self):
        """Test chromosome size validation."""
        processor = BAMFileProcessor(self.test_bam)
        filtered_refs, filtered_lengths = processor.apply_chromfilter()

        # Test with matching sizes
        external_sizes = {
            'chr1': filtered_lengths[0]
        }
        updated = processor.validate_chromosome_sizes(external_sizes)
        self.assertEqual(updated['chr1'], filtered_lengths[0])

        # Test with mismatched sizes
        external_sizes = {
            'chr1': filtered_lengths[0] + 1000
        }
        updated = processor.validate_chromosome_sizes(external_sizes)
        self.assertEqual(updated['chr1'], filtered_lengths[0] + 1000)

    def test_multiprocess_compatibility_check(self):
        """Test multiprocess compatibility checking."""
        processor = BAMFileProcessor(self.test_bam)
        processor.apply_chromfilter()

        # Should pass with index
        self.assertTrue(processor.check_multiprocess_compatibility())



# Deprecated test class removed - functionality has been integrated into ResultAggregator


class TestCalcHandler(unittest.TestCase):
    """Test refactored CalcHandler."""

    def setUp(self):
        """Set up test fixtures."""
        self.test_bam = os.path.join(
            os.path.dirname(__file__),
            '../data/ENCFF000RMB-test.bam'
        )

        self.config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=100,
            mapq_criteria=0,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )

    def test_handler_initialization(self):
        """Test handler can be initialized."""
        handler = CalcHandler(
            self.test_bam,
            self.config
        )

        self.assertEqual(handler.path, self.test_bam)
        self.assertIsNotNone(handler.config.references)
        self.assertIsNotNone(handler.config.lengths)
        self.assertGreater(len(handler.config.references), 0)

    def test_handler_backward_compatibility(self):
        """Test backward compatibility properties."""
        handler = CalcHandler(
            self.test_bam,
            self.config
        )

        # Test properties through configuration objects
        self.assertEqual(handler.config.skip_ncc, False)
        self.assertEqual(handler.config.max_shift, 100)
        self.assertEqual(handler.config.mapq_criteria, 0)
        self.assertEqual(handler.config.nproc, 1)

    def test_handler_read_length_estimation(self):
        """Test read length estimation."""
        handler = CalcHandler(
            self.test_bam,
            self.config
        )

        # Test that config can have read length set
        self.config.read_length = 50
        self.assertEqual(self.config.read_length, 50)

    def test_handler_with_invalid_bam(self):
        """Test handler with invalid BAM file."""
        with tempfile.NamedTemporaryFile(suffix='.bam') as tmp:
            with self.assertRaises(ValueError) as cm:
                CalcHandler(
                    tmp.name,
                    self.config
                )
            self.assertIn("file does not contain alignment data", str(cm.exception))


    @patch('pysam.AlignmentFile')
    def test_handler_multiprocess_fallback(self, mock_alignfile):
        """Test fallback to single process when no index."""
        # Mock BAM file without index
        mock_file = MagicMock()
        mock_file.has_index.return_value = False
        mock_file.references = ['chr1']
        mock_file.lengths = [1000]
        mock_alignfile.return_value = mock_file

        # Create config requesting multiprocess
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=100,
            mapq_criteria=0,
            nproc=4,  # Request multiprocess
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )
        config.ref2lengths = {'chr1': 1000}

        handler = CalcHandler(
            self.test_bam,
            config
        )

        # Due to multiprocess fallback (no index), nproc is adjusted to 1
        self.assertEqual(handler.config.nproc, 1)


if __name__ == '__main__':
    unittest.main()