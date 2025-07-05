"""Tests for handler simplification and new services."""

import unittest
from unittest.mock import Mock, MagicMock, patch, mock_open
import tempfile
import json
import numpy as np
from pathlib import Path

from PyMaSC.services.validation import (
    ValidationService, StandardValidationService, ValidationResult,
    create_validation_service
)
from PyMaSC.services.mappability import (
    MappabilityService, BigWigMappabilityService, JSONMappabilityService,
    MappabilityStats, create_mappability_service
)
from PyMaSC.handler.simplified import (
    SimplifiedCalcHandler, HandlerDependencies, HandlerBuilder,
    create_simplified_handler
)
from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, AlgorithmType
)


class TestValidationService(unittest.TestCase):
    """Test cases for validation service."""

    def setUp(self):
        """Set up test fixtures."""
        self.service = StandardValidationService()

    def test_validate_bam_file_not_found(self):
        """Test BAM file validation with missing file."""
        result = self.service.validate_bam_file("/nonexistent/file.bam")

        self.assertFalse(result.is_valid)
        self.assertIn("BAM file not found", result.errors[0])

    def test_validate_bam_file_valid(self):
        """Test BAM file validation with mock valid file."""
        with patch('PyMaSC.services.validation.Path') as mock_path:
            mock_path.return_value.exists.return_value = True

            with patch('PyMaSC.services.validation.AlignmentFile') as mock_bam:
                # Mock BAM file
                mock_instance = MagicMock()
                mock_instance.closed = False
                mock_instance.references = ['chr1', 'chr2']
                mock_instance.lengths = [1000, 2000]
                mock_instance.nreferences = 2
                mock_instance.header = {'HD': {'SO': 'coordinate'}}

                # Mock context manager
                mock_bam.return_value.__enter__.return_value = mock_instance
                mock_bam.return_value.__exit__.return_value = None

                result = self.service.validate_bam_file("/test/file.bam")

                self.assertTrue(result.is_valid)
                self.assertEqual(len(result.errors), 0)
                self.assertIsNotNone(result.metadata)
                self.assertEqual(result.metadata['n_references'], 2)

    def test_validate_calculation_config(self):
        """Test calculation configuration validation."""
        # Valid config
        config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=500,
            mapq_criteria=30,
            references=['chr1'],
            lengths=[1000]
        )

        result = self.service.validate_calculation_config(config)
        self.assertTrue(result.is_valid)

        # Invalid max_shift
        config.max_shift = -1
        result = self.service.validate_calculation_config(config)
        self.assertFalse(result.is_valid)
        self.assertIn("max_shift must be positive", result.errors[0])

    def test_validate_workflow_request(self):
        """Test complete workflow validation."""
        config = CalculationConfig(
            algorithm=AlgorithmType.MSCC,
            max_shift=500,
            mapq_criteria=30,
            references=['chr1'],
            lengths=[1000]
        )

        # Missing mappability for MSCC
        result = self.service.validate_workflow_request(
            bam_path="/test/file.bam",
            output_prefix="/test/output",
            calculation_config=config
        )

        self.assertFalse(result.is_valid)
        self.assertIn("mscc requires mappability configuration", result.errors[-1])


class TestMappabilityService(unittest.TestCase):
    """Test cases for mappability services."""

    def test_bigwig_mappability_service(self):
        """Test BigWig mappability service."""
        with patch('PyMaSC.services.mappability.BigWigReader') as mock_reader_class:
            # Mock reader instance
            mock_reader = MagicMock()
            mock_reader.get_as_array.return_value = [0.5, 0.7, 0.9, None, 1.0]
            mock_reader.chroms.return_value = {'chr1': 1000}
            mock_reader_class.return_value = mock_reader

            service = BigWigMappabilityService("/test/file.bw")

            # Test get values
            values = service.get_mappability_values('chr1', 100, 105)
            self.assertEqual(len(values), 5)
            self.assertEqual(values[3], 0.0)  # None converted to 0

            # Test is mappable
            mock_reader.get_as_array.return_value = [0.8]
            self.assertTrue(service.is_position_mappable('chr1', 100))

            # Test caching
            service.get_mappability_values('chr1', 100, 105)
            # Should use cache, not call reader again
            self.assertEqual(mock_reader.get_as_array.call_count, 2)

    def test_json_mappability_service(self):
        """Test JSON mappability service."""
        test_data = {
            'stat': {
                'total_mappable_length': 90000000,
                'total_genome_length': 100000000,
                'genome_mappable_fraction': 0.9
            },
            'chromosomes': {
                'chr1': {
                    'mappable_length': 900000,
                    'mappable_fraction': 0.9,
                    'mean_mappability': 0.85
                }
            }
        }

        with patch('builtins.open', mock_open(read_data=json.dumps(test_data))):
            service = JSONMappabilityService("/test/file.json")

            # Test stats
            stats = service.calculate_mappability_stats('chr1', 1000000)
            self.assertEqual(stats.mappable_length, 900000)
            self.assertEqual(stats.mappable_fraction, 0.9)

            # Test genome-wide stats
            genome_stats = service.calculate_genome_wide_stats({'chr1': 1000000})
            self.assertEqual(genome_stats.total_mappable_length, 90000000)

    def test_mappability_service_factory(self):
        """Test mappability service factory."""
        # BigWig file
        with patch('PyMaSC.services.mappability.BigWigReader'):
            service = create_mappability_service("/test/file.bw")
            self.assertIsInstance(service, BigWigMappabilityService)

        # JSON file
        with patch('builtins.open', mock_open(read_data='{}')):
            service = create_mappability_service("/test/file.json")
            self.assertIsInstance(service, JSONMappabilityService)

        # Unsupported format
        with self.assertRaises(ValueError):
            create_mappability_service("/test/file.txt")


class TestSimplifiedHandler(unittest.TestCase):
    """Test cases for simplified handler."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock dependencies
        self.deps = HandlerDependencies()
        self.deps.calculation_service = Mock()
        self.deps.io_service = Mock()
        self.deps.validation_service = Mock()
        self.deps.workflow_service = Mock()

        # Mock BAM info
        self.deps.io_service.get_bam_info.return_value = Mock(
            references=['chr1'],
            lengths=[1000],
            has_index=True
        )

        # Create config
        self.config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=100,
            mapq_criteria=20,
            references=['chr1'],
            lengths=[1000]
        )

    def test_handler_initialization(self):
        """Test handler initialization with dependencies."""
        handler = SimplifiedCalcHandler(
            bam_path="/test/file.bam",
            config=self.config,
            dependencies=self.deps
        )

        self.assertEqual(handler.bam_path, "/test/file.bam")
        self.assertEqual(handler.config, self.config)
        self.assertIsNotNone(handler.deps.calculation_service)

    def test_handler_validation(self):
        """Test handler validation."""
        # Mock validation result
        self.deps.validation_service.validate_workflow_request.return_value = \
            ValidationResult(is_valid=True, errors=[], warnings=['Test warning'])

        handler = SimplifiedCalcHandler(
            bam_path="/test/file.bam",
            config=self.config,
            dependencies=self.deps
        )

        self.assertTrue(handler.validate())
        self.deps.validation_service.validate_workflow_request.assert_called_once()

    def test_handler_calculation(self):
        """Test handler calculation execution."""
        # Mock successful workflow
        mock_result = Mock()
        mock_result.is_successful = True
        mock_result.calculation_result = Mock(total_forward_reads=100)
        self.deps.workflow_service.execute.return_value = mock_result

        # Mock validation
        self.deps.validation_service.validate_workflow_request.return_value = \
            ValidationResult(is_valid=True, errors=[], warnings=[])

        handler = SimplifiedCalcHandler(
            bam_path="/test/file.bam",
            config=self.config,
            dependencies=self.deps
        )

        result = handler.run_calculation()
        self.assertIsNotNone(result)
        self.assertEqual(result.total_forward_reads, 100)
        self.deps.workflow_service.execute.assert_called_once()

    def test_handler_set_mappability(self):
        """Test setting mappability file."""
        handler = SimplifiedCalcHandler(
            bam_path="/test/file.bam",
            config=self.config,
            dependencies=self.deps
        )

        with patch('PyMaSC.handler.simplified.create_mappability_service') as mock_create:
            mock_service = Mock()
            mock_create.return_value = mock_service

            handler.set_mappability_file("/test/mappability.bw", read_len=36)

            # Algorithm should stay as BITARRAY (it's already mappability-aware)
            self.assertEqual(handler.config.algorithm, AlgorithmType.BITARRAY)

            # Check service was created
            self.assertEqual(handler.deps.mappability_service, mock_service)
            mock_create.assert_called_once()

    def test_handler_set_mappability_naive_cc(self):
        """Test setting mappability with NAIVE_CC algorithm."""
        # Create config with NAIVE_CC
        config = CalculationConfig(
            algorithm=AlgorithmType.NAIVE_CC,
            max_shift=100,
            mapq_criteria=20,
            references=['chr1'],
            lengths=[1000]
        )

        handler = SimplifiedCalcHandler(
            bam_path="/test/file.bam",
            config=config,
            dependencies=self.deps
        )

        with patch('PyMaSC.handler.simplified.create_mappability_service') as mock_create:
            mock_service = Mock()
            mock_create.return_value = mock_service

            handler.set_mappability_file("/test/mappability.bw", read_len=36)

            # NAIVE_CC should be upgraded to MSCC
            self.assertEqual(handler.config.algorithm, AlgorithmType.MSCC)


class TestHandlerBuilder(unittest.TestCase):
    """Test cases for handler builder."""

    def test_builder_basic(self):
        """Test basic handler building."""
        with patch('PyMaSC.handler.simplified.create_io_service') as mock_io:
            mock_io_service = Mock()
            mock_io_service.get_bam_info.return_value = Mock(
                references=['chr1'],
                lengths=[1000]
            )
            mock_io.return_value = mock_io_service

            handler = HandlerBuilder() \
                .with_bam_file("/test/file.bam") \
                .with_algorithm('bitarray') \
                .with_max_shift(200) \
                .with_mapq_threshold(25) \
                .build()

            self.assertEqual(handler.bam_path, "/test/file.bam")
            self.assertEqual(handler.config.algorithm, AlgorithmType.BITARRAY)
            self.assertEqual(handler.config.max_shift, 200)
            self.assertEqual(handler.config.mapq_criteria, 25)

    def test_builder_with_mappability(self):
        """Test builder with mappability configuration."""
        with patch('PyMaSC.handler.simplified.create_io_service') as mock_io:
            mock_io_service = Mock()
            mock_io_service.get_bam_info.return_value = Mock(
                references=['chr1'],
                lengths=[1000]
            )
            mock_io.return_value = mock_io_service

            with patch('PyMaSC.handler.simplified.create_mappability_service') as mock_map:
                mock_map_service = Mock()
                mock_map.return_value = mock_map_service

                handler = HandlerBuilder() \
                    .with_bam_file("/test/file.bam") \
                    .with_algorithm('ncc') \
                    .with_mappability("/test/mappability.bw") \
                    .build()

                # Algorithm should be upgraded to MSCC
                self.assertEqual(handler.config.algorithm, AlgorithmType.MSCC)
                mock_map.assert_called_once()

    def test_builder_missing_bam(self):
        """Test builder with missing BAM file."""
        with self.assertRaises(ValueError) as ctx:
            HandlerBuilder().build()

        self.assertIn("BAM file path is required", str(ctx.exception))

    def test_create_simplified_handler(self):
        """Test convenience function."""
        with patch('PyMaSC.handler.simplified.create_io_service') as mock_io:
            mock_io_service = Mock()
            mock_io_service.get_bam_info.return_value = Mock(
                references=['chr1'],
                lengths=[1000]
            )
            mock_io.return_value = mock_io_service

            handler = create_simplified_handler(
                bam_path="/test/file.bam",
                algorithm='successive',
                max_shift=300,
                n_workers=4
            )

            self.assertEqual(handler.bam_path, "/test/file.bam")
            self.assertEqual(handler.config.algorithm, AlgorithmType.SUCCESSIVE)


if __name__ == '__main__':
    unittest.main()