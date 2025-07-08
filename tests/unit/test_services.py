"""Tests for PyMaSC service layer architecture."""

import unittest
from unittest.mock import Mock, MagicMock, patch
import numpy as np

from PyMaSC.services.calculation import (
    ChromosomeData, CalculationResult, GenomeWideResult,
    StandardCalculationService, create_calculation_service
)
from PyMaSC.services.io import (
    BAMFileInfo, InMemoryIOService, FileIOService, create_io_service
)
from PyMaSC.services.workflow import (
    WorkflowRequest, WorkflowResult, WorkflowStatus,
    StandardWorkflowService, create_workflow_service
)
from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig, 
    CalculationTarget, ImplementationAlgorithm
)


class TestCalculationService(unittest.TestCase):
    """Test cases for calculation service."""

    def setUp(self):
        """Set up test fixtures."""
        self.service = StandardCalculationService()

        # Create test chromosome data
        self.test_data = ChromosomeData(
            chromosome="chr1",
            forward_reads=[(100, 50), (200, 50), (300, 50)],
            reverse_reads=[(150, 50), (250, 50), (350, 50)],
            length=1000
        )

        # Create test configuration
        self.config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=100,
            mapq_criteria=20,
            references=["chr1"],
            lengths=[1000]
        )

    def test_chromosome_data_validation(self):
        """Test chromosome data validation."""
        # Valid data should work
        data = ChromosomeData("chr1", [], [], 1000)
        self.assertEqual(data.chromosome, "chr1")
        self.assertEqual(data.length, 1000)

        # Invalid length should raise error
        with self.assertRaises(ValueError):
            ChromosomeData("chr1", [], [], -1)

    def test_calculate_chromosome_basic(self):
        """Test basic chromosome calculation."""
        with patch('PyMaSC.services.calculation.CalculatorFactory.create_calculator') as mock_factory:
            # Mock calculator
            mock_calc = Mock()
            mock_calc.ref2forward_sum = {"chr1": 3}
            mock_calc.ref2reverse_sum = {"chr1": 3}
            mock_calc.ref2ccbins = {"chr1": np.array([0.1, 0.2, 0.3])}
            # Explicitly configure the mock to not have mappability attributes
            mock_calc.configure_mock(spec=['ref2forward_sum', 'ref2reverse_sum', 'ref2ccbins', 
                                           'feed_forward_read', 'feed_reverse_read', 
                                           'finishup_calculation'])
            mock_factory.return_value = mock_calc

            # Perform calculation
            result = self.service.calculate_chromosome(self.test_data, self.config)

            # Verify result
            self.assertEqual(result.chromosome, "chr1")
            self.assertEqual(result.forward_count, 3)
            self.assertEqual(result.reverse_count, 3)
            self.assertEqual(len(result.correlation_bins), 3)

            # Verify calculator was used correctly
            self.assertEqual(mock_calc.feed_forward_read.call_count, 3)
            self.assertEqual(mock_calc.feed_reverse_read.call_count, 3)
            mock_calc.finishup_calculation.assert_called_once()

    def test_genome_wide_calculation(self):
        """Test genome-wide calculation aggregation."""
        # Create multiple chromosome data
        data_list = [
            self.test_data,
            ChromosomeData(
                chromosome="chr2",
                forward_reads=[(100, 50), (200, 50)],
                reverse_reads=[(150, 50), (250, 50)],
                length=800
            )
        ]

        with patch('PyMaSC.services.calculation.CalculatorFactory.create_calculator') as mock_factory:
            # Mock calculator with different results per chromosome
            mock_calc = Mock()
            mock_calc.ref2forward_sum = {"chr1": 3, "chr2": 2}
            mock_calc.ref2reverse_sum = {"chr1": 3, "chr2": 2}
            mock_calc.ref2ccbins = {
                "chr1": np.array([0.1, 0.2, 0.3]),
                "chr2": np.array([0.2, 0.3, 0.4])
            }
            # Explicitly configure the mock to not have mappability attributes
            mock_calc.configure_mock(spec=['ref2forward_sum', 'ref2reverse_sum', 'ref2ccbins', 
                                           'feed_forward_read', 'feed_reverse_read', 
                                           'finishup_calculation'])
            mock_factory.return_value = mock_calc

            # Update config for multiple chromosomes
            self.config.references = ["chr1", "chr2"]
            self.config.lengths = [1000, 800]

            # Perform calculation
            result = self.service.calculate_genome_wide(data_list, self.config)

            # Verify aggregated results
            self.assertEqual(result.total_forward_reads, 5)
            self.assertEqual(result.total_reverse_reads, 5)
            self.assertEqual(len(result.chromosome_results), 2)
            self.assertIn("chr1", result.chromosome_results)
            self.assertIn("chr2", result.chromosome_results)

    def test_calculation_result_quality_metrics(self):
        """Test quality metrics calculation."""
        result = CalculationResult(
            chromosome="chr1",
            forward_count=100,
            reverse_count=100,
            correlation_bins=np.array([0.1, 0.2, 0.5, 0.3, 0.2, 0.1])
        )

        metrics = result.get_quality_metrics(read_length=2)

        self.assertIn('nsc', metrics)
        self.assertIn('rsc', metrics)
        self.assertIn('fragment_peak_pos', metrics)
        self.assertIn('phantom_peak_pos', metrics)


class TestIOService(unittest.TestCase):
    """Test cases for I/O service."""

    def setUp(self):
        """Set up test fixtures."""
        self.service = InMemoryIOService()

        # Add test data
        self.service.add_test_data(
            bam_path="/test/file.bam",
            chromosome="chr1",
            forward_reads=[(100, 50), (200, 50)],
            reverse_reads=[(150, 50)],
            length=1000
        )

    def test_get_bam_info(self):
        """Test BAM info retrieval."""
        info = self.service.get_bam_info("/test/file.bam")

        self.assertEqual(info.path, "/test/file.bam")
        self.assertTrue(info.has_index)
        self.assertIn("chr1", info.references)
        self.assertEqual(info.n_references, 1)

    def test_read_chromosome_data(self):
        """Test chromosome data reading."""
        data = self.service.read_chromosome_data(
            "/test/file.bam", "chr1", mapq_threshold=20
        )

        self.assertEqual(data.chromosome, "chr1")
        self.assertEqual(len(data.forward_reads), 2)
        self.assertEqual(len(data.reverse_reads), 1)
        self.assertEqual(data.length, 1000)

    def test_stream_reads(self):
        """Test read streaming."""
        reads = list(self.service.stream_reads(
            "/test/file.bam", "chr1", mapq_threshold=20
        ))

        self.assertEqual(len(reads), 3)

        # Check first forward read
        self.assertEqual(reads[0].chromosome, "chr1")
        self.assertEqual(reads[0].position, 100)
        self.assertEqual(reads[0].length, 50)
        self.assertFalse(reads[0].is_reverse)

    def test_write_results(self):
        """Test result writing."""
        # Create test result
        result = GenomeWideResult(
            chromosome_results={
                "chr1": CalculationResult(
                    chromosome="chr1",
                    forward_count=100,
                    reverse_count=100,
                    correlation_bins=np.array([0.1, 0.2, 0.3])
                )
            },
            total_forward_reads=100,
            total_reverse_reads=100,
            aggregated_correlation=np.array([0.1, 0.2, 0.3])
        )

        # Write results
        paths = self.service.write_results(
            result, "/test/output", ['table', 'stats']
        )

        self.assertIn('table', paths)
        self.assertIn('stats', paths)
        self.assertEqual(paths['table'], "/test/output.table")
        self.assertEqual(paths['stats'], "/test/output.stats")

        # Verify results were stored
        self.assertIn("/test/output.table", self.service.written_results)


class TestWorkflowService(unittest.TestCase):
    """Test cases for workflow service."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock services
        self.calc_service = Mock()
        self.io_service = Mock()

        # Create workflow service with mocks
        self.service = StandardWorkflowService(
            calculation_service=self.calc_service,
            io_service=self.io_service
        )

        # Create test request
        self.request = WorkflowRequest(
            bam_path="/test/file.bam",
            output_prefix="/test/output",
            calculation_config=CalculationConfig(
                target=CalculationTarget.NCC,
                implementation=ImplementationAlgorithm.SUCCESSIVE,
                max_shift=100,
                mapq_criteria=20
            ),
            output_formats=['table']
        )

    def test_validate_request_valid(self):
        """Test request validation with valid request."""
        errors = self.service.validate_request(self.request)
        self.assertEqual(len(errors), 0)

    def test_validate_request_invalid(self):
        """Test request validation with invalid request."""
        # Missing BAM path
        self.request.bam_path = ""
        errors = self.service.validate_request(self.request)
        self.assertIn("BAM file path is required", errors)

        # Invalid max shift
        self.request.bam_path = "/test/file.bam"
        self.request.calculation_config.max_shift = -1
        errors = self.service.validate_request(self.request)
        self.assertIn("Max shift must be positive", errors)

        # Target requires mappability
        self.request.calculation_config.max_shift = 100
        self.request.calculation_config.target = CalculationTarget.MSCC
        errors = self.service.validate_request(self.request)
        self.assertIn("mscc target requires mappability configuration", errors)

    def test_execute_workflow_success(self):
        """Test successful workflow execution."""
        # Mock BAM info
        self.io_service.get_bam_info.return_value = BAMFileInfo(
            path="/test/file.bam",
            has_index=True,
            references=["chr1"],
            lengths=[1000]
        )

        # Mock chromosome data
        chrom_data = ChromosomeData(
            chromosome="chr1",
            forward_reads=[(100, 50)],
            reverse_reads=[(150, 50)],
            length=1000
        )
        self.io_service.read_chromosome_data.return_value = chrom_data

        # Mock calculation result
        calc_result = GenomeWideResult(
            chromosome_results={
                "chr1": CalculationResult(
                    chromosome="chr1",
                    forward_count=1,
                    reverse_count=1,
                    correlation_bins=np.array([0.1])
                )
            },
            total_forward_reads=1,
            total_reverse_reads=1,
            aggregated_correlation=np.array([0.1])
        )
        self.calc_service.calculate_genome_wide.return_value = calc_result

        # Mock write results
        self.io_service.write_results.return_value = {
            'table': '/test/output.txt'
        }

        # Execute workflow
        result = self.service.execute(self.request)

        # Verify result
        self.assertEqual(result.status, WorkflowStatus.COMPLETED)
        self.assertTrue(result.is_successful)
        self.assertIsNotNone(result.calculation_result)
        self.assertIsNotNone(result.output_paths)
        self.assertIsNone(result.error)

        # Verify service calls
        self.io_service.get_bam_info.assert_called_once()
        self.io_service.read_chromosome_data.assert_called_once()
        self.calc_service.calculate_genome_wide.assert_called_once()
        self.io_service.write_results.assert_called_once()

    def test_execute_workflow_failure(self):
        """Test workflow execution with failure."""
        # Mock BAM info to raise error
        self.io_service.get_bam_info.side_effect = IOError("File not found")

        # Execute workflow
        result = self.service.execute(self.request)

        # Verify failure result
        self.assertEqual(result.status, WorkflowStatus.FAILED)
        self.assertFalse(result.is_successful)
        self.assertIsNone(result.calculation_result)
        self.assertIsNotNone(result.error)
        self.assertIn("File not found", result.error)

    def test_workflow_progress_notifications(self):
        """Test progress notifications during workflow."""
        # Track progress events
        progress_events = []

        def track_progress(event):
            progress_events.append(event)

        # Create a simple observer
        class TestObserver:
            def __init__(self, callback):
                self.callback = callback
            def update(self, event):
                self.callback(event)
            def should_handle(self, event):
                return True

        observer = TestObserver(track_progress)
        self.service.attach(observer)

        # Mock successful execution
        self.io_service.get_bam_info.return_value = BAMFileInfo(
            path="/test/file.bam",
            has_index=True,
            references=["chr1"],
            lengths=[1000]
        )

        # Execute workflow (will fail at chromosome read, but that's OK)
        try:
            self.service.execute(self.request)
        except:
            pass

        # Verify progress events were sent
        self.assertTrue(len(progress_events) > 0)

        # Check for expected progress messages
        messages = [e.message for e in progress_events]
        self.assertIn("workflow_started", messages)


class TestServiceIntegration(unittest.TestCase):
    """Integration tests for service layer."""

    def test_service_factory_functions(self):
        """Test service factory functions."""
        # Test calculation service creation
        calc_service = create_calculation_service(parallel=False)
        self.assertIsInstance(calc_service, StandardCalculationService)

        # Test I/O service creation
        io_service = create_io_service(test_mode=True)
        self.assertIsInstance(io_service, InMemoryIOService)

        # Test workflow service creation
        workflow_service = create_workflow_service(parallel=False)
        self.assertIsInstance(workflow_service, StandardWorkflowService)

    def test_end_to_end_with_test_services(self):
        """Test end-to-end workflow with test services."""
        # Create test I/O service with data
        io_service = InMemoryIOService()
        io_service.add_test_data(
            bam_path="/test/file.bam",
            chromosome="chr1",
            forward_reads=[(100, 50), (200, 50)],
            reverse_reads=[(150, 50)],
            length=1000
        )

        # Create workflow service with test I/O
        workflow_service = create_workflow_service(
            parallel=False,
            io_service=io_service
        )

        # Create request
        request = WorkflowRequest(
            bam_path="/test/file.bam",
            output_prefix="/test/output",
            calculation_config=CalculationConfig(
                target=CalculationTarget.NCC,
                implementation=ImplementationAlgorithm.SUCCESSIVE,
                max_shift=100,
                mapq_criteria=20,
                references=["chr1"],
                lengths=[1000]
            ),
            output_formats=['table']
        )

        # Execute workflow
        result = workflow_service.execute(request)

        # Verify success
        self.assertEqual(result.status, WorkflowStatus.COMPLETED)
        self.assertTrue(result.is_successful)
        self.assertIsNotNone(result.calculation_result)

        # Verify calculation result
        calc_result = result.calculation_result
        self.assertEqual(len(calc_result.chromosome_results), 1)
        self.assertIn("chr1", calc_result.chromosome_results)


if __name__ == '__main__':
    unittest.main()