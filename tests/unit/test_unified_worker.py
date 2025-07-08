"""Tests for unified worker architecture."""

import unittest
from unittest.mock import Mock, MagicMock, patch
from multiprocessing import Queue, Lock

from PyMaSC.core.worker import UnifiedWorker, StandardReadProcessor, DualReadProcessor
from PyMaSC.core.models import (
    WorkerConfig, CalculationConfig, MappabilityConfig, 
    CalculationTarget, ImplementationAlgorithm
)
from PyMaSC.core.worker_compat import (
    NaiveCCCalcWorker, MSCCCalcWorker, NCCandMSCCCalcWorker, BACalcWorker
)


class TestUnifiedWorker(unittest.TestCase):
    """Test cases for UnifiedWorker class."""

    def setUp(self):
        """Set up test fixtures."""
        self.order_queue = Queue()
        self.report_queue = Queue()
        self.logger_lock = Lock()
        self.bam_path = "/test/path.bam"

        # Create basic configuration
        self.calc_config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20,
            references=["chr1", "chr2"],
            lengths=[100000, 90000]
        )

        self.worker_config = WorkerConfig(
            calculation_config=self.calc_config,
            mappability_config=None,
            progress_enabled=False,
            logger_lock=self.logger_lock
        )

    def test_unified_worker_initialization(self):
        """Test UnifiedWorker initialization."""
        worker = UnifiedWorker(
            self.worker_config,
            self.order_queue,
            self.report_queue,
            self.logger_lock,
            self.bam_path
        )

        self.assertIsNotNone(worker)
        self.assertEqual(worker.config, self.worker_config)
        self.assertEqual(worker.bam_path, self.bam_path)
        self.assertIsNone(worker._calculator)
        self.assertIsNone(worker._processor)

    def test_unified_worker_with_mappability(self):
        """Test UnifiedWorker with mappability configuration."""
        map_config = MappabilityConfig(
            mappability_path="/test/mappability.bw",
            read_len=50
        )

        calc_config = CalculationConfig(
            target=CalculationTarget.MSCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20,
            references=["chr1", "chr2"],
            lengths=[100000, 90000],
            read_length=50
        )

        worker_config = WorkerConfig(
            calculation_config=calc_config,
            mappability_config=map_config,
            progress_enabled=False,
            logger_lock=self.logger_lock
        )

        worker = UnifiedWorker(
            worker_config,
            self.order_queue,
            self.report_queue,
            self.logger_lock,
            self.bam_path
        )

        self.assertIsNotNone(worker)
        self.assertEqual(worker.config.mappability_config, map_config)


class TestStandardReadProcessor(unittest.TestCase):
    """Test cases for StandardReadProcessor class."""

    def setUp(self):
        """Set up test fixtures."""
        self.calculator = Mock()
        self.processor = StandardReadProcessor(self.calculator, mapq_criteria=20)

    def test_should_skip_read_criteria(self):
        """Test read filtering criteria."""
        # Mock read that should be skipped
        read = Mock()
        read.is_read2 = True
        read.mapping_quality = 30
        read.is_unmapped = False
        read.is_duplicate = False

        self.assertTrue(self.processor.should_skip_read(read))

        # Mock read that should be processed
        read.is_read2 = False
        read.mapping_quality = 30
        read.is_unmapped = False
        read.is_duplicate = False

        self.assertFalse(self.processor.should_skip_read(read))

    def test_process_forward_read(self):
        """Test processing forward read."""
        read = Mock()
        read.is_reverse = False
        read.reference_name = "chr1"
        read.reference_start = 100
        read.infer_query_length.return_value = 50
        # Add required attributes for filtering
        read.is_read2 = False
        read.mapping_quality = 30
        read.is_unmapped = False
        read.is_duplicate = False

        self.processor.process_read(read)

        self.calculator.feed_forward_read.assert_called_once_with("chr1", 101, 50)

    def test_process_reverse_read(self):
        """Test processing reverse read."""
        read = Mock()
        read.is_reverse = True
        read.reference_name = "chr1"
        read.reference_start = 200
        read.infer_query_length.return_value = 50
        # Add required attributes for filtering
        read.is_read2 = False
        read.mapping_quality = 30
        read.is_unmapped = False
        read.is_duplicate = False

        self.processor.process_read(read)

        self.calculator.feed_reverse_read.assert_called_once_with("chr1", 201, 50)


class TestDualReadProcessor(unittest.TestCase):
    """Test cases for DualReadProcessor class."""

    def setUp(self):
        """Set up test fixtures."""
        self.primary_calculator = Mock()
        self.secondary_calculator = Mock()
        self.processor = DualReadProcessor(
            self.primary_calculator,
            self.secondary_calculator,
            mapq_criteria=20
        )

    def test_process_forward_read_dual(self):
        """Test processing forward read with dual calculators."""
        read = Mock()
        read.is_reverse = False
        read.reference_name = "chr1"
        read.reference_start = 100
        read.infer_query_length.return_value = 50
        # Add required attributes for filtering
        read.is_read2 = False
        read.mapping_quality = 30
        read.is_unmapped = False
        read.is_duplicate = False

        self.processor.process_read(read)

        self.primary_calculator.feed_forward_read.assert_called_once_with("chr1", 101, 50)
        self.secondary_calculator.feed_forward_read.assert_called_once_with("chr1", 101, 50)

    def test_process_reverse_read_dual(self):
        """Test processing reverse read with dual calculators."""
        read = Mock()
        read.is_reverse = True
        read.reference_name = "chr1"
        read.reference_start = 200
        read.infer_query_length.return_value = 50
        # Add required attributes for filtering
        read.is_read2 = False
        read.mapping_quality = 30
        read.is_unmapped = False
        read.is_duplicate = False

        self.processor.process_read(read)

        self.primary_calculator.feed_reverse_read.assert_called_once_with("chr1", 201, 50)
        self.secondary_calculator.feed_reverse_read.assert_called_once_with("chr1", 201, 50)


class TestWorkerCompatibility(unittest.TestCase):
    """Test cases for worker compatibility wrappers."""

    def setUp(self):
        """Set up test fixtures."""
        self.order_queue = Queue()
        self.report_queue = Queue()
        self.logger_lock = Lock()
        self.bam_path = "/test/path.bam"
        self.references = ["chr1", "chr2"]
        self.lengths = [100000, 90000]

    def test_naive_cc_worker_compatibility(self):
        """Test NaiveCCCalcWorker compatibility wrapper."""
        worker = NaiveCCCalcWorker(
            self.order_queue,
            self.report_queue,
            self.logger_lock,
            self.bam_path,
            mapq_criteria=20,
            max_shift=200,
            references=self.references,
            lengths=self.lengths
        )

        self.assertIsNotNone(worker)
        self.assertIsNotNone(worker._unified_worker)
        self.assertEqual(worker._unified_worker.config.calculation_config.target, CalculationTarget.NCC)
        self.assertEqual(worker._unified_worker.config.calculation_config.implementation, ImplementationAlgorithm.SUCCESSIVE)

    def test_mscc_worker_compatibility(self):
        """Test MSCCCalcWorker compatibility wrapper."""
        worker = MSCCCalcWorker(
            self.order_queue,
            self.report_queue,
            self.logger_lock,
            self.bam_path,
            mapq_criteria=20,
            max_shift=200,
            references=self.references,
            lengths=self.lengths,
            mappable_path="/test/mappability.bw",
            read_len=50,
            chrom2mappable_len={}
        )

        self.assertIsNotNone(worker)
        self.assertIsNotNone(worker._unified_worker)
        self.assertEqual(worker._unified_worker.config.calculation_config.target, CalculationTarget.MSCC)
        self.assertEqual(worker._unified_worker.config.calculation_config.implementation, ImplementationAlgorithm.SUCCESSIVE)
        self.assertTrue(worker._unified_worker.config.calculation_config.skip_ncc)

    def test_ncc_and_mscc_worker_compatibility(self):
        """Test NCCandMSCCCalcWorker compatibility wrapper."""
        worker = NCCandMSCCCalcWorker(
            self.order_queue,
            self.report_queue,
            self.logger_lock,
            self.bam_path,
            mapq_criteria=20,
            max_shift=200,
            references=self.references,
            lengths=self.lengths,
            mappable_path="/test/mappability.bw",
            read_len=50,
            chrom2mappable_len={}
        )

        self.assertIsNotNone(worker)
        self.assertIsNotNone(worker._unified_worker)
        self.assertEqual(worker._unified_worker.config.calculation_config.target, CalculationTarget.BOTH)
        self.assertEqual(worker._unified_worker.config.calculation_config.implementation, ImplementationAlgorithm.SUCCESSIVE)
        self.assertFalse(worker._unified_worker.config.calculation_config.skip_ncc)

    def test_ba_calc_worker_compatibility(self):
        """Test BACalcWorker compatibility wrapper."""
        worker = BACalcWorker(
            self.order_queue,
            self.report_queue,
            self.logger_lock,
            self.bam_path,
            mapq_criteria=20,
            max_shift=200,
            read_len=50,
            references=self.references,
            lengths=self.lengths,
            mappability_handler=None,
            skip_ncc=False
        )

        self.assertIsNotNone(worker)
        self.assertIsNotNone(worker._unified_worker)
        self.assertEqual(worker._unified_worker.config.calculation_config.target, CalculationTarget.BOTH)
        self.assertEqual(worker._unified_worker.config.calculation_config.implementation, ImplementationAlgorithm.BITARRAY)
        self.assertFalse(worker._unified_worker.config.calculation_config.skip_ncc)


if __name__ == '__main__':
    unittest.main()