"""Tests for unified worker architecture."""

import unittest
from unittest.mock import Mock
from multiprocessing import Queue, Lock

from PyMaSC.core.worker import CalcWorker
from PyMaSC.core.models import (
    WorkerConfig, CalculationConfig, MappabilityConfig,
    CalculationTarget, ImplementationAlgorithm
)
# Worker compatibility classes have been removed


class TestCalcWorker(unittest.TestCase):
    """Test cases for CalcWorker class."""

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
        """Test CalcWorker initialization."""
        worker = CalcWorker(
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
        """Test CalcWorker with mappability configuration."""
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

        worker = CalcWorker(
            worker_config,
            self.order_queue,
            self.report_queue,
            self.logger_lock,
            self.bam_path
        )

        self.assertIsNotNone(worker)
        self.assertEqual(worker.config.mappability_config, map_config)


# Removed TestStandardReadProcessor and TestDualReadProcessor
# These processor classes have been refactored into a unified ReadProcessor
# architecture in PyMaSC.handler.read module


if __name__ == '__main__':
    unittest.main()