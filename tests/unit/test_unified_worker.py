"""Tests for unified worker architecture."""

import unittest
from unittest.mock import Mock
from multiprocessing import Queue, Lock

from PyMaSC.handler.worker import CalcWorker
from PyMaSC.interfaces.config import (
    PyMaSCConfig, CalculationTarget, Algorithm, EstimationType
)
from pathlib import Path
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
        self.config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )
        self.config.ref2lengths = {"chr1": 100000, "chr2": 90000}

    def test_unified_worker_initialization(self):
        """Test CalcWorker initialization."""
        worker = CalcWorker(
            self.config,
            self.order_queue,
            self.report_queue,
            self.logger_lock,
            self.bam_path
        )

        self.assertIsNotNone(worker)
        self.assertEqual(worker.config, self.config)
        self.assertEqual(worker.bam_path, self.bam_path)

    def test_unified_worker_with_mappability(self):
        """Test CalcWorker with mappability configuration."""
        # Add mappability configuration to the base config
        self.config.mappability_path = Path("/test/mappability.bw")
        self.config.read_length = 50

        config = PyMaSCConfig(
            target=CalculationTarget.MSCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10,
            read_length=50,
            mappability_path=Path("/test/mappability.bw")
        )
        config.ref2lengths = {"chr1": 100000, "chr2": 90000}

        worker = CalcWorker(
            config,
            self.order_queue,
            self.report_queue,
            self.logger_lock,
            self.bam_path
        )

        self.assertIsNotNone(worker)
        self.assertEqual(worker.config.mappability_path, Path("/test/mappability.bw"))


# Removed TestStandardReadProcessor and TestDualReadProcessor
# These processor classes have been refactored into a unified ReadProcessor
# architecture in PyMaSC.handler.read module


if __name__ == '__main__':
    unittest.main()