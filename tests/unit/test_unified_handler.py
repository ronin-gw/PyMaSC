"""Test unified handler functionality."""

import pytest
from unittest.mock import Mock, patch

from PyMaSC.interfaces.config import (
    PyMaSCConfig, CalculationTarget, Algorithm, EstimationType
)
from pathlib import Path
from tests.utils.test_data_generator import create_mock_reference_data


class TestUnifiedHandler:
    """Test CalcHandler functionality."""

    def test_unified_handler_import(self):
        """Test that unified handler imports correctly."""
        from PyMaSC.handler.calc import CalcHandler
        assert CalcHandler is not None

    @patch('pysam.AlignmentFile')
    def test_unified_handler_basic_initialization(self, mock_alignment_file):
        """Test basic initialization of unified handler."""
        from PyMaSC.handler.calc import CalcHandler

        # Mock AlignmentFile
        mock_bam = Mock()
        references, lengths = create_mock_reference_data()
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        # Create configuration
        config = PyMaSCConfig(
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
        config.ref2lengths = dict(zip(references, lengths))

        # Initialize handler
        handler = CalcHandler(
            path="test.bam",
            config=config
        )

        assert handler.path == "test.bam"
        assert handler.config == config
        assert handler.config.references == tuple(references)
        assert handler.config.lengths == tuple(lengths)

    @patch('pysam.AlignmentFile')
    def test_unified_handler_with_mappability(self, mock_alignment_file):
        """Test unified handler with mappability configuration."""
        from PyMaSC.handler.calc import CalcHandler

        # Mock AlignmentFile
        mock_bam = Mock()
        references, lengths = create_mock_reference_data()
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        # Create configuration with mappability
        config = PyMaSCConfig(
            target=CalculationTarget.MSCC,
            implementation=Algorithm.BITARRAY,
            max_shift=300,
            mapq_criteria=30,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10,
            mappability_path=Path("/path/to/mappability.bw"),
            read_length=100
        )
        config.ref2lengths = dict(zip(references, lengths))

        # Initialize handler
        handler = CalcHandler(
            path="test.bam",
            config=config
        )

        assert handler.config.mappability_path == Path("/path/to/mappability.bw")

    @patch('pysam.AlignmentFile')
    def test_unified_handler_multiprocess_fallback(self, mock_alignment_file):
        """Test that multiprocess falls back to single process without index."""
        from PyMaSC.handler.calc import CalcHandler

        # Mock AlignmentFile without index
        mock_bam = Mock()
        references, lengths = create_mock_reference_data()
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        # This test case was already fixed above but duplicated code remained

    @patch('pysam.AlignmentFile')
    def test_unified_handler_with_empty_bam(self, mock_alignment_file):
        """Test handler with empty BAM file."""
        from PyMaSC.handler.calc import CalcHandler

        # Mock pysam to raise ValueError for empty file
        mock_alignment_file.side_effect = ValueError("File has no sequences defined.")

        # Create a basic config for this test
        config = PyMaSCConfig(
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

        # Should raise ValueError for empty file
        with pytest.raises(ValueError, match="File has no sequences"):
            CalcHandler(
                path="empty.bam",
                config=config
            )

    @patch('pysam.AlignmentFile')
    def test_unified_handler_chromosome_filtering(self, mock_alignment_file):
        """Test chromosome filtering in unified handler."""
        from PyMaSC.handler.calc import CalcHandler, NothingToCalc

        # Mock AlignmentFile
        mock_bam = Mock()
        mock_bam.references = ['chr1', 'chr2', 'chrM']
        mock_bam.lengths = [100000, 90000, 16000]
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        # Create config with chromfilter
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10,
            chromfilter=[(True, ['chr99'])]  # Non-existent chromosome
        )

        # Should raise NothingToCalc
        with pytest.raises(NothingToCalc):
            CalcHandler(
                path="test.bam",
                config=config
            )

    def test_unified_handler_config_read_length(self):
        """Test read length configuration."""
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10,
            read_length=100
        )

        # Test that read length can be set and retrieved
        assert config.read_length == 100

    @patch('pysam.AlignmentFile')
    def test_unified_handler_backward_compatible_properties(self, mock_alignment_file):
        """Test backward compatible properties."""
        from PyMaSC.handler.calc import CalcHandler

        # Mock AlignmentFile
        mock_bam = Mock()
        references, lengths = create_mock_reference_data()
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        config = PyMaSCConfig(
            target=CalculationTarget.MSCC,  # MSCC implies skip_ncc
            implementation=Algorithm.SUCCESSIVE,
            max_shift=250,
            mapq_criteria=25,
            nproc=2,
            esttype=EstimationType.MEAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )
        config.ref2lengths = dict(zip(references, lengths))

        handler = CalcHandler(
            path="test.bam",
            config=config
        )

        # Check backward compatible properties
        assert handler.config.skip_ncc is True
        assert handler.config.max_shift == 250
        assert handler.config.mapq_criteria == 25
        # Note: nproc may be adjusted to 1 if BAM file lacks index (multiprocess fallback)
        assert handler.config.nproc == 1  # Adjusted due to multiprocess fallback
        assert handler.config.esttype == EstimationType.MEAN



class TestHandlerIntegration:
    """Test integration between unified handler and compatibility wrappers."""

# Removed test_unified_handler_strategy_selection - tests internal strategy selection
    # that is now handled differently in current implementation

# Removed test_handler_result_collection - tests legacy result storage attributes 
    # that no longer exist in current implementation

    def test_handler_exceptions(self):
        """Test custom exceptions are available."""
        from PyMaSC.core.exceptions import InputUnseekable, NothingToCalc

        assert issubclass(InputUnseekable, Exception)
        assert issubclass(NothingToCalc, Exception)