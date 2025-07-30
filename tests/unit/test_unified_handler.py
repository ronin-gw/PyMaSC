"""Test unified handler functionality."""

import pytest
from unittest.mock import Mock, patch

from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig,
    CalculationTarget, ImplementationAlgorithm, ExecutionMode
)
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
        calc_config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20,
            skip_ncc=False
        )

        exec_config = ExecutionConfig(
            mode=ExecutionMode.SINGLE_PROCESS,
            worker_count=1
        )

        # Initialize handler
        handler = CalcHandler(
            path="test.bam",
            config=calc_config,
            execution_config=exec_config
        )

        assert handler.path == "test.bam"
        assert handler.config == calc_config
        assert handler.execution_config == exec_config
        assert handler.references == references
        assert handler.lengths == lengths

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
        calc_config = CalculationConfig(
            target=CalculationTarget.MSCC,
            implementation=ImplementationAlgorithm.BITARRAY,
            max_shift=300,
            mapq_criteria=30,
            skip_ncc=True
        )

        mappability_config = MappabilityConfig(
            mappability_path="/path/to/mappability.bw",
            read_len=100,
            shift_size=300
        )

        # Initialize handler
        handler = CalcHandler(
            path="test.bam",
            config=calc_config,
            execution_config=ExecutionConfig(),
            mappability_config=mappability_config
        )

        assert handler.mappability_config == mappability_config

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

        # Request multiprocess
        calc_config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20
        )

        exec_config = ExecutionConfig(
            mode=ExecutionMode.MULTI_PROCESS,
            worker_count=4
        )

        handler = CalcHandler(
            path="test.bam",
            config=calc_config,
            execution_config=exec_config
        )

        # Should fall back to single process
        assert handler.execution_config.mode == ExecutionMode.SINGLE_PROCESS
        assert handler.execution_config.worker_count == 1

    @patch('pysam.AlignmentFile')
    def test_unified_handler_with_empty_bam(self, mock_alignment_file):
        """Test handler with empty BAM file."""
        from PyMaSC.handler.calc import CalcHandler

        # Mock pysam to raise ValueError for empty file
        mock_alignment_file.side_effect = ValueError("File has no sequences defined.")

        calc_config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20
        )

        exec_config = ExecutionConfig(
            mode=ExecutionMode.SINGLE_PROCESS,
            worker_count=1
        )

        # Should raise ValueError for empty file
        with pytest.raises(ValueError, match="File has no sequences"):
            CalcHandler(
                path="empty.bam",
                config=calc_config,
                execution_config=exec_config
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
        calc_config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20
        )
        calc_config.chromfilter = [(True, ['chr99'])]  # Non-existent chromosome

        exec_config = ExecutionConfig(
            mode=ExecutionMode.SINGLE_PROCESS,
            worker_count=1
        )

        # Should raise NothingToCalc
        with pytest.raises(NothingToCalc):
            CalcHandler(
                path="test.bam",
                config=calc_config,
                execution_config=exec_config
            )

    def test_unified_handler_set_or_estimate_readlen(self):
        """Test read length setting."""
        from PyMaSC.handler.calc import CalcHandler

        # Create mock handler
        handler = Mock(spec=CalcHandler)
        handler.read_len = None
        handler.config = Mock(max_shift=500)

        # Manually set read length
        CalcHandler.set_or_estimate_readlen(handler, readlen=100)
        assert handler.read_len == 100

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

        calc_config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=250,
            mapq_criteria=25,
            skip_ncc=True
        )
        calc_config.esttype = 'mean'

        exec_config = ExecutionConfig(
            mode=ExecutionMode.SINGLE_PROCESS,
            worker_count=2
        )

        handler = CalcHandler(
            path="test.bam",
            config=calc_config,
            execution_config=exec_config
        )

        # Check backward compatible properties
        assert handler.config.skip_ncc is True
        assert handler.config.max_shift == 250
        assert handler.config.mapq_criteria == 25
        assert handler.execution_config.worker_count == 2
        assert handler.config.esttype == 'mean'



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