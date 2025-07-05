"""Test unified handler functionality."""

import pytest
from unittest.mock import Mock, patch, MagicMock

from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig,
    AlgorithmType, ExecutionMode
)
from tests.utils.test_data_generator import create_mock_reference_data


class TestUnifiedHandler:
    """Test UnifiedCalcHandler functionality."""

    def test_unified_handler_import(self):
        """Test that unified handler imports correctly."""
        from PyMaSC.handler.unified import UnifiedCalcHandler
        assert UnifiedCalcHandler is not None

    @patch('pysam.AlignmentFile')
    def test_unified_handler_basic_initialization(self, mock_alignment_file):
        """Test basic initialization of unified handler."""
        from PyMaSC.handler.unified import UnifiedCalcHandler

        # Mock AlignmentFile
        mock_bam = Mock()
        references, lengths = create_mock_reference_data()
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        # Create configuration
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20,
            skip_ncc=False
        )

        exec_config = ExecutionConfig(
            mode=ExecutionMode.SINGLE_PROCESS,
            worker_count=1
        )

        # Initialize handler
        handler = UnifiedCalcHandler(
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
        from PyMaSC.handler.unified import UnifiedCalcHandler

        # Mock AlignmentFile
        mock_bam = Mock()
        references, lengths = create_mock_reference_data()
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        # Create configuration with mappability
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
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
        handler = UnifiedCalcHandler(
            path="test.bam",
            config=calc_config,
            execution_config=ExecutionConfig(),
            mappability_config=mappability_config
        )

        assert handler.mappability_config == mappability_config

    @patch('pysam.AlignmentFile')
    def test_unified_handler_multiprocess_fallback(self, mock_alignment_file):
        """Test that multiprocess falls back to single process without index."""
        from PyMaSC.handler.unified import UnifiedCalcHandler

        # Mock AlignmentFile without index
        mock_bam = Mock()
        references, lengths = create_mock_reference_data()
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        # Request multiprocess
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20
        )

        exec_config = ExecutionConfig(
            mode=ExecutionMode.MULTI_PROCESS,
            worker_count=4
        )

        handler = UnifiedCalcHandler(
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
        from PyMaSC.handler.unified import UnifiedCalcHandler

        # Mock pysam to raise ValueError for empty file
        mock_alignment_file.side_effect = ValueError("File has no sequences defined.")

        calc_config = CalculationConfig(
            algorithm=AlgorithmType.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20
        )

        # Should raise ValueError for empty file
        with pytest.raises(ValueError, match="File has no sequences"):
            UnifiedCalcHandler(
                path="empty.bam",
                config=calc_config
            )

    @patch('pysam.AlignmentFile')
    def test_unified_handler_chromosome_filtering(self, mock_alignment_file):
        """Test chromosome filtering in unified handler."""
        from PyMaSC.handler.unified import UnifiedCalcHandler, NothingToCalc

        # Mock AlignmentFile
        mock_bam = Mock()
        mock_bam.references = ['chr1', 'chr2', 'chrM']
        mock_bam.lengths = [100000, 90000, 16000]
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        # Create config with chromfilter
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20
        )
        calc_config.chromfilter = [(True, ['chr99'])]  # Non-existent chromosome

        # Should raise NothingToCalc
        with pytest.raises(NothingToCalc):
            UnifiedCalcHandler(
                path="test.bam",
                config=calc_config
            )

    def test_unified_handler_set_readlen(self):
        """Test read length setting."""
        from PyMaSC.handler.unified import UnifiedCalcHandler

        # Create mock handler
        handler = Mock(spec=UnifiedCalcHandler)
        handler.read_len = None
        handler.config = Mock(max_shift=500)

        # Manually set read length
        UnifiedCalcHandler.set_readlen(handler, readlen=100)
        assert handler.read_len == 100

    @patch('pysam.AlignmentFile')
    def test_unified_handler_backward_compatible_properties(self, mock_alignment_file):
        """Test backward compatible properties."""
        from PyMaSC.handler.unified import UnifiedCalcHandler

        # Mock AlignmentFile
        mock_bam = Mock()
        references, lengths = create_mock_reference_data()
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        calc_config = CalculationConfig(
            algorithm=AlgorithmType.SUCCESSIVE,
            max_shift=250,
            mapq_criteria=25,
            skip_ncc=True
        )
        calc_config.esttype = 'mean'

        exec_config = ExecutionConfig(
            mode=ExecutionMode.SINGLE_PROCESS,
            worker_count=2
        )

        handler = UnifiedCalcHandler(
            path="test.bam",
            config=calc_config,
            execution_config=exec_config
        )

        # Check backward compatible properties
        assert handler.skip_ncc == True
        assert handler.max_shift == 250
        assert handler.mapq_criteria == 25
        assert handler.nworker == 2
        assert handler.esttype == 'mean'



class TestHandlerIntegration:
    """Test integration between unified handler and compatibility wrappers."""

    @patch('PyMaSC.handler.unified.CalculationContext')
    @patch('pysam.AlignmentFile')
    def test_unified_handler_strategy_selection(self, mock_alignment_file, mock_context):
        """Test that unified handler selects correct strategy."""
        from PyMaSC.handler.unified import UnifiedCalcHandler

        # Mock AlignmentFile
        mock_bam = Mock()
        references, lengths = create_mock_reference_data()
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        # Test different algorithms
        algorithms = [
            AlgorithmType.SUCCESSIVE,
            AlgorithmType.BITARRAY,
            AlgorithmType.MSCC
        ]

        for algo in algorithms:
            calc_config = CalculationConfig(
                algorithm=algo,
                max_shift=200,
                mapq_criteria=20
            )

            # Reset mock
            mock_context.create_from_config.reset_mock()

            handler = UnifiedCalcHandler(
                path="test.bam",
                config=calc_config
            )

            # Should have created context with config
            mock_context.create_from_config.assert_called_once_with(calc_config)

    @patch('pysam.AlignmentFile')
    def test_handler_result_collection(self, mock_alignment_file):
        """Test result collection in unified handler."""
        from PyMaSC.handler.unified import UnifiedCalcHandler

        # Mock AlignmentFile
        mock_bam = Mock()
        references = ['chr1']
        lengths = [100000]
        mock_bam.references = references
        mock_bam.lengths = lengths
        mock_bam.has_index.return_value = False
        mock_alignment_file.return_value = mock_bam

        calc_config = CalculationConfig(
            algorithm=AlgorithmType.SUCCESSIVE,
            max_shift=200,
            mapq_criteria=20
        )

        handler = UnifiedCalcHandler(
            path="test.bam",
            config=calc_config
        )

        # Test result storage initialization
        assert hasattr(handler, 'ref2forward_sum')
        assert hasattr(handler, 'ref2reverse_sum')
        assert hasattr(handler, 'ref2ccbins')
        assert hasattr(handler, 'mappable_ref2forward_sum')
        assert hasattr(handler, 'mappable_ref2reverse_sum')
        assert hasattr(handler, 'mappable_ref2ccbins')
        assert hasattr(handler, 'ref2mappable_len')

    def test_handler_exceptions(self):
        """Test custom exceptions are available."""
        from PyMaSC.handler.unified import InputUnseekable, NothingToCalc

        assert issubclass(InputUnseekable, Exception)
        assert issubclass(NothingToCalc, Exception)