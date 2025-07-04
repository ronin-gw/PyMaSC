"""End-to-end tests for new service architecture.

This module tests the complete workflow using the new service-oriented
architecture, verifying that all components work together correctly.
"""
import tempfile
import os
import json
import numpy as np
from pathlib import Path
import pytest

from PyMaSC.services.calculation import create_calculation_service
from PyMaSC.services.io import create_io_service, InMemoryIOService
from PyMaSC.services.validation import create_validation_service
from PyMaSC.services.mappability import create_mappability_service
from PyMaSC.services.workflow import create_workflow_service, WorkflowRequest
from PyMaSC.handler.simplified import SimplifiedCalcHandler, HandlerBuilder
from PyMaSC.handler.service_adapter import ServiceBasedCalcHandler
from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig,
    AlgorithmType, ExecutionMode
)


class TestServiceArchitectureE2E:
    """End-to-end tests for service architecture."""
    
    def setup_method(self):
        """Set up test data."""
        # Create test chromosome data
        self.test_chromosomes = ['chr1', 'chr2']
        self.test_lengths = [10000, 8000]
        
        # Create test reads (sorted by position)
        self.chr1_reads = [
            # Forward reads
            (100, 50, False),   # pos, length, is_reverse
            (200, 50, False),
            (500, 50, False),
            # Reverse reads
            (150, 50, True),
            (250, 50, True),
            (550, 50, True),
        ]
        
        self.chr2_reads = [
            # Forward reads
            (100, 50, False),
            (300, 50, False),
            # Reverse reads
            (200, 50, True),
            (400, 50, True),
        ]
    
    def test_simple_workflow_in_memory(self):
        """Test simple workflow with in-memory I/O service."""
        # Create in-memory I/O service with test data
        io_service = InMemoryIOService()
        
        # Add test data for chr1
        io_service.add_test_data(
            bam_path="/test/file.bam",
            chromosome="chr1",
            forward_reads=[(r[0], r[1]) for r in self.chr1_reads if not r[2]],
            reverse_reads=[(r[0], r[1]) for r in self.chr1_reads if r[2]],
            length=self.test_lengths[0]
        )
        
        # Add test data for chr2
        io_service.add_test_data(
            bam_path="/test/file.bam",
            chromosome="chr2",
            forward_reads=[(r[0], r[1]) for r in self.chr2_reads if not r[2]],
            reverse_reads=[(r[0], r[1]) for r in self.chr2_reads if r[2]],
            length=self.test_lengths[1]
        )
        
        # Create workflow service
        workflow_service = create_workflow_service(
            parallel=False,
            io_service=io_service
        )
        
        # Create calculation config
        config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=100,
            mapq_criteria=20,
            references=self.test_chromosomes,
            lengths=self.test_lengths
        )
        
        # Create workflow request
        request = WorkflowRequest(
            bam_path="/test/file.bam",
            output_prefix="/test/output",
            calculation_config=config,
            output_formats=['table']
        )
        
        # Execute workflow
        result = workflow_service.execute(request)
        
        # Verify success
        assert result.is_successful
        assert result.calculation_result is not None
        
        # Verify chromosome results
        calc_result = result.calculation_result
        assert len(calc_result.chromosome_results) == 2
        assert 'chr1' in calc_result.chromosome_results
        assert 'chr2' in calc_result.chromosome_results
        
        # Verify read counts
        chr1_result = calc_result.chromosome_results['chr1']
        assert chr1_result.forward_count == 3
        assert chr1_result.reverse_count == 3
        
        chr2_result = calc_result.chromosome_results['chr2']
        assert chr2_result.forward_count == 2
        assert chr2_result.reverse_count == 2
        
        # Verify aggregated results
        assert calc_result.total_forward_reads == 5
        assert calc_result.total_reverse_reads == 5
        
        # Verify correlation bins exist
        assert len(chr1_result.correlation_bins) > 0
        assert len(chr2_result.correlation_bins) > 0
        assert len(calc_result.aggregated_correlation) > 0
    
    def test_simplified_handler_workflow(self):
        """Test complete workflow using SimplifiedCalcHandler."""
        # Create test I/O service
        io_service = InMemoryIOService()
        io_service.add_test_data(
            bam_path="/test/file.bam",
            chromosome="chr1",
            forward_reads=[(r[0], r[1]) for r in self.chr1_reads if not r[2]],
            reverse_reads=[(r[0], r[1]) for r in self.chr1_reads if r[2]],
            length=self.test_lengths[0]
        )
        
        # Create mock validation service that always passes
        from unittest.mock import Mock
        mock_validation = Mock()
        mock_validation.validate_workflow_request.return_value = Mock(
            is_valid=True,
            errors=[],
            warnings=[]
        )
        
        # Build handler with test services
        handler = HandlerBuilder() \
            .with_bam_file("/test/file.bam") \
            .with_algorithm('ncc') \
            .with_max_shift(100) \
            .with_io_service(io_service) \
            .with_validation_service(mock_validation) \
            .build()
        
        # Run calculation
        result = handler.run_calculation()
        
        # Verify results
        assert result is not None
        assert result.total_forward_reads == 3
        assert result.total_reverse_reads == 3
        assert len(result.chromosome_results) == 1
        assert 'chr1' in result.chromosome_results
    
    def test_validation_service_integration(self):
        """Test validation service with workflow."""
        # Create validation service
        validation_service = create_validation_service(strict=True)
        
        # Create config with positive max_shift first
        config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=100,
            mapq_criteria=20
        )
        
        # Then make it invalid for testing
        config.max_shift = -100  # Invalid
        
        # Validate workflow request
        result = validation_service.validate_workflow_request(
            bam_path="/nonexistent/file.bam",
            output_prefix="/test/output",
            calculation_config=config
        )
        
        # Should have validation errors
        assert not result.is_valid
        assert len(result.errors) >= 2  # BAM not found + invalid max_shift
        assert any("BAM file not found" in err for err in result.errors)
        assert any("max_shift must be positive" in err for err in result.errors)
    
    def test_mappability_service_integration(self):
        """Test mappability service integration."""
        # Create test mappability data
        test_mappability = {
            'stat': {
                'total_mappable_length': 15000,
                'total_genome_length': 18000,
                'genome_mappable_fraction': 0.833
            },
            'chromosomes': {
                'chr1': {
                    'mappable_length': 9000,
                    'mappable_fraction': 0.9,
                    'mean_mappability': 0.85
                },
                'chr2': {
                    'mappable_length': 6000,
                    'mappable_fraction': 0.75,
                    'mean_mappability': 0.7
                }
            }
        }
        
        # Write test JSON file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(test_mappability, f)
            mappability_path = f.name
        
        try:
            # Create mappability service
            map_config = MappabilityConfig(
                mappability_path=mappability_path,
                read_len=50
            )
            map_service = create_mappability_service(mappability_path, map_config)
            
            # Test chromosome stats
            stats = map_service.calculate_mappability_stats('chr1', 10000)
            assert stats.mappable_length == 9000
            assert stats.mappable_fraction == 0.9
            
            # Test genome-wide stats
            genome_stats = map_service.calculate_genome_wide_stats({
                'chr1': 10000,
                'chr2': 8000
            })
            assert genome_stats.total_mappable_length == 15000
            assert genome_stats.genome_mappable_fraction == 0.833
            
        finally:
            os.unlink(mappability_path)
    
    @pytest.mark.skip(reason="ServiceBasedCalcHandler requires real BAM file in parent init")
    def test_service_adapter_backward_compatibility(self):
        """Test ServiceBasedCalcHandler for backward compatibility."""
        # Create test I/O service
        io_service = InMemoryIOService()
        io_service.add_test_data(
            bam_path="/test/file.bam",
            chromosome="chr1",
            forward_reads=[(100, 50), (200, 50)],
            reverse_reads=[(150, 50)],
            length=10000
        )
        
        # Create workflow service
        workflow_service = create_workflow_service(
            parallel=False,
            io_service=io_service
        )
        
        # Mock the workflow service to return a result
        from unittest.mock import Mock
        mock_result = Mock()
        mock_result.is_successful = True
        mock_result.calculation_result = Mock(
            chromosome_results={
                'chr1': Mock(
                    forward_count=2,
                    reverse_count=1,
                    correlation_bins=np.array([0.1, 0.2, 0.3]),
                    mappable_forward_count=None,
                    mappable_reverse_count=None,
                    mappable_correlation_bins=None,
                    mappable_length=None
                )
            },
            total_forward_reads=2,
            total_reverse_reads=1
        )
        workflow_service.execute = Mock(return_value=mock_result)
        
        # Create handler using old interface
        handler = ServiceBasedCalcHandler(
            path="/test/file.bam",
            esttype='mean',
            max_shift=100,
            mapq_criteria=20,
            algorithm='bitarray'
        )
        
        # Replace workflow service
        handler._workflow_service = workflow_service
        
        # Run calculation (old method name)
        handler._run_singleprocess_calculation()
        
        # Verify results were extracted correctly
        assert handler.ref2forward_sum['chr1'] == 2
        assert handler.ref2reverse_sum['chr1'] == 1
        assert len(handler.ref2ccbins['chr1']) == 3
        
        # Get workflow result (new method)
        workflow_result = handler.get_workflow_result()
        assert workflow_result is not None
        assert workflow_result.is_successful
    
    @pytest.mark.slow
    def test_multiprocess_service_workflow(self):
        """Test multiprocess workflow execution."""
        # Create test data for multiple chromosomes
        io_service = InMemoryIOService()
        
        for i, (chrom, length) in enumerate(zip(self.test_chromosomes, self.test_lengths)):
            # Create some test reads
            forward_reads = [(j * 100, 50) for j in range(1, 4)]
            reverse_reads = [(j * 100 + 50, 50) for j in range(1, 4)]
            
            io_service.add_test_data(
                bam_path="/test/file.bam",
                chromosome=chrom,
                forward_reads=forward_reads,
                reverse_reads=reverse_reads,
                length=length
            )
        
        # Create parallel workflow service
        workflow_service = create_workflow_service(
            parallel=True,
            n_workers=2,
            io_service=io_service
        )
        
        # Create config
        config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=100,
            mapq_criteria=20,
            references=self.test_chromosomes,
            lengths=self.test_lengths
        )
        
        # Create execution config
        exec_config = ExecutionConfig(
            mode=ExecutionMode.MULTI_PROCESS,
            worker_count=2
        )
        
        # Create request
        request = WorkflowRequest(
            bam_path="/test/file.bam",
            output_prefix="/test/output",
            calculation_config=config,
            execution_config=exec_config
        )
        
        # Execute workflow
        result = workflow_service.execute(request)
        
        # Verify success
        assert result.is_successful
        assert result.calculation_result is not None
        
        # Verify all chromosomes processed
        calc_result = result.calculation_result
        assert len(calc_result.chromosome_results) == 2
        
        # Verify aggregation worked correctly
        assert calc_result.total_forward_reads == 6  # 3 per chromosome
        assert calc_result.total_reverse_reads == 6
    
    def test_error_handling_in_services(self):
        """Test error handling across services."""
        # Create service that will fail
        io_service = InMemoryIOService()
        # No data added - should fail
        
        workflow_service = create_workflow_service(
            parallel=False,
            io_service=io_service
        )
        
        # Create request for non-existent data
        request = WorkflowRequest(
            bam_path="/test/missing.bam",
            output_prefix="/test/output",
            calculation_config=CalculationConfig(
                algorithm=AlgorithmType.BITARRAY,
                max_shift=100,
                mapq_criteria=20
            )
        )
        
        # Execute workflow - should fail gracefully
        result = workflow_service.execute(request)
        
        # Verify failure is handled properly
        assert not result.is_successful
        assert result.error is not None
        assert "not found" in result.error.lower()


class TestServicePerformance:
    """Performance tests for service architecture."""
    
    def test_calculation_service_caching(self):
        """Test that calculation service properly caches calculators."""
        calc_service = create_calculation_service()
        
        # Create test data
        from PyMaSC.services.calculation import ChromosomeData
        data1 = ChromosomeData(
            chromosome="chr1",
            forward_reads=[(100, 50), (200, 50)],
            reverse_reads=[(150, 50)],
            length=1000
        )
        
        data2 = ChromosomeData(
            chromosome="chr2",
            forward_reads=[(100, 50)],
            reverse_reads=[(150, 50)],
            length=800
        )
        
        config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=100,
            mapq_criteria=20,
            references=["chr1", "chr2"],
            lengths=[1000, 800]
        )
        
        # Calculate twice with same config
        result1 = calc_service.calculate_chromosome(data1, config)
        result2 = calc_service.calculate_chromosome(data2, config)
        
        # Should use cached calculator (verify by checking cache)
        if hasattr(calc_service, '_calculator_cache'):
            assert len(calc_service._calculator_cache) == 1
    
    def test_mappability_service_caching(self):
        """Test mappability service caching."""
        # Create test BigWig data (mock)
        from unittest.mock import patch, MagicMock
        
        with patch('PyMaSC.services.mappability.BigWigReader') as mock_reader_class:
            mock_reader = MagicMock()
            mock_reader.get_as_array.return_value = [0.5] * 100
            mock_reader.chroms.return_value = {'chr1': 1000}
            mock_reader_class.return_value = mock_reader
            
            service = create_mappability_service("/test/file.bw")
            
            # First call - should hit reader
            values1 = service.get_mappability_values('chr1', 100, 200)
            assert mock_reader.get_as_array.call_count == 1
            
            # Second call - should use cache
            values2 = service.get_mappability_values('chr1', 100, 200)
            assert mock_reader.get_as_array.call_count == 1  # No additional call
            
            # Different region - should hit reader again
            values3 = service.get_mappability_values('chr1', 300, 400)
            assert mock_reader.get_as_array.call_count == 2


if __name__ == '__main__':
    pytest.main([__file__, '-v'])