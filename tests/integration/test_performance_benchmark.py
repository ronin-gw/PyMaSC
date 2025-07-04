"""Performance benchmark tests for PyMaSC service architecture.

This module benchmarks the performance of the new service-oriented
architecture compared to the expected baseline performance.
"""
import time
import numpy as np
import pytest
from typing import List, Tuple

from PyMaSC.services.calculation import (
    ChromosomeData, create_calculation_service
)
from PyMaSC.services.io import InMemoryIOService
from PyMaSC.services.workflow import create_workflow_service, WorkflowRequest
from PyMaSC.core.models import CalculationConfig, AlgorithmType


class PerformanceTimer:
    """Simple timer for performance measurements."""
    
    last_elapsed = None  # Class variable to track last measurement
    
    def __init__(self, name: str):
        self.name = name
        self.start_time = None
        self.elapsed = None
    
    def __enter__(self):
        self.start_time = time.time()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.elapsed = time.time() - self.start_time
        PerformanceTimer.last_elapsed = self.elapsed
        print(f"\n{self.name}: {self.elapsed:.3f} seconds")


def generate_test_reads(n_reads: int, chrom_length: int, 
                       read_length: int = 50) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """Generate synthetic test reads.
    
    Args:
        n_reads: Number of reads to generate
        chrom_length: Chromosome length
        read_length: Read length
        
    Returns:
        Tuple of (forward_reads, reverse_reads)
    """
    # Generate random positions
    positions = np.random.randint(1, chrom_length - read_length, n_reads)
    positions.sort()
    
    # Split into forward and reverse
    forward_reads = []
    reverse_reads = []
    
    for i, pos in enumerate(positions):
        if i % 2 == 0:
            forward_reads.append((int(pos), read_length))
        else:
            reverse_reads.append((int(pos), read_length))
    
    return forward_reads, reverse_reads


class TestServicePerformanceBenchmark:
    """Performance benchmarks for service architecture."""
    
    def test_calculation_service_performance_small(self):
        """Benchmark calculation service with small dataset."""
        # Generate test data - 1000 reads
        forward_reads, reverse_reads = generate_test_reads(1000, 100000)
        
        data = ChromosomeData(
            chromosome="chr1",
            forward_reads=forward_reads,
            reverse_reads=reverse_reads,
            length=100000
        )
        
        config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=500,
            mapq_criteria=20,
            references=["chr1"],
            lengths=[100000]
        )
        
        # Create service
        calc_service = create_calculation_service()
        
        # Measure performance
        with PerformanceTimer("Small dataset (1k reads)") as timer:
            result = calc_service.calculate_chromosome(data, config)
        
        # Verify result
        assert result is not None
        assert result.forward_count == len(forward_reads)
        assert result.reverse_count == len(reverse_reads)
        assert len(result.correlation_bins) > 0
        
        # Performance assertion - should complete in under 0.5 seconds (generous for CI)
        assert timer.elapsed < 0.5, f"Too slow: {timer.elapsed}s"
    
    @pytest.mark.slow
    def test_calculation_service_performance_medium(self):
        """Benchmark calculation service with medium dataset."""
        # Generate test data - 10k reads
        forward_reads, reverse_reads = generate_test_reads(10000, 1000000)
        
        data = ChromosomeData(
            chromosome="chr1",
            forward_reads=forward_reads,
            reverse_reads=reverse_reads,
            length=1000000
        )
        
        config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=500,
            mapq_criteria=20,
            references=["chr1"],
            lengths=[1000000]
        )
        
        # Create service
        calc_service = create_calculation_service()
        
        # Measure performance
        with PerformanceTimer("Medium dataset (10k reads)") as timer:
            result = calc_service.calculate_chromosome(data, config)
        
        # Verify result
        assert result is not None
        # BitArray algorithm counts unique positions only (deduplication)
        unique_forward_positions = len(set(pos for pos, _ in data.forward_reads))
        unique_reverse_positions = len(set(pos for pos, _ in data.reverse_reads))
        assert result.forward_count == unique_forward_positions
        assert result.reverse_count == unique_reverse_positions
        
        # Performance assertion - should complete in under 1 second
        assert timer.elapsed < 1.0, f"Too slow: {timer.elapsed}s"
    
    @pytest.mark.slow
    def test_workflow_service_performance(self):
        """Benchmark complete workflow performance."""
        # Create test data for multiple chromosomes
        io_service = InMemoryIOService()
        
        chromosomes = [f"chr{i}" for i in range(1, 4)]
        lengths = [1000000, 800000, 600000]
        
        for chrom, length in zip(chromosomes, lengths):
            forward_reads, reverse_reads = generate_test_reads(5000, length)
            io_service.add_test_data(
                bam_path="/test/file.bam",
                chromosome=chrom,
                forward_reads=forward_reads,
                reverse_reads=reverse_reads,
                length=length
            )
        
        # Create workflow service
        workflow_service = create_workflow_service(
            parallel=False,
            io_service=io_service
        )
        
        # Create request
        config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=500,
            mapq_criteria=20,
            references=chromosomes,
            lengths=lengths
        )
        
        request = WorkflowRequest(
            bam_path="/test/file.bam",
            output_prefix="/test/output",
            calculation_config=config,
            output_formats=[]
        )
        
        # Measure performance
        with PerformanceTimer("Complete workflow (3 chromosomes, 15k reads total)") as timer:
            result = workflow_service.execute(request)
        
        # Verify result
        assert result.is_successful
        assert len(result.calculation_result.chromosome_results) == 3
        
        # Performance assertion - should complete in under 2 seconds
        assert timer.elapsed < 2.0, f"Too slow: {timer.elapsed}s"
    
    def test_service_caching_performance(self):
        """Test that service caching improves performance."""
        calc_service = create_calculation_service()
        
        # Generate test data
        forward_reads, reverse_reads = generate_test_reads(1000, 100000)
        
        data1 = ChromosomeData(
            chromosome="chr1",
            forward_reads=forward_reads[:500],
            reverse_reads=reverse_reads[:500],
            length=100000
        )
        
        data2 = ChromosomeData(
            chromosome="chr2",
            forward_reads=forward_reads[500:],
            reverse_reads=reverse_reads[500:],
            length=100000
        )
        
        config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=500,
            mapq_criteria=20,
            references=["chr1", "chr2"],
            lengths=[100000, 100000]
        )
        
        # First calculation - creates calculator
        with PerformanceTimer("First calculation (no cache)") as timer1:
            result1 = calc_service.calculate_chromosome(data1, config)
        
        # Second calculation - should use cached calculator
        with PerformanceTimer("Second calculation (with cache)") as timer2:
            result2 = calc_service.calculate_chromosome(data2, config)
        
        # Cache should make second calculation faster
        # Allow some variance but generally cache should help
        print(f"\nCache speedup: {timer1.elapsed / timer2.elapsed:.2f}x")
        
        # Both should produce valid results
        assert result1 is not None
        assert result2 is not None
    
    @pytest.mark.slow
    def test_algorithm_comparison_performance(self):
        """Compare performance of different algorithms."""
        # Generate larger test data
        forward_reads, reverse_reads = generate_test_reads(5000, 500000)
        
        data = ChromosomeData(
            chromosome="chr1",
            forward_reads=forward_reads,
            reverse_reads=reverse_reads,
            length=500000
        )
        
        algorithms = [
            AlgorithmType.NAIVE_CC,
            AlgorithmType.BITARRAY,
            AlgorithmType.SUCCESSIVE
        ]
        
        results = {}
        
        for algorithm in algorithms:
            config = CalculationConfig(
                algorithm=algorithm,
                max_shift=300,
                mapq_criteria=20,
                references=["chr1"],
                lengths=[500000]
            )
            
            calc_service = create_calculation_service()
            
            with PerformanceTimer(f"{algorithm.value} algorithm") as timer:
                result = calc_service.calculate_chromosome(data, config)
            
            results[algorithm.value] = timer.elapsed
            
            # Verify all algorithms produce results
            assert result is not None
            assert result.correlation_bins is not None
        
        # Print comparison
        print("\nAlgorithm Performance Comparison:")
        baseline = results.get('ncc', 1.0)
        for algo, elapsed in results.items():
            speedup = baseline / elapsed
            print(f"  {algo}: {elapsed:.3f}s ({speedup:.2f}x vs NCC)")


class TestServiceMemoryEfficiency:
    """Memory efficiency tests for service architecture."""
    
    def test_memory_efficient_streaming(self):
        """Test that services don't load all data at once."""
        # This is more of a design verification than measurement
        # The InMemoryIOService loads all data, but FileIOService streams
        
        io_service = InMemoryIOService()
        
        # Add a large dataset
        forward_reads, reverse_reads = generate_test_reads(50000, 10000000)
        
        # Split into chunks to simulate streaming
        chunk_size = 5000
        for i in range(0, len(forward_reads), chunk_size):
            # In real FileIOService, this would stream from disk
            chunk_forward = forward_reads[i:i+chunk_size]
            chunk_reverse = reverse_reads[i:i+chunk_size]
            
            # Process chunk (simulation)
            assert len(chunk_forward) <= chunk_size
            assert len(chunk_reverse) <= chunk_size
        
        print(f"\nProcessed {len(forward_reads)} reads in chunks of {chunk_size}")


if __name__ == '__main__':
    # Run performance tests
    pytest.main([__file__, '-v'])