"""Abstract interfaces and protocols for PyMaSC architecture.

This module defines the core interfaces and protocols that provide a unified
foundation for the PyMaSC refactored architecture. These interfaces enable
polymorphism, dependency injection, and clear separation of concerns.

Key interfaces:
- CrossCorrelationCalculator: Abstract base for all correlation calculators
- ReadProcessor: Protocol for read processing operations
- ProgressReporter: Protocol for progress reporting
- ResultCollector: Protocol for result collection and aggregation
"""
from abc import ABC, abstractmethod
from typing import Protocol, Dict, Any, Optional, Iterator, Tuple
import numpy as np


class CrossCorrelationCalculator(ABC):
    """Abstract base class for all cross-correlation calculators.
    
    This interface defines the common contract that all correlation calculation
    implementations must follow, enabling polymorphic usage across different
    algorithms (NCC, MSCC, BitArray, etc.).
    
    All implementations must provide methods for:
    - Processing forward and reverse strand reads
    - Finalizing calculations
    - Flushing intermediate results
    - Accessing result data
    """
    
    @abstractmethod
    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a forward strand read for cross-correlation calculation.
        
        Args:
            chrom: Chromosome name
            pos: Read position (1-based)
            readlen: Read length
        """
        pass
    
    @abstractmethod
    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a reverse strand read for cross-correlation calculation.
        
        Args:
            chrom: Chromosome name
            pos: Read position (1-based)
            readlen: Read length
        """
        pass
    
    @abstractmethod
    def finishup_calculation(self) -> None:
        """Complete cross-correlation calculation for all chromosomes.
        
        This method should be called after all reads have been processed
        to finalize the correlation calculations.
        """
        pass
    
    @abstractmethod
    def flush(self) -> None:
        """Flush any intermediate calculation results.
        
        Used in multiprocessing scenarios to ensure results are available
        for collection before worker process termination.
        """
        pass
    
    @property
    @abstractmethod
    def ref2forward_sum(self) -> Dict[str, int]:
        """Forward read counts by chromosome."""
        pass
    
    @property
    @abstractmethod
    def ref2reverse_sum(self) -> Dict[str, int]:
        """Reverse read counts by chromosome."""
        pass
    
    @property
    @abstractmethod
    def ref2ccbins(self) -> Dict[str, np.ndarray]:
        """Cross-correlation bins by chromosome."""
        pass


class ReadProcessor(Protocol):
    """Protocol for read processing and filtering operations.
    
    Defines the interface for components that process sequencing reads,
    including filtering, validation, and information extraction.
    """
    
    def should_skip_read(self, read: Any, mapq_criteria: int) -> bool:
        """Determine if a read should be skipped based on quality criteria.
        
        Args:
            read: Alignment read object (typically pysam.AlignedSegment)
            mapq_criteria: Minimum mapping quality threshold
            
        Returns:
            True if read should be skipped, False otherwise
        """
        ...
    
    def extract_read_info(self, read: Any) -> Tuple[str, int, int, bool]:
        """Extract relevant information from a read.
        
        Args:
            read: Alignment read object
            
        Returns:
            Tuple of (chromosome, position, read_length, is_reverse)
        """
        ...


class ProgressReporter(Protocol):
    """Protocol for progress reporting and monitoring.
    
    Defines the interface for components that report calculation progress,
    enabling different progress reporting implementations (console, GUI, etc.).
    """
    
    def set_context(self, context: str, total: int) -> None:
        """Set the progress context and total work units.
        
        Args:
            context: Description of the current operation
            total: Total number of work units to complete
        """
        ...
    
    def update(self, current: int) -> None:
        """Update progress with current completion status.
        
        Args:
            current: Current number of completed work units
        """
        ...
    
    def finish(self) -> None:
        """Mark progress as complete and clean up resources."""
        ...


class ResultCollector(Protocol):
    """Protocol for collecting and aggregating calculation results.
    
    Defines the interface for components that collect results from
    workers and aggregate them into final output data structures.
    """
    
    def collect_chromosome_result(self, 
                                  chrom: str, 
                                  forward_sum: int,
                                  reverse_sum: int, 
                                  ccbins: np.ndarray,
                                  mappable_data: Optional[Dict[str, Any]] = None) -> None:
        """Collect results for a single chromosome.
        
        Args:
            chrom: Chromosome name
            forward_sum: Forward read count
            reverse_sum: Reverse read count
            ccbins: Cross-correlation bins
            mappable_data: Optional mappability-related data
        """
        ...
    
    def finalize_results(self) -> Dict[str, Any]:
        """Finalize and return aggregated results.
        
        Returns:
            Dictionary containing aggregated results across all chromosomes
        """
        ...


class CalculatorFactory(Protocol):
    """Protocol for calculator creation factories.
    
    Defines the interface for factories that create calculator instances
    based on configuration parameters.
    """
    
    def create_calculator(self, 
                          algorithm: str,
                          config: 'CalculationConfig') -> CrossCorrelationCalculator:
        """Create a calculator instance for the specified algorithm.
        
        Args:
            algorithm: Algorithm name ('ncc', 'mscc', 'bitarray')
            config: Configuration parameters
            
        Returns:
            Calculator instance implementing CrossCorrelationCalculator
            
        Raises:
            ValueError: If algorithm is not supported
        """
        ...


class WorkerFactory(Protocol):
    """Protocol for worker creation factories.
    
    Defines the interface for factories that create worker instances
    for multiprocessing calculation scenarios.
    """
    
    def create_worker(self, config: 'WorkerConfig') -> 'CalculationWorker':
        """Create a worker instance with the specified configuration.
        
        Args:
            config: Worker configuration parameters
            
        Returns:
            Worker instance ready for multiprocessing execution
        """
        ...


# Type aliases for common data structures
ChromosomeResults = Dict[str, Tuple[int, int, np.ndarray]]  # chrom -> (forward, reverse, ccbins)
CalculationResults = Dict[str, Any]  # Final aggregated results
ProgressCallback = callable  # Function for progress notifications