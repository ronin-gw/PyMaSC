"""Calculation service for pure cross-correlation computation logic.

This module provides the CalculationService that encapsulates all
cross-correlation calculation logic, completely separated from I/O
operations and workflow management. This is the domain layer of
the application.

Key features:
- Pure calculation functions with no side effects
- Algorithm-agnostic interface
- Testable without external dependencies
- Immutable data processing
"""
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any

import numpy as np

from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, CalculationTarget, ImplementationAlgorithm
)
from PyMaSC.core.factory import CalculatorFactory
from PyMaSC.utils.stats_utils import (
    CrossCorrelationMetrics, ResultAggregator, calculate_quality_metrics
)

logger = logging.getLogger(__name__)


@dataclass
class ChromosomeData:
    """Input data for a single chromosome calculation.

    Immutable data structure containing all reads for one chromosome.
    """
    chromosome: str
    forward_reads: List[Tuple[int, int]]  # (position, length) pairs
    reverse_reads: List[Tuple[int, int]]  # (position, length) pairs
    length: int

    def __post_init__(self) -> None:
        """Validate chromosome data."""
        if self.length <= 0:
            raise ValueError(f"Invalid chromosome length: {self.length}")


@dataclass
class CalculationResult:
    """Result of cross-correlation calculation for one chromosome.

    Immutable result structure containing all calculation outputs.
    """
    chromosome: str
    forward_count: int
    reverse_count: int
    correlation_bins: np.ndarray
    mappable_forward_count: Optional[int] = None
    mappable_reverse_count: Optional[int] = None
    mappable_correlation_bins: Optional[np.ndarray] = None
    mappable_length: Optional[int] = None

    def get_quality_metrics(self, read_length: int) -> Dict[str, float]:
        """Calculate quality metrics for this result.

        Args:
            read_length: Read length for phantom peak detection

        Returns:
            Dictionary of quality metrics (NSC, RSC, etc.)
        """
        metrics = CrossCorrelationMetrics()

        # Find fragment peak with appropriate min_shift
        min_shift = min(50, len(self.correlation_bins) // 2) if len(self.correlation_bins) > 0 else 0
        fragment_pos, _ = metrics.find_fragment_peak(self.correlation_bins, min_shift=min_shift)

        # Calculate all metrics
        return calculate_quality_metrics(
            self.correlation_bins,
            fragment_pos,
            read_length
        )


@dataclass
class GenomeWideResult:
    """Aggregated result across all chromosomes.

    Contains both per-chromosome and genome-wide statistics.
    """
    chromosome_results: Dict[str, CalculationResult]
    total_forward_reads: int
    total_reverse_reads: int
    aggregated_correlation: np.ndarray
    mappable_fraction: Optional[float] = None

    @classmethod
    def from_chromosome_results(cls,
                               results: List[CalculationResult],
                               chromosome_lengths: Dict[str, int]) -> 'GenomeWideResult':
        """Create genome-wide result from chromosome results.

        Args:
            results: List of per-chromosome results
            chromosome_lengths: Total lengths of chromosomes

        Returns:
            Aggregated genome-wide result
        """
        # Build dictionaries for aggregation
        ref2forward = {}
        ref2reverse = {}
        ref2ccbins = {}
        ref2mappable_len = {}

        for result in results:
            ref2forward[result.chromosome] = result.forward_count
            ref2reverse[result.chromosome] = result.reverse_count
            ref2ccbins[result.chromosome] = result.correlation_bins
            if result.mappable_length is not None:
                ref2mappable_len[result.chromosome] = result.mappable_length

        # Use aggregator utilities
        aggregator = ResultAggregator()
        total_forward, total_reverse = aggregator.aggregate_read_counts(
            ref2forward, ref2reverse
        )
        aggregated_cc = aggregator.aggregate_ccbins(
            ref2ccbins, ref2forward, ref2reverse
        )

        # Calculate mappable fraction if available
        mappable_fraction = None
        if ref2mappable_len:
            references = list(chromosome_lengths.keys())
            lengths = list(chromosome_lengths.values())
            mappable_fraction = aggregator.calculate_mappable_fraction(
                ref2mappable_len, lengths, references
            )

        return cls(
            chromosome_results={r.chromosome: r for r in results},
            total_forward_reads=total_forward,
            total_reverse_reads=total_reverse,
            aggregated_correlation=aggregated_cc,
            mappable_fraction=mappable_fraction
        )


class CalculationService(ABC):
    """Abstract base for cross-correlation calculation services.

    Defines the interface for all calculation services, enabling
    different implementations while maintaining a consistent API.
    """

    @abstractmethod
    def calculate_chromosome(self,
                           data: ChromosomeData,
                           config: CalculationConfig,
                           mappability_config: Optional[MappabilityConfig] = None
                           ) -> CalculationResult:
        """Calculate cross-correlation for a single chromosome.

        Args:
            data: Chromosome read data
            config: Calculation configuration
            mappability_config: Optional mappability configuration

        Returns:
            Calculation result for the chromosome
        """
        pass

    @abstractmethod
    def calculate_genome_wide(self,
                            chromosome_data: List[ChromosomeData],
                            config: CalculationConfig,
                            mappability_config: Optional[MappabilityConfig] = None
                            ) -> GenomeWideResult:
        """Calculate cross-correlation for entire genome.

        Args:
            chromosome_data: List of chromosome data
            config: Calculation configuration
            mappability_config: Optional mappability configuration

        Returns:
            Genome-wide calculation result
        """
        pass


class StandardCalculationService(CalculationService):
    """Standard implementation of calculation service.

    Uses the existing PyMaSC calculators through the factory pattern,
    providing a clean service interface while leveraging proven
    calculation logic.
    """

    def __init__(self) -> None:
        """Initialize calculation service."""
        self._calculator_cache: Dict[Tuple[CalculationTarget, ImplementationAlgorithm, int, int, bool], Any] = {}

    def calculate_chromosome(self,
                           data: ChromosomeData,
                           config: CalculationConfig,
                           mappability_config: Optional[MappabilityConfig] = None
                           ) -> CalculationResult:
        """Calculate cross-correlation for a single chromosome.

        Performs pure calculation without any I/O operations.
        """
        # Create or reuse calculator
        cache_key = (config.target, config.implementation, config.max_shift,
                    config.read_length or 50, config.skip_ncc)

        if cache_key not in self._calculator_cache:
            # Create calculator for this configuration
            # Use all references from original config, not just current chromosome
            calc_config = CalculationConfig(
                target=config.target,
                implementation=config.implementation,
                max_shift=config.max_shift,
                mapq_criteria=config.mapq_criteria,
                references=config.references,  # Use all references
                lengths=config.lengths,        # Use all lengths
                read_length=config.read_length,
                skip_ncc=config.skip_ncc
            )

            calculator = CalculatorFactory.create_calculator(
                config.target,
                config.implementation,
                calc_config,
                mappability_config
            )
            self._calculator_cache[cache_key] = calculator
        else:
            calculator = self._calculator_cache[cache_key]

        # Feed reads to calculator in position order
        # Merge forward and reverse reads with strand information
        all_reads = []
        for pos, length in data.forward_reads:
            all_reads.append((pos, length, True))  # True = forward
        for pos, length in data.reverse_reads:
            all_reads.append((pos, length, False))  # False = reverse

        # Sort by position
        all_reads.sort(key=lambda x: x[0])

        # Feed sorted reads
        for pos, length, is_forward in all_reads:
            if is_forward:
                calculator.feed_forward_read(data.chromosome, pos, length)
            else:
                calculator.feed_reverse_read(data.chromosome, pos, length)

        # Finalize calculation
        calculator.finishup_calculation()

        # Extract results
        forward_count = calculator.ref2forward_sum.get(data.chromosome, 0)
        reverse_count = calculator.ref2reverse_sum.get(data.chromosome, 0)
        correlation_bins = calculator.ref2ccbins.get(data.chromosome, np.array([]))

        # Extract mappability results if available
        mappable_forward = None
        mappable_reverse = None
        mappable_bins = None
        mappable_length = None

        # Only extract mappability data if the attributes are real dictionaries
        if hasattr(calculator, 'ref2mappable_forward_sum') and isinstance(getattr(calculator, 'ref2mappable_forward_sum', None), dict):
            mappable_forward = calculator.ref2mappable_forward_sum.get(data.chromosome)
            if hasattr(calculator, 'ref2mappable_reverse_sum') and isinstance(getattr(calculator, 'ref2mappable_reverse_sum', None), dict):
                mappable_reverse = calculator.ref2mappable_reverse_sum.get(data.chromosome)

        if hasattr(calculator, 'ref2mascbins') and isinstance(getattr(calculator, 'ref2mascbins', None), dict):
            mappable_bins = calculator.ref2mascbins.get(data.chromosome)

        if hasattr(calculator, 'ref2mappable_len') and isinstance(getattr(calculator, 'ref2mappable_len', None), dict):
            mappable_length = calculator.ref2mappable_len.get(data.chromosome)

        return CalculationResult(
            chromosome=data.chromosome,
            forward_count=forward_count,
            reverse_count=reverse_count,
            correlation_bins=correlation_bins,
            mappable_forward_count=mappable_forward,
            mappable_reverse_count=mappable_reverse,
            mappable_correlation_bins=mappable_bins,
            mappable_length=mappable_length
        )

    def calculate_genome_wide(self,
                            chromosome_data: List[ChromosomeData],
                            config: CalculationConfig,
                            mappability_config: Optional[MappabilityConfig] = None
                            ) -> GenomeWideResult:
        """Calculate cross-correlation for entire genome.

        Processes each chromosome and aggregates results.
        """
        # Calculate each chromosome
        results = []
        for data in chromosome_data:
            result = self.calculate_chromosome(data, config, mappability_config)
            results.append(result)

        # Build chromosome lengths dictionary
        chromosome_lengths = {
            data.chromosome: data.length
            for data in chromosome_data
        }

        # Aggregate results
        return GenomeWideResult.from_chromosome_results(results, chromosome_lengths)


class ParallelCalculationService(CalculationService):
    """Parallel implementation of calculation service.

    Provides parallel processing capabilities while maintaining
    the same interface as the standard service.
    """

    def __init__(self, n_workers: int = 1):
        """Initialize parallel calculation service.

        Args:
            n_workers: Number of worker processes
        """
        self.n_workers = n_workers
        self._standard_service = StandardCalculationService()

    def calculate_chromosome(self,
                           data: ChromosomeData,
                           config: CalculationConfig,
                           mappability_config: Optional[MappabilityConfig] = None
                           ) -> CalculationResult:
        """Calculate using standard service for single chromosome."""
        return self._standard_service.calculate_chromosome(
            data, config, mappability_config
        )

    def calculate_genome_wide(self,
                            chromosome_data: List[ChromosomeData],
                            config: CalculationConfig,
                            mappability_config: Optional[MappabilityConfig] = None
                            ) -> GenomeWideResult:
        """Calculate using parallel processing for multiple chromosomes.

        Note: This is a placeholder for true parallel implementation.
        In production, this would use multiprocessing.
        """
        # For now, delegate to standard service
        # True parallel implementation would distribute work across processes
        return self._standard_service.calculate_genome_wide(
            chromosome_data, config, mappability_config
        )


# Factory function for service creation
def create_calculation_service(parallel: bool = False,
                             n_workers: int = 1) -> CalculationService:
    """Create appropriate calculation service.

    Args:
        parallel: Whether to use parallel processing
        n_workers: Number of worker processes

    Returns:
        Calculation service instance
    """
    if parallel and n_workers > 1:
        return ParallelCalculationService(n_workers)
    else:
        return StandardCalculationService()