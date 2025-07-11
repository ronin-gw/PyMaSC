"""Result aggregation system for PyMaSC cross-correlation analysis.

This module provides a unified system for aggregating calculation results
from various sources (single-process, multi-process, different algorithms).
It replaces the scattered aggregation logic in unified.py and result.py.

Key features:
- Unified aggregation interface for all calculation types
- Type-safe result handling with proper validation
- Support for both NCC and MSCC calculations
- Dual implementation system for gradual migration
- Comprehensive error handling and logging
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional, Union, Any, Tuple
from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from PyMaSC.core.result_container import (
    ResultContainer, ChromosomeResult, NCCData, MSCCData
)
from PyMaSC.core.interfaces import CrossCorrelationCalculator
from PyMaSC.utils.stats_utils import ArrayAggregator

# Import WorkerResult for the aggregate method
from PyMaSC.core.worker import WorkerResult

logger = logging.getLogger(__name__)


class AggregationMode(Enum):
    """Result aggregation modes."""
    LEGACY = "legacy"  # Direct attribute access (backward compatibility)
    CONTAINER = "container"  # Type-safe container-based aggregation
    HYBRID = "hybrid"  # Dual mode during migration


@dataclass
class AggregationConfig:
    """Configuration for result aggregation."""
    mode: AggregationMode = AggregationMode.HYBRID
    enable_validation: bool = True
    enable_logging: bool = True
    fallback_to_legacy: bool = True
    
    # Validation thresholds
    max_array_size: int = 10_000_000  # Maximum array elements
    max_chromosome_count: int = 1000  # Maximum chromosomes
    numerical_tolerance: float = 1e-12  # Numerical precision tolerance


@dataclass
class AggregationResult:
    """Result of aggregation process."""
    container: Optional[ResultContainer] = None
    legacy_attributes: Dict[str, Any] = field(default_factory=dict)
    aggregation_stats: Dict[str, Any] = field(default_factory=dict)
    validation_errors: List[str] = field(default_factory=list)
    aggregation_mode: AggregationMode = AggregationMode.HYBRID


class ResultAggregator:
    """Unified result aggregation system for PyMaSC calculations.
    
    This class provides a single interface for aggregating calculation results
    from various sources while maintaining backward compatibility with existing
    code. It supports both legacy direct attribute access and new type-safe
    container-based aggregation.
    
    Features:
    - Unified interface for all aggregation needs
    - Type-safe result handling
    - Validation and error reporting
    - Dual implementation during migration
    - Comprehensive logging and debugging
    
    Usage:
        # Single-process aggregation
        aggregator = ResultAggregator(config)
        result = aggregator.aggregate_from_calculator(calculator, references, lengths)
        
        # Multi-process aggregation (legacy format)
        result = aggregator.aggregate_from_results(worker_results)
        
        # Multi-process aggregation (WorkerResult objects)
        result = aggregator.aggregate(worker_results, references, lengths, read_length)
        
        # Legacy compatibility
        handler_attrs = result.legacy_attributes
        handler.ref2forward_sum = handler_attrs['ref2forward_sum']
    """
    
    def __init__(self, config: Optional[AggregationConfig] = None):
        """Initialize result aggregator.
        
        Args:
            config: Aggregation configuration (uses defaults if None)
        """
        self.config = config or AggregationConfig()
        self.logger = logger if self.config.enable_logging else logging.getLogger('null')
        
        # Internal state
        self._aggregation_stats = {
            'total_operations': 0,
            'successful_operations': 0,
            'failed_operations': 0,
            'validation_errors': 0,
            'fallback_to_legacy': 0
        }
    
    def aggregate_from_calculator(
        self,
        calculator: CrossCorrelationCalculator,
        references: List[str],
        lengths: List[int],
        read_length: int,
        skip_ncc: bool = False
    ) -> AggregationResult:
        """Aggregate results from a single calculator instance.
        
        Args:
            calculator: Calculator with completed calculations
            references: List of chromosome names
            lengths: List of chromosome lengths
            read_length: Read length in base pairs
            skip_ncc: Whether NCC calculation was skipped
            
        Returns:
            AggregationResult with both container and legacy attributes
        """
        self._aggregation_stats['total_operations'] += 1
        
        try:
            # Create container-based aggregation
            container = self._create_container_from_calculator(
                calculator, references, lengths, read_length, skip_ncc
            )
            
            # Create legacy attributes for backward compatibility
            legacy_attrs = self._create_legacy_attributes_from_calculator(
                calculator, references, lengths
            )
            
            # Validation
            if self.config.enable_validation:
                validation_errors = self._validate_calculator_results(
                    calculator, container, legacy_attrs
                )
            else:
                validation_errors = []
            
            self._aggregation_stats['successful_operations'] += 1
            
            return AggregationResult(
                container=container,
                legacy_attributes=legacy_attrs,
                aggregation_stats=self._aggregation_stats.copy(),
                validation_errors=validation_errors,
                aggregation_mode=self.config.mode
            )
            
        except Exception as e:
            self._aggregation_stats['failed_operations'] += 1
            self.logger.error(f"Calculator aggregation failed: {e}")
            
            # Fallback to legacy if enabled
            if self.config.fallback_to_legacy:
                return self._fallback_calculator_aggregation(
                    calculator, references, lengths, read_length, skip_ncc
                )
            else:
                raise
    
    def aggregate_from_results(
        self,
        worker_results: List[Tuple[str, Any]],
        references: List[str],
        lengths: List[int],
        read_length: int,
        skip_ncc: bool = False
    ) -> AggregationResult:
        """Aggregate results from multiple worker processes.
        
        Args:
            worker_results: List of (chromosome, result_tuple) pairs
            references: List of chromosome names
            lengths: List of chromosome lengths
            read_length: Read length in base pairs
            skip_ncc: Whether NCC calculation was skipped
            
        Returns:
            AggregationResult with aggregated data
        """
        self._aggregation_stats['total_operations'] += 1
        
        try:
            # Parse worker results
            parsed_results = self._parse_worker_results(worker_results)
            
            # Create container-based aggregation
            container = self._create_container_from_worker_results(
                parsed_results, references, lengths, read_length, skip_ncc
            )
            
            # Create legacy attributes
            legacy_attrs = self._create_legacy_attributes_from_worker_results(
                parsed_results, references, lengths
            )
            
            # Validation
            if self.config.enable_validation:
                validation_errors = self._validate_worker_results(
                    parsed_results, container, legacy_attrs
                )
            else:
                validation_errors = []
            
            self._aggregation_stats['successful_operations'] += 1
            
            return AggregationResult(
                container=container,
                legacy_attributes=legacy_attrs,
                aggregation_stats=self._aggregation_stats.copy(),
                validation_errors=validation_errors,
                aggregation_mode=self.config.mode
            )
            
        except Exception as e:
            self._aggregation_stats['failed_operations'] += 1
            self.logger.error(f"Worker results aggregation failed: {e}")
            
            # Fallback to legacy if enabled
            if self.config.fallback_to_legacy:
                return self._fallback_worker_aggregation(
                    worker_results, references, lengths, read_length, skip_ncc
                )
            else:
                raise
    
    def aggregate(
        self,
        worker_results: List[WorkerResult],
        references: List[str],
        lengths: List[int],
        read_length: int,
        skip_ncc: bool = False
    ) -> AggregationResult:
        """Aggregate results from WorkerResult objects.
        
        This method provides a modern interface for aggregating results from
        the new WorkerResult dataclass format, converting them internally to
        the legacy format for existing aggregation logic.
        
        Args:
            worker_results: List of WorkerResult objects
            references: List of chromosome names
            lengths: List of chromosome lengths
            read_length: Read length in base pairs
            skip_ncc: Whether NCC calculation was skipped
            
        Returns:
            AggregationResult with aggregated data
        """
        # Convert WorkerResult objects to legacy format
        legacy_worker_results = []
        
        for result in worker_results:
            # Convert to legacy tuple format: (chromosome, (mappable_len, cc_stats, masc_stats))
            cc_stats = None
            masc_stats = None
            
            # Build NCC stats if available
            if (result.ncc_forward_sum is not None and 
                result.ncc_reverse_sum is not None and 
                result.ncc_bins is not None):
                cc_stats = (result.ncc_forward_sum, result.ncc_reverse_sum, result.ncc_bins)
            
            # Build MSCC stats if available
            if (result.mscc_forward_sum is not None and 
                result.mscc_reverse_sum is not None and 
                result.mscc_bins is not None):
                masc_stats = (result.mscc_forward_sum, result.mscc_reverse_sum, result.mscc_bins)
            
            # Create legacy tuple
            legacy_tuple = (result.mappable_length, cc_stats, masc_stats)
            legacy_worker_results.append((result.chromosome, legacy_tuple))
        
        # Use existing aggregation logic
        return self.aggregate_from_results(
            legacy_worker_results, references, lengths, read_length, skip_ncc
        )
    
    def _create_container_from_calculator(
        self,
        calculator: CrossCorrelationCalculator,
        references: List[str],
        lengths: List[int],
        read_length: int,
        skip_ncc: bool
    ) -> ResultContainer:
        """Create type-safe container from calculator results."""
        container = ResultContainer(
            read_length=read_length,
            references=references,
            skip_ncc=skip_ncc,
            has_mappability=self._has_mappability_data(calculator)
        )
        
        # Process each chromosome
        for i, chrom in enumerate(references):
            ncc_data = None
            mscc_data = None
            
            # Extract NCC data if available
            if (not skip_ncc and 
                hasattr(calculator, 'ref2ccbins') and 
                chrom in calculator.ref2ccbins and
                calculator.ref2ccbins[chrom] is not None):
                ncc_data = NCCData(
                    forward_sum=calculator.ref2forward_sum.get(chrom, 0) if hasattr(calculator.ref2forward_sum, 'get') else calculator.ref2forward_sum[chrom],
                    reverse_sum=calculator.ref2reverse_sum.get(chrom, 0) if hasattr(calculator.ref2reverse_sum, 'get') else calculator.ref2reverse_sum[chrom],
                    ccbins=calculator.ref2ccbins[chrom],
                    genomelen=lengths[i]
                )
            
            # Extract MSCC data if available
            if (hasattr(calculator, 'ref2mappable_forward_sum') and 
                calculator.ref2mappable_forward_sum and
                chrom in calculator.ref2mappable_forward_sum and
                calculator.ref2mappable_forward_sum[chrom] is not None):
                mscc_data = MSCCData(
                    forward_sum=calculator.ref2mappable_forward_sum[chrom],
                    reverse_sum=calculator.ref2mappable_reverse_sum[chrom] if chrom in calculator.ref2mappable_reverse_sum else np.array([]),
                    ccbins=calculator.ref2mascbins[chrom] if chrom in calculator.ref2mascbins else np.array([]),
                    mappable_len=calculator.ref2mappable_len[chrom] if hasattr(calculator, 'ref2mappable_len') and chrom in calculator.ref2mappable_len else None
                )
            
            if ncc_data or mscc_data:
                result = ChromosomeResult(
                    chromosome=chrom,
                    length=lengths[i],
                    ncc_data=ncc_data,
                    mscc_data=mscc_data
                )
                container.add_chromosome_result(result)
        
        return container
    
    def _create_legacy_attributes_from_calculator(
        self,
        calculator: CrossCorrelationCalculator,
        references: List[str],
        lengths: List[int]
    ) -> Dict[str, Any]:
        """Create legacy handler attributes from calculator."""
        legacy_attrs = {}
        
        # Basic NCC attributes
        if hasattr(calculator, 'ref2forward_sum'):
            legacy_attrs['ref2forward_sum'] = calculator.ref2forward_sum
        if hasattr(calculator, 'ref2reverse_sum'):
            legacy_attrs['ref2reverse_sum'] = calculator.ref2reverse_sum
        if hasattr(calculator, 'ref2ccbins'):
            legacy_attrs['ref2ccbins'] = calculator.ref2ccbins
        
        # MSCC attributes from adapter
        if hasattr(calculator, 'ref2mappable_forward_sum'):
            legacy_attrs['mappable_ref2forward_sum'] = calculator.ref2mappable_forward_sum
        if hasattr(calculator, 'ref2mappable_reverse_sum'):
            legacy_attrs['mappable_ref2reverse_sum'] = calculator.ref2mappable_reverse_sum
        if hasattr(calculator, 'ref2mascbins'):
            legacy_attrs['mappable_ref2ccbins'] = calculator.ref2mascbins
        
        # Mappable length from BitArray calculator
        if hasattr(calculator, '_calculator') and hasattr(calculator._calculator, 'ref2mappable_len'):
            legacy_attrs['ref2mappable_len'] = calculator._calculator.ref2mappable_len
        
        # Aggregate mappable results (mimics unified.py behavior)
        legacy_attrs.update(self._aggregate_mappable_results_legacy(legacy_attrs))
        
        return legacy_attrs
    
    def _aggregate_mappable_results_legacy(self, legacy_attrs: Dict[str, Any]) -> Dict[str, Any]:
        """Aggregate mappable results to scalar attributes (legacy compatibility)."""
        aggregated = {}
        
        # Initialize as None by default
        aggregated['mappable_forward_sum'] = None
        aggregated['mappable_reverse_sum'] = None
        aggregated['mappable_ccbins'] = None
        
        # Aggregate if we have mappable data
        mappable_forward = legacy_attrs.get('mappable_ref2forward_sum', {})
        mappable_reverse = legacy_attrs.get('mappable_ref2reverse_sum', {})
        mappable_ccbins = legacy_attrs.get('mappable_ref2ccbins', {})
        
        if mappable_forward:
            valid_forward = [v for v in mappable_forward.values() if v is not None]
            if valid_forward:
                aggregated['mappable_forward_sum'] = np.sum(valid_forward, axis=0)
        
        if mappable_reverse:
            valid_reverse = [v for v in mappable_reverse.values() if v is not None]
            if valid_reverse:
                aggregated['mappable_reverse_sum'] = np.sum(valid_reverse, axis=0)
        
        if mappable_ccbins:
            valid_ccbins = [v for v in mappable_ccbins.values() if v is not None]
            if valid_ccbins:
                aggregated['mappable_ccbins'] = np.sum(valid_ccbins, axis=0)
        
        return aggregated
    
    def _parse_worker_results(self, worker_results: List[Tuple[str, Any]]) -> Dict[str, Any]:
        """Parse worker results into structured format."""
        parsed = {
            'ref2forward_sum': {},
            'ref2reverse_sum': {},
            'ref2ccbins': {},
            'mappable_ref2forward_sum': {},
            'mappable_ref2reverse_sum': {},
            'mappable_ref2ccbins': {},
            'ref2mappable_len': {}
        }
        
        for chrom, result_tuple in worker_results:
            if len(result_tuple) >= 3:
                mappable_len, cc_stats, masc_stats = result_tuple[:3]
                
                # Parse NCC stats
                if cc_stats and len(cc_stats) >= 3:
                    f_sum, r_sum, ccbins = cc_stats
                    if None not in (f_sum, r_sum, ccbins):
                        parsed['ref2forward_sum'][chrom] = f_sum
                        parsed['ref2reverse_sum'][chrom] = r_sum
                        parsed['ref2ccbins'][chrom] = ccbins
                
                # Parse MSCC stats
                if masc_stats and len(masc_stats) >= 3:
                    mf_sum, mr_sum, mccbins = masc_stats
                    if None not in (mf_sum, mr_sum, mccbins):
                        parsed['mappable_ref2forward_sum'][chrom] = mf_sum
                        parsed['mappable_ref2reverse_sum'][chrom] = mr_sum
                        parsed['mappable_ref2ccbins'][chrom] = mccbins
                
                # Parse mappable length
                if mappable_len is not None:
                    parsed['ref2mappable_len'][chrom] = mappable_len
        
        return parsed
    
    def _create_container_from_worker_results(
        self,
        parsed_results: Dict[str, Any],
        references: List[str],
        lengths: List[int],
        read_length: int,
        skip_ncc: bool
    ) -> ResultContainer:
        """Create container from parsed worker results."""
        # Use the existing factory function with proper data structure
        handler_data = {
            'skip_ncc': skip_ncc,
            'ref2forward_sum': parsed_results['ref2forward_sum'],
            'ref2reverse_sum': parsed_results['ref2reverse_sum'],
            'ref2ccbins': parsed_results['ref2ccbins'],
            'mappable_ref2forward_sum': parsed_results['mappable_ref2forward_sum'],
            'mappable_ref2reverse_sum': parsed_results['mappable_ref2reverse_sum'],
            'mappable_ref2ccbins': parsed_results['mappable_ref2ccbins'],
            'ref2mappable_len': parsed_results['ref2mappable_len']
        }
        
        from PyMaSC.core.result_container import create_from_handler_data
        return create_from_handler_data(handler_data, read_length, references, lengths)
    
    def _create_legacy_attributes_from_worker_results(
        self,
        parsed_results: Dict[str, Any],
        references: List[str],
        lengths: List[int]
    ) -> Dict[str, Any]:
        """Create legacy attributes from parsed worker results."""
        legacy_attrs = parsed_results.copy()
        
        # Add aggregated mappable results
        legacy_attrs.update(self._aggregate_mappable_results_legacy(legacy_attrs))
        
        return legacy_attrs
    
    def _has_mappability_data(self, calculator: CrossCorrelationCalculator) -> bool:
        """Check if calculator has mappability data."""
        return (hasattr(calculator, 'ref2mappable_forward_sum') and 
                calculator.ref2mappable_forward_sum is not None and
                len(calculator.ref2mappable_forward_sum) > 0)
    
    def _validate_calculator_results(
        self,
        calculator: CrossCorrelationCalculator,
        container: ResultContainer,
        legacy_attrs: Dict[str, Any]
    ) -> List[str]:
        """Validate calculator results for consistency."""
        errors = []
        
        # Check basic consistency
        if not container.chromosome_results:
            errors.append("No chromosome results found")
        
        # Check array sizes
        for result in container.chromosome_results.values():
            if result.ncc_data and len(result.ncc_data.ccbins) > self.config.max_array_size:
                errors.append(f"NCC array too large for {result.chromosome}")
            if result.mscc_data and len(result.mscc_data.ccbins) > self.config.max_array_size:
                errors.append(f"MSCC array too large for {result.chromosome}")
        
        # Check numerical consistency (NCC vs MSCC read counts)
        # This is a simplified check - more sophisticated validation can be added
        
        if errors:
            self._aggregation_stats['validation_errors'] += len(errors)
        
        return errors
    
    def _validate_worker_results(
        self,
        parsed_results: Dict[str, Any],
        container: ResultContainer,
        legacy_attrs: Dict[str, Any]
    ) -> List[str]:
        """Validate worker results for consistency."""
        errors = []
        
        # Check that we have some results
        if not any(parsed_results.values()):
            errors.append("No valid worker results found")
        
        # Check chromosome consistency
        ncc_chroms = set(parsed_results['ref2ccbins'].keys())
        mscc_chroms = set(parsed_results['mappable_ref2ccbins'].keys())
        
        if ncc_chroms and mscc_chroms:
            if not ncc_chroms.issubset(mscc_chroms) and not mscc_chroms.issubset(ncc_chroms):
                errors.append("Chromosome sets differ between NCC and MSCC")
        
        if errors:
            self._aggregation_stats['validation_errors'] += len(errors)
        
        return errors
    
    def _fallback_calculator_aggregation(
        self,
        calculator: CrossCorrelationCalculator,
        references: List[str],
        lengths: List[int],
        read_length: int,
        skip_ncc: bool
    ) -> AggregationResult:
        """Fallback aggregation using legacy methods."""
        self._aggregation_stats['fallback_to_legacy'] += 1
        
        legacy_attrs = self._create_legacy_attributes_from_calculator(
            calculator, references, lengths
        )
        
        return AggregationResult(
            container=None,
            legacy_attributes=legacy_attrs,
            aggregation_stats=self._aggregation_stats.copy(),
            validation_errors=["Fallback to legacy aggregation"],
            aggregation_mode=AggregationMode.LEGACY
        )
    
    def _fallback_worker_aggregation(
        self,
        worker_results: List[Tuple[str, Any]],
        references: List[str],
        lengths: List[int],
        read_length: int,
        skip_ncc: bool
    ) -> AggregationResult:
        """Fallback aggregation for worker results."""
        self._aggregation_stats['fallback_to_legacy'] += 1
        
        parsed_results = self._parse_worker_results(worker_results)
        legacy_attrs = self._create_legacy_attributes_from_worker_results(
            parsed_results, references, lengths
        )
        
        return AggregationResult(
            container=None,
            legacy_attributes=legacy_attrs,
            aggregation_stats=self._aggregation_stats.copy(),
            validation_errors=["Fallback to legacy aggregation"],
            aggregation_mode=AggregationMode.LEGACY
        )
    
    def get_aggregation_stats(self) -> Dict[str, Any]:
        """Get aggregation statistics."""
        return self._aggregation_stats.copy()
    
    def reset_stats(self) -> None:
        """Reset aggregation statistics."""
        self._aggregation_stats = {
            'total_operations': 0,
            'successful_operations': 0,
            'failed_operations': 0,
            'validation_errors': 0,
            'fallback_to_legacy': 0
        }


# Convenience functions for backward compatibility
def aggregate_calculator_results(
    calculator: CrossCorrelationCalculator,
    references: List[str],
    lengths: List[int],
    read_length: int,
    skip_ncc: bool = False,
    config: Optional[AggregationConfig] = None
) -> AggregationResult:
    """Convenience function to aggregate calculator results."""
    aggregator = ResultAggregator(config)
    return aggregator.aggregate_from_calculator(
        calculator, references, lengths, read_length, skip_ncc
    )


def aggregate_worker_results(
    worker_results: List[Tuple[str, Any]],
    references: List[str],
    lengths: List[int],
    read_length: int,
    skip_ncc: bool = False,
    config: Optional[AggregationConfig] = None
) -> AggregationResult:
    """Convenience function to aggregate worker results."""
    aggregator = ResultAggregator(config)
    return aggregator.aggregate_from_results(
        worker_results, references, lengths, read_length, skip_ncc
    )


def aggregate_from_worker_results(
    worker_results: List[WorkerResult],
    references: List[str],
    lengths: List[int],
    read_length: int,
    skip_ncc: bool = False,
    config: Optional[AggregationConfig] = None
) -> AggregationResult:
    """Convenience function to aggregate from WorkerResult objects."""
    aggregator = ResultAggregator(config)
    return aggregator.aggregate(
        worker_results, references, lengths, read_length, skip_ncc
    )