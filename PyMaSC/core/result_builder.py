"""Builder pattern for constructing PyMaSC analysis results.

This module provides a clean builder interface for constructing
complex PyMaSC result objects, replacing the unwieldy initialization
logic in the original PyMaSCStats class.

Key improvements:
- Step-by-step construction with validation
- Clear separation of data sources (handler vs. external data)
- Type-safe parameter setting
- Flexible configuration of analysis options
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Union, Tuple
import logging

import numpy as np

from PyMaSC.core.result_container import (
    ResultContainer, ChromosomeResult, NCCData, MSCCData, create_from_handler_data
)
from PyMaSC.core.statistics_calculator import (
    StatisticsCalculator, StatisticsConfig, StatisticsResult
)
from PyMaSC.utils.stats_utils import ResultAggregator

logger = logging.getLogger(__name__)


class ResultBuilder:
    """Builder for constructing PyMaSC analysis results.
    
    Provides a fluent interface for building complex result objects
    with proper validation and error handling. Supports both handler-based
    initialization and external data loading.
    
    Usage:
        # From handler
        result = (ResultBuilder()
                 .from_handler(handler)
                 .with_statistics_config(mv_avr_filter_len=15)
                 .with_fragment_length(expected=200)
                 .build())
        
        # From external data
        result = (ResultBuilder()
                 .with_basic_params(read_len=50, references=['chr1'])
                 .with_ncc_data(ref2forward={...}, ref2reverse={...})
                 .build())
    """
    
    def __init__(self):
        """Initialize empty builder."""
        # Basic parameters
        self._read_length: Optional[int] = None
        self._references: Optional[List[str]] = None
        self._chromosome_lengths: Optional[List[int]] = None
        
        # Data containers
        self._ncc_data: Dict[str, Dict[str, Any]] = {}
        self._mscc_data: Dict[str, Dict[str, Any]] = {}
        
        # Configuration
        self._statistics_config = StatisticsConfig()
        self._expected_fragment_length: Optional[int] = None
        self._chi2_pval: float = 0.01
        
        # Flags
        self._skip_ncc: bool = False
        self._has_mappability: bool = False
        
        # Built container
        self._result_container: Optional[ResultContainer] = None
    
    def from_handler(self, handler: Any) -> 'ResultBuilder':
        """Initialize builder from calculation handler.
        
        Args:
            handler: Calculation handler with results
            
        Returns:
            Self for method chaining
        """
        self._read_length = handler.read_len
        self._references = handler.references
        # Use complete chromosome lengths for genome length calculation if available
        self._chromosome_lengths = getattr(handler, 'all_lengths', handler.lengths)
        self._analysis_lengths = handler.lengths  # Store analysis lengths separately
        
        # Extract data based on what's available
        handler_data = {
            'skip_ncc': getattr(handler, 'skip_ncc', False),
            'ref2forward_sum': getattr(handler, 'ref2forward_sum', {}),
            'ref2reverse_sum': getattr(handler, 'ref2reverse_sum', {}),
            'ref2ccbins': getattr(handler, 'ref2ccbins', {}),
            'mappable_ref2forward_sum': getattr(handler, 'mappable_ref2forward_sum', {}),
            'mappable_ref2reverse_sum': getattr(handler, 'mappable_ref2reverse_sum', {}),
            'mappable_ref2ccbins': getattr(handler, 'mappable_ref2ccbins', {}),
            'ref2mappable_len': getattr(handler, 'ref2mappable_len', {})
        }
        
        # Create container from handler data (use analysis lengths for per-chromosome results)
        self._result_container = create_from_handler_data(
            handler_data, self._read_length, self._references, self._analysis_lengths
        )
        
        # Set flags
        self._skip_ncc = handler_data['skip_ncc']
        self._has_mappability = bool(handler_data['mappable_ref2ccbins'])
        
        return self
    
    def with_basic_params(self, 
                         read_length: int, 
                         references: List[str],
                         chromosome_lengths: Optional[List[int]] = None) -> 'ResultBuilder':
        """Set basic analysis parameters.
        
        Args:
            read_length: Read length in base pairs
            references: List of chromosome names
            chromosome_lengths: Optional chromosome lengths
            
        Returns:
            Self for method chaining
        """
        self._read_length = read_length
        self._references = references
        self._chromosome_lengths = chromosome_lengths
        return self
    
    def with_ncc_data(self,
                     ref2forward_sum: Dict[str, int],
                     ref2reverse_sum: Dict[str, int],
                     ref2ccbins: Dict[str, np.ndarray],
                     ref2genomelen: Optional[Dict[str, int]] = None) -> 'ResultBuilder':
        """Add NCC (Naive Cross-Correlation) data.
        
        Args:
            ref2forward_sum: Forward read counts by chromosome
            ref2reverse_sum: Reverse read counts by chromosome
            ref2ccbins: Cross-correlation bins by chromosome
            ref2genomelen: Optional genome lengths by chromosome
            
        Returns:
            Self for method chaining
        """
        self._ncc_data = {
            'forward_sum': ref2forward_sum,
            'reverse_sum': ref2reverse_sum,
            'ccbins': ref2ccbins,
            'genomelen': ref2genomelen or {}
        }
        return self
    
    def with_mscc_data(self,
                      mappable_ref2forward_sum: Dict[str, np.ndarray],
                      mappable_ref2reverse_sum: Dict[str, np.ndarray],
                      mappable_ref2ccbins: Dict[str, np.ndarray],
                      ref2mappable_len: Dict[str, Union[int, np.ndarray]]) -> 'ResultBuilder':
        """Add MSCC (Mappability-Sensitive Cross-Correlation) data.
        
        Args:
            mappable_ref2forward_sum: Mappable forward counts by chromosome
            mappable_ref2reverse_sum: Mappable reverse counts by chromosome
            mappable_ref2ccbins: Mappable CC bins by chromosome
            ref2mappable_len: Mappable lengths by chromosome
            
        Returns:
            Self for method chaining
        """
        self._mscc_data = {
            'forward_sum': mappable_ref2forward_sum,
            'reverse_sum': mappable_ref2reverse_sum,
            'ccbins': mappable_ref2ccbins,
            'mappable_len': ref2mappable_len
        }
        self._has_mappability = True
        return self
    
    def with_statistics_config(self,
                             mv_avr_filter_len: int = 15,
                             filter_mask_len: int = 5,
                             min_calc_width: int = 10,
                             output_warnings: bool = True,
                             do_fragment_estimation: bool = False) -> 'ResultBuilder':
        """Configure statistics calculation parameters.
        
        Args:
            mv_avr_filter_len: Moving average filter length
            filter_mask_len: Mask length around read length
            min_calc_width: Width for minimum calculation
            output_warnings: Whether to output warnings
            do_fragment_estimation: Whether to estimate fragment length
            
        Returns:
            Self for method chaining
        """
        self._statistics_config = StatisticsConfig(
            mv_avr_filter_len=mv_avr_filter_len,
            filter_mask_len=filter_mask_len,
            min_calc_width=min_calc_width,
            output_warnings=output_warnings,
            do_fragment_estimation=do_fragment_estimation,
            expected_fragment_length=self._expected_fragment_length
        )
        return self
    
    def with_fragment_length(self, expected: Optional[int] = None) -> 'ResultBuilder':
        """Set expected fragment length.
        
        Args:
            expected: Expected fragment length in base pairs
            
        Returns:
            Self for method chaining
        """
        self._expected_fragment_length = expected
        if self._statistics_config:
            self._statistics_config.expected_fragment_length = expected
        return self
    
    def with_chi2_threshold(self, pvalue: float = 0.01) -> 'ResultBuilder':
        """Set chi-squared test p-value threshold.
        
        Args:
            pvalue: P-value threshold for strand bias testing
            
        Returns:
            Self for method chaining
        """
        self._chi2_pval = pvalue
        return self
    
    def skip_ncc_calculation(self, skip: bool = True) -> 'ResultBuilder':
        """Configure whether to skip NCC calculation.
        
        Args:
            skip: Whether to skip NCC calculation
            
        Returns:
            Self for method chaining
        """
        self._skip_ncc = skip
        return self
    
    def build(self) -> 'BuiltResult':
        """Build the final result object.
        
        Returns:
            Complete result object with calculated statistics
            
        Raises:
            ValueError: If required parameters are missing
        """
        # Validate required parameters
        if self._read_length is None:
            raise ValueError("Read length must be specified")
        if self._references is None:
            raise ValueError("References must be specified")
        
        # Create result container if not already created
        if self._result_container is None:
            self._result_container = self._build_container_from_data()
        
        # Calculate statistics for each chromosome
        chromosome_statistics = {}
        for chrom, result in self._result_container.chromosome_results.items():
            chrom_stats = {}
            
            # Calculate NCC statistics
            if result.ncc_data and not self._skip_ncc:
                ncc_calculator = StatisticsCalculator(
                    result.ncc_data, self._read_length, self._statistics_config
                )
                chrom_stats['ncc'] = ncc_calculator.calculate_statistics()
            
            # Calculate MSCC statistics
            if result.mscc_data:
                mscc_config = StatisticsConfig(
                    **{**self._statistics_config.__dict__, 'do_fragment_estimation': True}
                )
                mscc_calculator = StatisticsCalculator(
                    result.mscc_data, self._read_length, mscc_config
                )
                chrom_stats['mscc'] = mscc_calculator.calculate_statistics()
            
            chromosome_statistics[chrom] = chrom_stats
        
        # Calculate aggregate statistics
        aggregate_stats = self._calculate_aggregate_statistics()
        
        return BuiltResult(
            container=self._result_container,
            chromosome_statistics=chromosome_statistics,
            aggregate_statistics=aggregate_stats,
            config=self._statistics_config,
            chi2_pval=self._chi2_pval
        )
    
    def _build_container_from_data(self) -> ResultContainer:
        """Build result container from manually provided data."""
        container = ResultContainer(
            read_length=self._read_length,
            references=self._references,
            skip_ncc=self._skip_ncc,
            has_mappability=self._has_mappability
        )
        
        # Create chromosome results
        for i, chrom in enumerate(self._references):
            ncc_data = None
            mscc_data = None
            
            # Create NCC data if available
            if (not self._skip_ncc and 
                chrom in self._ncc_data.get('ccbins', {})):
                length = (self._chromosome_lengths[i] if self._chromosome_lengths 
                         else self._ncc_data.get('genomelen', {}).get(chrom, 1))
                ncc_data = NCCData(
                    forward_sum=self._ncc_data['forward_sum'].get(chrom, 0),
                    reverse_sum=self._ncc_data['reverse_sum'].get(chrom, 0),
                    ccbins=self._ncc_data['ccbins'][chrom],
                    genomelen=length
                )
            
            # Create MSCC data if available
            if chrom in self._mscc_data.get('ccbins', {}):
                mscc_data = MSCCData(
                    forward_sum=self._mscc_data['forward_sum'].get(chrom, np.array([])),
                    reverse_sum=self._mscc_data['reverse_sum'].get(chrom, np.array([])),
                    ccbins=self._mscc_data['ccbins'][chrom],
                    mappable_len=self._mscc_data['mappable_len'].get(chrom, 0)
                )
            
            if ncc_data or mscc_data:
                length = (self._chromosome_lengths[i] if self._chromosome_lengths 
                         else (ncc_data.genomelen if ncc_data else 1))
                result = ChromosomeResult(
                    chromosome=chrom,
                    length=length,
                    ncc_data=ncc_data,
                    mscc_data=mscc_data
                )
                container.add_chromosome_result(result)
        
        return container
    
    def _calculate_aggregate_statistics(self) -> Dict[str, Any]:
        """Calculate genome-wide aggregate statistics."""
        aggregator = ResultAggregator()
        
        # Get total read counts
        read_counts = self._result_container.get_total_reads()
        
        # Calculate aggregate cross-correlation if NCC data available
        ncc_data = self._result_container.get_ncc_data()
        if ncc_data and not self._skip_ncc:
            ref2ccbins = {chrom: data.ccbins for chrom, data in ncc_data.items()}
            ref2forward = {chrom: data.forward_sum for chrom, data in ncc_data.items()}
            ref2reverse = {chrom: data.reverse_sum for chrom, data in ncc_data.items()}
            
            aggregated_cc = aggregator.aggregate_ccbins(
                ref2ccbins, ref2forward, ref2reverse
            )
        else:
            aggregated_cc = np.array([])
        
        # Calculate mappable fraction if MSCC data available
        mscc_data = self._result_container.get_mscc_data()
        mappable_fraction = 0.0
        if mscc_data:
            ref2mappable_len = {}
            for chrom, data in mscc_data.items():
                if data.mappable_len is None:
                    # Skip chromosomes with no mappable length data
                    continue
                elif isinstance(data.mappable_len, (int, np.integer)):
                    ref2mappable_len[chrom] = int(data.mappable_len)
                else:
                    ref2mappable_len[chrom] = int(np.sum(data.mappable_len))
            
            lengths = [result.length for result in self._result_container.chromosome_results.values()]
            mappable_fraction = aggregator.calculate_mappable_fraction(
                ref2mappable_len, lengths, self._references
            )
        
        return {
            'total_reads': read_counts,
            'genome_length': sum(self._chromosome_lengths),  # Use complete chromosome lengths
            'aggregated_cc': aggregated_cc,
            'mappable_fraction': mappable_fraction,
            'num_chromosomes': len(self._result_container.chromosome_results)
        }


class BuiltResult:
    """Container for complete PyMaSC analysis results.
    
    This class represents the final result of the builder process,
    containing all calculated statistics and aggregate data.
    
    Attributes:
        container: Structured data container
        chromosome_statistics: Per-chromosome statistical results
        aggregate_statistics: Genome-wide aggregate statistics
        config: Statistics calculation configuration
        chi2_pval: Chi-squared test p-value threshold
    """
    
    def __init__(self,
                 container: ResultContainer,
                 chromosome_statistics: Dict[str, Dict[str, StatisticsResult]],
                 aggregate_statistics: Dict[str, Any],
                 config: StatisticsConfig,
                 chi2_pval: float):
        """Initialize built result.
        
        Args:
            container: Data container
            chromosome_statistics: Per-chromosome stats
            aggregate_statistics: Aggregate stats
            config: Statistics configuration
            chi2_pval: Chi-squared threshold
        """
        self.container = container
        self.chromosome_statistics = chromosome_statistics
        self.aggregate_statistics = aggregate_statistics
        self.config = config
        self.chi2_pval = chi2_pval
    
    @property
    def read_length(self) -> int:
        """Get read length."""
        return self.container.read_length
    
    @property
    def references(self) -> List[str]:
        """Get analyzed chromosome names."""
        return self.container.references
    
    @property
    def skip_ncc(self) -> bool:
        """Check if NCC was skipped."""
        return self.container.skip_ncc
    
    @property
    def has_mappability(self) -> bool:
        """Check if mappability analysis was performed."""
        return self.container.has_mappability
    
    def get_ncc_statistics(self, chromosome: Optional[str] = None) -> Union[StatisticsResult, Dict[str, StatisticsResult]]:
        """Get NCC statistics for chromosome(s).
        
        Args:
            chromosome: Specific chromosome, or None for all
            
        Returns:
            Statistics for specified chromosome or all chromosomes
        """
        if chromosome:
            return self.chromosome_statistics.get(chromosome, {}).get('ncc')
        
        return {
            chrom: stats.get('ncc')
            for chrom, stats in self.chromosome_statistics.items()
            if 'ncc' in stats
        }
    
    def get_mscc_statistics(self, chromosome: Optional[str] = None) -> Union[StatisticsResult, Dict[str, StatisticsResult]]:
        """Get MSCC statistics for chromosome(s).
        
        Args:
            chromosome: Specific chromosome, or None for all
            
        Returns:
            Statistics for specified chromosome or all chromosomes
        """
        if chromosome:
            return self.chromosome_statistics.get(chromosome, {}).get('mscc')
        
        return {
            chrom: stats.get('mscc')
            for chrom, stats in self.chromosome_statistics.items()
            if 'mscc' in stats
        }
    
    def get_estimated_fragment_length(self) -> Optional[int]:
        """Get estimated fragment length from MSCC analysis.
        
        Returns:
            Estimated fragment length or None if not available
        """
        for stats in self.chromosome_statistics.values():
            if 'mscc' in stats and stats['mscc'].estimated_fragment_length:
                return stats['mscc'].estimated_fragment_length
        return None


# Convenience function for legacy compatibility
def build_from_handler(handler: Any,
                      mv_avr_filter_len: int = 15,
                      expected_fragment_length: Optional[int] = None,
                      chi2_pval: float = 0.01,
                      filter_mask_len: int = 5,
                      min_calc_width: int = 10) -> BuiltResult:
    """Build result from handler using legacy parameter interface.
    
    Args:
        handler: Calculation handler
        mv_avr_filter_len: Moving average filter length
        expected_fragment_length: Expected fragment length
        chi2_pval: Chi-squared p-value threshold
        filter_mask_len: Filter mask length
        min_calc_width: Minimum calculation width
        
    Returns:
        Complete built result
    """
    return (ResultBuilder()
            .from_handler(handler)
            .with_statistics_config(
                mv_avr_filter_len=mv_avr_filter_len,
                filter_mask_len=filter_mask_len,
                min_calc_width=min_calc_width,
                output_warnings=True,
                do_fragment_estimation=True
            )
            .with_fragment_length(expected=expected_fragment_length)
            .with_chi2_threshold(chi2_pval)
            .build())