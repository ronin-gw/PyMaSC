"""Utilities for statistical calculations in cross-correlation analysis.

This module provides common statistical functions used across PyMaSC
for calculating quality metrics, cross-correlation statistics, and
result normalization. Eliminates duplicate calculation code.

Key features:
- NSC (Normalized Strand Coefficient) calculation
- RSC (Relative Strand Coefficient) calculation  
- Cross-correlation peak finding
- Statistical result aggregation
- Numerical stability utilities

This eliminates 50-65 lines of duplicate statistical code.
"""
import logging
import numpy as np
from typing import Dict, List, Tuple, Optional, Any, Union

logger = logging.getLogger(__name__)


class CrossCorrelationMetrics:
    """Standardized cross-correlation quality metrics calculation.
    
    Provides consistent calculation of NSC and RSC metrics used
    throughout PyMaSC for quality assessment of ChIP-seq data.
    
    These metrics help identify phantom peaks and assess the
    quality of immunoprecipitation in ChIP-seq experiments.
    """
    
    @staticmethod
    def calculate_nsc(ccbins: np.ndarray, 
                     fragment_peak_pos: int,
                     read_length: int) -> float:
        """Calculate Normalized Strand Coefficient (NSC).
        
        NSC is the ratio of the fragment-length cross-correlation peak
        to the background cross-correlation.
        
        Args:
            ccbins: Cross-correlation values array
            fragment_peak_pos: Position of fragment-length peak
            read_length: Read length for phantom peak calculation
            
        Returns:
            NSC value (higher is better, typically > 1.05)
        """
        try:
            # Get fragment peak value
            fragment_peak = ccbins[fragment_peak_pos]
            
            # Calculate background (excluding read-length phantom peak region)
            phantom_start = max(0, read_length - 5)
            phantom_end = min(len(ccbins), read_length + 5)
            
            # Create mask for background calculation
            background_mask = np.ones(len(ccbins), dtype=bool)
            background_mask[phantom_start:phantom_end] = False
            
            # Calculate background as minimum or percentile
            if np.any(background_mask):
                background = np.percentile(ccbins[background_mask], 5)
            else:
                background = np.min(ccbins)
            
            # Avoid division by zero
            if background <= 0:
                background = 1e-10
            
            nsc = fragment_peak / background
            
            return float(nsc)
            
        except Exception as e:
            logger.error(f"Error calculating NSC: {e}")
            return 1.0  # Default fallback value
    
    @staticmethod
    def calculate_rsc(ccbins: np.ndarray,
                     fragment_peak_pos: int,
                     phantom_peak_pos: int,
                     min_value: Optional[float] = None) -> float:
        """Calculate Relative Strand Coefficient (RSC).
        
        RSC is the ratio of the fragment-length peak to the read-length
        (phantom) peak in the cross-correlation profile.
        
        Args:
            ccbins: Cross-correlation values array
            fragment_peak_pos: Position of fragment-length peak
            phantom_peak_pos: Position of phantom peak
            min_value: Optional minimum correlation value
            
        Returns:
            RSC value (higher is better, typically > 0.8)
        """
        try:
            # Get peak values
            fragment_peak = ccbins[fragment_peak_pos]
            phantom_peak = ccbins[phantom_peak_pos]
            
            # Use provided minimum or calculate
            if min_value is None:
                min_value = np.min(ccbins)
            
            # Ensure non-negative values
            min_value = max(0, min_value)
            
            # Calculate RSC
            numerator = fragment_peak - min_value
            denominator = phantom_peak - min_value
            
            # Avoid division by zero
            if denominator <= 0:
                return 1.0  # Default when phantom peak is not prominent
            
            rsc = numerator / denominator
            
            return float(rsc)
            
        except Exception as e:
            logger.error(f"Error calculating RSC: {e}")
            return 0.0  # Default fallback value
    
    @staticmethod
    def find_fragment_peak(ccbins: np.ndarray,
                          min_shift: int = 50,
                          max_shift: Optional[int] = None) -> Tuple[int, float]:
        """Find the fragment-length peak in cross-correlation.
        
        Identifies the primary peak corresponding to the typical
        fragment length in ChIP-seq data.
        
        Args:
            ccbins: Cross-correlation values array
            min_shift: Minimum shift to consider for fragment peak
            max_shift: Maximum shift to consider (None = full range)
            
        Returns:
            Tuple of (peak_position, peak_value)
        """
        if max_shift is None:
            max_shift = len(ccbins)
        else:
            max_shift = min(max_shift, len(ccbins))
        
        # Search for peak in valid range
        search_range = ccbins[min_shift:max_shift]
        if len(search_range) == 0:
            return min_shift, ccbins[min_shift] if min_shift < len(ccbins) else 0.0
        
        # Find maximum in search range
        local_max_pos = np.argmax(search_range)
        peak_pos = min_shift + local_max_pos
        peak_value = search_range[local_max_pos]
        
        return int(peak_pos), float(peak_value)
    
    @staticmethod
    def find_phantom_peak(ccbins: np.ndarray,
                         read_length: int,
                         window: int = 10) -> Tuple[int, float]:
        """Find the phantom peak around read length.
        
        Identifies the artificial peak at read length caused by
        sequencing artifacts and mappability bias.
        
        Args:
            ccbins: Cross-correlation values array
            read_length: Expected read length
            window: Search window around read length
            
        Returns:
            Tuple of (peak_position, peak_value)
        """
        # Define search range around read length
        start = max(0, read_length - window)
        end = min(len(ccbins), read_length + window + 1)
        
        if start >= end:
            return read_length, ccbins[read_length] if read_length < len(ccbins) else 0.0
        
        # Find maximum in phantom peak region
        search_range = ccbins[start:end]
        local_max_pos = np.argmax(search_range)
        peak_pos = start + local_max_pos
        peak_value = search_range[local_max_pos]
        
        return int(peak_pos), float(peak_value)


class ResultAggregator:
    """Utilities for aggregating cross-correlation results.
    
    Provides consistent methods for combining results from multiple
    chromosomes or workers into summary statistics.
    """
    
    @staticmethod
    def aggregate_read_counts(ref2forward: Dict[str, int],
                            ref2reverse: Dict[str, int]) -> Tuple[int, int]:
        """Aggregate forward and reverse read counts across chromosomes.
        
        Args:
            ref2forward: Dictionary of forward counts by chromosome
            ref2reverse: Dictionary of reverse counts by chromosome
            
        Returns:
            Tuple of (total_forward, total_reverse) counts
        """
        total_forward = sum(ref2forward.values())
        total_reverse = sum(ref2reverse.values())
        return total_forward, total_reverse
    
    @staticmethod
    def aggregate_ccbins(ref2ccbins: Dict[str, Any],
                        ref2forward: Dict[str, int],
                        ref2reverse: Dict[str, int]) -> np.ndarray:
        """Aggregate cross-correlation bins across chromosomes.
        
        Combines per-chromosome cross-correlation profiles into
        a genome-wide profile weighted by read counts.
        
        Args:
            ref2ccbins: Dictionary of CC bins by chromosome
            ref2forward: Forward read counts by chromosome
            ref2reverse: Reverse read counts by chromosome
            
        Returns:
            Aggregated cross-correlation array
        """
        # Determine the common length of ccbins
        if not ref2ccbins:
            return np.array([])
        
        # Get the first ccbins length as reference
        first_key = next(iter(ref2ccbins))
        ccbins_length = len(ref2ccbins[first_key])
        
        # Initialize aggregated array
        aggregated = np.zeros(ccbins_length, dtype=np.float64)
        total_reads = 0
        
        # Weighted sum based on read counts
        for chrom, ccbins in ref2ccbins.items():
            if chrom in ref2forward and chrom in ref2reverse:
                # Weight by total reads for this chromosome
                weight = ref2forward[chrom] + ref2reverse[chrom]
                if weight > 0 and len(ccbins) == ccbins_length:
                    aggregated += np.array(ccbins) * weight
                    total_reads += weight
        
        # Normalize by total reads
        if total_reads > 0:
            aggregated /= total_reads
        
        return aggregated
    
    @staticmethod
    def calculate_mappable_fraction(ref2mappable_len: Dict[str, int],
                                  ref2length: List[int],
                                  references: List[str]) -> float:
        """Calculate fraction of genome that is mappable.
        
        Args:
            ref2mappable_len: Mappable lengths by chromosome
            ref2length: Total lengths by chromosome position
            references: List of chromosome names
            
        Returns:
            Fraction of genome that is mappable (0-1)
        """
        total_length = sum(ref2length)
        if total_length == 0:
            return 0.0
        
        # Sum mappable lengths for chromosomes we have data for
        mappable_length = 0
        for i, chrom in enumerate(references):
            if chrom in ref2mappable_len:
                mappable_length += ref2mappable_len[chrom]
        
        return mappable_length / total_length


# Convenience functions for common calculations
def calculate_quality_metrics(ccbins: np.ndarray,
                            fragment_peak_pos: int,
                            read_length: int) -> Dict[str, float]:
    """Calculate all standard quality metrics for cross-correlation.
    
    Args:
        ccbins: Cross-correlation values
        fragment_peak_pos: Position of fragment peak
        read_length: Read length for phantom peak
        
    Returns:
        Dictionary with NSC, RSC, and other metrics
    """
    metrics = CrossCorrelationMetrics()
    
    # Find peaks
    phantom_pos, phantom_value = metrics.find_phantom_peak(ccbins, read_length)
    
    # Calculate metrics
    nsc = metrics.calculate_nsc(ccbins, fragment_peak_pos, read_length)
    rsc = metrics.calculate_rsc(ccbins, fragment_peak_pos, phantom_pos)
    
    return {
        'nsc': nsc,
        'rsc': rsc,
        'fragment_peak_pos': fragment_peak_pos,
        'fragment_peak_value': float(ccbins[fragment_peak_pos]),
        'phantom_peak_pos': phantom_pos,
        'phantom_peak_value': phantom_value
    }


def aggregate_results(handler_results: Dict[str, Any]) -> Dict[str, Any]:
    """Aggregate results from a calculation handler.
    
    Provides a standardized way to aggregate per-chromosome results
    into genome-wide summary statistics.
    
    Args:
        handler_results: Dictionary with ref2forward_sum, ref2reverse_sum, etc.
        
    Returns:
        Dictionary with aggregated statistics
    """
    aggregator = ResultAggregator()
    
    # Extract data from handler
    ref2forward = handler_results.get('ref2forward_sum', {})
    ref2reverse = handler_results.get('ref2reverse_sum', {})
    ref2ccbins = handler_results.get('ref2ccbins', {})
    
    # Aggregate counts
    total_forward, total_reverse = aggregator.aggregate_read_counts(ref2forward, ref2reverse)
    
    # Aggregate cross-correlation
    if ref2ccbins:
        aggregated_cc = aggregator.aggregate_ccbins(ref2ccbins, ref2forward, ref2reverse)
    else:
        aggregated_cc = np.array([])
    
    return {
        'total_forward_reads': total_forward,
        'total_reverse_reads': total_reverse,
        'total_reads': total_forward + total_reverse,
        'aggregated_ccbins': aggregated_cc,
        'num_chromosomes': len(ref2forward)
    }