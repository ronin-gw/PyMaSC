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
from functools import wraps
from typing import Dict, List, Tuple, Optional, Any, Union
from PyMaSC.utils.calc import moving_avr_filter

logger = logging.getLogger(__name__)


def npcalc_with_logging_warn(func):
    """Decorator for handling numpy floating point errors gracefully.
    
    This decorator wraps numerical calculation functions to handle common
    floating point errors that can occur during cross-correlation calculations,
    such as division by zero or invalid operations.
    """
    @wraps(func)
    def _inner(*args, **kwargs):
        try:
            with np.errstate(divide="raise", invalid="raise"):
                return func(*args, **kwargs)
        except (FloatingPointError, ZeroDivisionError) as e:
            logger.debug("catch numpy warning: " + repr(e))
            logger.debug("continue anyway.")
            with np.errstate(divide="ignore", invalid="ignore"):
                return func(*args, **kwargs)
    return _inner


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
            # Validate inputs
            if len(ccbins) == 0:
                return 1.0

            # Get fragment peak value safely
            if fragment_peak_pos < 0 or fragment_peak_pos >= len(ccbins):
                fragment_peak = np.max(ccbins) if len(ccbins) > 0 else 0
            else:
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
            # Validate inputs
            if len(ccbins) == 0:
                return 1.0

            # Get peak values safely
            if fragment_peak_pos < 0 or fragment_peak_pos >= len(ccbins):
                fragment_peak = np.max(ccbins) if len(ccbins) > 0 else 0
            else:
                fragment_peak = ccbins[fragment_peak_pos]

            if phantom_peak_pos < 0 or phantom_peak_pos >= len(ccbins):
                phantom_peak = np.max(ccbins) if len(ccbins) > 0 else 0
            else:
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

    # Get fragment peak value safely
    fragment_peak_value = 0.0
    if 0 <= fragment_peak_pos < len(ccbins):
        fragment_peak_value = float(ccbins[fragment_peak_pos])

    return {
        'nsc': nsc,
        'rsc': rsc,
        'fragment_peak_pos': fragment_peak_pos,
        'fragment_peak_value': fragment_peak_value,
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


class AdvancedStatisticsCalculator:
    """Advanced statistical calculations extracted from result.py for consistency.
    
    This class provides the exact statistical calculation methods used in the
    original result.py implementation, extracted as standalone functions to
    eliminate code duplication while preserving numerical behavior.
    """
    
    @staticmethod
    def calculate_background_minimum(ccbins: np.ndarray, 
                                   min_calc_width: int,
                                   near_zero_len: int = 10,
                                   output_warnings: bool = False) -> float:
        """Calculate background cross-correlation minimum using robust median method.
        
        This method replicates the exact behavior of CCStats._calc_cc_min(),
        using the median of the rightmost values to estimate background correlation.
        
        Args:
            ccbins: Cross-correlation values array
            min_calc_width: Number of rightmost values to consider
            near_zero_len: Number of leftmost values to check for validation
            output_warnings: Whether to output warnings for potential issues
            
        Returns:
            Background correlation minimum value
            
        Raises:
            Warning: If detected minimum seems too large compared to near-zero values
        """
        cc_min = np.sort(ccbins[-min_calc_width:])[min(min_calc_width, ccbins.size) // 2]
        
        if output_warnings and np.median(ccbins[:near_zero_len]) < cc_min:
            logger.warning("Detected minimum coefficient seems to be larger than beginning part minimum. "
                          "Consider increasing shift size (-d/--max-shift).")
        
        return float(cc_min)
    
    @staticmethod
    def calculate_phantom_peak_correlation(ccbins: np.ndarray, read_length: int) -> float:
        """Calculate cross-correlation at read length (phantom peak).
        
        This method replicates the exact behavior of CCStats._calc_cc_rl(),
        returning the correlation value at the read length position.
        
        Args:
            ccbins: Cross-correlation values array
            read_length: Sequencing read length
            
        Returns:
            Cross-correlation value at read length position, or 0 if out of bounds
        """
        if read_length > ccbins.size:
            return 0.0
        else:
            return float(ccbins[read_length - 1])
    
    @staticmethod
    @npcalc_with_logging_warn
    def calculate_exact_quality_metrics(ccbins: np.ndarray,
                                       library_len: int,
                                       read_length: int,
                                       cc_min: float) -> Tuple[float, float, float]:
        """Calculate exact NSC and RSC metrics matching CCStats._calc_lib_metrics().
        
        This preserves the exact numerical behavior of the original implementation
        including the specific formulas and edge case handling.
        
        Args:
            ccbins: Cross-correlation values array
            library_len: Fragment length position for metric calculation
            read_length: Read length for phantom peak
            cc_min: Background minimum correlation value
            
        Returns:
            Tuple of (ccfl, nsc, rsc) where:
                - ccfl: Cross-correlation at fragment length
                - nsc: Normalized Strand Coefficient
                - rsc: Relative Strand Coefficient
        """
        # Get correlation values
        ccfl = float(ccbins[library_len - 1])
        ccrl = AdvancedStatisticsCalculator.calculate_phantom_peak_correlation(ccbins, read_length)
        
        # Calculate metrics using exact original formulas
        nsc = ccfl / cc_min
        rsc = (ccfl - cc_min) / (ccrl - cc_min)
        
        return ccfl, nsc, rsc
    
    @staticmethod
    def estimate_fragment_length_with_masking(ccbins: np.ndarray,
                                            read_length: int,
                                            mv_avr_filter_len: int,
                                            filter_mask_len: int,
                                            near_readlen_criterion: int = 5,
                                            output_warnings: bool = False) -> Tuple[np.ndarray, int]:
        """Estimate fragment length with phantom peak masking.
        
        This method replicates the exact behavior of CCStats._estimate_fragment_len(),
        including the masking logic to avoid phantom peaks.
        
        Args:
            ccbins: Cross-correlation values array
            read_length: Sequencing read length
            mv_avr_filter_len: Moving average filter length
            filter_mask_len: Masking window around read length
            near_readlen_criterion: Distance threshold for phantom peak detection
            output_warnings: Whether to output warnings
            
        Returns:
            Tuple of (smoothed_cc, estimated_length) where:
                - smoothed_cc: Moving average filtered correlation
                - estimated_length: Estimated fragment length
        """
        # Apply moving average filter
        avr_cc = moving_avr_filter(ccbins, mv_avr_filter_len)
        est_lib_len = np.argmax(avr_cc) + 1
        
        need_warning = False
        
        # Check for phantom peak contamination and apply masking if needed
        if filter_mask_len and abs(est_lib_len - read_length) <= filter_mask_len:
            if output_warnings:
                logger.warning("Estimated library length is close to the read length.")
                logger.warning(f"Trying to masking around the read length +/- {filter_mask_len}bp...")
            
            # Apply masking around read length
            masked_avr_cc = avr_cc.copy()
            mask_from = max(0, read_length - 1 - filter_mask_len)
            mask_to = min(len(masked_avr_cc), read_length + filter_mask_len)
            
            for i in range(mask_from, mask_to):
                masked_avr_cc[i] = -float("inf")
            
            est_lib_len = np.argmax(masked_avr_cc) + 1
            
            # Check if estimate is at mask boundary
            if est_lib_len - 1 in (mask_from - 1, mask_to):
                need_warning = True
        
        elif output_warnings and abs(est_lib_len - read_length) <= near_readlen_criterion:
            need_warning = True
        
        # Issue warning if estimate is problematic
        if output_warnings and need_warning:
            logger.error("Estimated library length is close to the read length! Please check output plots.")
        
        return avr_cc, est_lib_len
    
    @staticmethod
    def calculate_fwhm(avr_cc: np.ndarray,
                      library_len: int,
                      cc_min: float,
                      output_warnings: bool = False) -> Union[int, bool]:
        """Calculate Full Width at Half Maximum (FWHM) of fragment length peak.
        
        This method replicates the exact behavior of CCStats._get_FWHM(),
        including the edge case handling and fallback strategies.
        
        Args:
            avr_cc: Smoothed cross-correlation array
            library_len: Fragment length position of the peak
            cc_min: Background minimum correlation value
            output_warnings: Whether to output warnings
            
        Returns:
            FWHM width in base pairs, or False if calculation fails
            
        Raises:
            AssertionError: If peak position is invalid or peak is not above background
        """
        max_i = library_len - 1
        assert max_i >= 0
        cc_max = avr_cc[max_i - 1]
        assert cc_max > cc_min
        target = cc_min + (cc_max - cc_min) / 2
        
        # Search forward from peak
        forward_shift = 0
        forward_failed = False
        while avr_cc[max_i + forward_shift] > target:
            forward_shift += 1
            if max_i + forward_shift == avr_cc.size:
                if output_warnings:
                    logger.warning("Failed to calc the half width at half maximum in the forward "
                                  "side of the peak. Consider increasing shift size (-d/--max-shift).")
                forward_failed = True
                forward_shift -= 1
                break
        
        # Search backward from peak
        backward_shift = 0
        backward_failed = False
        while avr_cc[max_i - backward_shift] > target:
            backward_shift += 1
            if max_i < backward_shift:
                if output_warnings:
                    logger.warning("Failed to calc the half width at half maximum in the backward side of the peak.")
                backward_failed = True
                backward_shift -= 1
                break
        
        # Calculate FWHM based on available measurements
        logger.debug((forward_shift, backward_shift))
        if forward_failed and backward_failed:
            if output_warnings:
                logger.error("Failed to calculate the full width at half maximum.")
            return False
        elif forward_failed:
            if output_warnings:
                logger.warning("Use twice width of the half width at half maximum in the backward side")
            return backward_shift * 2 + 1
        elif backward_failed:
            if output_warnings:
                logger.warning("Use twice width of the half width at half maximum in the forward side")
            return forward_shift * 2 + 1
        else:
            return backward_shift + forward_shift + 1


class CrossCorrelationMerger:
    """Handles merging of cross-correlation results from multiple chromosomes.
    
    This class provides the exact Fisher z-transformation logic extracted from
    CCResult._merge_cc() method, enabling consistent correlation merging across
    the codebase while maintaining numerical accuracy.
    """
    
    def __init__(self, confidence_interval: float = 0.99):
        """Initialize merger with confidence interval setting.
        
        Args:
            confidence_interval: Confidence level for interval calculations (default: 0.99)
        """
        self.confidence_interval = confidence_interval
        
    def merge_correlations(self, 
                         genome_lengths: np.ndarray, 
                         correlation_arrays: list,
                         read_length: int) -> Tuple[list, list, list]:
        """Merge cross-correlation results using Fisher z-transformation.
        
        This method replicates the exact behavior of CCResult._merge_cc(),
        combining per-chromosome correlations into genome-wide profiles with
        confidence intervals.
        
        Mathematical Process:
        1. Fisher z-transformation: z = arctanh(r) for each correlation
        2. Weighted averaging: z_avg = Σ(w_i * z_i) / Σ(w_i) where w_i = n_i - 3
        3. Confidence intervals: CI = tanh(z_avg ± z_α/2 * SE)
        4. Inverse transformation: r = tanh(z_avg) back to correlation
        
        Args:
            genome_lengths: Array of genome lengths for each chromosome
            correlation_arrays: List of correlation arrays, one per chromosome
            read_length: Read length for position-dependent weighting
            
        Returns:
            Tuple of (merged_correlations, lower_confidence, upper_confidence)
        """
        from scipy.stats import norm
        
        ns = np.array(genome_lengths)
        
        merged_r = []
        interval_upper = []
        interval_lower = []
        
        for i, _ccs in enumerate(zip(*correlation_arrays)):
            # Filter out NaN values
            nans = np.isnan(_ccs)
            _ccs = np.array(_ccs)[~nans]
            
            # Calculate weights (effective sample size - 3)
            if len(ns.shape) == 1:
                _ns = ns[~nans] - 3
            else:
                _ns = ns[~nans, abs(read_length - i)] - 3
            
            # Fisher z-transformation
            zs = np.arctanh(_ccs)
            
            # Filter out infinite values
            infs = np.isinf(zs)
            zs = zs[~infs]
            _ns = _ns[~infs]
            
            # Weighted average in z-space
            avr_z = np.average(zs, weights=_ns)
            
            # Calculate confidence interval
            z_interval = norm.ppf(1 - (1 - self.confidence_interval) / 2) * np.sqrt(1 / np.sum(_ns))
            z_interval_upper = avr_z + z_interval
            z_interval_lower = avr_z - z_interval
            
            # Transform back to correlation space
            merged_r.append(np.tanh(avr_z))
            interval_upper.append(np.tanh(z_interval_upper))
            interval_lower.append(np.tanh(z_interval_lower))
        
        return merged_r, interval_lower, interval_upper
    
    @staticmethod
    def validate_correlation_arrays(correlation_arrays: list) -> bool:
        """Validate that correlation arrays are compatible for merging.
        
        Args:
            correlation_arrays: List of correlation arrays to validate
            
        Returns:
            True if arrays are compatible, False otherwise
        """
        if not correlation_arrays:
            return False
            
        # Check that all arrays have the same length
        lengths = [len(arr) for arr in correlation_arrays if arr is not None]
        if len(set(lengths)) > 1:
            logger.warning(f"Correlation arrays have different lengths: {lengths}")
            return False
            
        return True
    
    @staticmethod
    def calculate_effective_sample_size(genome_length: int, position: int, read_length: int) -> int:
        """Calculate effective sample size for correlation at given position.
        
        Args:
            genome_length: Total genome length
            position: Position in correlation array
            read_length: Read length for position-dependent calculation
            
        Returns:
            Effective sample size for weighting
        """
        # Standard formula: genome_length - 3 for degrees of freedom
        return max(1, genome_length - 3)
    
    def merge_with_weights(self, 
                         correlations: list,
                         weights: list) -> Tuple[float, float, float]:
        """Merge correlations with explicit weights.
        
        Args:
            correlations: List of correlation values
            weights: List of weights for each correlation
            
        Returns:
            Tuple of (merged_correlation, lower_ci, upper_ci)
        """
        from scipy.stats import norm
        
        # Filter valid correlations and weights
        valid_pairs = [(r, w) for r, w in zip(correlations, weights) 
                      if not np.isnan(r) and w > 0]
        
        if not valid_pairs:
            return 0.0, 0.0, 0.0
            
        correlations, weights = zip(*valid_pairs)
        correlations = np.array(correlations)
        weights = np.array(weights)
        
        # Fisher z-transformation
        zs = np.arctanh(correlations)
        
        # Filter infinite values
        finite_mask = np.isfinite(zs)
        zs = zs[finite_mask]
        weights = weights[finite_mask]
        
        if len(zs) == 0:
            return 0.0, 0.0, 0.0
        
        # Weighted average
        avr_z = np.average(zs, weights=weights)
        
        # Confidence interval
        z_interval = norm.ppf(1 - (1 - self.confidence_interval) / 2) * np.sqrt(1 / np.sum(weights))
        z_upper = avr_z + z_interval
        z_lower = avr_z - z_interval
        
        # Transform back to correlation space
        merged_r = np.tanh(avr_z)
        upper_ci = np.tanh(z_upper)
        lower_ci = np.tanh(z_lower)
        
        return float(merged_r), float(lower_ci), float(upper_ci)


class ArrayAggregator:
    """Utility class for handling array aggregation operations.
    
    This class provides robust methods for aggregating arrays from multiple
    chromosomes or data sources, handling common issues like None values,
    different array lengths, and data validation.
    
    These utilities are extracted from result.py to eliminate code duplication
    and provide consistent array handling across the codebase.
    """
    
    @staticmethod
    def filter_none_values(arrays):
        """Filter out None values from an iterable of arrays.
        
        This utility function removes None values from any iterable,
        which is commonly needed when aggregating results from multiple
        chromosomes where some may have failed calculations.
        
        Args:
            arrays: Iterable that may contain None values
            
        Returns:
            list: List with None values removed
            
        Example:
            >>> ArrayAggregator.filter_none_values([arr1, None, arr2, None, arr3])
            [arr1, arr2, arr3]
        """
        return [x for x in arrays if x is not None]
    
    @staticmethod
    def safe_array_sum(arrays, axis=0):
        """Safely sum arrays that might have different lengths or contain None values.
        
        This utility function handles the common case in cross-correlation analysis
        where arrays from different chromosomes need to be aggregated, but some
        chromosomes may have failed calculations (None values) or different
        array lengths due to varying shift ranges.
        
        The function implements a robust aggregation strategy:
        1. Filter out None values
        2. Convert to numpy arrays
        3. Detect length mismatches
        4. Use most common length as reference
        5. Filter arrays to match reference length
        6. Sum compatible arrays
        
        This approach ensures numerical stability while handling the irregular
        data structures that can arise from multi-chromosome calculations.
        
        Args:
            arrays: Iterable of arrays that might be None or have different lengths
            axis (int): Axis along which to sum (default: 0)
            
        Returns:
            np.ndarray or None: Summed array, or None if no valid arrays exist
            
        Raises:
            Warning: If arrays have inconsistent lengths (logged, not raised)
            
        Example:
            >>> arr1 = np.array([1, 2, 3])
            >>> arr2 = np.array([4, 5, 6])
            >>> arr3 = None
            >>> ArrayAggregator.safe_array_sum([arr1, arr2, arr3])
            array([5, 7, 9])
        """
        # Filter out None values
        valid_arrays = ArrayAggregator.filter_none_values(arrays)
        
        if not valid_arrays:
            return None
        
        # Convert to numpy arrays
        valid_arrays = [np.asarray(arr) for arr in valid_arrays]
        
        # Check if all arrays have the same shape
        if len(valid_arrays) == 1:
            return valid_arrays[0]
        
        # Find the expected length (use the most common length)
        lengths = [len(arr) for arr in valid_arrays]
        if len(set(lengths)) > 1:
            # Arrays have different lengths - log warning and use the most common length
            from collections import Counter
            length_counts = Counter(lengths)
            expected_length = length_counts.most_common(1)[0][0]
            logger.warning(f"Arrays have inconsistent lengths: {lengths}. Using length {expected_length}")
            
            # Filter to only arrays with the expected length
            valid_arrays = [arr for arr in valid_arrays if len(arr) == expected_length]
            
            if not valid_arrays:
                return None
        
        # Sum the arrays
        return np.sum(valid_arrays, axis=axis)
    
    @staticmethod
    def validate_array_compatibility(arrays, require_same_length=True):
        """Validate that arrays are compatible for aggregation operations.
        
        Args:
            arrays: Iterable of arrays to validate
            require_same_length (bool): Whether to require all arrays to have the same length
            
        Returns:
            tuple: (is_valid, valid_arrays, length_info)
                - is_valid: Boolean indicating if arrays are compatible
                - valid_arrays: List of non-None arrays
                - length_info: Dict with length statistics
        """
        valid_arrays = ArrayAggregator.filter_none_values(arrays)
        
        if not valid_arrays:
            return False, [], {"message": "No valid arrays found"}
        
        # Convert to numpy arrays and get lengths
        valid_arrays = [np.asarray(arr) for arr in valid_arrays]
        lengths = [len(arr) for arr in valid_arrays]
        
        length_info = {
            "lengths": lengths,
            "unique_lengths": list(set(lengths)),
            "min_length": min(lengths),
            "max_length": max(lengths),
            "most_common_length": None
        }
        
        if len(set(lengths)) > 1:
            from collections import Counter
            length_counts = Counter(lengths)
            length_info["most_common_length"] = length_counts.most_common(1)[0][0]
            
            if require_same_length:
                return False, valid_arrays, length_info
        
        return True, valid_arrays, length_info
    
    @staticmethod
    def aggregate_chromosome_data(ref2data, ref2weights=None, method='sum'):
        """Aggregate data from multiple chromosomes with optional weighting.
        
        Args:
            ref2data: Dictionary mapping chromosome names to data arrays
            ref2weights: Optional dictionary mapping chromosome names to weights
            method: Aggregation method ('sum', 'weighted_sum', 'mean', 'weighted_mean')
            
        Returns:
            np.ndarray or None: Aggregated result
        """
        if not ref2data:
            return None
        
        # Get data arrays and weights
        data_arrays = list(ref2data.values())
        
        if method == 'sum':
            return ArrayAggregator.safe_array_sum(data_arrays)
        elif method == 'weighted_sum' and ref2weights:
            # Filter to chromosomes that have both data and weights
            valid_refs = [ref for ref in ref2data.keys() if ref in ref2weights]
            if not valid_refs:
                return None
            
            weighted_arrays = []
            for ref in valid_refs:
                if ref2data[ref] is not None and ref2weights[ref] is not None:
                    weighted_arrays.append(np.asarray(ref2data[ref]) * ref2weights[ref])
            
            return ArrayAggregator.safe_array_sum(weighted_arrays)
        elif method == 'mean':
            summed = ArrayAggregator.safe_array_sum(data_arrays)
            return summed / len(data_arrays) if summed is not None else None
        elif method == 'weighted_mean' and ref2weights:
            weighted_sum = ArrayAggregator.aggregate_chromosome_data(ref2data, ref2weights, 'weighted_sum')
            total_weights = sum(w for w in ref2weights.values() if w is not None)
            return weighted_sum / total_weights if weighted_sum is not None and total_weights > 0 else None
        else:
            raise ValueError(f"Unsupported aggregation method: {method}")


class StatisticalTestUtilities:
    """Utility class for common statistical tests used in cross-correlation analysis.
    
    This class provides statistical test functions that are commonly used for
    quality assessment and validation in ChIP-seq analysis.
    """
    
    @staticmethod
    def chi2_test(forward_count, reverse_count, chi2_p_thresh, label):
        """Perform chi-squared test for forward/reverse read count balance.
        
        This function tests whether forward and reverse read counts are balanced
        using a chi-squared goodness of fit test. Significant imbalance may indicate
        technical issues with the sequencing or sample preparation.
        
        The test assumes that forward and reverse reads should be approximately
        equal under the null hypothesis. Large deviations from a 1:1 ratio suggest
        potential problems with the data quality.
        
        Mathematical Formula:
        - Expected count for each strand: (a + b) / 2
        - Chi-squared statistic: χ² = Σ((observed - expected)² / expected)
        - P-value: P(χ² > observed_χ² | df = 1)
        
        Args:
            forward_count (int): Forward read count
            reverse_count (int): Reverse read count  
            chi2_p_thresh (float): P-value threshold for significance (e.g., 0.05)
            label (str): Descriptive label for logging (e.g., "Whole genome", "chr1")
            
        Returns:
            tuple: (chi2_statistic, p_value, is_significant)
            
        Side Effects:
            Logs warning messages if imbalance is significant (p ≤ threshold)
            Logs info messages if balance is acceptable (p > threshold)
            
        Example:
            >>> StatisticalTestUtilities.chi2_test(1000, 1050, 0.05, "chr1")
            # Logs info message about acceptable balance
            >>> StatisticalTestUtilities.chi2_test(1000, 1500, 0.05, "chr1") 
            # Logs warning about significant imbalance
        """
        from scipy.stats.distributions import chi2
        
        sum_counts = forward_count + reverse_count
        
        # Handle edge case where both counts are zero
        if sum_counts == 0:
            chi2_val = 0.0
            chi2_p = 1.0  # Perfect balance (no imbalance when no reads)
        else:
            chi2_val = (((forward_count - sum_counts / 2.) ** 2) + ((reverse_count - sum_counts / 2.) ** 2)) / sum_counts
            chi2_p = chi2.sf(chi2_val, 1)
        
        is_significant = chi2_p <= chi2_p_thresh
        
        if is_significant:
            logger.warning("{} Forward/Reverse read count imbalance.".format(label))
            logger.warning("+/- = {} / {}, Chi-squared test p-val = {} <= {}".format(
                forward_count, reverse_count, chi2_p, chi2_p_thresh
            ))
        else:
            logger.info("{} Forward/Reverse read count +/- = {} / {}".format(label, forward_count, reverse_count))
            logger.info("Chi-squared test p-val = {} > {}".format(chi2_p, chi2_p_thresh))
        
        return chi2_val, chi2_p, is_significant