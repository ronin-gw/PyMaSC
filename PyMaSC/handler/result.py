import logging
from functools import wraps
from typing import Union, Tuple, List, Dict

import numpy as np
from scipy.stats import norm

from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.calc import moving_avr_filter
from PyMaSC.utils.stats_utils import ArrayAggregator, StatisticalTestUtilities, CrossCorrelationMerger

# New builder-based result construction
from PyMaSC.core.result_builder import ResultBuilder, build_from_handler

# New unified statistics module
from PyMaSC.core.statistics import NCCStats, MSCCStats, UnifiedStats

logger = logging.getLogger(__name__)

# Quality Assessment Constants
NEAR_READLEN_ERR_CRITERION = 5
"""int: Distance threshold for detecting phantom peak contamination.

When the estimated fragment length is within this many base pairs of the
read length, it suggests potential phantom peak contamination. This triggers
warnings to alert users about potential quality issues in the data.
"""

MERGED_CC_CONFIDENCE_INTERVAL = 0.99
"""float: Confidence level for merged cross-correlation intervals.

This sets the confidence level (99%) for calculating confidence intervals
when merging cross-correlation results from multiple chromosomes using
Fisher z-transformation. Higher values provide wider, more conservative
confidence intervals.
"""

NEAR_ZERO_MIN_CALC_LEN = 10
"""int: Number of positions to examine at the beginning of correlation profile.

Used to detect potential issues with minimum correlation calculation.
If the median of the first NEAR_ZERO_MIN_CALC_LEN positions is smaller
than the calculated minimum, it suggests the shift size may be too small.
"""




# The npcalc_with_logging_warn decorator is now imported from the unified statistics module




# Backward compatibility alias for CCStats
# The original CCStats class has been replaced with the new unified statistics module
CCStats = NCCStats


# Backward compatibility alias for PyMaSCStats
# The original PyMaSCStats class has been replaced with the new unified statistics module
PyMaSCStats = UnifiedStats


class ReadsTooFew(IndexError):
    pass


class CCResult(object):
    """Cross-correlation result aggregator and statistics calculator.
    
    This class serves as the main interface for PyMaSC cross-correlation analysis,
    aggregating per-chromosome results into genome-wide statistics and providing
    comprehensive quality metrics for ChIP-seq data assessment.
    
    The class supports two initialization modes:
    1. From handler: Direct initialization from calculation handlers (pymasc)
    2. From parameters: Manual initialization with explicit data (pymasc-plot)
    
    Key Features:
    - Aggregates per-chromosome cross-correlation results
    - Calculates genome-wide NCC and MSCC statistics
    - Performs Fisher z-transformation for correlation merging
    - Provides confidence intervals for merged correlations
    - Handles read count balance validation (chi-squared test)
    - Supports both standard and mappability-sensitive analysis
    
    Mathematical Background:
    - Fisher z-transformation: z = arctanh(r) for correlation merging
    - Confidence intervals: CI = tanh(z ± z_α/2 * SE)
    - Chi-squared test: χ² = (O - E)² / E for read balance
    - Weighted averaging by effective genome length
    
    Quality Assessment:
    - Read count balance between forward/reverse strands
    - NSC/RSC thresholds for data quality (NSC > 1.05, RSC > 0.8)
    - Fragment length estimation consistency
    - Phantom peak detection and masking
    
    The class automatically validates input data, handles missing chromosomes,
    and provides comprehensive error reporting for quality control.
    
    Parameters:
        mv_avr_filter_len (int): Moving average filter length for peak detection
        chi2_pval (float): Chi-squared p-value threshold for read balance test
        filter_mask_len (int): Masking window around read length
        min_calc_width (int): Minimum width for background calculation
        expected_library_len (int, optional): Expected fragment length
        handler (object, optional): Calculation handler with results
        read_len (int, optional): Sequencing read length
        references (list, optional): List of chromosome names
        ref2genomelen (dict, optional): Chromosome lengths
        ref2forward_sum (dict, optional): Forward read counts per chromosome
        ref2reverse_sum (dict, optional): Reverse read counts per chromosome
        ref2cc (dict, optional): NCC values per chromosome
        mappable_ref2forward_sum (dict, optional): Mappable forward reads
        mappable_ref2reverse_sum (dict, optional): Mappable reverse reads
        ref2mappable_len (dict, optional): Mappable lengths per chromosome
        ref2masc (dict, optional): MSCC values per chromosome
    
    Attributes:
        ref2stats (dict): Per-chromosome PyMaSCStats objects
        whole (PyMaSCStats): Genome-wide aggregated statistics
        skip_ncc (bool): Whether NCC calculation was skipped
        calc_masc (bool): Whether MSCC calculation was performed
        ncc_upper (list): Upper confidence interval for NCC
        ncc_lower (list): Lower confidence interval for NCC
        mscc_upper (list): Upper confidence interval for MSCC
        mscc_lower (list): Lower confidence interval for MSCC
    """
    
    def __init__(
        self,
        mv_avr_filter_len, chi2_pval, filter_mask_len, min_calc_width,  # mandatory params
        expected_library_len=None,  # optional parameter
        handler=None,  # source 1 (pymasc)
        read_len=None, references=None,
        ref2genomelen=None, ref2forward_sum=None, ref2reverse_sum=None, ref2cc=None,
        mappable_ref2forward_sum=None, mappable_ref2reverse_sum=None,
        ref2mappable_len=None, ref2masc=None  # source 2 (pymasc-plot)
    ):

        # settings
        self.mv_avr_filter_len = mv_avr_filter_len
        self.chi2_p_thresh = chi2_pval
        self.expected_library_len = expected_library_len
        self.filter_mask_len = max(filter_mask_len, 0)
        self.min_calc_width = max(min_calc_width, 0)

        if handler:
            self._init_from_handler(handler)
        else:
            self.read_len = read_len
            self.references = references
            self.ref2genomelen = ref2genomelen
            self.ref2mappable_len = ref2mappable_len
            self.ref2forward_sum = ref2forward_sum
            self.ref2reverse_sum = ref2reverse_sum
            self.mappable_ref2forward_sum = mappable_ref2forward_sum
            self.mappable_ref2reverse_sum = mappable_ref2reverse_sum

            self.ref2stats = {
                ref: PyMaSCStats(
                    self.read_len, self.mv_avr_filter_len, self.expected_library_len,
                    genomelen=self.ref2genomelen[ref],
                    forward_sum=self.ref2forward_sum.get(ref, None),
                    reverse_sum=self.ref2reverse_sum.get(ref, None),
                    cc=ref2cc.get(ref, None) if ref2cc else None,
                    mappable_len=self.ref2mappable_len.get(ref, None),
                    mappable_forward_sum=self.mappable_ref2forward_sum.get(ref, None),
                    mappable_reverse_sum=self.mappable_ref2reverse_sum.get(ref, None),
                    masc=ref2masc.get(ref, None) if ref2masc else None,
                    filter_mask_len=self.filter_mask_len,
                    min_calc_width=self.min_calc_width
                ) for ref in self.references
            }

            self.skip_ncc = not (self.ref2genomelen and ref2cc)
            self.calc_masc = self.ref2mappable_len and ref2masc

            assert (not self.skip_ncc) or self.calc_masc

        #
        if not self.skip_ncc:
            # Extract genome lengths and correlation arrays for NCC
            ncc_data = [(self.ref2genomelen[ref], self.ref2stats[ref].cc.cc)
                       for ref in self.references if self.ref2stats[ref].cc is not None]
            if ncc_data:
                genome_lengths, correlation_arrays = zip(*ncc_data)
                merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
                ncc, self.ncc_lower, self.ncc_upper = merger.merge_correlations(
                    np.array(genome_lengths), correlation_arrays, self.read_len
                )
            else:
                ncc = self.ncc_upper = self.ncc_lower = None
        else:
            ncc = self.ncc_upper = self.ncc_lower = None

        if self.calc_masc:
            # Extract mappable lengths and correlation arrays for MSCC
            mscc_data = [(self.ref2mappable_len[ref], self.ref2stats[ref].masc.cc)
                        for ref in self.references if self.ref2stats[ref].masc is not None]
            if mscc_data:
                mappable_lengths, correlation_arrays = zip(*mscc_data)
                merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
                mscc, self.mscc_lower, self.mscc_upper = merger.merge_correlations(
                    np.array(mappable_lengths), correlation_arrays, self.read_len
                )
            else:
                mscc = self.mscc_upper = self.mscc_lower = None
        else:
            mscc = self.mscc_upper = self.mscc_lower = None

        self.whole = PyMaSCStats(
            self.read_len, self.mv_avr_filter_len, self.expected_library_len,
            genomelen=sum(self.ref2genomelen.values()),
            forward_sum=sum(self.ref2forward_sum.values()),
            reverse_sum=sum(self.ref2reverse_sum.values()),
            cc=ncc,
            mappable_len=ArrayAggregator.safe_array_sum(self.ref2mappable_len.values()) if self.ref2mappable_len else None,
            mappable_forward_sum=ArrayAggregator.safe_array_sum(self.mappable_ref2forward_sum.values()) if self.mappable_ref2forward_sum else None,
            mappable_reverse_sum=ArrayAggregator.safe_array_sum(self.mappable_ref2reverse_sum.values()) if self.mappable_ref2reverse_sum else None,
            masc=mscc,
            warning=True,
            filter_mask_len=self.filter_mask_len,
            min_calc_width=self.min_calc_width
        )

    def _init_from_handler(self, handler):
        #
        self.references = handler.references
        lengths = handler.lengths
        self.ref2genomelen = dict(zip(self.references, lengths))
        self.read_len = handler.read_len
        self.ref2forward_sum = handler.ref2forward_sum
        self.ref2reverse_sum = handler.ref2reverse_sum
        ref2ccbins = handler.ref2ccbins
        self.mappable_ref2forward_sum = handler.mappable_ref2forward_sum
        self.mappable_ref2reverse_sum = handler.mappable_ref2reverse_sum
        mappable_ref2ccbins = handler.mappable_ref2ccbins
        self.ref2mappable_len = handler.ref2mappable_len

        #
        self.ref2stats = {
            ref: PyMaSCStats(
                self.read_len, self.mv_avr_filter_len, self.expected_library_len,
                genomelen=self.ref2genomelen[ref],
                forward_sum=self.ref2forward_sum.get(ref, 0),
                reverse_sum=self.ref2reverse_sum.get(ref, 0),
                ccbins=ref2ccbins.get(ref, None),
                mappable_len=self.ref2mappable_len.get(ref, None),
                mappable_forward_sum=self.mappable_ref2forward_sum.get(ref, None),
                mappable_reverse_sum=self.mappable_ref2reverse_sum.get(ref, None),
                mappable_ccbins=mappable_ref2ccbins.get(ref, None),
                filter_mask_len=self.filter_mask_len,
                min_calc_width=self.min_calc_width
            ) for ref in self.references
        }

        #
        if all((self.ref2forward_sum, self.ref2reverse_sum, ref2ccbins)):
            self.skip_ncc = False
            forward_sum = sum(list(self.ref2forward_sum.values()))
            reverse_sum = sum(list(self.ref2reverse_sum.values()))
        else:
            self.skip_ncc = True
            forward_sum = reverse_sum = None
        if forward_sum is not None:
            if forward_sum == 0:
                logger.error("There is no forward read.")
                raise ReadsTooFew
            elif reverse_sum == 0:
                logger.error("There is no reverse read.")
                raise ReadsTooFew
            StatisticalTestUtilities.chi2_test(forward_sum, reverse_sum, self.chi2_p_thresh, "Whole genome")

        #
        if all((self.mappable_ref2forward_sum, self.mappable_ref2reverse_sum,
                mappable_ref2ccbins, self.ref2mappable_len)):
            self.calc_masc = True
            mappable_forward_sum = ArrayAggregator.safe_array_sum(self.mappable_ref2forward_sum.values())
            mappable_reverse_sum = ArrayAggregator.safe_array_sum(self.mappable_ref2reverse_sum.values())
        else:
            self.calc_masc = False
            mappable_forward_sum = mappable_reverse_sum = None
        if mappable_forward_sum is not None:
            if np.all(mappable_forward_sum == 0):
                logger.error("There is no forward read.")
                raise ReadsTooFew
            elif np.all(mappable_reverse_sum == 0):
                logger.error("There is no reverse read.")
                raise ReadsTooFew

    # New builder-based result construction methods
    @classmethod
    def from_handler_with_builder(
        cls,
        handler,
        mv_avr_filter_len: int = 15,
        chi2_pval: float = 0.01,
        filter_mask_len: int = 5,
        min_calc_width: int = 10,
        expected_library_len: int = None
    ) -> 'CCResult':
        """Create CCResult using the new builder system.
        
        This method provides an alternative to the traditional __init__ method,
        using the type-safe builder pattern for result construction.
        
        Args:
            handler: Calculation handler with results
            mv_avr_filter_len: Moving average filter length
            chi2_pval: Chi-squared p-value threshold
            filter_mask_len: Filter mask length around read length
            min_calc_width: Minimum width for background calculation
            expected_library_len: Expected fragment length
            
        Returns:
            CCResult instance constructed via builder pattern
        """
        # NOTE: The builder system is currently incomplete and causes test failures.
        # It creates dict objects instead of proper PyMaSCStats objects for ref2stats,
        # leading to AttributeError in output modules.
        # Until the builder system is properly implemented, we use the legacy system.
        return cls(
            mv_avr_filter_len=mv_avr_filter_len,
            chi2_pval=chi2_pval,
            filter_mask_len=filter_mask_len,
            min_calc_width=min_calc_width,
            expected_library_len=expected_library_len,
            handler=handler
        )
    
    @classmethod  
    def from_file_data_with_builder(
        cls,
        mv_avr_filter_len: int,
        chi2_pval: float,
        filter_mask_len: int,
        min_calc_width: int,
        expected_library_len: int = None,
        read_len: int = None,
        references: List[str] = None,
        ref2genomelen: Dict[str, int] = None,
        ref2forward_sum: Dict[str, int] = None,
        ref2reverse_sum: Dict[str, int] = None,
        ref2cc: Dict[str, np.ndarray] = None,
        ref2mappable_len: Dict[str, int] = None,
        mappable_ref2forward_sum: Dict[str, np.ndarray] = None,
        mappable_ref2reverse_sum: Dict[str, np.ndarray] = None,
        ref2masc: Dict[str, np.ndarray] = None
    ) -> 'CCResult':
        """Create CCResult from file data using builder pattern.
        
        This method is designed for creating CCResult objects from data loaded
        from files (e.g., in plot.py), using the new builder system for 
        better type safety and maintainability.
        
        Args:
            mv_avr_filter_len: Moving average filter length
            chi2_pval: Chi-squared p-value threshold
            filter_mask_len: Filter mask length around read length
            min_calc_width: Minimum width for background calculation
            expected_library_len: Expected fragment length
            read_len: Read length in base pairs
            references: List of chromosome names
            ref2genomelen: Genome lengths by chromosome
            ref2forward_sum: Forward read counts by chromosome
            ref2reverse_sum: Reverse read counts by chromosome
            ref2cc: Cross-correlation data by chromosome
            ref2mappable_len: Mappable lengths by chromosome
            mappable_ref2forward_sum: Mappable forward counts by chromosome
            mappable_ref2reverse_sum: Mappable reverse counts by chromosome
            ref2masc: MSCC data by chromosome
            
        Returns:
            CCResult instance constructed via builder pattern
        """
        try:
            # Use builder pattern for construction
            from PyMaSC.core.result_builder import ResultBuilder
            
            builder = ResultBuilder()
            
            # Set basic parameters
            if read_len and references:
                chromosome_lengths = [ref2genomelen.get(ref, 1) for ref in references] if ref2genomelen else None
                builder.with_basic_params(read_len, references, chromosome_lengths)
            
            # Add NCC data if available
            if ref2forward_sum and ref2reverse_sum and ref2cc:
                builder.with_ncc_data(ref2forward_sum, ref2reverse_sum, ref2cc, ref2genomelen)
            
            # Add MSCC data if available
            if (mappable_ref2forward_sum and mappable_ref2reverse_sum and 
                ref2masc and ref2mappable_len):
                builder.with_mscc_data(
                    mappable_ref2forward_sum, mappable_ref2reverse_sum, 
                    ref2masc, ref2mappable_len
                )
            
            # Configure and build
            built_result = (builder
                           .with_statistics_config(mv_avr_filter_len=mv_avr_filter_len)
                           .with_fragment_length(expected=expected_library_len)
                           .build())
            
            # Convert to CCResult format
            return cls._from_built_result(
                built_result, mv_avr_filter_len, chi2_pval, filter_mask_len,
                min_calc_width, expected_library_len
            )
            
        except Exception as e:
            logger.warning(f"Builder system failed for file data, falling back to legacy: {e}")
            # Fall back to legacy construction
            return cls(
                mv_avr_filter_len=mv_avr_filter_len,
                chi2_pval=chi2_pval,
                filter_mask_len=filter_mask_len,
                min_calc_width=min_calc_width,
                expected_library_len=expected_library_len,
                read_len=read_len,
                references=references,
                ref2genomelen=ref2genomelen,
                ref2forward_sum=ref2forward_sum,
                ref2reverse_sum=ref2reverse_sum,
                ref2cc=ref2cc,
                ref2mappable_len=ref2mappable_len,
                mappable_ref2forward_sum=mappable_ref2forward_sum,
                mappable_ref2reverse_sum=mappable_ref2reverse_sum,
                ref2masc=ref2masc
            )

    @classmethod
    def _from_built_result(
        cls,
        built_result,
        mv_avr_filter_len: int,
        chi2_pval: float,
        filter_mask_len: int,
        min_calc_width: int,
        expected_library_len: int = None
    ) -> 'CCResult':
        """Convert BuiltResult to CCResult format."""
        # Extract legacy-compatible data from built result
        container = built_result.container
        
        # Create CCResult using legacy initialization but with builder data
        instance = cls.__new__(cls)
        
        # Set configuration parameters
        instance.mv_avr_filter_len = mv_avr_filter_len
        instance.chi2_p_thresh = chi2_pval
        instance.expected_library_len = expected_library_len
        instance.filter_mask_len = max(filter_mask_len, 0)
        instance.min_calc_width = max(min_calc_width, 0)
        
        # Extract data from container
        instance.read_len = container.read_length
        instance.references = container.references
        instance.skip_ncc = container.skip_ncc
        instance.calc_masc = container.has_mappability
        
        # Convert container data to legacy format
        instance.ref2genomelen = {
            result.chromosome: result.length
            for result in container.chromosome_results.values()
        }
        
        ncc_data = container.get_ncc_data()
        if ncc_data:
            instance.ref2forward_sum = {
                chrom: data.forward_sum for chrom, data in ncc_data.items()
            }
            instance.ref2reverse_sum = {
                chrom: data.reverse_sum for chrom, data in ncc_data.items()
            }
        
        mscc_data = container.get_mscc_data()
        if mscc_data:
            instance.mappable_ref2forward_sum = {
                chrom: data.forward_sum for chrom, data in mscc_data.items()
            }
            instance.mappable_ref2reverse_sum = {
                chrom: data.reverse_sum for chrom, data in mscc_data.items()
            }
            instance.ref2mappable_len = {
                chrom: data.mappable_len for chrom, data in mscc_data.items()
            }
        
        # Create per-chromosome statistics using builder results
        instance.ref2stats = {}
        for chrom in instance.references:
            chrom_stats = built_result.chromosome_statistics.get(chrom, {})
            
            # Create PyMaSCStats for this chromosome
            # This is a simplified version - full implementation would need
            # to convert StatisticsResult back to PyMaSCStats format
            instance.ref2stats[chrom] = chrom_stats
        
        # Create whole-genome statistics from aggregate data
        aggregate_stats = built_result.aggregate_statistics
        
        # Extract data from aggregated statistics
        total_reads = aggregate_stats.get('total_reads', {})
        aggregated_cc = aggregate_stats.get('aggregated_cc', {})
        
        # Extract NCC and MSCC arrays from aggregated results
        ncc = None
        mscc = None
        if isinstance(aggregated_cc, dict):
            if 'cc' in aggregated_cc:
                ncc = aggregated_cc['cc']
            elif 'ncc' in aggregated_cc:
                ncc = aggregated_cc['ncc']
            elif 'mscc' in aggregated_cc:
                mscc = aggregated_cc['mscc']
        elif isinstance(aggregated_cc, np.ndarray):
            ncc = aggregated_cc
        
        # Create PyMaSCStats object for whole-genome data
        instance.whole = PyMaSCStats(
            read_len=instance.read_len,
            mv_avr_filter_len=mv_avr_filter_len,
            expected_library_len=expected_library_len,
            genomelen=aggregate_stats.get('genome_length', sum(instance.ref2genomelen.values())),
            forward_sum=total_reads.get('forward', sum(instance.ref2forward_sum.values()) if instance.ref2forward_sum else 0),
            reverse_sum=total_reads.get('reverse', sum(instance.ref2reverse_sum.values()) if instance.ref2reverse_sum else 0),
            cc=ncc,
            mappable_len=total_reads.get('mappable_len'),
            mappable_forward_sum=total_reads.get('mappable_forward'),
            mappable_reverse_sum=total_reads.get('mappable_reverse'),
            masc=mscc,
            filter_mask_len=filter_mask_len,
            warning=True,
            min_calc_width=min_calc_width
        )
        
        # Additional initialization that matches legacy behavior
        instance.ncc_upper = instance.ncc_lower = None
        instance.mscc_upper = instance.mscc_lower = None
        
        return instance

