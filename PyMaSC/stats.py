"""Statistics calculation and quality metrics for cross-correlation analysis.

Provides calculation of quality metrics, fragment length estimation,
and statistical analysis for NCC and MSCC results.
"""
import logging
from dataclasses import dataclass, field
from typing import Union, Optional, Type, Dict, List, Tuple, TypeVar, Protocol, Generic, overload, cast

import numpy as np
import numpy.typing as npt

from .interfaces.config import StatConfig
from .interfaces.result import (
    NCCResultModel, MSCCResultModel,
    GenomeWideResultModel,
    NCCGenomeWideResultModel, MSCCGenomeWideResultModel, BothGenomeWideResultModel
)
from .result import EmptyNCCResult, EmptyMSCCResult, EmptyBothChromResult
from .interfaces.stats import (
    CCQualityMetrics as CCQualityMetricsModel,
    CCStats as CCStatsModel,
    ChromosomeStats,
    WholeGenomeStats,
    GenomeWideStats,
    TCount
)
from .core.exceptions import ReadsTooFew

from .utils.calc import moving_avr_filter, merge_correlations

logger = logging.getLogger(__name__)

# Constants
NEAR_READLEN_ERR_CRITERION = 5
NEAR_ZERO_MIN_CALC_LEN = 10


@dataclass
class CCQualityMetrics(CCQualityMetricsModel):
    """Quality metrics for cross-correlation analysis with calculation methods.

    Concrete implementation that adds calculation logic to the base model.
    """
    fragment_length: Optional[int] = None
    ccfl: Optional[float] = None
    fwhm: Optional[int] = None

    #
    nsc: Optional[float] = None
    rsc: Optional[float] = None
    vsn: Optional[float] = None

    def calc_metrics(self, stats: 'CCStats') -> None:
        """Calculate NSC, RSC, and VSN quality metrics from statistics.

        Args:
            stats: CCStats object containing required values for calculation
        """
        if self.fragment_length is None:
            return
        else:
            assert isinstance(self.ccfl, float), "ccfl must be a float."
        self.nsc = self.ccfl / stats.cc_min
        self.rsc = (self.ccfl - stats.cc_min) / (stats.ccrl - stats.cc_min)
        if self.fwhm is not None:
            self.vsn = 2 * self.ccfl * self.fwhm / (stats.forward_reads_repr + stats.reverse_reads_repr)


@dataclass
class EmptyChromosomeStats:
    """Represents statistics for a chromosome with no data.

    This class is used to explicitly mark chromosomes that had no reads
    but need to be tracked for genome length consistency in parallel processing.
    """
    # Representative genome length for this chromosome
    genomelen_repr: int
    # These attributes maintain interface compatibility with ChromosomeStats
    stats: None = None
    cc: None = None
    avr_cc: None = None
    est_lib_len: int = 0
    mv_avr_filter_len: int = 0


@dataclass
class CCStats(CCStatsModel[TCount]):
    """Cross-correlation statistics with quality metrics calculation.

    Base implementation of CCStatsModel that automatically computes
    quality metrics after initialization.
    """
    read_len: int
    cc_min: float
    ccrl: float
    genomelen: TCount
    forward_reads: TCount
    reverse_reads: TCount
    metrics_at_expected_length: CCQualityMetrics
    metrics_at_estimated_length: CCQualityMetrics

    def __post_init__(self) -> None:
        self.metrics_at_expected_length.calc_metrics(self)
        self.metrics_at_estimated_length.calc_metrics(self)


@dataclass
class NCCStats(CCStats[int]):
    """Naive Cross-Correlation statistics with integer read counts."""
    @property
    def genomelen_repr(self) -> int:
        """Representative genome length value (returns scalar value)."""
        return self.genomelen

    @property
    def forward_reads_repr(self) -> int:
        """Representative forward read count (returns scalar value)."""
        return self.forward_reads

    @property
    def reverse_reads_repr(self) -> int:
        """Representative reverse read count (returns scalar value)."""
        return self.reverse_reads


@dataclass
class MSCCStats(CCStats[npt.NDArray[np.int64]]):
    """Mappability-Sensitive Cross-Correlation statistics with array data."""
    @property
    def genomelen_repr(self) -> int:
        """Representative genome length value (returns value at read length position)."""
        return int(self.genomelen[self.read_len - 1])

    @property
    def forward_reads_repr(self) -> int:
        """Representative forward read count (returns value at read length position)."""
        return int(self.forward_reads[self.read_len - 1])

    @property
    def reverse_reads_repr(self) -> int:
        """Representative reverse read count (returns value at read length position)."""
        return int(self.reverse_reads[self.read_len - 1])


@dataclass
class CCContainer:
    """Container for cross-correlation processing and analysis.

    Manages correlation data and provides methods for smoothing,
    background estimation, and fragment length estimation.
    """
    cc: npt.NDArray[np.float64]
    output_warnings: bool
    window_size: int
    min_calc_width: int
    read_len: int
    filter_mask_len: int

    #
    avr_cc: npt.NDArray[np.float64] = field(init=False)
    cc_min: float = field(init=False)
    est_lib_len: int = field(init=False)

    def __post_init__(self) -> None:
        self.calc_avr_cc()
        self.calc_cc_min()
        self.estimate_fragment_length()

    def calc_avr_cc(self) -> None:
        """Apply moving average filter to smooth cross-correlation values."""
        self.avr_cc = moving_avr_filter(self.cc, self.window_size)

    def calc_cc_min(self) -> None:
        """Calculate minimum cross-correlation value for background estimation."""
        cc_min = np.sort(self.cc[-self.min_calc_width:])[
            min(self.min_calc_width, self.cc.size) // 2
        ]

        if (np.median(self.cc[:NEAR_ZERO_MIN_CALC_LEN]) < cc_min and
                self.output_warnings):
            logger.warning(
                "Detected minimum coefficient seems to be larger than beginning part minimum. "
                "Consider increasing shift size (-d/--max-shift)."
            )
        self.cc_min = cc_min

    def estimate_fragment_length(self) -> None:
        """Estimate fragment length from cross-correlation peak."""
        self.est_lib_len = int(np.argmax(self.avr_cc)) + 1

        need_warning = False

        # Handle phantom peak masking
        if self.filter_mask_len and abs(self.est_lib_len - self.read_len) <= self.filter_mask_len:
            logger.warning("Estimated library length is close to the read length.")
            logger.warning(
                f"Trying to masking around the read length +/- {self.filter_mask_len}bp..."
            )

            _avr_cc = self.avr_cc.copy()
            mask_from = max(0, self.read_len - 1 - self.filter_mask_len)
            mask_to = min(len(_avr_cc), self.read_len + self.filter_mask_len)
            for i in range(mask_from, mask_to):
                _avr_cc[i] = -float("inf")
            self.est_lib_len = int(np.argmax(_avr_cc)) + 1

            if self.est_lib_len - 1 in (mask_from - 1, mask_to):
                need_warning = True

        elif (self.output_warnings and
              abs(self.est_lib_len - self.read_len) <= NEAR_READLEN_ERR_CRITERION):
            need_warning = True

        if self.output_warnings and need_warning:
            logger.error(
                "Estimated library length is close to the read length! "
                "Please check output plots."
            )

    def calc_FWHM(self, library_len: int) -> Union[int, bool]:
        """Calculate Full Width at Half Maximum (FWHM) of fragment length peak."""
        assert hasattr(self, 'avr_cc'), "Call calc_avr_cc() before calc_FWHM()."
        assert hasattr(self, 'cc_min'), "Call calc_cc_min() before calc_FWHM()."

        if np.isnan(self.cc_min):
            #
            return False

        max_i = library_len - 1
        assert max_i >= 0, max_i  # Original implementation check
        cc_max = self.avr_cc[max_i - 1]  # Original implementation uses max_i - 1
        assert cc_max > self.cc_min, (cc_max, self.cc_min)  # Original implementation check

        target = self.cc_min + (cc_max - self.cc_min) / 2

        # Forward direction
        forward_shift = 0
        forward_failed = False
        while self.avr_cc[max_i + forward_shift] > target:
            forward_shift += 1
            if max_i + forward_shift == self.avr_cc.size:
                logger.warning(
                    "Failed to calc the half width at half maximum in the forward "
                    "side of the peak. Consider increasing shift size (-d/--max-shift)."
                )
                forward_failed = True
                forward_shift -= 1
                break

        # Backward direction
        backward_shift = 0
        backward_failed = False
        while self.avr_cc[max_i - backward_shift] > target:
            backward_shift += 1
            if max_i < backward_shift:
                logger.warning(
                    "Failed to calc the half width at half maximum in the backward side of the peak."
                )
                backward_failed = True
                backward_shift -= 1
                break

        if forward_failed and backward_failed:
            logger.error("Failed to calcurate the full width at half maximum.")
            return False
        elif forward_failed:
            logger.warning("Use twice width of the half width at half maximum in the backward side")
            return backward_shift * 2 + 1
        elif backward_failed:
            logger.warning("Use twice width of the half width at half maximum in the forward side")
            return forward_shift * 2 + 1
        else:
            return backward_shift + forward_shift + 1


TLen = TypeVar("TLen", int, npt.NDArray[np.int64])


class CorrLike(Protocol, Generic[TLen, TCount]):
    """Protocol for objects that behave like CorrelationResult."""
    cc: npt.NDArray[np.float64]
    genomelen: TLen
    forward_sum: TCount
    reverse_sum: TCount


TStats = TypeVar('TStats', NCCStats, MSCCStats)


@dataclass
class CorrParams(CorrLike):
    """Parameter container for correlation statistics.

    Used to pass correlation data and read statistics to
    statistics calculation functions.
    """
    cc: npt.NDArray[np.float64]
    genomelen: Union[int, npt.NDArray[np.int64]]
    forward_sum: Union[int, npt.NDArray[np.int64]]
    reverse_sum: Union[int, npt.NDArray[np.int64]]


@overload
def _prepare_chromosome_stat(
    result: NCCResultModel,
    config: StatConfig,
    stats_type: None = None,
    output_warnings: bool = ...,
    estimated_library_len: Optional[int] = ...
) -> Tuple[NCCStats, CCContainer]: ...


@overload
def _prepare_chromosome_stat(
    result: MSCCResultModel,
    config: StatConfig,
    stats_type: None = None,
    output_warnings: bool = ...,
    estimated_library_len: Optional[int] = ...
) -> Tuple[MSCCStats, CCContainer]: ...


@overload
def _prepare_chromosome_stat(
    result: CorrParams,
    config: StatConfig,
    stats_type: Type[TStats],
    output_warnings: bool = ...,
    estimated_library_len: Optional[int] = ...
) -> Tuple[TStats, CCContainer]: ...


def _prepare_chromosome_stat(
    result: CorrLike,
    config: StatConfig,
    stats_type: Optional[Type[TStats]] = None,
    output_warnings: bool = True,
    estimated_library_len: Optional[int] = None
) -> Tuple[TStats, CCContainer]:
    """Prepare chromosome statistics from correlation results.

    Args:
        result: Correlation result data
        config: Statistical configuration
        stats_type: Target statistics type (inferred if None)
        output_warnings: Enable warning output
        estimated_library_len: Override estimated fragment length

    Returns:
        Tuple of statistics object and correlation container
    """
    #
    cc_container = CCContainer(
        cc=result.cc,
        output_warnings=output_warnings,
        window_size=config.mv_avr_filter_len,
        min_calc_width=config.min_calc_width,
        read_len=config.read_length,
        filter_mask_len=config.filter_mask_len
    )

    #
    if config.expected_library_length is not None:
        metrics_at_expected_length = CCQualityMetrics(
            fragment_length=config.expected_library_length,
            ccfl=cc_container.cc[config.expected_library_length - 1],
            fwhm=cc_container.calc_FWHM(config.expected_library_length)
        )
    else:
        metrics_at_expected_length = CCQualityMetrics()

    #
    if estimated_library_len is None:
        estimated_library_len = cc_container.est_lib_len

    metrics_at_estimated_length = CCQualityMetrics(
        fragment_length=estimated_library_len,
        ccfl=cc_container.cc[estimated_library_len - 1],
        fwhm=cc_container.calc_FWHM(estimated_library_len)
    )

    stats: Union[NCCStats, MSCCStats, TStats]
    if isinstance(result, NCCResultModel):
        stats = NCCStats(
            read_len=config.read_length,
            genomelen=result.genomelen,
            forward_reads=result.forward_sum,
            reverse_reads=result.reverse_sum,
            cc_min=cc_container.cc_min,
            ccrl=result.cc[config.read_length - 1],
            metrics_at_expected_length=metrics_at_expected_length,
            metrics_at_estimated_length=metrics_at_estimated_length
        )
    elif isinstance(result, MSCCResultModel):
        stats = MSCCStats(
            read_len=config.read_length,
            genomelen=np.array(result.mappable_len, dtype=np.int64),
            forward_reads=result.forward_sum,
            reverse_reads=result.reverse_sum,
            cc_min=cc_container.cc_min,
            ccrl=result.cc[config.read_length - 1],
            metrics_at_expected_length=metrics_at_expected_length,
            metrics_at_estimated_length=metrics_at_estimated_length
        )
    elif stats_type is not None:
        stats = stats_type(
            read_len=config.read_length,
            genomelen=result.genomelen,
            forward_reads=result.forward_sum,
            reverse_reads=result.reverse_sum,
            cc_min=cc_container.cc_min,
            ccrl=result.cc[config.read_length - 1],
            metrics_at_expected_length=metrics_at_expected_length,
            metrics_at_estimated_length=metrics_at_estimated_length
        )
    else:
        raise TypeError("Unsupported CorrelationResult type.")

    #
    return cast(TStats, stats), cc_container


def make_chromosome_stat(
    result: Union[NCCResultModel, MSCCResultModel],
    config: StatConfig,
    output_warnings: bool = False,
    estimated_library_len: Optional[int] = None
) -> Union[ChromosomeStats, EmptyChromosomeStats]:
    """Create chromosome statistics, returning EmptyChromosomeStats for empty results."""

    # Use class-based detection for consistency with aggregate_chromosome_stats
    if isinstance(result, (EmptyNCCResult, EmptyMSCCResult, EmptyBothChromResult)):
        return EmptyChromosomeStats(genomelen_repr=result.genomelen)

    # For regular results, proceed with normal statistical processing
    stats, cc_container = _prepare_chromosome_stat(
        result,
        config,
        None,
        output_warnings,
        estimated_library_len
    )

    return ChromosomeStats(
        stats=stats,
        cc=cc_container.cc,
        avr_cc=cc_container.avr_cc,
        est_lib_len=cc_container.est_lib_len,
        mv_avr_filter_len=config.mv_avr_filter_len
    )


def aggregate_chromosome_stats(
    chrom_stats: Optional[Dict[str, Union[ChromosomeStats[TStats], EmptyChromosomeStats]]],
    config: StatConfig,
    output_warnings: bool,
    estimated_library_len: Optional[int] = None
) -> Optional[WholeGenomeStats[TStats]]:
    """Aggregate chromosome statistics into genome-wide statistics.

    Args:
        chrom_stats: Dictionary of per-chromosome statistics
        config: Statistical configuration
        output_warnings: Enable warning output
        estimated_library_len: Fragment length estimate for processing

    Returns:
        Aggregated genome-wide statistics or None if no data
    """
    if chrom_stats is None:
        return None

    # Separate empty results from regular results using clear class-based identification
    regular_stats: Dict[str, ChromosomeStats[TStats]] = {}
    empty_genomelen_contributions: List[int] = []  # Store only for tracking, not for calculation

    for chrom, stats_obj in chrom_stats.items():
        if isinstance(stats_obj, EmptyChromosomeStats):
            # EmptyChromosomeStats: exclude from statistical processing
            empty_genomelen_contributions.append(stats_obj.genomelen_repr)
        else:
            # ChromosomeStats: participate in full statistical aggregation
            assert stats_obj.stats is not None, f"Stats for chromosome {chrom} is None."
            regular_stats[chrom] = stats_obj

    # If no regular stats, return None (should not happen in normal cases)
    if not regular_stats:
        return None

    #
    first_stats = next(iter(regular_stats.values())).stats
    assert first_stats is not None, "No stats available for aggregation."
    stats_type = type(first_stats)

    # Prepare data for aggregation (only regular results)
    genome_lengths: List[Union[int, npt.NDArray[np.int64]]] = []
    forward_reads: List[Union[int, npt.NDArray[np.int64]]] = []
    reverse_reads: List[Union[int, npt.NDArray[np.int64]]] = []
    cc_arrays: List[npt.NDArray[np.float64]] = []
    representative_genome_lengths: List[int] = []  # For cc merging

    for chrom, stats_obj in regular_stats.items():
        assert stats_obj.stats is not None, f"Stats for chromosome {chrom} is None."

        genome_lengths.append(stats_obj.stats.genomelen)
        forward_reads.append(stats_obj.stats.forward_reads)
        reverse_reads.append(stats_obj.stats.reverse_reads)
        representative_genome_lengths.append(stats_obj.stats.genomelen_repr)

        assert stats_obj.cc is not None, f"cross-correlation for chromosome {chrom} is None."
        cc_arrays.append(stats_obj.cc)

    # Aggregate basic statistics from regular results
    regular_genomelen: Union[int, npt.NDArray[np.int64]] = np.sum(np.asarray(genome_lengths, dtype=np.int64), axis=0)
    total_forward_reads: Union[int, npt.NDArray[np.int64]] = np.sum(np.asarray(forward_reads, dtype=np.int64), axis=0)
    total_reverse_reads: Union[int, npt.NDArray[np.int64]] = np.sum(np.asarray(reverse_reads, dtype=np.int64), axis=0)

    # Add EmptyChromosomeStats genome lengths to maintain consistency with single-process mode
    # EmptyChromosomeStats represent chromosomes with no reads but valid genome lengths
    total_genomelen: Union[int, npt.NDArray[np.int64]]
    if isinstance(regular_genomelen, np.ndarray):
        # For MSCC: EmptyChromosomeStats don't contribute to mappable length calculation
        # because they represent chromosomes with no actual mappable data
        # Only regular results with actual mappability data should be included
        total_genomelen = regular_genomelen
    else:
        # For NCC: add empty contributions to total genome length
        empty_genomelen_total = sum(empty_genomelen_contributions)
        total_genomelen = regular_genomelen + empty_genomelen_total

    # Aggregate raw cross-correlation using Fisher z-transformation (only regular results)
    aggregated_cc, interval_lower, interval_upper = merge_correlations(
        np.array(representative_genome_lengths, dtype=np.int64),
        cc_arrays,
        first_stats.read_len
    )
    # Convert to numpy array for consistency
    aggregated_cc = np.array(aggregated_cc, dtype=np.float64)

    #
    return make_whole_genome_stat(
        CorrParams(
            cc=aggregated_cc,
            genomelen=total_genomelen,
            forward_sum=total_forward_reads,
            reverse_sum=total_reverse_reads
        ),
        config,
        interval_upper=interval_upper,
        interval_lower=interval_lower,
        stats_type=stats_type,
        output_warnings=output_warnings,
        estimated_library_len=estimated_library_len
    )


def make_whole_genome_stat(
    result: CorrParams,
    config: StatConfig,
    interval_upper: npt.NDArray[np.float64],
    interval_lower: npt.NDArray[np.float64],
    stats_type: Type[TStats],
    output_warnings: bool = True,
    estimated_library_len: Optional[int] = None
) -> WholeGenomeStats[TStats]:
    """Create whole-genome statistics with confidence intervals.

    Args:
        result: Correlation parameters
        config: Statistical configuration
        interval_upper: Upper confidence bounds
        interval_lower: Lower confidence bounds
        stats_type: Target statistics type
        output_warnings: Enable warning output
        estimated_library_len: Fragment length estimate

    Returns:
        Complete whole-genome statistics with intervals
    """
    stat, cc_container = _prepare_chromosome_stat(
        result,
        config,
        stats_type,
        output_warnings,
        estimated_library_len
    )

    return WholeGenomeStats(
        stats=stat,
        cc=cc_container.cc,
        avr_cc=cc_container.avr_cc,
        est_lib_len=cc_container.est_lib_len,
        cc_upper=interval_upper,
        cc_lower=interval_lower,
        mv_avr_filter_len=config.mv_avr_filter_len
    )


def make_genome_wide_stat(
    result: GenomeWideResultModel,
    config: StatConfig,
    output_warnings: bool
) -> GenomeWideStats:
    """Create a GenomeWideStats object from a GenomeWideResult."""
    #
    ncc_stats = mscc_stats = None

    #
    if isinstance(result, MSCCGenomeWideResultModel):
        mscc_stats = {
            chrom: make_chromosome_stat(chromres, config)
            for chrom, chromres in result.chroms.items()
        }
    elif isinstance(result, BothGenomeWideResultModel):
        mscc_stats = {
            chrom: make_chromosome_stat(chromres, config)
            for chrom, chromres in result.mappable_chroms.items()
        }

    #
    if isinstance(result, NCCGenomeWideResultModel):
        ncc_stats = {
            chrom: make_chromosome_stat(chromres, config)
            for chrom, chromres in result.chroms.items()
        }
    elif isinstance(result, BothGenomeWideResultModel):
        assert mscc_stats is not None, "MSCC stats must be available for BothGenomeWideResultModel."
        ncc_stats = {}
        for chrom, chromres in result.chroms.items():
            if chrom in mscc_stats:
                chromstat = mscc_stats[chrom]
                assert not isinstance(chromstat, EmptyChromosomeStats)
                estimated_library_len = chromstat.stats.metrics_at_estimated_length.fragment_length
            else:
                estimated_library_len = None

            ncc_stats[chrom] = make_chromosome_stat(
                chromres,
                config,
                estimated_library_len=estimated_library_len
            )

    if ncc_stats is None and mscc_stats is None:
        raise TypeError("Unsupported GenomeWideResult type.")

    #
    whole_mscc_stats: Optional[WholeGenomeStats[MSCCStats]] = aggregate_chromosome_stats(
        mscc_stats,
        config,
        output_warnings
    )

    if whole_mscc_stats is None:
        estimated_library_len = None
    else:
        estimated_library_len = whole_mscc_stats.est_lib_len

    whole_ncc_stats = aggregate_chromosome_stats(
        ncc_stats,
        config,
        output_warnings,
        estimated_library_len=estimated_library_len
    )

    #
    if whole_ncc_stats is not None:
        if whole_ncc_stats.stats.forward_reads == 0:
            logger.error("There is no forward read.")
            raise ReadsTooFew
        if whole_ncc_stats.stats.reverse_reads == 0:
            logger.error("There is no reverse read.")
            raise ReadsTooFew

        # Check strand balance for NCC reads
        whole_ncc_stats.stats.check_strand_balance(config.chi2_pval, "NCC")

    if whole_mscc_stats is not None:
        errormsg = "There is no forward read in mappable regions."
        if whole_mscc_stats.stats.forward_reads.sum() == 0:
            if whole_ncc_stats is not None:
                logger.warning(errormsg)
            else:
                logger.error(errormsg)
                raise ReadsTooFew

        errormsg = "There is no reverse read in mappable regions."
        if whole_mscc_stats.stats.reverse_reads.sum() == 0:
            if whole_ncc_stats is not None:
                logger.warning(errormsg)
            else:
                logger.error(errormsg)
                raise ReadsTooFew

        # Check strand balance for MSCC mappable reads
        whole_mscc_stats.stats.check_strand_balance(config.chi2_pval, "MSCC")

    #
    cleaned_ncc_stats = cleaned_mscc_stats = None
    if ncc_stats is not None:
        cleaned_ncc_stats = {chrom: stats for chrom, stats in ncc_stats.items()
                             if not isinstance(stats, EmptyChromosomeStats)}
    if mscc_stats is not None:
        cleaned_mscc_stats = {chrom: stats for chrom, stats in mscc_stats.items()
                              if not isinstance(stats, EmptyChromosomeStats)}

    return GenomeWideStats(
        whole_ncc_stats=whole_ncc_stats,
        whole_mscc_stats=whole_mscc_stats,
        ncc_stats=cleaned_ncc_stats,
        mscc_stats=cleaned_mscc_stats
    )
