import logging
from abc import abstractmethod
from dataclasses import dataclass, field
from typing import Union, Optional, Type, Dict, List, Tuple, TypeVar, Protocol, Generic, runtime_checkable
from functools import cached_property

import numpy as np
import numpy.typing as npt
from scipy.stats import norm

from .interfaces.result import (
    CorrelationResult, NCCResult, MSCCResult,
    GenomeWideResult,
    NCCGenomeWideResult, MSCCGenomeWideResult, BothGenomeWideResult
)
from .interfaces.stats import (
    CCQualityMetrics as CCQualityMetricsModel,
    CCStats as CCStatsModel,
    ChromosomeStats,
    GenomeWideStats
)

from PyMaSC.utils.calc import moving_avr_filter
from PyMaSC.utils.stats_utils import npcalc_with_logging_warn, CrossCorrelationMerger, ArrayAggregator

logger = logging.getLogger(__name__)

# Constants
NEAR_READLEN_ERR_CRITERION = 5
NEAR_ZERO_MIN_CALC_LEN = 10


@dataclass
class CCQualityMetrics(CCQualityMetricsModel):
    fragment_length: Optional[int] = None
    ccfl: Optional[float] = None
    fwhm: Optional[int] = None

    #
    nsc: Optional[float] = None
    rsc: Optional[float] = None
    vsn: Optional[float] = None

    def calc_metrics(self, stats: 'CCStats') -> None:
        if self.fragment_length is None:
            return
        else:
            assert isinstance(self.ccfl, float), "ccfl must be a float."
        self.nsc = self.ccfl / stats.cc_min
        self.rsc = (self.ccfl - stats.cc_min) / (stats.ccrl - stats.cc_min)
        if self.fwhm is not None:
            self.vsn = 2 * self.ccfl * self.fwhm / (stats.forward_reads + stats.reverse_reads)


@dataclass
class CCStats(CCStatsModel):
    read_len: int
    raw_genomelen: Union[int, npt.NDArray[np.int64]]
    raw_forward_reads: Union[int, npt.NDArray[np.int64]]
    raw_reverse_reads: Union[int, npt.NDArray[np.int64]]

    def calc_qc_metrics(self) -> None:
        self.metrics_at_expected_length.calc_metrics(self)
        self.metrics_at_estimated_length.calc_metrics(self)


@dataclass
class NCCStats(CCStats):
    @property
    def genomelen(self) -> int:
        """Return the representative value of genomelen.
        i.e., for MSCCStats, genomelen is an array and return the value at read_len index.
        """
        if isinstance(self.raw_genomelen, np.integer):
            self.raw_genomelen = int(self.raw_genomelen)
        assert isinstance(self.raw_genomelen, int), "genomelen must be an integer."
        return self.raw_genomelen

    @property
    def forward_reads(self) -> int:
        if isinstance(self.raw_forward_reads, np.integer):
            self.raw_forward_reads = int(self.raw_forward_reads)
        assert isinstance(self.raw_forward_reads, int), "forward_reads must be an integer."
        return self.raw_forward_reads

    @property
    def reverse_reads(self) -> int:
        if isinstance(self.raw_reverse_reads, np.integer):
            self.raw_reverse_reads = int(self.raw_reverse_reads)
        assert isinstance(self.raw_reverse_reads, int), "reverse_reads must be an integer."
        return self.raw_reverse_reads


@dataclass
class MSCCStats(CCStats):
    @property
    def genomelen(self) -> int:
        assert isinstance(self.raw_genomelen, np.ndarray), "genomelen must be an array."
        return int(self.raw_genomelen[self.read_len - 1])

    @property
    def forward_reads(self) -> int:
        assert isinstance(self.raw_forward_reads, np.ndarray), "forward_reads must be an array."
        return int(self.raw_forward_reads[self.read_len - 1])

    @property
    def reverse_reads(self) -> int:
        assert isinstance(self.raw_reverse_reads, np.ndarray), "reverse_reads must be an array."
        return int(self.raw_reverse_reads[self.read_len - 1])


@dataclass
class CCContainer:
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
TCount = TypeVar("TCount", int, npt.NDArray[np.int64])


class CorrLike(Protocol, Generic[TLen, TCount]):
    """Protocol for objects that behave like CorrelationResult."""
    cc: npt.NDArray[np.float64]
    genomelen: TLen
    forward_sum: TCount
    reverse_sum: TCount


def make_chromosome_stat(
    result: CorrLike,
    read_len: int,
    mv_avr_filter_len: int = 15,
    expected_library_len: Optional[int] = None,
    filter_mask_len: int = 5,
    min_calc_width: int = 10,
    output_warnings: bool = True,
    stats_type: Optional[Type[CCStats]] = None
) -> ChromosomeStats:
    #
    cc_container = CCContainer(
        cc=result.cc,
        output_warnings=output_warnings,
        window_size=mv_avr_filter_len,
        min_calc_width=min_calc_width,
        read_len=read_len,
        filter_mask_len=filter_mask_len
    )

    #
    if expected_library_len is not None:
        metrics_at_expected_length = CCQualityMetrics(
            fragment_length=expected_library_len,
            ccfl=cc_container.cc[expected_library_len - 1],
            fwhm=cc_container.calc_FWHM(expected_library_len)
        )
    else:
        metrics_at_expected_length = CCQualityMetrics()

    #
    metrics_at_estimated_length = CCQualityMetrics(
        fragment_length=cc_container.est_lib_len,
        ccfl=cc_container.cc[cc_container.est_lib_len - 1],
        fwhm=cc_container.calc_FWHM(cc_container.est_lib_len)
    )

    statclass: Type[CCStats]
    genomelen: Union[int, npt.NDArray[np.int64]]
    if isinstance(result, NCCResult):
        statclass = NCCStats
        genomelen = result.genomelen
    elif isinstance(result, MSCCResult):
        statclass = MSCCStats
        assert result.mappable_len is not None, "mappable_len must be set for MSCCResult."
        genomelen = np.array(result.mappable_len, dtype=np.int64)
    elif stats_type is not None:
        statclass = stats_type
        genomelen = result.genomelen
    else:
        raise TypeError("Unsupported CorrelationResult type.")

    stats = statclass(
        read_len=read_len,
        raw_genomelen=genomelen,
        raw_forward_reads=result.forward_sum,
        raw_reverse_reads=result.reverse_sum,
        cc_min=cc_container.cc_min,
        ccrl=result.cc[read_len - 1],
        metrics_at_expected_length=metrics_at_expected_length,
        metrics_at_estimated_length=metrics_at_estimated_length
    )

    #
    return ChromosomeStats(
        stats=stats,
        cc=cc_container.cc,
        avr_cc=cc_container.avr_cc,
        est_lib_len=cc_container.est_lib_len
    )


@dataclass
class CorrParams(CorrLike):
    cc: npt.NDArray[np.float64]
    genomelen: Union[int, npt.NDArray[np.int64]]
    forward_sum: Union[int, npt.NDArray[np.int64]]
    reverse_sum: Union[int, npt.NDArray[np.int64]]


TStats = TypeVar('TStats', NCCStats, MSCCStats)


def aggregate_chromosome_stats(
    chrom_stats: Optional[Dict[str, ChromosomeStats[TStats]]],
    read_len: int,
    mv_avr_filter_len: int,
    filter_mask_len: int,
    min_calc_width: int,
    output_warnings: bool,
    expected_library_len: Optional[int] = None
) -> Optional[ChromosomeStats[TStats]]:
    if chrom_stats is None:
        return None

    #
    first_stats = next(iter(chrom_stats.values())).stats
    assert first_stats is not None, "No stats available for aggregation."
    stats_type = type(first_stats)

    # Prepare data for aggregation
    genome_lengths: List[Union[int, npt.NDArray[np.int64]]] = []
    forward_reads: List[Union[int, npt.NDArray[np.int64]]] = []
    reverse_reads: List[Union[int, npt.NDArray[np.int64]]] = []
    cc_arrays: List[npt.NDArray[np.float64]] = []
    representative_genome_lengths: List[int] = []  # For cc merging

    for chrom, stats_obj in chrom_stats.items():
        assert stats_obj.stats is not None, f"Stats for chromosome {chrom} is None."

        # Collect data
        genome_lengths.append(stats_obj.stats.raw_genomelen)
        forward_reads.append(stats_obj.stats.raw_forward_reads)
        reverse_reads.append(stats_obj.stats.raw_reverse_reads)
        representative_genome_lengths.append(stats_obj.stats.genomelen)

        assert stats_obj.cc is not None, f"Raw cross-correlation for chromosome {chrom} is None."
        cc_arrays.append(stats_obj.cc)

    # Aggregate basic statistics
    total_genomelen: Union[int, npt.NDArray[np.int64]] = np.sum(np.asarray(genome_lengths, dtype=np.int64), axis=0)
    total_forward_reads: Union[int, npt.NDArray[np.int64]] = np.sum(np.asarray(forward_reads, dtype=np.int64), axis=0)
    total_reverse_reads: Union[int, npt.NDArray[np.int64]] = np.sum(np.asarray(reverse_reads, dtype=np.int64), axis=0)

    # Aggregate raw cross-correlation using Fisher z-transformation
    from PyMaSC.core.constants import MERGED_CC_CONFIDENCE_INTERVAL
    merger = CrossCorrelationMerger(confidence_interval=MERGED_CC_CONFIDENCE_INTERVAL)
    _, aggregated_cc, _ = merger.merge_correlations(
        np.array(representative_genome_lengths, dtype=np.int64),
        cc_arrays,
        first_stats.read_len
    )
    # Convert to numpy array for consistency
    aggregated_cc = np.array(aggregated_cc, dtype=np.float64)

    #
    return make_chromosome_stat(
        CorrParams(
            cc=aggregated_cc,
            genomelen=total_genomelen,
            forward_sum=total_forward_reads,
            reverse_sum=total_reverse_reads
        ),
        read_len=read_len,
        mv_avr_filter_len=mv_avr_filter_len,
        expected_library_len=expected_library_len,
        filter_mask_len=filter_mask_len,
        min_calc_width=min_calc_width,
        output_warnings=output_warnings,
        stats_type=stats_type
    )


def make_genome_wide_stat(
    result: GenomeWideResult,
    read_len: int,
    mv_avr_filter_len: int,
    expected_library_len: Optional[int],
    filter_mask_len: int,
    min_calc_width: int,
    output_warnings: bool
) -> GenomeWideStats:
    """Create a GenomeWideStats object from a GenomeWideResult."""
    #
    ncc_stats = mscc_stats = None

    if isinstance(result, (NCCGenomeWideResult, BothGenomeWideResult)):
        ncc_stats = {
            chrom: make_chromosome_stat(
                chromres,
                read_len,
                mv_avr_filter_len,
                expected_library_len,
                filter_mask_len,
                min_calc_width
            )
            for chrom, chromres in result.chroms.items()
        }
    if isinstance(result, MSCCGenomeWideResult):
        mscc_stats = {
            chrom: make_chromosome_stat(
                chromres,
                read_len,
                mv_avr_filter_len,
                expected_library_len,
                filter_mask_len,
                min_calc_width
            )
            for chrom, chromres in result.chroms.items()
        }
    elif isinstance(result, BothGenomeWideResult):
        mscc_stats = {
            chrom: make_chromosome_stat(
                chromres,
                read_len,
                mv_avr_filter_len,
                expected_library_len,
                filter_mask_len,
                min_calc_width
            )
            for chrom, chromres in result.mappable_chroms.items()
        }

    if ncc_stats is None and mscc_stats is None:
        raise TypeError("Unsupported GenomeWideResult type.")

    return GenomeWideStats(
        whole_ncc_stats=aggregate_chromosome_stats(
            ncc_stats,
            read_len,
            mv_avr_filter_len,
            filter_mask_len,
            min_calc_width,
            output_warnings
        ),
        whole_mscc_stats=aggregate_chromosome_stats(
            mscc_stats,
            read_len,
            mv_avr_filter_len,
            filter_mask_len,
            min_calc_width,
            output_warnings
        ),
        ncc_stats=ncc_stats,
        mscc_stats=mscc_stats
    )


# TODO: Remove this duplicate implementation after verifying CrossCorrelationMerger works correctly
# def merge_correlations(
#     genome_lengths: npt.NDArray[np.int64],
#     correlation_arrays: List[npt.NDArray[np.float64]],
#     read_length: int,
#     confidence_interval: float = 0.99
# ) -> Tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], npt.NDArray[np.float64]]:
#     """Merge cross-correlation results using Fisher z-transformation.
# 
#     This method replicates the exact behavior of CCResult._merge_cc(),
#     combining per-chromosome correlations into genome-wide profiles with
#     confidence intervals.
# 
#     Mathematical Process:
#     1. Fisher z-transformation: z = arctanh(r) for each correlation
#     2. Weighted averaging: z_avg = Σ(w_i * z_i) / Σ(w_i) where w_i = n_i - 3
#     3. Confidence intervals: CI = tanh(z_avg ± z_α/2 * SE)
#     4. Inverse transformation: r = tanh(z_avg) back to correlation
# 
#     Args:
#         genome_lengths: Array of genome lengths for each chromosome
#         correlation_arrays: List of correlation arrays, one per chromosome
#         read_length: Read length for position-dependent weighting
# 
#     Returns:
#         Tuple of (merged_correlations, lower_confidence, upper_confidence)
#     """
#     ns = genome_lengths
# 
#     merged_r = []
#     interval_upper = []
#     interval_lower = []
# 
#     for i, __ccs in enumerate(zip(*correlation_arrays)):
#         _ccs: Tuple[np.float64, ...] = __ccs
#         # Filter out NaN values
#         nans = np.isnan(_ccs)
#         ccs = np.array(_ccs)[~nans]
# 
#         # Calculate weights (effective sample size - 3)
#         if len(ns.shape) == 1:
#             _ns = ns[~nans] - 3
#         else:
#             _ns = ns[~nans, abs(read_length - i)] - 3
# 
#         # Fisher z-transformation
#         zs = np.arctanh(ccs)
# 
#         # Filter out infinite values
#         infs = np.isinf(zs)
#         zs = zs[~infs]
#         _ns = _ns[~infs]
# 
#         # Weighted average in z-space
#         avr_z = np.average(zs, weights=_ns)
# 
#         # Calculate confidence interval
#         z_interval = norm.ppf(1 - (1 - confidence_interval) / 2) * np.sqrt(1 / np.sum(_ns))
#         z_interval_upper = avr_z + z_interval
#         z_interval_lower = avr_z - z_interval
# 
#         # Transform back to correlation space
#         merged_r.append(np.tanh(avr_z))
#         interval_upper.append(np.tanh(z_interval_upper))
#         interval_lower.append(np.tanh(z_interval_lower))
# 
#     return (
#         np.array(merged_r, dtype=np.float64),
#         np.array(interval_lower, dtype=np.float64),
#         np.array(interval_upper, dtype=np.float64)
#     )
