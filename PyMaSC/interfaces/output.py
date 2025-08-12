"""Output data models and statistics structures.

Defines dataclasses for formatted output and statistics reporting.
"""
from dataclasses import dataclass, asdict, make_dataclass, field
from typing import Generator, Tuple, Dict, Any, Union, Literal

from .result import AbstractDataclass
from .result import (
    CorrelationResultModel, NCCResultModel, MSCCResultModel, BothChromResultModel,
    NCCGenomeWideResultModel, MSCCGenomeWideResultModel, BothGenomeWideResultModel
)


@dataclass
class StatItems(AbstractDataclass):
    """Base class for statistical data items."""
    pass


@dataclass
class StatLabels(AbstractDataclass):
    """Base class for statistical label definitions."""
    pass


@dataclass
class StatWithLabels(AbstractDataclass):
    """Container combining statistical data with descriptive labels.

    Pairs statistical values with human-readable labels for output formatting.
    """
    stats: StatItems
    labels: StatLabels

    def get_with_labels(self) -> Generator[Tuple[str, Any], None, None]:
        """Generate (label, value) pairs for formatted output.

        Yields:
            Tuple of (label_string, value) for each statistic
        """
        for attr, value in asdict(self.stats).items():
            yield getattr(self.labels, attr), value


nanint = Union[int, Literal["nan"]]
nanfloat = Union[float, Literal["nan"]]


@dataclass
class SummaryItems(StatItems):
    """Summary statistics for PyMaSC analysis results.

    Contains basic analysis parameters and library length information.
    """
    name: str
    read_len: int
    expected_lib_len: nanint
    est_lib_len: nanint


@dataclass
class CorrelationItems(StatItems):
    """Cross-correlation analysis statistics and quality metrics.

    Contains detailed statistics from NCC/MSCC analysis including
    read counts, correlation values, and quality indicators.
    """
    genomelen: nanint
    forward_reads: nanint
    reverse_reads: nanint
    min_cc: nanfloat
    ccrl: nanfloat
    ccfl: nanfloat
    ccfl_est: nanfloat
    nsc: nanfloat
    rsc: nanfloat
    nsc_est: nanfloat
    rsc_est: nanfloat
    fwhm: nanfloat
    vsn: nanfloat
    fwhm_est: nanfloat
    vsn_est: nanfloat


def _make_labeled_dataclass(name: str, dc: type) -> type:
    """Create a dataclass with string labels for statistical data.

    Args:
        name: Name of the new dataclass
        dc: Source dataclass to create labels for

    Returns:
        New dataclass type with string labels matching source fields
    """
    return make_dataclass(
        name,
        [(attr, str) for attr in dc.__annotations__.keys()],
        bases=(StatLabels, )
    )


SUMMARY_LABELS = _make_labeled_dataclass('SUMMARY_LABELS', SummaryItems)(
    "Name",
    "Read length",
    "Expected library length",
    "Estimated library length"
)


NCC_LABELS = _make_labeled_dataclass('NCC_LABELS', CorrelationItems)(
    "Genome length",
    "Forward reads",
    "Reverse reads",
    "Minimum NCC",
    "NCC at read length",
    "NCC at expected library length",
    "NCC at estimated library length",
    "NSC",
    "RSC",
    "Estimated NSC",
    "Estimated RSC",
    "FWHM",
    "VSN",
    "Estimated FWHM",
    "Estimated VSN"
)


MSCC_LABELS = _make_labeled_dataclass('MSCC_LABELS', CorrelationItems)(
    "DMP length",
    "Forward reads in DMP",
    "Reverse reads in DMP",
    "Minimum MSCC",
    "MSCC at read length",
    "MSCC at expected library length",
    "MSCC at estimated library length",
    "MSCC NSC",
    "MSCC RSC",
    "Estimated MSCC NSC",
    "Estimated MSCC RSC",
    "MSCC FWHM",
    "MSCC VSN",
    "Estimated MSCC FWHM",
    "Estimated MSCC VSN"
)


@dataclass
class SummaryStats(StatWithLabels):
    """Summary statistics with predefined labels.

    Combines SummaryItems with standard summary labels for output formatting.
    """
    stats: SummaryItems
    labels: StatLabels = field(init=False)

    def __post_init__(self) -> None:
        """Initialize labels with predefined summary labels."""
        self.labels = SUMMARY_LABELS


@dataclass
class CorrelationStats(StatWithLabels):
    """Cross-correlation statistics with configurable labels.

    Combines CorrelationItems with specified labels for NCC or MSCC output.
    """
    stats: CorrelationItems
    labels: StatLabels


@dataclass
class OutputStats(SummaryStats):
    """Complete output statistics combining summary and correlation data.

    Contains all statistics needed for PyMaSC output files including
    both NCC and MSCC correlation statistics.
    """
    ncc_stats: CorrelationStats
    mscc_stats: CorrelationStats


@dataclass
class CorrelationResult(CorrelationResultModel):
    """Base class for cross-correlation results with automatic max_shift calculation.

    Computes max_shift from the correlation array length during initialization.
    """
    max_shift: int = field(init=False)

    def __post_init__(self) -> None:
        """Calculate max_shift from correlation array length.

        Raises:
            ValueError: If correlation array is empty
        """
        self.max_shift = len(self.cc) - 1
        if self.max_shift < 0:
            raise ValueError("Cross-correlation array must have at least one element.")


@dataclass
class NCCResult(CorrelationResult, NCCResultModel):
    """Naive Cross-Correlation (NCC) analysis result."""
    pass


@dataclass
class MSCCResult(CorrelationResult, MSCCResultModel):
    """Mappability-Sensitive Cross-Correlation (MSCC) analysis result."""
    pass


@dataclass
class BothChromResult(BothChromResultModel):
    """Combined NCC and MSCC results for a single chromosome."""
    chrom: NCCResult
    mappable_chrom: MSCCResult


@dataclass
class NCCGenomeWideResult(NCCGenomeWideResultModel):
    """Genome-wide NCC analysis results for all chromosomes."""
    chroms: Dict[str, NCCResult]


@dataclass
class MSCCGenomeWideResult(MSCCGenomeWideResultModel):
    """Genome-wide MSCC analysis results for all chromosomes."""
    chroms: Dict[str, MSCCResult]


@dataclass
class BothGenomeWideResult(BothGenomeWideResultModel):
    """Combined genome-wide NCC and MSCC analysis results."""
    chroms: Dict[str, NCCResult]
    mappable_chroms: Dict[str, MSCCResult]
