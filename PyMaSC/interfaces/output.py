from dataclasses import dataclass, asdict, make_dataclass, field
from typing import Generator, Tuple, Dict, Any, Union, Literal

from .result import AbstractDataclass
from .result import (
    CorrelationResultModel, NCCResultModel, MSCCResultModel, BothChromResultModel,
    NCCGenomeWideResultModel, MSCCGenomeWideResultModel, BothGenomeWideResultModel
)


@dataclass
class StatItems(AbstractDataclass):
    pass


@dataclass
class StatLabels(AbstractDataclass):
    pass


@dataclass
class StatWithLabels(AbstractDataclass):
    stats: StatItems
    labels: StatLabels

    def get_with_labels(self) -> Generator[Tuple[str, Any], None, None]:
        for attr, value in asdict(self.stats).items():
            yield getattr(self.labels, attr), value


nanint = Union[int, Literal["nan"]]
nanfloat = Union[float, Literal["nan"]]


@dataclass
class SummaryItems(StatItems):
    name: str
    read_len: int
    expected_lib_len: nanint
    est_lib_len: nanint


@dataclass
class CorrelationItems(StatItems):
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
    return make_dataclass(
        name,
        [(attr, str) for attr in dc.__annotations__.keys()],
        bases=(StatLabels, ),
        frozen=True
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
    stats: SummaryItems
    labels: StatLabels = field(init=False)

    def __post_init__(self) -> None:
        self.labels = SUMMARY_LABELS


@dataclass
class CorrelationStats(StatWithLabels):
    stats: CorrelationItems
    labels: StatLabels


@dataclass
class OutputStats(SummaryStats):
    ncc_stats: CorrelationStats
    mscc_stats: CorrelationStats


# Concrete dataclass models for results read from files
@dataclass
class CorrelationResult(CorrelationResultModel):
    max_shift: int = field(init=False)

    def __post_init__(self) -> None:
        self.max_shift = len(self.cc) - 1
        if self.max_shift < 0:
            raise ValueError("Cross-correlation array must have at least one element.")


@dataclass
class NCCResult(CorrelationResult, NCCResultModel):
    pass


@dataclass
class MSCCResult(CorrelationResult, MSCCResultModel):
    pass


@dataclass
class BothChromResult(BothChromResultModel):
    chrom: NCCResult
    mappable_chrom: MSCCResult


@dataclass
class NCCGenomeWideResult(NCCGenomeWideResultModel):
    chroms: Dict[str, NCCResult]


@dataclass
class MSCCGenomeWideResult(MSCCGenomeWideResultModel):
    chroms: Dict[str, MSCCResult]


@dataclass
class BothGenomeWideResult(BothGenomeWideResultModel):
    chroms: Dict[str, NCCResult]
    mappable_chroms: Dict[str, MSCCResult]
