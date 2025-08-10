from dataclasses import dataclass
from abc import abstractmethod
from functools import cached_property
from typing import Optional, Dict, Tuple, TypeVar, Generic, Protocol
import logging

import numpy as np
import numpy.typing as npt

from scipy.stats import chi2

from .result import AbstractDataclass

logger = logging.getLogger(__name__)


@dataclass
class CCQualityMetrics:
    fragment_length: Optional[int] = None
    ccfl: Optional[float] = None
    fwhm: Optional[int] = None
    nsc: Optional[float] = None
    rsc: Optional[float] = None
    vsn: Optional[float] = None


TCount = TypeVar("TCount", int, npt.NDArray[np.int64])


class CCStats(Protocol, Generic[TCount]):
    read_len: int
    cc_min: float
    ccrl: float

    metrics_at_expected_length: CCQualityMetrics
    metrics_at_estimated_length: CCQualityMetrics
    genomelen: TCount
    forward_reads: TCount
    reverse_reads: TCount

    # Attributes for the representative values
    # i.e. return itself in NCC, return 1st item of array in MSCC
    @property
    @abstractmethod
    def genomelen_repr(self) -> int:
        pass

    @property
    @abstractmethod
    def forward_reads_repr(self) -> int:
        pass

    @property
    @abstractmethod
    def reverse_reads_repr(self) -> int:
        pass

    def check_strand_balance(self, chi2_p_thresh: float, label: str) -> None:
        a = self.forward_reads_repr
        b = self.reverse_reads_repr
        if a == 0 and b == 0:
            return

        sum_ = a + b
        chi2_val = (((a - sum_ / 2.) ** 2) + ((b - sum_ / 2.) ** 2)) / sum_
        chi2_p = chi2.sf(chi2_val, 1)

        if chi2_p <= chi2_p_thresh:
            logger.warning(f"{label} Forward/Reverse read count imbalance.")
            logger.warning(f"+/- = {a} / {b}, Chi-squared test p-val = {chi2_p:.5g} <= {chi2_p_thresh}")
        else:
            logger.info(f"{label} Forward/Reverse read count +/- = {a} / {b}")
            logger.info(f"Chi-squared test p-val = {chi2_p:.5g} > {chi2_p_thresh}")


NCCStats = CCStats[int]
MSCCStats = CCStats[npt.NDArray[np.int64]]


@dataclass
class CCArrays(AbstractDataclass):
    cc: npt.NDArray[np.float64]
    avr_cc: npt.NDArray[np.float64]
    est_lib_len: Optional[int]
    mv_avr_filter_len: int


TStats = TypeVar("TStats", bound=CCStats, covariant=True)


@dataclass
class ChromosomeStats(Generic[TStats], CCArrays):
    stats: TStats


@dataclass
class WholeGenomeStats(ChromosomeStats[TStats]):
    cc_upper: npt.NDArray[np.float64]
    cc_lower: npt.NDArray[np.float64]


@dataclass
class GenomeWideStats:
    whole_ncc_stats: Optional[WholeGenomeStats[NCCStats]] = None
    whole_mscc_stats: Optional[WholeGenomeStats[MSCCStats]] = None
    ncc_stats: Optional[Dict[str, ChromosomeStats[NCCStats]]] = None
    mscc_stats: Optional[Dict[str, ChromosomeStats[MSCCStats]]] = None

    @property
    def has_ncc(self) -> bool:
        return self.whole_ncc_stats is not None

    @property
    def has_mscc(self) -> bool:
        return self.whole_mscc_stats is not None

    @property
    def read_len(self) -> int:
        if self.whole_ncc_stats is not None:
            return self.whole_ncc_stats.stats.read_len
        elif self.whole_mscc_stats is not None:
            return self.whole_mscc_stats.stats.read_len
        else:
            raise ValueError("No read length available in GenomeWideStats.")

    @property
    def expected_lib_len(self) -> Optional[int]:
        if self.whole_ncc_stats is not None:
            return self.whole_ncc_stats.stats.metrics_at_expected_length.fragment_length
        elif self.whole_mscc_stats is not None:
            return self.whole_mscc_stats.stats.metrics_at_expected_length.fragment_length
        else:
            raise ValueError("No expected library length available in GenomeWideStats.")

    @property
    def est_lib_len(self) -> Optional[int]:
        if self.whole_mscc_stats is not None:
            return self.whole_mscc_stats.est_lib_len
        elif self.whole_ncc_stats is not None:
            return self.whole_ncc_stats.est_lib_len
        else:
            raise ValueError("No estimated library length available in GenomeWideStats.")

    @cached_property
    def references(self) -> Tuple[str, ...]:
        if self.ncc_stats is not None:
            return tuple(self.ncc_stats.keys())
        elif self.mscc_stats is not None:
            return tuple(self.mscc_stats.keys())
        else:
            raise ValueError("No chromosome stats available in GenomeWideStats.")
