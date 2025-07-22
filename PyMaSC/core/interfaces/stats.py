from dataclasses import dataclass, field
from abc import abstractmethod
from typing import Optional, Dict, TypeVar, Generic

import numpy as np
import numpy.typing as npt

from .result import AbstractDataclass


@dataclass
class CCQualityMetrics:
    fragment_length: Optional[int] = None
    ccfl: Optional[float] = None
    fwhm: Optional[int] = None
    nsc: Optional[float] = None
    rsc: Optional[float] = None
    vsn: Optional[float] = None


@dataclass
class CCStats(AbstractDataclass):
    read_len: int
    cc_min: float
    ccrl: float

    metrics_at_expected_length: CCQualityMetrics
    metrics_at_estimated_length: CCQualityMetrics

    @property
    @abstractmethod
    def genomelen(self) -> int:
        pass

    @property
    @abstractmethod
    def forward_reads(self) -> int:
        pass

    @property
    @abstractmethod
    def reverse_reads(self) -> int:
        pass


@dataclass
class CCArrays(AbstractDataclass):
    cc: npt.NDArray[np.float64]
    avr_cc: npt.NDArray[np.float64]
    est_lib_len: Optional[int]


TStats = TypeVar("TStats", bound=CCStats)


@dataclass
class ChromosomeStats(Generic[TStats], CCArrays):
    stats: TStats


@dataclass
class GenomeWideStats:
    whole_ncc_stats: Optional[ChromosomeStats] = None
    whole_mscc_stats: Optional[ChromosomeStats] = None
    ncc_stats: Optional[Dict[str, ChromosomeStats]] = None
    mscc_stats: Optional[Dict[str, ChromosomeStats]] = None

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
