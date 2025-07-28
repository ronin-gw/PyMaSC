from abc import ABC, abstractmethod
from typing import Tuple, List, Dict, Any, Optional, Union
from dataclasses import dataclass, field

import sys
if sys.version_info >= (3, 9):
    from typing import Self
else:
    from typing_extensions import Self

import numpy as np
import numpy.typing as npt

from PyMaSC.utils.calc import npcalc_with_logging_warn


@dataclass
class AbstractDataclass(ABC):
    def __new__(cls, *args: Any, **kwargs: Any) -> Self:
        if cls == AbstractDataclass or (cls.__bases__ and cls.__bases__[0] == AbstractDataclass):
            raise TypeError("Cannot instantiate abstract class.")
        return super().__new__(cls)


@dataclass
class CorrelationResult(AbstractDataclass):
    max_shift: int
    read_len: int
    genomelen: int
    forward_sum: Any
    reverse_sum: Any
    ccbins: Union[List[float], npt.NDArray[np.int64]]

    cc: npt.NDArray[np.float64] = field(init=False)

    @abstractmethod
    def calc_cc(self) -> None:
        pass

    @staticmethod
    @npcalc_with_logging_warn
    def _calc_cc(
        forward_sum: Union[float, npt.NDArray[np.float64]],
        reverse_sum: Union[float, npt.NDArray[np.float64]],
        ccbins: Union[List[float], npt.NDArray[np.int64]],
        totlen: Union[int, npt.NDArray[np.float64]],
        denom: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Calculate normalized cross-correlation using binomial variance model."""
        ccbins = np.array(ccbins, dtype=np.int64)
        #
        if ccbins.sum() == 0:
            return np.full_like(ccbins, np.nan, dtype=np.float64)

        forward_mean = forward_sum / totlen
        reverse_mean = reverse_sum / totlen

        forward_var = forward_mean * (1 - forward_mean)
        reverse_var = reverse_mean * (1 - reverse_mean)

        sum_prod = forward_mean * reverse_mean
        var_geomean: Union[float, npt.NDArray[np.float64]] = (forward_var * reverse_var) ** 0.5
        return (ccbins / denom - sum_prod) / var_geomean


class ChromResult:
    pass


@dataclass
class NCCResult(ChromResult, CorrelationResult):
    forward_sum: int
    reverse_sum: int
    forward_read_len_sum: int
    reverse_read_len_sum: int
    ccbins: Union[List[float], npt.NDArray[np.int64]]

    def calc_cc(self) -> None:
        denom = self.genomelen - np.array(range(self.max_shift + 1), dtype=np.float64)
        self.cc = self._calc_cc(
            float(self.forward_sum),
            float(self.reverse_sum),
            self.ccbins[:self.max_shift + 1],
            self.genomelen,
            denom
        )


@dataclass
class MSCCResult(ChromResult, CorrelationResult):
    forward_sum: npt.NDArray[np.int64]
    reverse_sum: npt.NDArray[np.int64]
    forward_read_len_sum: Optional[int]
    reverse_read_len_sum: Optional[int]
    ccbins: Union[List[float], npt.NDArray[np.int64]]

    mappable_len: Optional[Tuple[int]]

    def calc_cc(self) -> None:
        assert self.mappable_len is not None, "mappable_len must be set before calculating CC."
        totlen = np.array(self.mappable_len, dtype=np.float64)
        totlen = np.concatenate((
            totlen[:self.read_len][::-1], totlen[1:]
        ))[:self.max_shift + 1]

        self.cc = self._calc_cc(
            np.array(self.forward_sum[:self.max_shift + 1], dtype=np.float64),
            np.array(self.reverse_sum[:self.max_shift + 1], dtype=np.float64),
            self.ccbins[:self.max_shift + 1],
            totlen,
            totlen
        )


@dataclass
class BothChromResult(ChromResult):
    """For multiprocessing workers"""

    chrom: NCCResult
    mappable_chrom: MSCCResult


@dataclass
class GenomeWideResult(AbstractDataclass):
    genomelen: int
    forward_read_len_sum: int
    reverse_read_len_sum: int


@dataclass
class NCCGenomeWideResult(GenomeWideResult):
    forward_sum: int
    reverse_sum: int
    chroms: Dict[str, NCCResult]


@dataclass
class MSCCGenomeWideResult(GenomeWideResult):
    chroms: Dict[str, MSCCResult]


@dataclass
class BothGenomeWideResult(NCCGenomeWideResult):
    mappable_chroms: Dict[str, MSCCResult]
