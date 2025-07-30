from typing import List, Dict, Any, Union, Optional, cast
from dataclasses import dataclass, field, asdict
from abc import abstractmethod

import numpy as np
import numpy.typing as npt

from .interfaces.result import (
    CorrelationResultModel,
    ChromResult, NCCResultModel, MSCCResultModel, BothChromResultModel,
    GenomeWideResultModel, NCCGenomeWideResultModel, MSCCGenomeWideResultModel, BothGenomeWideResultModel
)
from PyMaSC.utils.calc import npcalc_with_logging_warn


@dataclass
class CorrelationResult(CorrelationResultModel):
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


@dataclass
class NCCResult(NCCResultModel, CorrelationResult):
    forward_read_len_sum: int
    reverse_read_len_sum: int
    ccbins: Union[List[float], npt.NDArray[np.int64]]
    cc: npt.NDArray[np.float64] = field(init=False)

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
class MSCCResult(MSCCResultModel, CorrelationResult):
    forward_read_len_sum: Optional[int]
    reverse_read_len_sum: Optional[int]
    ccbins: Union[List[float], npt.NDArray[np.int64]]
    cc: npt.NDArray[np.float64] = field(init=False)

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
class BothChromResult(BothChromResultModel):
    """For multiprocessing workers"""

    chrom: NCCResult
    mappable_chrom: MSCCResult


@dataclass
class GenomeWideResult(GenomeWideResultModel):
    forward_read_len_sum: int
    reverse_read_len_sum: int


@dataclass
class NCCGenomeWideResult(NCCGenomeWideResultModel, GenomeWideResult):
    forward_sum: int
    reverse_sum: int
    chroms: Dict[str, NCCResult]


@dataclass
class MSCCGenomeWideResult(MSCCGenomeWideResultModel, GenomeWideResult):
    chroms: Dict[str, MSCCResult]


@dataclass
class BothGenomeWideResult(BothGenomeWideResultModel, GenomeWideResult):
    mappable_chroms: Dict[str, MSCCResult]


def aggregate_results(results: Dict[str, ChromResult]) -> GenomeWideResult:
    """Aggregate results from multiple chromosomes into a genome-wide result."""
    first_item = next(iter(results.values()))

    if isinstance(first_item, NCCResult):
        assert all(isinstance(res, NCCResult) for res in results.values())
        return _aggregate_ncc_results(cast(Dict[str, NCCResult], results))
    elif isinstance(first_item, MSCCResult):
        assert all(isinstance(res, MSCCResult) for res in results.values())
        return _aggregate_mscc_results(cast(Dict[str, MSCCResult], results))
    elif isinstance(first_item, BothChromResult):
        assert all(isinstance(res, BothChromResult) for res in results.values())
        return _aggregate_both_results(cast(Dict[str, BothChromResult], results))
    else:
        raise TypeError("Unknown result type")


@dataclass
class AggregateItems:
    genomelen: int = 0
    forward_read_len_sum: int = 0
    reverse_read_len_sum: int = 0

    @abstractmethod
    def __add__(self, other: Any) -> "AggregateItems":
        pass


@dataclass
class NCCAggregateItems(AggregateItems):
    forward_sum: int = 0
    reverse_sum: int = 0

    def __add__(self, other: NCCResult) -> "NCCAggregateItems":
        return NCCAggregateItems(
            genomelen=self.genomelen + other.genomelen,
            forward_sum=self.forward_sum + other.forward_sum,
            reverse_sum=self.reverse_sum + other.reverse_sum,
            forward_read_len_sum=self.forward_read_len_sum + other.forward_read_len_sum,
            reverse_read_len_sum=self.reverse_read_len_sum + other.reverse_read_len_sum,
        )


def _aggregate_ncc_results(results: Dict[str, NCCResult]) -> NCCGenomeWideResult:
    total = NCCAggregateItems()

    for res in results.values():
        total += res

    return NCCGenomeWideResult(chroms=results, **asdict(total))


@dataclass
class MSCCAggregateItems(AggregateItems):
    def __add__(self, other: MSCCResult) -> "MSCCAggregateItems":
        assert other.forward_read_len_sum is not None
        assert other.reverse_read_len_sum is not None
        return MSCCAggregateItems(
            genomelen=self.genomelen + other.genomelen,
            forward_read_len_sum=self.forward_read_len_sum + other.forward_read_len_sum,
            reverse_read_len_sum=self.reverse_read_len_sum + other.reverse_read_len_sum,
        )


def _aggregate_mscc_results(results: Dict[str, MSCCResult]) -> MSCCGenomeWideResult:
    total = MSCCAggregateItems()

    for res in results.values():
        total += res

    return MSCCGenomeWideResult(chroms=results, **asdict(total))


def _aggregate_both_results(results: Dict[str, BothChromResult]) -> BothGenomeWideResult:
    ncc = _aggregate_ncc_results({chrom: res.chrom for chrom, res in results.items()})
    mscc = _aggregate_mscc_results({chrom: res.mappable_chrom for chrom, res in results.items()})
    return BothGenomeWideResult(
        genomelen=ncc.genomelen,
        forward_sum=ncc.forward_sum,
        reverse_sum=ncc.reverse_sum,
        forward_read_len_sum=ncc.forward_read_len_sum,
        reverse_read_len_sum=ncc.reverse_read_len_sum,
        chroms=ncc.chroms,
        mappable_chroms=mscc.chroms
    )
