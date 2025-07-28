from typing import Dict, cast, Any
from dataclasses import dataclass, asdict
from abc import abstractmethod

from .interfaces.result import (
    ChromResult, NCCResult, MSCCResult, BothChromResult,
    GenomeWideResult, NCCGenomeWideResult, MSCCGenomeWideResult, BothGenomeWideResult
)


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
