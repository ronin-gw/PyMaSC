"""Result data structures and aggregation for cross-correlation calculations.

Provides concrete implementations of result models and functions for
combining cross-correlation results across chromosomes.
"""
from typing import List, Dict, Any, Union, Optional, cast
from dataclasses import dataclass, field, asdict
from abc import ABC, abstractmethod

import sys
if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

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
    """Base class for cross-correlation calculation results.

    Provides common functionality for calculating normalized cross-correlation
    from read count data using binomial variance model.
    """
    ccbins: Union[List[float], npt.NDArray[np.int64]]
    cc: npt.NDArray[np.float64] = field(init=False)

    @abstractmethod
    def calc_cc(self) -> None:
        """Calculate cross-correlation values from count data."""
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
    """Naive Cross-Correlation (NCC) analysis result implementation.

    Contains NCC calculation results with integer read counts and
    standard cross-correlation computation.
    """
    forward_read_len_sum: int
    reverse_read_len_sum: int
    ccbins: Union[List[float], npt.NDArray[np.int64]]
    cc: npt.NDArray[np.float64] = field(init=False)

    def calc_cc(self) -> None:
        """Calculate NCC values using genome length denominator."""
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
    """Mappability-Sensitive Cross-Correlation (MSCC) analysis result implementation.

    Contains MSCC calculation results with mappability-corrected denominators
    for regions with variable mappability.
    """
    forward_read_len_sum: Optional[int]
    reverse_read_len_sum: Optional[int]
    ccbins: Union[List[float], npt.NDArray[np.int64]]
    cc: npt.NDArray[np.float64] = field(init=False)

    def calc_cc(self) -> None:
        """Calculate MSCC values using mappable length denominators."""
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


class EmptyResult(ChromResult, ABC):
    """Base class for empty chromosome results.

    Used to create results for chromosomes with no reads, ensuring
    consistent genome length totals between single-process and
    parallel processing modes.
    """

    @classmethod
    @abstractmethod
    def create_empty(cls, genome_length: int, max_shift: int, read_len: int) -> Self:
        """Create an empty result with zero/NaN values but correct chromosome info."""
        pass


@dataclass
class EmptyNCCResult(EmptyResult, NCCResult):
    """Empty NCC result for chromosomes with no data.

    Used in parallel processing to ensure genome length consistency
    when chromosomes have no reads but need to be included in statistics.
    """

    @classmethod
    def create_empty(cls, genome_length: int, max_shift: int, read_len: int) -> Self:
        """Create an empty NCC result for a chromosome with no data.

        Args:
            chromosome: Chromosome name
            genome_length: Length of the chromosome
            max_shift: Maximum shift value
            read_len: Read length

        Returns:
            EmptyNCCResult with zero/NaN values but correct chromosome info
        """
        # Create zero ccbins array
        ccbins = [0.0] * (max_shift + 1)

        result = cls(
            genomelen=genome_length,
            max_shift=max_shift,
            read_len=read_len,
            forward_sum=0,
            reverse_sum=0,
            forward_read_len_sum=0,
            reverse_read_len_sum=0,
            ccbins=ccbins
        )

        # Calculate CC (will be all NaN due to zero sums)
        result.calc_cc()
        return result


@dataclass
class EmptyMSCCResult(EmptyResult, MSCCResult):
    """Empty MSCC result for chromosomes with no data.

    Used in parallel processing to ensure genome length consistency
    when chromosomes have no reads but need to be included in statistics.
    """

    @classmethod
    def create_empty(cls, genome_length: int, max_shift: int, read_len: int) -> Self:
        """Create an empty MSCC result for a chromosome with no data.

        Args:
            chromosome: Chromosome name
            genome_length: Length of the chromosome
            max_shift: Maximum shift value
            read_len: Read length

        Returns:
            EmptyMSCCResult with zero/NaN values but correct chromosome info
        """
        # Create zero arrays - MSCC requires numpy arrays for proper aggregation
        ccbins = [0.0] * (max_shift + 1)
        forward_sum = np.zeros(max_shift + 1, dtype=np.int64)
        reverse_sum = np.zeros(max_shift + 1, dtype=np.int64)
        # For empty results, mappable_len should represent the full chromosome length
        # at all shift positions (no mappability filtering applied)
        # Must be tuple to match MSCCResultModel.mappable_len type definition
        mappable_len = tuple([0] * (max_shift + 1))

        result = cls(
            genomelen=genome_length,
            max_shift=max_shift,
            read_len=read_len,
            forward_sum=forward_sum,
            reverse_sum=reverse_sum,
            forward_read_len_sum=0,
            reverse_read_len_sum=0,
            ccbins=ccbins,
            mappable_len=mappable_len
        )

        # Calculate CC (will be all NaN due to zero sums)
        result.calc_cc()
        return result


@dataclass
class EmptyBothChromResult(EmptyResult, BothChromResult):
    """Empty Both result for chromosomes with no data.

    Contains both empty NCC and MSCC results for complete coverage.
    """

    @classmethod
    def create_empty(cls, genome_length: int, max_shift: int, read_len: int) -> Self:
        """Create an empty Both result for a chromosome with no data.

        Args:
            chromosome: Chromosome name
            genome_length: Length of the chromosome
            max_shift: Maximum shift value
            read_len: Read length

        Returns:
            EmptyBothChromResult with empty NCC and MSCC components
        """
        empty_ncc = EmptyNCCResult.create_empty(genome_length, max_shift, read_len)
        empty_mscc = EmptyMSCCResult.create_empty(genome_length, max_shift, read_len)

        return cls(
            chrom=empty_ncc,
            mappable_chrom=empty_mscc
        )


@dataclass
class GenomeWideResult(GenomeWideResultModel):
    """Base class for genome-wide analysis result implementations.

    Contains common fields for genome-wide aggregated results.
    """
    forward_read_len_sum: int
    reverse_read_len_sum: int


@dataclass
class NCCGenomeWideResult(NCCGenomeWideResultModel, GenomeWideResult):
    """Genome-wide NCC analysis result implementation.

    Contains aggregated NCC results from all chromosomes with
    total read counts and per-chromosome result dictionary.
    """
    forward_sum: int
    reverse_sum: int
    chroms: Dict[str, NCCResult]


@dataclass
class MSCCGenomeWideResult(MSCCGenomeWideResultModel, GenomeWideResult):
    """Genome-wide MSCC analysis result implementation.

    Contains aggregated MSCC results from all chromosomes with
    per-chromosome result dictionary.
    """
    chroms: Dict[str, MSCCResult]


@dataclass
class BothGenomeWideResult(BothGenomeWideResultModel, GenomeWideResult):
    """Combined genome-wide NCC and MSCC analysis result implementation.

    Extends NCC genome-wide results with additional MSCC chromosome data.
    """
    mappable_chroms: Dict[str, MSCCResult]


def aggregate_results(results: Dict[str, ChromResult]) -> GenomeWideResult:
    """Aggregate results from multiple chromosomes into a genome-wide result.

    This function handles both regular and empty results seamlessly. Empty results
    (EmptyNCCResult, EmptyMSCCResult, EmptyBothChromResult) are treated as their
    parent types for aggregation, ensuring genome length consistency in parallel
    processing scenarios where some chromosomes have no data.

    Args:
        results: Dictionary mapping chromosome names to their calculation results

    Returns:
        GenomeWideResult containing aggregated statistics across all chromosomes

    Note:
        Empty results contribute their genome length to totals but have zero/NaN
        statistical values, which is the desired behavior for maintaining
        consistency between single-process and parallel processing modes.
    """
    if not results:
        raise ValueError("Cannot aggregate empty results dictionary")

    first_item = next(iter(results.values()))

    # Check result type hierarchy - Empty classes inherit from base classes
    # so isinstance checks work correctly for mixed regular/empty results

    #
    if isinstance(first_item, BothChromResult):
        # This includes both BothChromResult and EmptyBothChromResult instances
        assert all(isinstance(res, BothChromResult) for res in results.values()), \
            "All results must be BothChromResult instances (including EmptyBothChromResult)"

        _results = cast(Dict[str, BothChromResult], results)
        if all(res.chrom is None for res in _results.values() if not isinstance(res, EmptyResult)):
            results = {chrom: res.mappable_chrom for chrom, res in _results.items()}
            return _aggregate_mscc_results(cast(Dict[str, MSCCResult], results))
        elif all(res.mappable_chrom is None for res in _results.values() if not isinstance(res, EmptyResult)):
            results = {chrom: res.chrom for chrom, res in _results.items()}
            return _aggregate_ncc_results(cast(Dict[str, NCCResult], results))
        else:
            return _aggregate_both_results(cast(Dict[str, BothChromResult], results))

    elif isinstance(first_item, NCCResult):
        # This includes both NCCResult and EmptyNCCResult instances
        assert all(isinstance(res, NCCResult) for res in results.values()), \
            "All results must be NCCResult instances (including EmptyNCCResult)"
        return _aggregate_ncc_results(cast(Dict[str, NCCResult], results))
    elif isinstance(first_item, MSCCResult):
        # This includes both MSCCResult and EmptyMSCCResult instances
        assert all(isinstance(res, MSCCResult) for res in results.values()), \
            "All results must be MSCCResult instances (including EmptyMSCCResult)"
        return _aggregate_mscc_results(cast(Dict[str, MSCCResult], results))

    else:
        raise TypeError(f"Unknown result type: {type(first_item)}")


@dataclass
class AggregateItems:
    """Base class for aggregating chromosome results into genome-wide totals.

    Contains common fields for accumulating statistics across chromosomes.
    """
    genomelen: int = 0
    forward_read_len_sum: int = 0
    reverse_read_len_sum: int = 0

    @abstractmethod
    def __add__(self, other: Any) -> "AggregateItems":
        pass


@dataclass
class NCCAggregateItems(AggregateItems):
    """Aggregation container for NCC results.

    Accumulates integer read counts from NCC chromosome results
    for genome-wide NCC statistics.
    """
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
    """Aggregate NCC results from all chromosomes.

    Args:
        results: Dictionary mapping chromosome names to NCC results

    Returns:
        NCCGenomeWideResult containing aggregated statistics
    """
    total = NCCAggregateItems()

    for res in results.values():
        total += res

    return NCCGenomeWideResult(chroms=results, **asdict(total))


@dataclass
class MSCCAggregateItems(AggregateItems):
    """Aggregation container for MSCC results.

    Accumulates read length statistics from MSCC chromosome results
    for genome-wide MSCC statistics.
    """
    def __add__(self, other: MSCCResult) -> "MSCCAggregateItems":
        assert other.forward_read_len_sum is not None
        assert other.reverse_read_len_sum is not None
        return MSCCAggregateItems(
            genomelen=self.genomelen + other.genomelen,
            forward_read_len_sum=self.forward_read_len_sum + other.forward_read_len_sum,
            reverse_read_len_sum=self.reverse_read_len_sum + other.reverse_read_len_sum,
        )


def _aggregate_mscc_results(results: Dict[str, MSCCResult]) -> MSCCGenomeWideResult:
    """Aggregate MSCC results from all chromosomes.

    Args:
        results: Dictionary mapping chromosome names to MSCC results

    Returns:
        MSCCGenomeWideResult containing aggregated statistics
    """
    total = MSCCAggregateItems()

    for res in results.values():
        total += res

    return MSCCGenomeWideResult(chroms=results, **asdict(total))


def _aggregate_both_results(results: Dict[str, BothChromResult]) -> BothGenomeWideResult:
    """Aggregate combined NCC and MSCC results from all chromosomes.

    Args:
        results: Dictionary mapping chromosome names to BothChromResult objects

    Returns:
        BothGenomeWideResult containing aggregated NCC and MSCC statistics
    """
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
