"""Result data models for cross-correlation calculations.

Defines abstract and concrete result structures for chromosome and genome-wide results.
"""
from abc import ABC
from typing import Tuple, Mapping, Any, Optional
from dataclasses import dataclass

import sys
if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

import numpy as np
import numpy.typing as npt


@dataclass
class AbstractDataclass(ABC):
    """Abstract base class for dataclass-based result objects.

    Prevents direct instantiation of abstract dataclasses.
    """
    def __new__(cls, *args: Any, **kwargs: Any) -> Self:
        if cls == AbstractDataclass or (cls.__bases__ and cls.__bases__[0] == AbstractDataclass):
            raise TypeError("Cannot instantiate abstract class.")
        return super().__new__(cls)


@dataclass
class CorrelationResultModel(AbstractDataclass):
    """Base model for cross-correlation calculation results.

    Contains common fields for all correlation analysis results including
    read counts, genomic parameters, and correlation arrays.
    """
    max_shift: int
    read_len: int
    genomelen: int
    forward_sum: Any
    reverse_sum: Any
    cc: npt.NDArray[np.float64]


class ChromResult(ABC):
    """Abstract base class for chromosome-level analysis results."""
    pass


@dataclass
class NCCResultModel(ChromResult, CorrelationResultModel):
    """Model for Naive Cross-Correlation (NCC) results.

    Extends base correlation results with integer read counts
    for standard cross-correlation analysis.
    """
    forward_sum: int
    reverse_sum: int


@dataclass
class MSCCResultModel(ChromResult, CorrelationResultModel):
    """Model for Mappability-Sensitive Cross-Correlation (MSCC) results.

    Extends base correlation results with mappability arrays and
    per-position read counts for mappability-corrected analysis.
    """
    forward_sum: npt.NDArray[np.int64]
    reverse_sum: npt.NDArray[np.int64]
    mappable_len: Optional[Tuple[int, ...]]


@dataclass
class BothChromResultModel(ChromResult):
    """Combined NCC and MSCC results for a chromosome."""

    chrom: NCCResultModel
    mappable_chrom: MSCCResultModel


@dataclass
class GenomeWideResultModel(AbstractDataclass):
    """Base model for genome-wide analysis results.

    Contains genome-level summary information common to all
    genome-wide correlation analyses.
    """
    genomelen: int


@dataclass
class NCCGenomeWideResultModel(GenomeWideResultModel):
    """Model for genome-wide NCC analysis results.

    Contains total read counts and per-chromosome NCC results
    for complete genome-wide naive cross-correlation analysis.
    """
    forward_sum: int
    reverse_sum: int
    chroms: Mapping[str, NCCResultModel]


@dataclass
class MSCCGenomeWideResultModel(GenomeWideResultModel):
    """Model for genome-wide MSCC analysis results.

    Contains per-chromosome MSCC results for complete
    genome-wide mappability-sensitive analysis.
    """
    chroms: Mapping[str, MSCCResultModel]


@dataclass
class BothGenomeWideResultModel(NCCGenomeWideResultModel):
    """Model for combined genome-wide NCC and MSCC results.

    Extends NCC genome-wide results with additional MSCC data
    for comprehensive correlation analysis.
    """
    mappable_chroms: Mapping[str, MSCCResultModel]
