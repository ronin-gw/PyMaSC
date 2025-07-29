from abc import ABC
from typing import Tuple, Mapping, Any, Optional
from dataclasses import dataclass

import sys
if sys.version_info >= (3, 9):
    from typing import Self
else:
    from typing_extensions import Self

import numpy as np
import numpy.typing as npt


@dataclass
class AbstractDataclass(ABC):
    def __new__(cls, *args: Any, **kwargs: Any) -> Self:
        if cls == AbstractDataclass or (cls.__bases__ and cls.__bases__[0] == AbstractDataclass):
            raise TypeError("Cannot instantiate abstract class.")
        return super().__new__(cls)


@dataclass
class CorrelationResultModel(AbstractDataclass):
    max_shift: int
    read_len: int
    genomelen: int
    forward_sum: Any
    reverse_sum: Any
    cc: npt.NDArray[np.float64]


class ChromResult:
    pass


@dataclass
class NCCResultModel(ChromResult, CorrelationResultModel):
    forward_sum: int
    reverse_sum: int


@dataclass
class MSCCResultModel(ChromResult, CorrelationResultModel):
    forward_sum: npt.NDArray[np.int64]
    reverse_sum: npt.NDArray[np.int64]
    mappable_len: Optional[Tuple[int, ...]]


@dataclass
class BothChromResultModel(ChromResult):
    """For multiprocessing workers"""

    chrom: NCCResultModel
    mappable_chrom: MSCCResultModel


@dataclass
class GenomeWideResultModel(AbstractDataclass):
    genomelen: int


@dataclass
class NCCGenomeWideResultModel(GenomeWideResultModel):
    forward_sum: int
    reverse_sum: int
    chroms: Mapping[str, NCCResultModel]


@dataclass
class MSCCGenomeWideResultModel(GenomeWideResultModel):
    chroms: Mapping[str, MSCCResultModel]


@dataclass
class BothGenomeWideResultModel(NCCGenomeWideResultModel):
    mappable_chroms: Mapping[str, MSCCResultModel]
