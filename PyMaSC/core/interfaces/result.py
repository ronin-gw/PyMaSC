from abc import ABC
from typing import Tuple, Dict, Any, Optional, Generic, TypeVar
from dataclasses import dataclass

import sys
if sys.version_info >= (3, 9):
    from typing import Self
else:
    from typing_extensions import Self

import numpy as np


@dataclass
class AbstractDataclass(ABC):
    def __new__(cls, *args: Any, **kwargs: Any) -> Self:
        if cls == AbstractDataclass or (cls.__bases__ and cls.__bases__[0] == AbstractDataclass):
            raise TypeError("Cannot instantiate abstract class.")
        return super().__new__(cls)


@dataclass
class CorrelationResult(AbstractDataclass):
    genomelen: int
    forward_sum: Any
    reverse_sum: Any
    ccbins: np.ndarray


@dataclass
class NCCResult(CorrelationResult):
    forward_sum: int
    reverse_sum: int
    ccbins: np.ndarray


@dataclass
class MSCCResult(CorrelationResult):
    forward_sum: np.ndarray
    reverse_sum: np.ndarray
    ccbins: np.ndarray

    mappable_len: Optional[Tuple[int]]


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
