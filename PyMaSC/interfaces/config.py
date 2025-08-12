"""Configuration models and type definitions for PyMaSC.

Defines configuration protocols, enums, and the main PyMaSCConfig dataclass
for cross-correlation calculation parameters.
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Protocol
from pathlib import Path
from enum import Enum
from argparse import Namespace

import sys

if sys.version_info >= (3, 10):
    from typing import TypeGuard
else:
    from typing_extensions import TypeGuard

if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self


class CalculationTarget(Enum):
    """What type of cross-correlation to calculate."""
    NCC = "ncc"
    MSCC = "mscc"
    BOTH = "both"


class Algorithm(Enum):
    """How to implement the cross-correlation calculation."""
    BITARRAY = "bitarray"
    SUCCESSIVE = "successive"


class EstimationType(Enum):
    """Type of estimation to use for read length."""
    MEAN = "MEAN"
    MEDIAN = "MEDIAN"
    MODE = "MODE"
    MIN = "MIN"
    MAX = "MAX"


class CalculationConfig(Protocol):
    """Protocol for calculation-specific configuration parameters."""

    max_shift: int
    mapq_criteria: int
    target: CalculationTarget
    implementation: Algorithm
    nproc: int
    skip_ncc: bool
    esttype: EstimationType
    read_length: Optional[int] = None
    chromfilter: Optional[List[Tuple[bool, List[str]]]] = None


class MappabilityConfig(Protocol):
    """Protocol for mappability correction configuration."""
    max_shift: int
    read_length: int
    nproc: int
    mappability_path: Path
    mappability_stats_path: Optional[Path] = None


class StatConfig(Protocol):
    """Protocol for statistics calculation configuration."""
    read_length: int
    chi2_pval: float
    mv_avr_filter_len: int
    filter_mask_len: int
    min_calc_width: int
    expected_library_length: Optional[int] = None


@dataclass
class PyMaSCConfig:
    """Configuration for PyMaSC cross-correlation analysis.

    Central configuration object containing all parameters needed
    for cross-correlation calculations including algorithm selection,
    quality thresholds, and statistical parameters.
    """
    max_shift: int
    mapq_criteria: int
    target: CalculationTarget
    implementation: Algorithm
    nproc: int
    esttype: EstimationType

    chi2_pval: float
    mv_avr_filter_len: int
    filter_mask_len: int
    min_calc_width: int

    read_length: Optional[int] = None
    chromfilter: Optional[List[Tuple[bool, List[str]]]] = None
    ref2lengths: Dict[str, int] = field(default_factory=dict)

    mappability_path: Optional[Path] = None
    mappability_stats_path: Optional[Path] = None

    expected_library_length: Optional[int] = None

    @property
    def skip_ncc(self) -> bool:
        """State that indicates mappability info is supplied but NCC is not needed."""
        return self.target is CalculationTarget.MSCC

    @property
    def multiprocess(self) -> bool:
        """Check if the configuration is set for multiprocess execution."""
        return self.nproc > 1

    @property
    def references(self) -> Tuple[str, ...]:
        """Get the list of reference names."""
        return tuple(self.ref2lengths.keys())

    @property
    def lengths(self) -> Tuple[int, ...]:
        """Get the list of reference lengths."""
        return tuple(self.ref2lengths.values())

    @classmethod
    def from_args(cls, args: Namespace) -> Self:
        """Create configuration from parsed command-line arguments."""
        if args.mappability:
            target = CalculationTarget.MSCC if args.skip_ncc else CalculationTarget.BOTH
        else:
            target = CalculationTarget.NCC

        implementation = Algorithm.SUCCESSIVE if args.successive else Algorithm.BITARRAY

        return cls(
            max_shift=args.max_shift,
            mapq_criteria=args.mapq,
            target=target,
            implementation=implementation,
            nproc=args.process,
            esttype=EstimationType[args.readlen_estimator],
            chi2_pval=args.chi2_pval,
            mv_avr_filter_len=args.smooth_window,
            filter_mask_len=args.mask_size,
            min_calc_width=args.bg_avr_width,
            chromfilter=args.chromfilter,
            mappability_path=args.mappability,
            mappability_stats_path=args.mappability_stats,
            expected_library_length=args.library_length
        )


def mappability_configured(config: PyMaSCConfig) -> TypeGuard[MappabilityConfig]:
    """Check if mappability path is set."""
    return (
        config.read_length is not None and
        config.mappability_path is not None
    )


def stats_configured(config: PyMaSCConfig) -> TypeGuard[StatConfig]:
    """Check if statistics configuration is set."""
    return config.read_length is not None
