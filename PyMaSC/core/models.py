"""Data models and configuration classes for PyMaSC architecture.

This module defines the data models, configuration classes, and result
structures used throughout the refactored PyMaSC architecture. These
classes provide structured data containers with validation and type safety.

Key model categories:
- Configuration models: Algorithm, execution, and I/O settings
- Result models: Calculation results and statistics
- Request/Response models: API-like data transfer objects
- Validation utilities: Configuration validation and error handling
"""
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
import numpy as np
from enum import Enum


class CalculationTarget(Enum):
    """What type of cross-correlation to calculate."""
    NCC = "ncc"           # Naive Cross-Correlation only
    MSCC = "mscc"         # Mappability-Sensitive Cross-Correlation only
    BOTH = "both"         # Both NCC and MSCC


class ImplementationAlgorithm(Enum):
    """How to implement the cross-correlation calculation."""
    BITARRAY = "bitarray"     # Memory-intensive but fast BitArray algorithm
    SUCCESSIVE = "successive"  # Memory-efficient successive algorithm


class ExecutionMode(Enum):
    """Execution mode for calculations."""
    SINGLE_PROCESS = "single"
    MULTI_PROCESS = "multi"


@dataclass
class CalculationConfig:
    """Configuration for cross-correlation calculations.

    This class encapsulates all parameters needed for cross-correlation
    calculation, providing validation and default values.

    Attributes:
        target: What type of cross-correlation to calculate
        implementation: How to implement the calculation
        max_shift: Maximum shift distance for correlation
        mapq_criteria: Minimum mapping quality threshold
        skip_ncc: Whether to skip naive cross-correlation calculation
        read_length: Read length (None for auto-detection)
        references: List of chromosome names to process
        lengths: List of chromosome lengths
        esttype: Read length estimation type (mean, median, mode, min, max)
        chromfilter: Chromosome filtering patterns (include/exclude rules)
    """
    max_shift: int
    mapq_criteria: int
    target: CalculationTarget
    implementation: ImplementationAlgorithm
    skip_ncc: bool = False
    read_length: Optional[int] = None
    references: List[str] = field(default_factory=list)
    lengths: List[int] = field(default_factory=list)
    esttype: str = "MEDIAN"
    chromfilter: Optional[List[Tuple[bool, List[str]]]] = None

    def __post_init__(self) -> None:
        """Validate configuration after initialization."""
        if self.max_shift <= 0:
            raise ValueError("max_shift must be positive")
        if self.mapq_criteria < 0:
            raise ValueError("mapq_criteria must be non-negative")
        if self.read_length is not None and self.read_length <= 0:
            raise ValueError("read_length must be positive if specified")
        if len(self.references) != len(self.lengths):
            raise ValueError("references and lengths must have same length")


@dataclass
class MappabilityConfig:
    """Configuration for mappability correction.

    Attributes:
        mappability_path: Path to BigWig mappability file
        mappability_stats_path: Optional path to pre-computed stats
        read_len: Read length for mappability calculation
        shift_size: Shift size for mappability calculation
    """
    mappability_path: Optional[Path] = None
    mappability_stats_path: Optional[Path] = None
    read_len: Optional[int] = None
    shift_size: Optional[int] = None

    def is_enabled(self) -> bool:
        """Check if mappability correction is enabled."""
        return self.mappability_path is not None


@dataclass
class ExecutionConfig:
    """Configuration for execution parameters.

    Attributes:
        mode: Execution mode (single or multi-process)
        worker_count: Number of worker processes
        chrom_filter: Chromosome filtering patterns
    """
    mode: ExecutionMode = ExecutionMode.SINGLE_PROCESS
    worker_count: int = 1
    chrom_filter: Optional[Dict[str, Any]] = None

    def __post_init__(self) -> None:
        """Validate execution configuration."""
        if self.worker_count <= 0:
            raise ValueError("worker_count must be positive")


@dataclass
class WorkerConfig:
    """Configuration for individual worker processes.

    Attributes:
        calculation_config: Core calculation parameters
        mappability_config: Mappability correction settings
        progress_enabled: Whether to enable progress reporting
        logger_lock: Thread lock for safe logging
    """
    calculation_config: CalculationConfig
    mappability_config: MappabilityConfig
    progress_enabled: bool = True
    logger_lock: Optional[Any] = None


# Type aliases for convenience
ConfigDict = Dict[str, Any]
ResultDict = Dict[str, Any]
ChromosomeData = Tuple[str, int, int, np.ndarray]  # (chrom, forward, reverse, bins)