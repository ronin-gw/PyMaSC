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
        timeout: Worker timeout in seconds
        chrom_filter: Chromosome filtering patterns
    """
    mode: ExecutionMode = ExecutionMode.SINGLE_PROCESS
    worker_count: int = 1
    timeout: Optional[float] = None
    chrom_filter: Optional[Dict[str, Any]] = None

    def __post_init__(self) -> None:
        """Validate execution configuration."""
        if self.worker_count <= 0:
            raise ValueError("worker_count must be positive")


@dataclass
class IOConfig:
    """Configuration for input/output operations.

    Attributes:
        input_paths: List of input BAM file paths
        output_dir: Output directory path
        output_names: Custom output names (optional)
        expected_suffixes: Expected output file suffixes
    """
    input_paths: List[Path]
    output_dir: Path
    output_names: Optional[List[str]] = None
    expected_suffixes: List[str] = field(default_factory=lambda: ['.pdf', '.tab', '.json'])

    def __post_init__(self) -> None:
        """Validate I/O configuration."""
        if not self.input_paths:
            raise ValueError("At least one input path must be specified")
        for path in self.input_paths:
            if not isinstance(path, Path):
                raise TypeError("input_paths must contain Path objects")


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


@dataclass
class ChromosomeResult:
    """Results for a single chromosome calculation.

    Attributes:
        chromosome: Chromosome name
        forward_count: Number of forward reads processed
        reverse_count: Number of reverse reads processed
        correlation_bins: Cross-correlation bins array
        mappable_length: Mappable length for the chromosome (if applicable)
        execution_time: Time taken for calculation
    """
    chromosome: str
    forward_count: int
    reverse_count: int
    correlation_bins: np.ndarray
    mappable_length: Optional[int] = None
    execution_time: Optional[float] = None


@dataclass
class CalculationResult:
    """Aggregated results from cross-correlation calculation.

    Attributes:
        chromosome_results: Results for each chromosome
        total_forward_reads: Total forward reads across all chromosomes
        total_reverse_reads: Total reverse reads across all chromosomes
        calculation_target: What type of cross-correlation was calculated
        implementation_algorithm: How the calculation was implemented
        execution_metadata: Additional execution information
    """
    chromosome_results: List[ChromosomeResult]
    total_forward_reads: int
    total_reverse_reads: int
    calculation_target: CalculationTarget
    implementation_algorithm: ImplementationAlgorithm
    execution_metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def processed_chromosomes(self) -> List[str]:
        """Get list of processed chromosome names."""
        return [result.chromosome for result in self.chromosome_results]

    @property
    def total_reads(self) -> int:
        """Get total number of reads processed."""
        return self.total_forward_reads + self.total_reverse_reads


@dataclass
class CalculationRequest:
    """Request object for cross-correlation calculation.

    This class encapsulates all configuration and parameters needed
    to execute a complete cross-correlation analysis.

    Attributes:
        calculation_config: Core calculation parameters
        mappability_config: Mappability correction settings
        execution_config: Execution mode and worker settings
        io_config: Input/output file settings
    """
    calculation_config: CalculationConfig
    mappability_config: MappabilityConfig
    execution_config: ExecutionConfig
    io_config: IOConfig

    def validate(self) -> List[str]:
        """Validate the entire request configuration.

        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []

        # Check algorithm-specific requirements using new conceptual approach
        if (self.calculation_config.target in [CalculationTarget.MSCC, CalculationTarget.BOTH]
            and not self.mappability_config.is_enabled()):
            errors.append(f"{self.calculation_config.target.value} target requires mappability configuration")

        # Check execution mode compatibility
        if (self.execution_config.mode == ExecutionMode.MULTI_PROCESS
            and self.execution_config.worker_count == 1):
            errors.append("Multi-process mode requires worker_count > 1")

        # Check I/O compatibility
        if (len(self.io_config.input_paths) > 1
            and self.io_config.output_names
            and len(self.io_config.output_names) != len(self.io_config.input_paths)):
            errors.append("output_names length must match input_paths length")

        return errors


@dataclass
class ValidationResult:
    """Result of configuration validation.

    Attributes:
        is_valid: Whether the configuration is valid
        errors: List of validation error messages
        warnings: List of validation warning messages
    """
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# Type aliases for convenience
ConfigDict = Dict[str, Any]
ResultDict = Dict[str, Any]
ChromosomeData = Tuple[str, int, int, np.ndarray]  # (chrom, forward, reverse, bins)