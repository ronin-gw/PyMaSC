"""Configuration management system for PyMaSC.

This module provides a centralized configuration management system using
the Builder pattern for complex configuration construction and validation.
It bridges the existing argparse-based configuration with the new structured
data models.

Key components:
- ConfigurationBuilder: Builder pattern for step-by-step configuration
- ConfigurationService: Centralized validation and configuration management
- ArgumentConverter: Bridge between argparse and new configuration models
- ConfigurationValidator: Comprehensive validation with detailed error reporting
"""
import argparse
import logging
from pathlib import Path
from typing import List, Optional, Dict, Any, Union
from dataclasses import replace

from .models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig, IOConfig,
    WorkerConfig, CalculationRequest, ValidationResult,
    AlgorithmType, ExecutionMode
)

logger = logging.getLogger(__name__)


class ConfigurationBuilder:
    """Builder pattern implementation for PyMaSC configuration.

    Provides a fluent interface for constructing complex configuration
    objects with validation at each step. Enables method chaining and
    provides clear error messages for invalid configurations.

    Example:
        config = (ConfigurationBuilder()
                  .with_algorithm('bitarray')
                  .with_max_shift(500)
                  .with_mappability('/path/to/bigwig')
                  .with_multiprocessing(workers=4)
                  .build())
    """

    def __init__(self):
        """Initialize builder with default values."""
        self._calculation_config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=300,
            mapq_criteria=20
        )
        self._mappability_config = MappabilityConfig()
        self._execution_config = ExecutionConfig()
        self._io_config = None
        self._errors = []
        self._warnings = []

    def with_algorithm(self, algorithm: Union[str, AlgorithmType]) -> 'ConfigurationBuilder':
        """Set the cross-correlation algorithm.

        Args:
            algorithm: Algorithm name ('ncc', 'mscc', 'bitarray', 'successive')
                      or AlgorithmType enum value

        Returns:
            Self for method chaining

        Raises:
            ValueError: If algorithm is not supported
        """
        if isinstance(algorithm, str):
            try:
                algorithm = AlgorithmType(algorithm.lower())
            except ValueError:
                valid_algorithms = [alg.value for alg in AlgorithmType]
                raise ValueError(f"Unsupported algorithm '{algorithm}'. "
                                 f"Valid options: {valid_algorithms}")

        self._calculation_config = replace(self._calculation_config, algorithm=algorithm)
        return self

    def with_max_shift(self, max_shift: int) -> 'ConfigurationBuilder':
        """Set maximum shift distance for correlation calculation.

        Args:
            max_shift: Maximum shift distance in base pairs

        Returns:
            Self for method chaining
        """
        if max_shift <= 0:
            raise ValueError("max_shift must be positive")

        self._calculation_config = replace(self._calculation_config, max_shift=max_shift)
        return self

    def with_mapq_criteria(self, mapq: int) -> 'ConfigurationBuilder':
        """Set minimum mapping quality threshold.

        Args:
            mapq: Minimum mapping quality (0-255)

        Returns:
            Self for method chaining
        """
        if mapq < 0 or mapq > 255:
            raise ValueError("mapq_criteria must be between 0 and 255")

        self._calculation_config = replace(self._calculation_config, mapq_criteria=mapq)
        return self

    def with_read_length(self, read_length: Optional[int]) -> 'ConfigurationBuilder':
        """Set read length (None for auto-detection).

        Args:
            read_length: Read length in base pairs, or None for auto-detection

        Returns:
            Self for method chaining
        """
        if read_length is not None and read_length <= 0:
            raise ValueError("read_length must be positive if specified")

        self._calculation_config = replace(self._calculation_config, read_length=read_length)
        return self

    def with_references(self, references: List[str], lengths: List[int]) -> 'ConfigurationBuilder':
        """Set chromosome references and lengths.

        Args:
            references: List of chromosome names
            lengths: List of corresponding chromosome lengths

        Returns:
            Self for method chaining
        """
        if len(references) != len(lengths):
            raise ValueError("references and lengths must have same length")
        if any(length <= 0 for length in lengths):
            raise ValueError("all chromosome lengths must be positive")

        self._calculation_config = replace(
            self._calculation_config,
            references=references.copy(),
            lengths=lengths.copy()
        )
        return self

    def with_mappability(self,
                         bigwig_path: Union[str, Path],
                         stats_path: Optional[Union[str, Path]] = None) -> 'ConfigurationBuilder':
        """Configure mappability correction.

        Args:
            bigwig_path: Path to BigWig mappability file
            stats_path: Optional path to pre-computed mappability statistics

        Returns:
            Self for method chaining
        """
        bigwig_path = Path(bigwig_path)
        if not bigwig_path.exists():
            self._warnings.append(f"Mappability file does not exist: {bigwig_path}")

        stats_path = Path(stats_path) if stats_path else None
        if stats_path and not stats_path.exists():
            self._warnings.append(f"Mappability stats file does not exist: {stats_path}")

        self._mappability_config = MappabilityConfig(
            mappability_path=bigwig_path,
            mappability_stats_path=stats_path
        )
        return self

    def with_multiprocessing(self, workers: int) -> 'ConfigurationBuilder':
        """Configure multiprocessing execution.

        Args:
            workers: Number of worker processes

        Returns:
            Self for method chaining
        """
        if workers <= 0:
            raise ValueError("worker count must be positive")

        mode = ExecutionMode.MULTI_PROCESS if workers > 1 else ExecutionMode.SINGLE_PROCESS
        self._execution_config = ExecutionConfig(
            mode=mode,
            worker_count=workers
        )
        return self

    def with_input_output(self,
                          input_paths: List[Union[str, Path]],
                          output_dir: Union[str, Path],
                          output_names: Optional[List[str]] = None) -> 'ConfigurationBuilder':
        """Configure input and output settings.

        Args:
            input_paths: List of input BAM file paths
            output_dir: Output directory path
            output_names: Optional custom output names

        Returns:
            Self for method chaining
        """
        input_paths = [Path(p) for p in input_paths]
        output_dir = Path(output_dir)

        # Validate input files exist
        for path in input_paths:
            if not path.exists():
                self._errors.append(f"Input file does not exist: {path}")

        # Validate output names if provided
        if output_names and len(output_names) != len(input_paths):
            self._errors.append("output_names length must match input_paths length")

        self._io_config = IOConfig(
            input_paths=input_paths,
            output_dir=output_dir,
            output_names=output_names
        )
        return self

    def skip_ncc(self, skip: bool = True) -> 'ConfigurationBuilder':
        """Configure whether to skip naive cross-correlation calculation.

        Args:
            skip: Whether to skip NCC calculation

        Returns:
            Self for method chaining
        """
        self._calculation_config = replace(self._calculation_config, skip_ncc=skip)
        return self

    def build(self) -> CalculationRequest:
        """Build and validate the complete configuration.

        Returns:
            Complete CalculationRequest object

        Raises:
            ValueError: If configuration is invalid
        """
        if self._io_config is None:
            raise ValueError("Input/output configuration is required")

        # Create request object
        request = CalculationRequest(
            calculation_config=self._calculation_config,
            mappability_config=self._mappability_config,
            execution_config=self._execution_config,
            io_config=self._io_config
        )

        # Validate the complete configuration
        validation_errors = request.validate()
        all_errors = self._errors + validation_errors

        if all_errors:
            raise ValueError("Configuration validation failed:\n" +
                             "\n".join(f"  - {error}" for error in all_errors))

        # Log warnings if any
        for warning in self._warnings:
            logger.warning(warning)

        return request


class ConfigurationService:
    """Service for configuration management and validation.

    Provides centralized configuration validation, transformation between
    different configuration formats, and integration with existing
    argparse-based systems.
    """

    @staticmethod
    def from_argparse(args: argparse.Namespace) -> CalculationRequest:
        """Convert argparse Namespace to structured configuration.

        Args:
            args: Parsed command-line arguments from argparse

        Returns:
            CalculationRequest object

        Raises:
            ValueError: If argument conversion fails
        """
        builder = ConfigurationBuilder()

        # Algorithm selection
        if hasattr(args, 'successive') and args.successive:
            builder.with_algorithm(AlgorithmType.SUCCESSIVE)
        else:
            builder.with_algorithm(AlgorithmType.BITARRAY)

        # Core calculation parameters
        if hasattr(args, 'max_shift'):
            builder.with_max_shift(args.max_shift)
        if hasattr(args, 'mapq'):
            builder.with_mapq_criteria(args.mapq)
        if hasattr(args, 'read_length'):
            builder.with_read_length(args.read_length)
        if hasattr(args, 'skip_ncc') and args.skip_ncc:
            builder.skip_ncc()

        # Mappability configuration
        if hasattr(args, 'mappability') and args.mappability:
            stats_path = getattr(args, 'mappability_stats', None)
            builder.with_mappability(args.mappability, stats_path)

        # Execution configuration
        worker_count = getattr(args, 'process', 1)
        builder.with_multiprocessing(worker_count)

        # I/O configuration
        input_paths = args.reads if hasattr(args, 'reads') else []
        output_dir = getattr(args, 'outdir', '.')
        output_names = getattr(args, 'name', None)

        builder.with_input_output(input_paths, output_dir, output_names)

        return builder.build()

    @staticmethod
    def validate_configuration(config: CalculationRequest) -> ValidationResult:
        """Perform comprehensive validation of configuration.

        Args:
            config: Configuration to validate

        Returns:
            ValidationResult with errors and warnings
        """
        errors = []
        warnings = []

        # Validate algorithm-specific requirements
        calc_config = config.calculation_config
        map_config = config.mappability_config

        if calc_config.algorithm in [AlgorithmType.MSCC, AlgorithmType.BITARRAY]:
            if not map_config.is_enabled():
                errors.append(f"{calc_config.algorithm.value} algorithm requires mappability data")

        if calc_config.skip_ncc and not map_config.is_enabled():
            errors.append("Cannot skip NCC without mappability data for MSCC")

        # Validate execution compatibility
        exec_config = config.execution_config
        if exec_config.mode == ExecutionMode.MULTI_PROCESS and exec_config.worker_count == 1:
            warnings.append("Multi-process mode with 1 worker, consider single-process mode")

        # Validate file accessibility
        io_config = config.io_config
        for input_path in io_config.input_paths:
            if not input_path.exists():
                errors.append(f"Input file not found: {input_path}")
            elif input_path.suffix.lower() not in ['.bam', '.sam']:
                warnings.append(f"Input file may not be BAM/SAM: {input_path}")

        if not io_config.output_dir.exists():
            try:
                io_config.output_dir.mkdir(parents=True)
            except OSError as e:
                errors.append(f"Cannot create output directory: {e}")

        return ValidationResult(
            is_valid=len(errors) == 0,
            errors=errors,
            warnings=warnings
        )

    @staticmethod
    def create_worker_config(request: CalculationRequest,
                             logger_lock: Optional[Any] = None) -> WorkerConfig:
        """Create worker configuration from calculation request.

        Args:
            request: Main calculation request
            logger_lock: Optional thread lock for logging

        Returns:
            WorkerConfig object for worker processes
        """
        return WorkerConfig(
            calculation_config=request.calculation_config,
            mappability_config=request.mappability_config,
            progress_enabled=True,
            logger_lock=logger_lock
        )