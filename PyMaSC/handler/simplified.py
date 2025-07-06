"""Simplified calculation handler using dependency injection.

This module provides SimplifiedCalcHandler that demonstrates
how handlers can be simplified through dependency injection
and delegation to services. This is the target architecture
for all PyMaSC handlers.

Key features:
- Constructor injection of services
- Minimal initialization logic
- All operations delegated to services
- Clean separation of concerns
"""
import logging
from dataclasses import dataclass
from typing import Optional, List, Dict, Any
from pathlib import Path

from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig,
    AlgorithmType, ExecutionMode
)
from PyMaSC.services.calculation import (
    CalculationService, GenomeWideResult, create_calculation_service
)
from PyMaSC.services.io import IOService, create_io_service
from PyMaSC.services.validation import ValidationService, create_validation_service
from PyMaSC.services.mappability import MappabilityService, create_mappability_service
from PyMaSC.services.workflow import (
    WorkflowService, WorkflowRequest, create_workflow_service
)

logger = logging.getLogger(__name__)


@dataclass
class HandlerDependencies:
    """Container for handler dependencies.

    Groups all services that a handler might need,
    enabling easy dependency injection.
    """
    calculation_service: Optional[CalculationService] = None
    io_service: Optional[IOService] = None
    validation_service: Optional[ValidationService] = None
    mappability_service: Optional[MappabilityService] = None
    workflow_service: Optional[WorkflowService] = None

    def ensure_services(self) -> None:
        """Ensure all services are initialized with defaults."""
        if self.calculation_service is None:
            self.calculation_service = create_calculation_service()
        if self.io_service is None:
            self.io_service = create_io_service()
        if self.validation_service is None:
            self.validation_service = create_validation_service()
        if self.workflow_service is None:
            self.workflow_service = create_workflow_service(
                calculation_service=self.calculation_service,
                io_service=self.io_service
            )


class SimplifiedCalcHandler:
    """Simplified calculation handler using dependency injection.

    This handler demonstrates the target architecture where:
    - All dependencies are injected
    - Initialization is minimal
    - All operations are delegated to services
    - Handler acts as a thin coordination layer
    """

    def __init__(self,
                 bam_path: str,
                 config: CalculationConfig,
                 dependencies: Optional[HandlerDependencies] = None):
        """Initialize handler with injected dependencies.

        Args:
            bam_path: Path to BAM file
            config: Calculation configuration
            dependencies: Injected service dependencies
        """
        self.bam_path = bam_path
        self.config = config

        # Use provided dependencies or create defaults
        self.deps = dependencies or HandlerDependencies()
        self.deps.ensure_services()

        # Handler state (minimal)
        self._result: Optional[GenomeWideResult] = None
        self._validated = False

    def validate(self) -> bool:
        """Validate inputs using validation service.

        Returns:
            True if validation passed
        """
        if self.deps.validation_service is None:
            raise RuntimeError("Validation service not initialized")
        result = self.deps.validation_service.validate_workflow_request(
            bam_path=self.bam_path,
            output_prefix="",  # Not writing files in handler mode
            calculation_config=self.config,
            mappability_config=self._get_mappability_config(),
            execution_config=self._get_execution_config()
        )

        # Log errors and warnings
        for error in result.errors:
            logger.error(f"Validation error: {error}")
        for warning in result.warnings:
            logger.warning(f"Validation warning: {warning}")

        self._validated = result.is_valid
        return result.is_valid

    def run_calculation(self) -> GenomeWideResult:
        """Run calculation using workflow service.

        Returns:
            Calculation result
        """
        if not self._validated:
            if not self.validate():
                raise ValueError("Validation failed")

        # Create workflow request
        request = WorkflowRequest(
            bam_path=self.bam_path,
            output_prefix="",  # Not writing files
            calculation_config=self.config,
            mappability_config=self._get_mappability_config(),
            execution_config=self._get_execution_config(),
            output_formats=[]  # No output in handler mode
        )

        # Execute workflow
        if self.deps.workflow_service is None:
            raise RuntimeError("Workflow service not initialized")
        workflow_result = self.deps.workflow_service.execute(request)

        if not workflow_result.is_successful:
            raise RuntimeError(f"Calculation failed: {workflow_result.error}")

        if workflow_result.calculation_result is None:
            raise RuntimeError("Workflow result has no calculation data")
        self._result = workflow_result.calculation_result
        return self._result

    def get_result(self) -> Optional[GenomeWideResult]:
        """Get calculation result if available.

        Returns:
            Calculation result or None
        """
        return self._result

    def set_mappability_file(self, path: str, read_len: int = 50) -> None:
        """Set mappability file and create service.

        Args:
            path: Path to mappability file
            read_len: Read length for mappability
        """
        # Create mappability configuration
        map_config = MappabilityConfig(
            mappability_path=Path(path),
            read_len=read_len
        )

        # Update configuration
        if self.config.algorithm == AlgorithmType.NAIVE_CC:
            self.config.algorithm = AlgorithmType.MSCC

        # Create mappability service
        self.deps.mappability_service = create_mappability_service(path, map_config)

    def set_multiprocess(self, n_workers: int) -> None:
        """Enable multiprocessing with specified workers.

        Args:
            n_workers: Number of worker processes
        """
        # Update execution mode in config
        if not hasattr(self, '_exec_config'):
            self._exec_config = ExecutionConfig()

        self._exec_config.mode = ExecutionMode.MULTI_PROCESS
        self._exec_config.worker_count = n_workers

        # Recreate services with parallel support
        self.deps.calculation_service = create_calculation_service(
            parallel=True,
            n_workers=n_workers
        )
        self.deps.workflow_service = create_workflow_service(
            parallel=True,
            n_workers=n_workers,
            calculation_service=self.deps.calculation_service,
            io_service=self.deps.io_service
        )

    def _get_mappability_config(self) -> Optional[MappabilityConfig]:
        """Get mappability configuration if available."""
        if self.deps.mappability_service and hasattr(self.deps.mappability_service, 'config'):
            return self.deps.mappability_service.config  # type: ignore[no-any-return]
        return None

    def _get_execution_config(self) -> ExecutionConfig:
        """Get execution configuration."""
        if hasattr(self, '_exec_config'):
            return self._exec_config
        return ExecutionConfig()

    def cleanup(self) -> None:
        """Clean up resources."""
        # Close mappability service if open
        if self.deps.mappability_service:
            if hasattr(self.deps.mappability_service, 'close'):
                self.deps.mappability_service.close()

        # Close I/O service if open
        if self.deps.io_service:
            if hasattr(self.deps.io_service, 'close'):
                self.deps.io_service.close()


class HandlerBuilder:
    """Builder for creating handlers with complex configurations.

    Provides a fluent interface for constructing handlers
    with all necessary services and configurations.
    """

    def __init__(self) -> None:
        """Initialize builder."""
        self._bam_path: Optional[str] = None
        self._algorithm = AlgorithmType.BITARRAY
        self._max_shift = 500
        self._mapq_criteria = 30
        self._read_length: Optional[int] = None
        self._skip_ncc = False
        self._mappability_path: Optional[str] = None
        self._mappability_read_len = 50
        self._n_workers = 1
        self._dependencies = HandlerDependencies()

    def with_bam_file(self, path: str) -> 'HandlerBuilder':
        """Set BAM file path."""
        self._bam_path = path
        return self

    def with_algorithm(self, algorithm: str) -> 'HandlerBuilder':
        """Set algorithm type."""
        self._algorithm = AlgorithmType(algorithm)
        return self

    def with_max_shift(self, max_shift: int) -> 'HandlerBuilder':
        """Set maximum shift distance."""
        self._max_shift = max_shift
        return self

    def with_mapq_threshold(self, threshold: int) -> 'HandlerBuilder':
        """Set mapping quality threshold."""
        self._mapq_criteria = threshold
        return self

    def with_read_length(self, length: int) -> 'HandlerBuilder':
        """Set read length."""
        self._read_length = length
        return self

    def with_mappability(self, path: str, read_len: int = 50) -> 'HandlerBuilder':
        """Set mappability file."""
        self._mappability_path = path
        self._mappability_read_len = read_len

        # Adjust algorithm if needed
        if self._algorithm == AlgorithmType.NAIVE_CC:
            self._algorithm = AlgorithmType.MSCC

        return self

    def with_multiprocessing(self, n_workers: int) -> 'HandlerBuilder':
        """Enable multiprocessing."""
        self._n_workers = n_workers
        return self

    def with_validation_service(self, service: ValidationService) -> 'HandlerBuilder':
        """Set custom validation service."""
        self._dependencies.validation_service = service
        return self

    def with_io_service(self, service: IOService) -> 'HandlerBuilder':
        """Set custom I/O service."""
        self._dependencies.io_service = service
        return self

    def with_skip_ncc(self, skip: bool = True) -> 'HandlerBuilder':
        """Set whether to skip NCC calculation."""
        self._skip_ncc = skip
        return self

    def build(self) -> SimplifiedCalcHandler:
        """Build the handler with all configurations.

        Returns:
            Configured handler

        Raises:
            ValueError: If required parameters are missing
        """
        if not self._bam_path:
            raise ValueError("BAM file path is required")

        # Create services based on configuration
        parallel = self._n_workers > 1

        # Create calculation service
        self._dependencies.calculation_service = create_calculation_service(
            parallel=parallel,
            n_workers=self._n_workers
        )

        # Create I/O service if not provided
        if self._dependencies.io_service is None:
            self._dependencies.io_service = create_io_service()

        # Create validation service if not provided
        if self._dependencies.validation_service is None:
            self._dependencies.validation_service = create_validation_service()

        # Create mappability service if needed
        if self._mappability_path:
            map_config = MappabilityConfig(
                mappability_path=Path(self._mappability_path),
                read_len=self._mappability_read_len
            )
            self._dependencies.mappability_service = create_mappability_service(
                self._mappability_path,
                map_config
            )

        # Create workflow service
        self._dependencies.workflow_service = create_workflow_service(
            parallel=parallel,
            n_workers=self._n_workers,
            calculation_service=self._dependencies.calculation_service,
            io_service=self._dependencies.io_service
        )

        # Get BAM info for configuration
        bam_info = self._dependencies.io_service.get_bam_info(self._bam_path)

        # Create calculation configuration
        calc_config = CalculationConfig(
            algorithm=self._algorithm,
            max_shift=self._max_shift,
            mapq_criteria=self._mapq_criteria,
            references=bam_info.references,
            lengths=bam_info.lengths,
            read_length=self._read_length,
            skip_ncc=self._skip_ncc
        )

        # Create handler
        handler = SimplifiedCalcHandler(
            bam_path=self._bam_path,
            config=calc_config,
            dependencies=self._dependencies
        )

        # Set multiprocessing if needed
        if parallel:
            handler.set_multiprocess(self._n_workers)

        return handler


# Convenience function for simple cases
def create_simplified_handler(bam_path: str,
                            algorithm: str = 'bitarray',
                            max_shift: int = 500,
                            mapq_criteria: int = 30,
                            mappability_path: Optional[str] = None,
                            n_workers: int = 1) -> SimplifiedCalcHandler:
    """Create a simplified handler with common options.

    Args:
        bam_path: Path to BAM file
        algorithm: Algorithm to use
        max_shift: Maximum shift distance
        mapq_criteria: Mapping quality threshold
        mappability_path: Optional mappability file
        n_workers: Number of workers for multiprocessing

    Returns:
        Configured handler
    """
    builder = HandlerBuilder() \
        .with_bam_file(bam_path) \
        .with_algorithm(algorithm) \
        .with_max_shift(max_shift) \
        .with_mapq_threshold(mapq_criteria) \
        .with_multiprocessing(n_workers)

    if mappability_path:
        builder.with_mappability(mappability_path)

    return builder.build()