"""Service adapter for integrating service layer with existing handlers.

This module provides adapter classes that allow existing PyMaSC handlers
to use the new service layer architecture. This enables gradual migration
to the service-oriented architecture while maintaining backward compatibility.

Key features:
- Seamless integration with existing code
- Service-based execution with handler interface
- Progressive migration path
- Full backward compatibility
"""
import logging
from typing import Optional, Dict, List, Any

from PyMaSC.handler.base import BaseCalcHandler
from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig, 
    AlgorithmType, WorkerConfig
)
from PyMaSC.services.workflow import (
    WorkflowService, WorkflowRequest, WorkflowResult,
    create_workflow_service
)
from PyMaSC.services.calculation import create_calculation_service
from PyMaSC.services.io import create_io_service

logger = logging.getLogger(__name__)


class ServiceBasedCalcHandler(BaseCalcHandler):
    """Calculation handler that uses service layer internally.

    This adapter allows existing code to use the new service-oriented
    architecture while maintaining the same handler interface. It
    demonstrates how the service layer can be integrated progressively.

    The handler delegates all actual work to the service layer while
    maintaining compatibility with existing PyMaSC code.
    """

    def __init__(self, 
                 path: str,
                 esttype: str,
                 max_shift: int,
                 mapq_criteria: int,
                 nworker: int = 1,
                 skip_ncc: bool = False,
                 chromfilter: Optional[List[str]] = None,
                 algorithm: str = 'bitarray'):
        """Initialize service-based handler.

        Args:
            path: Path to BAM file
            esttype: Read length estimation type
            max_shift: Maximum shift distance
            mapq_criteria: Minimum mapping quality
            nworker: Number of workers
            skip_ncc: Whether to skip NCC calculation
            chromfilter: Chromosome filter list
            algorithm: Algorithm to use
        """
        # Store parameters before parent init
        self.esttype = esttype
        self.algorithm = AlgorithmType(algorithm)
        self._chromfilter = chromfilter

        # Initialize parent (handles BAM file validation)
        super().__init__(
            path=path,
            mapq_criteria=mapq_criteria,
            max_shift=max_shift,
            nworker=nworker,
            mappability_handler=None,
            skip_ncc=skip_ncc,
            chroms=chromfilter
        )

        # Create services
        self._workflow_service = None
        self._calculation_service = None
        self._io_service = None

        # Workflow result storage
        self._workflow_result = None

    def set_mappability_handler(self, handler: Any) -> None:
        """Set mappability handler for compatibility.

        Args:
            handler: Mappability handler instance
        """
        self.mappability_handler = handler

    def _run_singleprocess_calculation(self) -> None:
        """Run calculation using service layer in single process mode."""
        self._execute_workflow(multiprocess=False)

    def _run_multiprocess_calculation(self) -> None:
        """Run calculation using service layer in multiprocess mode."""
        self._execute_workflow(multiprocess=True)

    def _execute_workflow(self, multiprocess: bool) -> None:
        """Execute workflow using service layer.

        Args:
            multiprocess: Whether to use multiprocessing
        """
        # Create calculation configuration
        calc_config = CalculationConfig(
            algorithm=self.algorithm,
            max_shift=self.max_shift,
            mapq_criteria=self.mapq_criteria,
            references=self.references,
            lengths=self.lengths,
            read_length=self.read_len,
            skip_ncc=self.skip_ncc
        )

        # Create mappability configuration if handler exists
        map_config = None
        if self.mappability_handler:
            map_config = MappabilityConfig(
                mappability_path=self.mappability_handler.path,
                read_len=self.read_len or 50
            )

        # Create execution configuration
        exec_config = ExecutionConfig(
            multiprocess=multiprocess,
            n_workers=self.nworker if multiprocess else 1
        )

        # Create workflow request
        request = WorkflowRequest(
            bam_path=self.path,
            output_prefix="",  # Not writing files in handler mode
            calculation_config=calc_config,
            mappability_config=map_config,
            execution_config=exec_config,
            output_formats=[],  # No output in handler mode
            chromosome_filter=self._chromfilter
        )

        # Create services if not already created
        if self._workflow_service is None:
            self._calculation_service = create_calculation_service(
                parallel=multiprocess,
                n_workers=self.nworker
            )
            self._io_service = create_io_service()
            self._workflow_service = create_workflow_service(
                parallel=multiprocess,
                n_workers=self.nworker,
                calculation_service=self._calculation_service,
                io_service=self._io_service
            )

        # Execute workflow
        self._workflow_result = self._workflow_service.execute(request)

        # Check for errors
        if not self._workflow_result.is_successful:
            raise RuntimeError(f"Workflow failed: {self._workflow_result.error}")

        # Extract results to handler attributes
        self._extract_results()

    def _extract_results(self) -> None:
        """Extract workflow results to handler attributes."""
        if not self._workflow_result or not self._workflow_result.calculation_result:
            return

        result = self._workflow_result.calculation_result

        # Clear existing results
        self.ref2forward_sum.clear()
        self.ref2reverse_sum.clear()
        self.ref2ccbins.clear()
        self.mappable_ref2forward_sum.clear()
        self.mappable_ref2reverse_sum.clear()
        self.mappable_ref2ccbins.clear()
        self.ref2mappable_len.clear()

        # Extract per-chromosome results
        for chrom, chrom_result in result.chromosome_results.items():
            # NCC results
            if not self.skip_ncc:
                self.ref2forward_sum[chrom] = chrom_result.forward_count
                self.ref2reverse_sum[chrom] = chrom_result.reverse_count
                self.ref2ccbins[chrom] = chrom_result.correlation_bins

            # MSCC results
            if chrom_result.mappable_forward_count is not None:
                self.mappable_ref2forward_sum[chrom] = chrom_result.mappable_forward_count
                self.mappable_ref2reverse_sum[chrom] = chrom_result.mappable_reverse_count
                self.mappable_ref2ccbins[chrom] = chrom_result.mappable_correlation_bins

            # Mappable length
            if chrom_result.mappable_length is not None:
                self.ref2mappable_len[chrom] = chrom_result.mappable_length

        # Update mappability handler if present
        if self.mappability_handler and self.ref2mappable_len:
            self.mappability_handler.chrom2mappable_len.update(self.ref2mappable_len)
            self.mappability_handler.mappable_len = list(
                self.ref2mappable_len.values()
            )

    def get_workflow_result(self) -> Optional[WorkflowResult]:
        """Get the complete workflow result.

        Returns:
            Workflow result if available
        """
        return self._workflow_result

    def cleanup(self) -> None:
        """Clean up resources."""
        super()._cleanup()

        # Close services
        if self._io_service and hasattr(self._io_service, 'close'):
            self._io_service.close()


def create_service_based_handler(path: str,
                               esttype: str = 'mean',
                               max_shift: int = 500,
                               mapq_criteria: int = 30,
                               nworker: int = 1,
                               skip_ncc: bool = False,
                               chromfilter: Optional[List[str]] = None,
                               algorithm: str = 'bitarray',
                               mappability_path: Optional[str] = None,
                               read_len: Optional[int] = None) -> ServiceBasedCalcHandler:
    """Create a service-based calculation handler.

    Factory function that creates a handler using the service layer
    architecture with all necessary configuration.

    Args:
        path: Path to BAM file
        esttype: Read length estimation type
        max_shift: Maximum shift distance
        mapq_criteria: Minimum mapping quality threshold
        nworker: Number of worker processes
        skip_ncc: Whether to skip NCC calculation
        chromfilter: Chromosome filter list
        algorithm: Algorithm to use
        mappability_path: Optional mappability file path
        read_len: Optional read length

    Returns:
        Configured service-based handler
    """
    # Create handler
    handler = ServiceBasedCalcHandler(
        path=path,
        esttype=esttype,
        max_shift=max_shift,
        mapq_criteria=mapq_criteria,
        nworker=nworker,
        skip_ncc=skip_ncc,
        chromfilter=chromfilter,
        algorithm=algorithm
    )

    # Set read length if provided
    if read_len:
        handler.set_readlen(read_len)

    # Set mappability if provided
    if mappability_path:
        from PyMaSC.handler.mappability import MappabilityHandler
        mappability_handler = MappabilityHandler(
            path=mappability_path,
            max_shift=max_shift,
            read_len=read_len or 50,
            n_worker=nworker
        )
        handler.set_mappability_handler(mappability_handler)

    return handler