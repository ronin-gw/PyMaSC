"""Workflow service for orchestrating PyMaSC calculations.

This module provides the WorkflowService that coordinates the overall
calculation process, managing the interaction between calculation and
I/O services. This is the application layer of the architecture.

Key features:
- Process orchestration
- Error handling and recovery
- Progress tracking
- Resource management
- Workflow customization
"""
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Dict, Any, Callable
from concurrent.futures import ProcessPoolExecutor, as_completed

from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, ExecutionConfig, AlgorithmType,
    ExecutionMode
)
from PyMaSC.core.observer import ProgressEvent, ProgressEventType, ProgressSubject
from PyMaSC.services.calculation import (
    CalculationService, ChromosomeData, GenomeWideResult,
    create_calculation_service
)
from PyMaSC.services.io import IOService, create_io_service
from PyMaSC.utils.calc import filter_chroms

logger = logging.getLogger(__name__)


class WorkflowStatus(Enum):
    """Status of workflow execution."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


@dataclass
class WorkflowRequest:
    """Request for workflow execution.
    
    Contains all parameters needed to execute a complete
    cross-correlation analysis workflow.
    """
    bam_path: str
    output_prefix: str
    calculation_config: CalculationConfig
    mappability_config: Optional[MappabilityConfig] = None
    execution_config: Optional[ExecutionConfig] = None
    output_formats: List[str] = None
    chromosome_filter: Optional[List[str]] = None
    
    def __post_init__(self):
        """Set defaults after initialization."""
        if self.execution_config is None:
            self.execution_config = ExecutionConfig()
        if self.output_formats is None:
            self.output_formats = ['table', 'stats']


@dataclass
class WorkflowResult:
    """Result of workflow execution.
    
    Contains the calculation results and execution metadata.
    """
    status: WorkflowStatus
    calculation_result: Optional[GenomeWideResult] = None
    output_paths: Optional[Dict[str, str]] = None
    error: Optional[str] = None
    execution_time: Optional[float] = None
    
    @property
    def is_successful(self) -> bool:
        """Check if workflow completed successfully."""
        return self.status == WorkflowStatus.COMPLETED


class WorkflowService(ABC):
    """Abstract base for workflow services.
    
    Defines the interface for workflow orchestration services.
    """
    
    @abstractmethod
    def execute(self, request: WorkflowRequest) -> WorkflowResult:
        """Execute a complete workflow.
        
        Args:
            request: Workflow execution request
            
        Returns:
            Workflow execution result
        """
        pass
    
    @abstractmethod
    def validate_request(self, request: WorkflowRequest) -> List[str]:
        """Validate a workflow request.
        
        Args:
            request: Request to validate
            
        Returns:
            List of validation errors (empty if valid)
        """
        pass


class StandardWorkflowService(WorkflowService, ProgressSubject):
    """Standard implementation of workflow service.
    
    Coordinates calculation and I/O services to execute the
    complete cross-correlation analysis workflow.
    """
    
    def __init__(self,
                 calculation_service: Optional[CalculationService] = None,
                 io_service: Optional[IOService] = None):
        """Initialize workflow service.
        
        Args:
            calculation_service: Service for calculations
            io_service: Service for I/O operations
        """
        # Initialize both parent classes explicitly
        WorkflowService.__init__(self)
        ProgressSubject.__init__(self)
        
        self._calc_service = calculation_service
        self._io_service = io_service
        self._status = WorkflowStatus.PENDING
    
    def execute(self, request: WorkflowRequest) -> WorkflowResult:
        """Execute a complete workflow."""
        import time
        start_time = time.time()
        
        # Update status
        self._status = WorkflowStatus.RUNNING
        self._notify_progress("workflow_started", 0, 100)
        
        try:
            # Validate request
            errors = self.validate_request(request)
            if errors:
                raise ValueError(f"Invalid request: {', '.join(errors)}")
            
            # Create services if not provided
            if self._calc_service is None:
                self._calc_service = create_calculation_service(
                    parallel=request.execution_config.mode == ExecutionMode.MULTI_PROCESS,
                    n_workers=request.execution_config.worker_count
                )
            
            if self._io_service is None:
                self._io_service = create_io_service()
            
            # Execute workflow steps
            self._notify_progress("reading_bam_info", 10, 100)
            bam_info = self._io_service.get_bam_info(request.bam_path)
            
            # Apply chromosome filtering
            target_chromosomes = self._filter_chromosomes(
                bam_info.references,
                request.chromosome_filter
            )
            
            if not target_chromosomes:
                raise ValueError("No chromosomes match the filter criteria")
            
            # Update configuration with actual references
            request.calculation_config.references = target_chromosomes
            request.calculation_config.lengths = [
                bam_info.lengths[bam_info.references.index(chrom)]
                for chrom in target_chromosomes
            ]
            
            # Read chromosome data
            self._notify_progress("reading_chromosomes", 20, 100)
            chromosome_data = self._read_chromosome_data(
                request, target_chromosomes, bam_info
            )
            
            # Perform calculation
            self._notify_progress("calculating", 50, 100)
            calculation_result = self._calc_service.calculate_genome_wide(
                chromosome_data,
                request.calculation_config,
                request.mappability_config
            )
            
            # Write results
            self._notify_progress("writing_results", 90, 100)
            output_paths = self._io_service.write_results(
                calculation_result,
                request.output_prefix,
                request.output_formats
            )
            
            # Complete workflow
            self._status = WorkflowStatus.COMPLETED
            self._notify_progress("workflow_completed", 100, 100)
            
            execution_time = time.time() - start_time
            
            return WorkflowResult(
                status=WorkflowStatus.COMPLETED,
                calculation_result=calculation_result,
                output_paths=output_paths,
                execution_time=execution_time
            )
            
        except Exception as e:
            # Handle errors
            logger.error(f"Workflow failed: {e}", exc_info=True)
            self._status = WorkflowStatus.FAILED
            self._notify_progress("workflow_failed", -1, -1)
            
            return WorkflowResult(
                status=WorkflowStatus.FAILED,
                error=str(e),
                execution_time=time.time() - start_time
            )
        
        finally:
            # Cleanup resources
            if hasattr(self._io_service, 'close'):
                self._io_service.close()
    
    def validate_request(self, request: WorkflowRequest) -> List[str]:
        """Validate a workflow request."""
        errors = []
        
        # Validate BAM file path
        if not request.bam_path:
            errors.append("BAM file path is required")
        
        # Validate output prefix (only required if output formats specified)
        if request.output_formats and not request.output_prefix:
            errors.append("Output prefix is required when output formats are specified")
        
        # Validate calculation config
        if request.calculation_config.max_shift <= 0:
            errors.append("Max shift must be positive")
        
        if request.calculation_config.mapq_criteria < 0:
            errors.append("MAPQ criteria must be non-negative")
        
        # Validate algorithm-specific requirements
        algorithm = request.calculation_config.algorithm
        if algorithm == AlgorithmType.MSCC:
            if not request.mappability_config or not request.mappability_config.is_enabled():
                errors.append(f"{algorithm.value} requires mappability configuration")
        
        # Validate output formats
        valid_formats = {'table', 'stats', 'figure'}
        for fmt in request.output_formats:
            if fmt not in valid_formats:
                errors.append(f"Invalid output format: {fmt}")
        
        return errors
    
    def _filter_chromosomes(self, 
                          all_chromosomes: List[str],
                          filter_list: Optional[List[str]]) -> List[str]:
        """Apply chromosome filtering."""
        if filter_list is None:
            return all_chromosomes
        
        return filter_chroms(all_chromosomes, filter_list)
    
    def _read_chromosome_data(self,
                            request: WorkflowRequest,
                            chromosomes: List[str],
                            bam_info: Any) -> List[ChromosomeData]:
        """Read data for all target chromosomes."""
        chromosome_data = []
        total_chroms = len(chromosomes)
        
        for i, chrom in enumerate(chromosomes):
            # Update progress
            progress = 20 + int(30 * (i / total_chroms))
            self._notify_progress(f"reading_chromosome_{chrom}", progress, 100)
            
            # Read chromosome data
            data = self._io_service.read_chromosome_data(
                request.bam_path,
                chrom,
                request.calculation_config.mapq_criteria
            )
            chromosome_data.append(data)
        
        return chromosome_data
    
    def _notify_progress(self, message: str, current: int, total: int) -> None:
        """Notify observers of progress."""
        # Determine event type based on message
        if "started" in message:
            event_type = ProgressEventType.STARTED
        elif "completed" in message or "finished" in message:
            event_type = ProgressEventType.COMPLETED
        elif "failed" in message or "error" in message:
            event_type = ProgressEventType.FAILED
        else:
            event_type = ProgressEventType.UPDATED
            
        event = ProgressEvent(
            event_type=event_type,
            source="workflow",
            message=message,
            current=current,
            total=total
        )
        self.notify(event)


class ParallelWorkflowService(StandardWorkflowService):
    """Parallel implementation of workflow service.
    
    Extends standard workflow with true parallel processing
    capabilities for chromosome-level parallelization.
    """
    
    def __init__(self,
                 n_workers: int = 4,
                 calculation_service: Optional[CalculationService] = None,
                 io_service: Optional[IOService] = None):
        """Initialize parallel workflow service.
        
        Args:
            n_workers: Number of parallel workers
            calculation_service: Service for calculations
            io_service: Service for I/O operations
        """
        super().__init__(calculation_service, io_service)
        self.n_workers = n_workers
    
    def _read_chromosome_data(self,
                            request: WorkflowRequest,
                            chromosomes: List[str],
                            bam_info: Any) -> List[ChromosomeData]:
        """Read chromosome data in parallel."""
        chromosome_data = []
        
        # Use process pool for parallel reading
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            # Submit all read tasks
            future_to_chrom = {
                executor.submit(
                    self._read_single_chromosome,
                    request.bam_path,
                    chrom,
                    request.calculation_config.mapq_criteria
                ): chrom
                for chrom in chromosomes
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_chrom):
                chrom = future_to_chrom[future]
                try:
                    data = future.result()
                    chromosome_data.append(data)
                    
                    # Update progress
                    completed = len(chromosome_data)
                    progress = 20 + int(30 * (completed / len(chromosomes)))
                    self._notify_progress(f"read_chromosome_{chrom}", progress, 100)
                    
                except Exception as e:
                    logger.error(f"Error reading chromosome {chrom}: {e}")
                    raise
        
        # Sort by original chromosome order
        chromosome_data.sort(key=lambda d: chromosomes.index(d.chromosome))
        
        return chromosome_data
    
    def _read_single_chromosome(self,
                              bam_path: str,
                              chromosome: str,
                              mapq_threshold: int) -> ChromosomeData:
        """Read single chromosome data (for parallel execution)."""
        # Create new I/O service for this process
        io_service = create_io_service()
        try:
            return io_service.read_chromosome_data(
                bam_path, chromosome, mapq_threshold
            )
        finally:
            if hasattr(io_service, 'close'):
                io_service.close()


# Factory function for workflow service creation
def create_workflow_service(parallel: bool = False,
                          n_workers: int = 1,
                          calculation_service: Optional[CalculationService] = None,
                          io_service: Optional[IOService] = None) -> WorkflowService:
    """Create appropriate workflow service.
    
    Args:
        parallel: Whether to use parallel processing
        n_workers: Number of worker processes
        calculation_service: Optional calculation service
        io_service: Optional I/O service
        
    Returns:
        Workflow service instance
    """
    if parallel and n_workers > 1:
        return ParallelWorkflowService(n_workers, calculation_service, io_service)
    else:
        return StandardWorkflowService(calculation_service, io_service)


# Convenience function for common use case
def execute_workflow(bam_path: str,
                    output_prefix: str,
                    algorithm: str = 'bitarray',
                    max_shift: int = 500,
                    mappability_path: Optional[str] = None,
                    n_workers: int = 1,
                    **kwargs) -> WorkflowResult:
    """Execute a workflow with common parameters.
    
    Convenience function that creates all necessary configurations
    and services to execute a complete workflow.
    
    Args:
        bam_path: Path to input BAM file
        output_prefix: Base path for output files
        algorithm: Algorithm to use ('bitarray', 'successive', 'naive')
        max_shift: Maximum shift for correlation
        mappability_path: Optional path to mappability BigWig
        n_workers: Number of parallel workers
        **kwargs: Additional configuration parameters
        
    Returns:
        Workflow execution result
    """
    # Create configurations
    calc_config = CalculationConfig(
        algorithm=AlgorithmType(algorithm),
        max_shift=max_shift,
        mapq_criteria=kwargs.get('mapq_criteria', 30),
        skip_ncc=kwargs.get('skip_ncc', False)
    )
    
    map_config = None
    if mappability_path:
        map_config = MappabilityConfig(
            mappability_path=mappability_path,
            read_len=kwargs.get('read_len', 50)
        )
    
    exec_config = ExecutionConfig(
        multiprocess=n_workers > 1,
        n_workers=n_workers
    )
    
    # Create request
    request = WorkflowRequest(
        bam_path=bam_path,
        output_prefix=output_prefix,
        calculation_config=calc_config,
        mappability_config=map_config,
        execution_config=exec_config,
        output_formats=kwargs.get('output_formats', ['table', 'stats']),
        chromosome_filter=kwargs.get('chromosome_filter')
    )
    
    # Create and execute workflow
    service = create_workflow_service(parallel=n_workers > 1, n_workers=n_workers)
    return service.execute(request)