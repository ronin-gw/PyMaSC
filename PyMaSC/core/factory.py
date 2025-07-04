"""Factory pattern implementations for PyMaSC calculators and workers.

This module provides factory classes that encapsulate the creation logic
for various calculator and worker implementations. The factories handle
the complexity of different constructor signatures and dependencies,
providing a unified interface for object creation.

Key factories:
- CalculatorFactory: Creates cross-correlation calculator instances
- WorkerFactory: Creates worker instances for multiprocessing
- CalculatorAdapter: Adapts existing calculators to new interface
"""
import logging
from typing import Optional, Any, Dict, List, Union
from multiprocessing import Queue, Lock

from .interfaces import CrossCorrelationCalculator, CalculatorFactory as ICalculatorFactory
from .models import (
    CalculationConfig, MappabilityConfig, WorkerConfig, 
    AlgorithmType, ExecutionMode
)

# Import existing calculator implementations
# These will be wrapped to conform to the new interface
from PyMaSC.core.ncc import NaiveCCCalculator as _NativeCCCalculator
from PyMaSC.core.mscc import MSCCCalculator as _NativeMSCCCalculator
from PyMaSC.bacore.mscc import CCBitArrayCalculator as _NativeBitArrayCalculator

# Import BigWig reader for mappability
from PyMaSC.reader.bigwig import BigWigReader

logger = logging.getLogger(__name__)


class CalculatorAdapter(CrossCorrelationCalculator):
    """Adapter to make existing calculators conform to new interface.
    
    This adapter wraps existing Cython calculator implementations to
    provide the standardized CrossCorrelationCalculator interface.
    It handles the differences in method signatures and ensures
    consistent behavior across all calculator types.
    """
    
    def __init__(self, native_calculator: Any):
        """Initialize adapter with native calculator instance.
        
        Args:
            native_calculator: Existing calculator implementation
        """
        self._calculator = native_calculator
        self._is_bitarray = hasattr(native_calculator, 'ref2mascbins')
    
    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a forward strand read."""
        self._calculator.feed_forward_read(chrom, pos, readlen)
    
    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a reverse strand read."""
        self._calculator.feed_reverse_read(chrom, pos, readlen)
    
    def finishup_calculation(self) -> None:
        """Complete cross-correlation calculation."""
        self._calculator.finishup_calculation()
    
    def flush(self) -> None:
        """Flush intermediate results."""
        if hasattr(self._calculator, 'flush'):
            self._calculator.flush()
        else:
            # For calculators without flush, finishup serves the same purpose
            self._calculator.finishup_calculation()
    
    @property
    def ref2forward_sum(self) -> Dict[str, int]:
        """Forward read counts by chromosome."""
        return self._calculator.ref2forward_sum
    
    @property
    def ref2reverse_sum(self) -> Dict[str, int]:
        """Reverse read counts by chromosome."""
        return self._calculator.ref2reverse_sum
    
    @property
    def ref2ccbins(self) -> Dict[str, Any]:
        """Cross-correlation bins by chromosome."""
        return self._calculator.ref2ccbins
    
    @property
    def ref2mappable_forward_sum(self) -> Optional[Dict[str, int]]:
        """Mappable forward read counts (BitArray only)."""
        if self._is_bitarray:
            return self._calculator.ref2mappable_forward_sum
        return None
    
    @property
    def ref2mappable_reverse_sum(self) -> Optional[Dict[str, int]]:
        """Mappable reverse read counts (BitArray only)."""
        if self._is_bitarray:
            return self._calculator.ref2mappable_reverse_sum
        return None
    
    @property
    def ref2mascbins(self) -> Optional[Dict[str, Any]]:
        """MSCC bins (BitArray only)."""
        if self._is_bitarray:
            return self._calculator.ref2mascbins
        return None


class CalculatorFactory:
    """Factory for creating cross-correlation calculator instances.
    
    This factory encapsulates the complex logic of creating different
    calculator types with their varying constructor signatures and
    dependencies. It provides a unified interface for calculator creation
    regardless of the underlying implementation.
    """
    
    @staticmethod
    def create_calculator(
        algorithm: AlgorithmType,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig] = None,
        logger_lock: Optional[Lock] = None,
        progress_hook: Optional[Any] = None
    ) -> CrossCorrelationCalculator:
        """Create a calculator instance for the specified algorithm.
        
        Args:
            algorithm: Algorithm type to create
            config: Calculation configuration
            mappability_config: Optional mappability configuration
            logger_lock: Optional thread lock for logging
            progress_hook: Optional progress reporting hook
            
        Returns:
            Calculator instance wrapped in adapter
            
        Raises:
            ValueError: If algorithm is not supported or configuration is invalid
        """
        # Validate algorithm-specific requirements
        if algorithm == AlgorithmType.MSCC:
            if not mappability_config or not mappability_config.is_enabled():
                raise ValueError(f"{algorithm.value} algorithm requires mappability configuration")
        
        # Create appropriate calculator based on algorithm
        if algorithm == AlgorithmType.NAIVE_CC or algorithm == AlgorithmType.SUCCESSIVE:
            return CalculatorFactory._create_ncc_calculator(config, logger_lock)
        
        elif algorithm == AlgorithmType.MSCC:
            return CalculatorFactory._create_mscc_calculator(
                config, mappability_config, logger_lock
            )
        
        elif algorithm == AlgorithmType.BITARRAY:
            return CalculatorFactory._create_bitarray_calculator(
                config, mappability_config, logger_lock, progress_hook
            )
        
        else:
            raise ValueError(f"Unsupported algorithm: {algorithm}")
    
    @staticmethod
    def _create_ncc_calculator(
        config: CalculationConfig,
        logger_lock: Optional[Lock] = None
    ) -> CrossCorrelationCalculator:
        """Create Naive CC calculator instance."""
        native_calc = _NativeCCCalculator(
            config.max_shift,
            config.references,
            config.lengths,
            logger_lock
        )
        return CalculatorAdapter(native_calc)
    
    @staticmethod
    def _create_mscc_calculator(
        config: CalculationConfig,
        mappability_config: MappabilityConfig,
        logger_lock: Optional[Lock] = None
    ) -> CrossCorrelationCalculator:
        """Create MSCC calculator instance."""
        # Create BigWig reader for mappability data
        bwfeeder = BigWigReader(str(mappability_config.mappability_path))
        
        # Use read_length from config or mappability config
        read_len = config.read_length
        if read_len is None:
            read_len = mappability_config.read_len
        if read_len is None:
            raise ValueError("Read length must be specified for MSCC algorithm")
        
        native_calc = _NativeMSCCCalculator(
            config.max_shift,
            read_len,
            config.references,
            config.lengths,
            bwfeeder,
            logger_lock
        )
        return CalculatorAdapter(native_calc)
    
    @staticmethod
    def _create_bitarray_calculator(
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock] = None,
        progress_hook: Optional[Any] = None
    ) -> CrossCorrelationCalculator:
        """Create BitArray calculator instance."""
        # Create BigWig reader if mappability is enabled
        bwfeeder = None
        if mappability_config and mappability_config.is_enabled():
            bwfeeder = BigWigReader(str(mappability_config.mappability_path))
        
        # Use read_length from config
        read_len = config.read_length or 50  # Default if not specified
        if mappability_config and mappability_config.read_len:
            read_len = mappability_config.read_len
        
        native_calc = _NativeBitArrayCalculator(
            config.max_shift,
            read_len,
            config.references,
            config.lengths,
            bwfeeder,
            config.skip_ncc,
            logger_lock,
            progress_hook
        )
        return CalculatorAdapter(native_calc)


class WorkerFactory:
    """Factory for creating worker instances for multiprocessing.
    
    This factory handles the complex logic of creating different worker
    types based on the calculation configuration, including the selection
    of appropriate worker classes and initialization of their dependencies.
    """
    
    @staticmethod
    def create_worker(
        worker_config: WorkerConfig,
        order_queue: Queue,
        report_queue: Queue,
        logger_lock: Lock,
        bam_path: str,
        use_unified: bool = True
    ) -> Any:  # Returns Process-based worker
        """Create a worker instance with the specified configuration.
        
        Args:
            worker_config: Worker configuration parameters
            order_queue: Queue for receiving work orders
            report_queue: Queue for reporting results
            logger_lock: Lock for thread-safe logging
            bam_path: Path to BAM file to process
            use_unified: Whether to use new unified worker (default True)
            
        Returns:
            Worker instance ready for multiprocessing execution
            
        Raises:
            ImportError: If required worker class is not available
        """
        # Use unified worker by default
        if use_unified:
            from PyMaSC.core.worker import UnifiedWorker
            return UnifiedWorker(
                worker_config,
                order_queue,
                report_queue,
                logger_lock,
                bam_path
            )
        
        # Legacy worker creation for backward compatibility
        calc_config = worker_config.calculation_config
        map_config = worker_config.mappability_config
        
        # Determine which worker class to use based on algorithm
        if calc_config.algorithm == AlgorithmType.BITARRAY:
            # Use BitArray worker
            from PyMaSC.handler.bamasc import BACalcWorker
            
            # Create handler with mappability if configured
            handler = None
            if map_config and map_config.is_enabled():
                from PyMaSC.handler.mappability import MappabilityHandler
                handler = MappabilityHandler(
                    str(map_config.mappability_path),
                    calc_config.max_shift,
                    calc_config.read_length or 50,
                    str(map_config.mappability_stats_path) if map_config.mappability_stats_path else None,
                    1  # Single worker in this context
                )
            
            return BACalcWorker(
                order_queue, report_queue, logger_lock,
                bam_path, calc_config.mapq_criteria, calc_config.max_shift,
                calc_config.read_length or 50,
                calc_config.references, calc_config.lengths,
                handler, calc_config.skip_ncc
            )
        
        else:
            # Use standard workers for NCC/MSCC
            if calc_config.skip_ncc and map_config and map_config.is_enabled():
                # MSCC only
                from PyMaSC.handler.masc_worker import MSCCCalcWorker
                return MSCCCalcWorker(
                    order_queue, report_queue, logger_lock,
                    bam_path, calc_config.mapq_criteria, calc_config.max_shift,
                    calc_config.references, calc_config.lengths,
                    str(map_config.mappability_path),
                    calc_config.read_length or 50,
                    {}  # chrom2mappable_len - will be populated during execution
                )
            
            elif map_config and map_config.is_enabled():
                # Both NCC and MSCC
                from PyMaSC.handler.masc_worker import NCCandMSCCCalcWorker
                return NCCandMSCCCalcWorker(
                    order_queue, report_queue, logger_lock,
                    bam_path, calc_config.mapq_criteria, calc_config.max_shift,
                    calc_config.references, calc_config.lengths,
                    str(map_config.mappability_path),
                    calc_config.read_length or 50,
                    {}  # chrom2mappable_len
                )
            
            else:
                # NCC only
                from PyMaSC.handler.masc_worker import NaiveCCCalcWorker
                return NaiveCCCalcWorker(
                    order_queue, report_queue, logger_lock,
                    bam_path, calc_config.mapq_criteria, calc_config.max_shift,
                    calc_config.references, calc_config.lengths
                )