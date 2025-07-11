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
from abc import ABC, abstractmethod
from typing import Optional, Any, Dict
from multiprocessing import Queue
from multiprocessing.synchronize import Lock

from .interfaces import CrossCorrelationCalculator
from .models import (
    CalculationConfig, MappabilityConfig, WorkerConfig,
    CalculationTarget, ImplementationAlgorithm
)

# Import existing calculator implementations
# These will be wrapped to conform to the new interface
from PyMaSC.core.ncc import NaiveCCCalculator as _NativeCCCalculator
from PyMaSC.core.mscc import MSCCCalculator as _NativeMSCCCalculator
from PyMaSC.bacore.mscc import CCBitArrayCalculator as _NativeBitArrayCalculator

# Import BigWig reader for mappability
from PyMaSC.reader.bigwig import BigWigReader
from PyMaSC.core.mappability import BWFeederWithMappableRegionSum

logger = logging.getLogger(__name__)


class ImplementationStrategy(ABC):
    """Abstract base class for implementation strategies.

    This class defines the interface for different implementation algorithms
    (SUCCESSIVE, BITARRAY) to create appropriate calculators based on the
    calculation target and configuration.

    Uses Template Method pattern to handle common validation and routing
    while delegating specific calculator creation to subclasses.
    """

    def create_calculator(
        self,
        target: CalculationTarget,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig] = None,
        logger_lock: Optional[Lock] = None,
        progress_hook: Optional[Any] = None
    ) -> CrossCorrelationCalculator:
        """Create a calculator instance for the specified target.

        This is the template method that handles common validation and
        routes to appropriate abstract methods based on target.

        Args:
            target: What type of cross-correlation to calculate
            config: Calculation configuration
            mappability_config: Optional mappability configuration
            logger_lock: Optional thread lock for logging
            progress_hook: Optional progress reporting hook

        Returns:
            Calculator instance wrapped in adapter

        Raises:
            ValueError: If target is not supported or requirements not met
        """
        # Common validation: MSCC and BOTH targets require mappability
        if target in [CalculationTarget.MSCC, CalculationTarget.BOTH]:
            if not mappability_config or not mappability_config.is_enabled():
                raise ValueError(
                    f"Target {target.value} requires mappability configuration"
                )

        # Route to specific implementation based on target
        if target == CalculationTarget.NCC:
            return self.create_ncc_calculator(
                config, logger_lock, progress_hook
            )
        elif target == CalculationTarget.MSCC:
            return self.create_mscc_calculator(
                config, mappability_config, logger_lock, progress_hook
            )
        elif target == CalculationTarget.BOTH:
            return self.create_both_calculator(
                config, mappability_config, logger_lock, progress_hook
            )
        else:
            raise ValueError(f"Unsupported target: {target.value}")

    @abstractmethod
    def create_ncc_calculator(
        self,
        config: CalculationConfig,
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create a calculator for NCC (Naive Cross-Correlation) only.

        Args:
            config: Calculation configuration
            logger_lock: Optional thread lock for logging
            progress_hook: Optional progress reporting hook

        Returns:
            Calculator instance for NCC calculation
        """
        pass

    @abstractmethod
    def create_mscc_calculator(
        self,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create a calculator for MSCC (Mappability-Sensitive CC) only.

        Args:
            config: Calculation configuration
            mappability_config: Mappability configuration (required)
            logger_lock: Optional thread lock for logging
            progress_hook: Optional progress reporting hook

        Returns:
            Calculator instance for MSCC calculation
        """
        pass

    @abstractmethod
    def create_both_calculator(
        self,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create a calculator for both NCC and MSCC calculations.

        Args:
            config: Calculation configuration
            mappability_config: Mappability configuration (required)
            logger_lock: Optional thread lock for logging
            progress_hook: Optional progress reporting hook

        Returns:
            Calculator instance for both NCC and MSCC calculations
        """
        pass


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
        self._is_mscc = hasattr(native_calculator, 'ref2ccbins') and not self._is_bitarray

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
        return self._calculator.ref2forward_sum  # type: ignore[no-any-return]

    @property
    def ref2reverse_sum(self) -> Dict[str, int]:
        """Reverse read counts by chromosome."""
        return self._calculator.ref2reverse_sum  # type: ignore[no-any-return]

    @property
    def ref2ccbins(self) -> Dict[str, Any]:
        """Cross-correlation bins by chromosome."""
        return self._calculator.ref2ccbins  # type: ignore[no-any-return]

    @property
    def ref2mappable_forward_sum(self) -> Optional[Dict[str, int]]:
        """Mappable forward read counts (BitArray and MSCC)."""
        if self._is_bitarray:
            return self._calculator.ref2mappable_forward_sum  # type: ignore[no-any-return]
        elif self._is_mscc:
            # For pure MSCC, forward sums are mappability-weighted
            return self._calculator.ref2forward_sum  # type: ignore[no-any-return]
        return None

    @property
    def ref2mappable_reverse_sum(self) -> Optional[Dict[str, int]]:
        """Mappable reverse read counts (BitArray and MSCC)."""
        if self._is_bitarray:
            return self._calculator.ref2mappable_reverse_sum  # type: ignore[no-any-return]
        elif self._is_mscc:
            # For pure MSCC, reverse sums are mappability-weighted
            return self._calculator.ref2reverse_sum  # type: ignore[no-any-return]
        return None

    @property
    def ref2mascbins(self) -> Optional[Dict[str, Any]]:
        """MSCC bins (BitArray and MSCC)."""
        if self._is_bitarray:
            return self._calculator.ref2mascbins  # type: ignore[no-any-return]
        elif self._is_mscc:
            # For pure MSCC, ccbins are mappability-weighted
            return self._calculator.ref2ccbins  # type: ignore[no-any-return]
        return None


class CompositeCalculator(CrossCorrelationCalculator):
    """Composite calculator that runs both NCC and MSCC calculations.

    This calculator is used when SUCCESSIVE algorithm is combined with
    mappability data, matching the behavior of the original implementation.
    """

    def __init__(self, ncc_calculator: Any, mscc_calculator: Any):
        """Initialize composite calculator.

        Args:
            ncc_calculator: NCC calculator instance
            mscc_calculator: MSCC calculator instance
        """
        self._ncc_calculator = ncc_calculator
        self._mscc_calculator = mscc_calculator

    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a forward strand read in both calculators."""
        self._ncc_calculator.feed_forward_read(chrom, pos, readlen)
        self._mscc_calculator.feed_forward_read(chrom, pos, readlen)

    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a reverse strand read in both calculators."""
        self._ncc_calculator.feed_reverse_read(chrom, pos, readlen)
        self._mscc_calculator.feed_reverse_read(chrom, pos, readlen)

    def finishup_calculation(self) -> None:
        """Complete calculation in both calculators."""
        self._ncc_calculator.finishup_calculation()
        self._mscc_calculator.finishup_calculation()

    def flush(self) -> None:
        """Flush both calculators."""
        self._ncc_calculator.finishup_calculation()
        self._mscc_calculator.finishup_calculation()

    @property
    def ref2forward_sum(self) -> Dict[str, int]:
        """Forward read counts from NCC calculator."""
        return self._ncc_calculator.ref2forward_sum  # type: ignore[no-any-return]

    @property
    def ref2reverse_sum(self) -> Dict[str, int]:
        """Reverse read counts from NCC calculator."""
        return self._ncc_calculator.ref2reverse_sum  # type: ignore[no-any-return]

    @property
    def ref2ccbins(self) -> Dict[str, Any]:
        """NCC bins from NCC calculator."""
        return self._ncc_calculator.ref2ccbins  # type: ignore[no-any-return]

    @property
    def ref2mappable_forward_sum(self) -> Optional[Dict[str, int]]:
        """Mappable forward read counts from MSCC calculator."""
        return self._mscc_calculator.ref2forward_sum  # type: ignore[no-any-return]

    @property
    def ref2mappable_reverse_sum(self) -> Optional[Dict[str, int]]:
        """Mappable reverse read counts from MSCC calculator."""
        return self._mscc_calculator.ref2reverse_sum  # type: ignore[no-any-return]

    @property
    def ref2mascbins(self) -> Optional[Dict[str, Any]]:
        """MSCC bins from MSCC calculator."""
        return self._mscc_calculator.ref2ccbins  # type: ignore[no-any-return]


class SuccessiveImplementationStrategy(ImplementationStrategy):
    """Implementation strategy for the SUCCESSIVE algorithm.

    This strategy handles the creation of calculators using the memory-efficient
    successive algorithm, which can run NCC, MSCC, or both calculations.
    """

    def create_ncc_calculator(
        self,
        config: CalculationConfig,
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create a calculator for NCC only using SUCCESSIVE algorithm."""
        return self._create_ncc_calculator(config, logger_lock)

    def create_mscc_calculator(
        self,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create a calculator for MSCC only using SUCCESSIVE algorithm."""
        return self._create_mscc_calculator(config, mappability_config, logger_lock)

    def create_both_calculator(
        self,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create a calculator for both NCC and MSCC using SUCCESSIVE algorithm."""
        # SUCCESSIVE with mappability runs both NCC and MSCC
        if mappability_config and mappability_config.is_enabled():
            # Reuse existing create methods to avoid duplication
            ncc_wrapped = self.create_ncc_calculator(config, logger_lock, progress_hook)
            mscc_wrapped = self.create_mscc_calculator(config, mappability_config, logger_lock, progress_hook)
            
            # Extract native calculators from adapters
            # CompositeCalculator needs the raw calculator instances
            if isinstance(ncc_wrapped, CalculatorAdapter):
                ncc_calc = ncc_wrapped._calculator
            else:
                ncc_calc = ncc_wrapped
                
            if isinstance(mscc_wrapped, CalculatorAdapter):
                mscc_calc = mscc_wrapped._calculator
            else:
                mscc_calc = mscc_wrapped
            
            return CompositeCalculator(ncc_calc, mscc_calc)
        else:
            # If no mappability, just return NCC calculator
            return self.create_ncc_calculator(config, logger_lock, progress_hook)

    def _create_ncc_calculator(
        self,
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

    def _create_mscc_calculator(
        self,
        config: CalculationConfig,
        mappability_config: MappabilityConfig,
        logger_lock: Optional[Lock] = None
    ) -> CrossCorrelationCalculator:
        """Create MSCC calculator instance."""
        # Create BWFeeder with mappability calculation for MSCC data
        bwfeeder = BWFeederWithMappableRegionSum(
            str(mappability_config.mappability_path),
            config.max_shift
        )

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


class BitArrayImplementationStrategy(ImplementationStrategy):
    """Implementation strategy for the BITARRAY algorithm.

    This strategy handles the creation of calculators using the memory-intensive
    but fast BitArray algorithm, which can run NCC, MSCC, or both calculations.
    """

    def create_ncc_calculator(
        self,
        config: CalculationConfig,
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create a calculator for NCC only using BITARRAY algorithm."""
        # NCC only: BitArray without mappability requirement
        return self._create_bitarray_calculator(
            config, None, logger_lock, progress_hook, skip_ncc=False
        )

    def create_mscc_calculator(
        self,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create a calculator for MSCC only using BITARRAY algorithm."""
        # MSCC only: BitArray with skip_ncc=True
        modified_config = self._create_mscc_config(config)
        return self._create_bitarray_calculator(
            modified_config, mappability_config, logger_lock, progress_hook, skip_ncc=True
        )

    def create_both_calculator(
        self,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock],
        progress_hook: Optional[Any]
    ) -> CrossCorrelationCalculator:
        """Create a calculator for both NCC and MSCC using BITARRAY algorithm."""
        # BOTH: Use original config, BitArray will handle both
        return self._create_bitarray_calculator(
            config, mappability_config, logger_lock, progress_hook, skip_ncc=False
        )

    def _create_mscc_config(self, config: CalculationConfig) -> CalculationConfig:
        """Create a modified config for MSCC-only calculation."""
        return CalculationConfig(
            target=config.target,
            implementation=config.implementation,
            max_shift=config.max_shift,
            mapq_criteria=config.mapq_criteria,
            skip_ncc=True,  # Force skip NCC for MSCC-only target
            read_length=config.read_length,
            references=config.references,
            lengths=config.lengths,
            esttype=config.esttype,
            chromfilter=config.chromfilter
        )

    def _create_bitarray_calculator(
        self,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig],
        logger_lock: Optional[Lock] = None,
        progress_hook: Optional[Any] = None,
        skip_ncc: Optional[bool] = None
    ) -> CrossCorrelationCalculator:
        """Create BitArray calculator instance."""
        # Create BigWig reader if mappability is enabled
        bwfeeder = None
        if mappability_config and mappability_config.is_enabled():
            bwfeeder = BigWigReader(str(mappability_config.mappability_path))

        # Use read_length from config
        read_len = config.read_length
        if mappability_config and mappability_config.read_len:
            read_len = mappability_config.read_len
        
        if read_len is None:
            raise ValueError("Read length must be specified for BitArray calculator")

        # Use skip_ncc from parameter if provided, otherwise from config
        effective_skip_ncc = skip_ncc if skip_ncc is not None else config.skip_ncc

        native_calc = _NativeBitArrayCalculator(
            config.max_shift,
            read_len,
            config.references,
            config.lengths,
            bwfeeder,
            effective_skip_ncc,
            logger_lock,
            progress_hook
        )
        return CalculatorAdapter(native_calc)


class CalculatorFactory:
    """Factory for creating cross-correlation calculator instances.

    This factory encapsulates the complex logic of creating different
    calculator types with their varying constructor signatures and
    dependencies. It provides a unified interface for calculator creation
    regardless of the underlying implementation.
    """

    @staticmethod
    def create_calculator(
        target: CalculationTarget,
        implementation: ImplementationAlgorithm,
        config: CalculationConfig,
        mappability_config: Optional[MappabilityConfig] = None,
        logger_lock: Optional[Lock] = None,
        progress_hook: Optional[Any] = None
    ) -> CrossCorrelationCalculator:
        """Create a calculator instance using strategy pattern.

        Args:
            target: What type of cross-correlation to calculate
            implementation: How to implement the calculation
            config: Calculation configuration
            mappability_config: Optional mappability configuration
            logger_lock: Optional thread lock for logging
            progress_hook: Optional progress reporting hook

        Returns:
            Calculator instance wrapped in adapter

        Raises:
            ValueError: If combination is not supported or configuration is invalid
        """
        # Select appropriate strategy based on implementation algorithm
        strategy = CalculatorFactory._get_strategy(implementation)

        # Delegate calculator creation to the strategy
        return strategy.create_calculator(
            target, config, mappability_config, logger_lock, progress_hook
        )

    @staticmethod
    def _get_strategy(implementation: ImplementationAlgorithm) -> ImplementationStrategy:
        """Get the appropriate implementation strategy.

        Args:
            implementation: Implementation algorithm enum

        Returns:
            Strategy instance for the specified implementation

        Raises:
            ValueError: If implementation is not supported
        """
        if implementation == ImplementationAlgorithm.SUCCESSIVE:
            return SuccessiveImplementationStrategy()
        elif implementation == ImplementationAlgorithm.BITARRAY:
            return BitArrayImplementationStrategy()
        else:
            raise ValueError(f"Unsupported implementation: {implementation.value}")


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
        # Use unified worker for all target-implementation combinations
        from PyMaSC.core.worker import UnifiedWorker
        return UnifiedWorker(
            worker_config,
            order_queue,
            report_queue,
            logger_lock,
            bam_path
        )