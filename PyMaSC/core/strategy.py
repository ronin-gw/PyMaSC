"""Strategy pattern implementation for algorithm selection.

This module implements the Strategy pattern to decouple algorithm selection
from the calculation logic. It enables runtime algorithm switching and
provides a clean abstraction for different cross-correlation algorithms.

Key components:
- CalculationStrategy: Abstract base for all algorithm strategies
- BitArrayStrategy: Strategy for BitArray algorithm
- SuccessiveStrategy: Strategy for Successive/NCC algorithm
- MSCCStrategy: Strategy for pure MSCC algorithm
- CalculationContext: Context class for strategy management
"""
from abc import ABC, abstractmethod
from typing import Optional, Any, Dict, List
import logging

from .interfaces import CrossCorrelationCalculator
from .factory import CalculatorFactory
from .models import (
    CalculationConfig, MappabilityConfig, WorkerConfig,
    AlgorithmType, ExecutionMode
)

logger = logging.getLogger(__name__)


class CalculationStrategy(ABC):
    """Abstract base class for calculation strategies.

    Defines the interface that all algorithm strategies must implement.
    Each strategy encapsulates the specific logic for creating and
    configuring calculators for a particular algorithm.
    """

    @abstractmethod
    def create_calculator(self,
                          config: CalculationConfig,
                          mappability_config: Optional[MappabilityConfig] = None,
                          logger_lock: Optional[Any] = None,
                          progress_hook: Optional[Any] = None) -> CrossCorrelationCalculator:
        """Create a calculator instance for this strategy.

        Args:
            config: Calculation configuration
            mappability_config: Optional mappability configuration
            logger_lock: Optional thread lock for logging
            progress_hook: Optional progress reporting hook

        Returns:
            Calculator instance implementing CrossCorrelationCalculator
        """
        pass

    @abstractmethod
    def get_algorithm_type(self) -> AlgorithmType:
        """Get the algorithm type for this strategy.

        Returns:
            AlgorithmType enum value
        """
        pass

    @abstractmethod
    def requires_mappability(self) -> bool:
        """Check if this strategy requires mappability data.

        Returns:
            True if mappability data is required
        """
        pass

    @abstractmethod
    def supports_skip_ncc(self) -> bool:
        """Check if this strategy supports skipping NCC calculation.

        Returns:
            True if skip_ncc option is supported
        """
        pass

    def validate_config(self, 
                        config: CalculationConfig,
                        mappability_config: Optional[MappabilityConfig]) -> List[str]:
        """Validate configuration for this strategy.

        Args:
            config: Calculation configuration to validate
            mappability_config: Optional mappability configuration

        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []

        if self.requires_mappability() and (not mappability_config or not mappability_config.is_enabled()):
            errors.append(f"{self.get_algorithm_type().value} algorithm requires mappability configuration")

        if config.skip_ncc and not self.supports_skip_ncc():
            errors.append(f"{self.get_algorithm_type().value} algorithm does not support skip_ncc option")

        return errors


class BitArrayStrategy(CalculationStrategy):
    """Strategy for BitArray algorithm.

    Implements the creation and configuration logic specific to the
    BitArray algorithm, which uses bit arrays for efficient memory usage
    and supports both NCC and MSCC calculations.
    """

    def create_calculator(self,
                          config: CalculationConfig,
                          mappability_config: Optional[MappabilityConfig] = None,
                          logger_lock: Optional[Any] = None,
                          progress_hook: Optional[Any] = None) -> CrossCorrelationCalculator:
        """Create BitArray calculator instance."""
        return CalculatorFactory.create_calculator(
            AlgorithmType.BITARRAY,
            config,
            mappability_config,
            logger_lock,
            progress_hook
        )

    def get_algorithm_type(self) -> AlgorithmType:
        """Return BitArray algorithm type."""
        return AlgorithmType.BITARRAY

    def requires_mappability(self) -> bool:
        """BitArray requires mappability data."""
        return True

    def supports_skip_ncc(self) -> bool:
        """BitArray supports skipping NCC."""
        return True


class SuccessiveStrategy(CalculationStrategy):
    """Strategy for Successive/NCC algorithm.

    Implements the creation and configuration logic specific to the
    Successive algorithm, which is memory-efficient and suitable for
    smaller datasets or when mappability is not required.
    """

    def create_calculator(self,
                          config: CalculationConfig,
                          mappability_config: Optional[MappabilityConfig] = None,
                          logger_lock: Optional[Any] = None,
                          progress_hook: Optional[Any] = None) -> CrossCorrelationCalculator:
        """Create Successive/NCC calculator instance."""
        return CalculatorFactory.create_calculator(
            AlgorithmType.SUCCESSIVE,
            config,
            mappability_config,
            logger_lock,
            progress_hook
        )

    def get_algorithm_type(self) -> AlgorithmType:
        """Return Successive algorithm type."""
        return AlgorithmType.SUCCESSIVE

    def requires_mappability(self) -> bool:
        """Successive does not require mappability data."""
        return False

    def supports_skip_ncc(self) -> bool:
        """Successive is NCC-only, cannot skip."""
        return False


class MSCCStrategy(CalculationStrategy):
    """Strategy for pure MSCC algorithm.

    Implements the creation and configuration logic specific to the
    MSCC algorithm, which requires mappability data and calculates
    mappability-sensitive cross-correlation.
    """

    def create_calculator(self,
                          config: CalculationConfig,
                          mappability_config: Optional[MappabilityConfig] = None,
                          logger_lock: Optional[Any] = None,
                          progress_hook: Optional[Any] = None) -> CrossCorrelationCalculator:
        """Create MSCC calculator instance."""
        return CalculatorFactory.create_calculator(
            AlgorithmType.MSCC,
            config,
            mappability_config,
            logger_lock,
            progress_hook
        )

    def get_algorithm_type(self) -> AlgorithmType:
        """Return MSCC algorithm type."""
        return AlgorithmType.MSCC

    def requires_mappability(self) -> bool:
        """MSCC requires mappability data."""
        return True

    def supports_skip_ncc(self) -> bool:
        """MSCC is mappability-only, NCC is always skipped."""
        return False


class CalculationContext:
    """Context class for managing calculation strategies.

    This class provides the context in which different calculation
    strategies operate. It manages strategy selection and provides
    a unified interface for calculator creation regardless of the
    underlying algorithm.

    The context enables runtime strategy switching, making it easy
    to change algorithms based on configuration or user preferences.
    """

    def __init__(self, strategy: Optional[CalculationStrategy] = None):
        """Initialize context with optional initial strategy.

        Args:
            strategy: Initial calculation strategy (can be None)
        """
        self._strategy = strategy
        self._available_strategies = {
            AlgorithmType.BITARRAY: BitArrayStrategy(),
            AlgorithmType.SUCCESSIVE: SuccessiveStrategy(),
            AlgorithmType.MSCC: MSCCStrategy(),
            AlgorithmType.NAIVE_CC: SuccessiveStrategy()  # NCC uses Successive implementation
        }

    def set_strategy(self, strategy: CalculationStrategy) -> None:
        """Set the calculation strategy.

        Args:
            strategy: New calculation strategy to use
        """
        self._strategy = strategy
        logger.debug(f"Strategy changed to {strategy.get_algorithm_type().value}")

    def set_strategy_by_algorithm(self, algorithm: AlgorithmType) -> None:
        """Set strategy based on algorithm type.

        Args:
            algorithm: Algorithm type to select strategy for

        Raises:
            ValueError: If algorithm is not supported
        """
        if algorithm not in self._available_strategies:
            raise ValueError(f"Unsupported algorithm: {algorithm}")

        self._strategy = self._available_strategies[algorithm]
        logger.debug(f"Strategy set to {algorithm.value}")

    def get_current_strategy(self) -> Optional[CalculationStrategy]:
        """Get the current calculation strategy.

        Returns:
            Current strategy or None if not set
        """
        return self._strategy

    def create_calculator(self,
                          config: CalculationConfig,
                          mappability_config: Optional[MappabilityConfig] = None,
                          logger_lock: Optional[Any] = None,
                          progress_hook: Optional[Any] = None) -> CrossCorrelationCalculator:
        """Create calculator using current strategy.

        Args:
            config: Calculation configuration
            mappability_config: Optional mappability configuration
            logger_lock: Optional thread lock for logging
            progress_hook: Optional progress reporting hook

        Returns:
            Calculator instance

        Raises:
            RuntimeError: If no strategy is set
        """
        if self._strategy is None:
            raise RuntimeError("No calculation strategy set")

        # Validate configuration for current strategy
        errors = self._strategy.validate_config(config, mappability_config)
        if errors:
            raise ValueError("Configuration validation failed:\n" + 
                           "\n".join(f"  - {error}" for error in errors))

        return self._strategy.create_calculator(
            config, mappability_config, logger_lock, progress_hook
        )

    def requires_mappability(self) -> bool:
        """Check if current strategy requires mappability.

        Returns:
            True if mappability is required

        Raises:
            RuntimeError: If no strategy is set
        """
        if self._strategy is None:
            raise RuntimeError("No calculation strategy set")

        return self._strategy.requires_mappability()

    def supports_skip_ncc(self) -> bool:
        """Check if current strategy supports skipping NCC.

        Returns:
            True if skip_ncc is supported

        Raises:
            RuntimeError: If no strategy is set
        """
        if self._strategy is None:
            raise RuntimeError("No calculation strategy set")

        return self._strategy.supports_skip_ncc()

    @staticmethod
    def create_from_config(config: CalculationConfig) -> 'CalculationContext':
        """Create context with strategy based on configuration.

        Args:
            config: Calculation configuration

        Returns:
            CalculationContext with appropriate strategy set
        """
        context = CalculationContext()
        context.set_strategy_by_algorithm(config.algorithm)
        return context