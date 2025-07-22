"""Strategy pattern implementation with mixin-based architecture.

This module implements the Strategy pattern using a sophisticated mixin architecture
that cleanly separates "what to calculate" (targets) from "how to calculate"
(implementations). This design eliminates code duplication and enables easy
extension with new targets or implementations.

Architecture Design:
- Target Mixins: NCCMixin, MSCCMixin, BothMixin (what to calculate)
- Implementation Mixins: BitArrayMixin, SuccessiveMixin (how to calculate)
- Strategy Classes: Combinations of target + implementation mixins using ClassVar

Technical Features:
- ClassVar approach: Properties defined as class variables for zero overhead
- Multiple inheritance: Clean composition of orthogonal concerns
- Type safety: Compile-time validation of strategy combinations
- Zero instantiation overhead: Properties accessed directly from class

Key components:
- CalculationTargetMixin: Abstract base for target-specific logic
- CalculationStrategy: Unified strategy interface with ClassVar properties
- Target Mixins: Encapsulate calculation target logic (NCC/MSCC/BOTH)
- Implementation Mixins: Encapsulate algorithm logic (BitArray/Successive)
- Strategy Classes: Minimal combinations using multiple inheritance
- CalculationContext: Context class for strategy management
"""
from abc import ABC
from typing import Optional, Any, List, ClassVar
import logging

from .interfaces.calculator import CrossCorrelationCalculator
from .factory import CalculatorFactory
from .models import (
    CalculationConfig, MappabilityConfig,
    CalculationTarget, ImplementationAlgorithm,
)

logger = logging.getLogger(__name__)


# === Target Mixins ===
# These mixins define "what to calculate" (NCC/MSCC/BOTH) and their properties
# using ClassVar for zero-overhead property access.

class CalculationTargetMixin(ABC):
    """Abstract base class for calculation target mixins.

    Defines the interface that target mixins must implement using ClassVar.
    This approach provides compile-time type safety and zero runtime overhead.

    Class Variables:
        calculation_target: The specific target this mixin handles
        requires_mappability: Whether this target requires mappability data
        supports_skip_ncc: Whether this target supports skipping NCC calculation
    """

    calculation_target: ClassVar[CalculationTarget]
    requires_mappability: ClassVar[bool]
    supports_skip_ncc: ClassVar[bool]


class NCCMixin(CalculationTargetMixin):
    """Mixin for NCC (Naive Cross-Correlation) target-specific logic.

    NCC calculates basic cross-correlation without mappability correction.
    This is the fastest calculation but may be affected by mappability bias.
    """

    calculation_target = CalculationTarget.NCC
    requires_mappability = False  # NCC doesn't need mappability data
    supports_skip_ncc = False     # Can't skip what we're calculating


class MSCCMixin(CalculationTargetMixin):
    """Mixin for MSCC (Mappability-Sensitive Cross-Correlation) target-specific logic.

    MSCC calculates cross-correlation with mappability correction to address
    the "phantom peak" problem caused by non-uniform mappability.
    """

    calculation_target = CalculationTarget.MSCC
    requires_mappability = True   # MSCC requires mappability data
    supports_skip_ncc = False     # MSCC-only mode, NCC already skipped


class BothMixin(CalculationTargetMixin):
    """Mixin for BOTH (NCC and MSCC) target-specific logic.

    Calculates both NCC and MSCC for comparison and analysis.
    This is the most comprehensive option but requires mappability data.
    """

    calculation_target = CalculationTarget.BOTH
    requires_mappability = True   # MSCC component requires mappability
    supports_skip_ncc = True      # Can skip NCC to get MSCC-only


# === Strategy Base Class ===
# This class combines target properties with implementation interface to create
# the unified strategy abstraction.

class CalculationStrategy(CalculationTargetMixin):
    """Abstract base class for calculation strategies.

    Combines target-specific properties (from CalculationTargetMixin) with
    implementation-specific behavior to create a unified strategy interface.
    Each strategy encapsulates the logic for creating and configuring
    calculators for a particular target-implementation combination.

    Class Variables:
        implementation_algorithm: Which implementation algorithm to use
    """

    implementation_algorithm: ClassVar[ImplementationAlgorithm]

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
        return CalculatorFactory.create_calculator(
            self.calculation_target,
            self.implementation_algorithm,
            config,
            mappability_config,
            logger_lock,
            progress_hook
        )

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

        if self.requires_mappability and (not mappability_config or not mappability_config.is_enabled()):
            errors.append(f"{self.calculation_target.value} target requires mappability configuration")

        if config.skip_ncc and not self.supports_skip_ncc:
            errors.append(f"{self.calculation_target.value} target does not support skip_ncc option")

        return errors


# === Implementation Mixins ===
# These mixins define "how to calculate" by specifying the implementation algorithm.
# They inherit from CalculationStrategy to get the unified interface.

class BitArrayMixin(CalculationStrategy):
    """Mixin for BitArray implementation-specific logic.

    BitArray implementation uses memory-intensive but fast binary arrays.
    It's optimized for speed but requires more memory (~250MB for human genome).
    """

    implementation_algorithm = ImplementationAlgorithm.BITARRAY


class SuccessiveMixin(CalculationStrategy):
    """Mixin for Successive implementation-specific logic.

    Successive implementation uses memory-efficient algorithms.
    It processes data sequentially with lower memory usage.
    """

    implementation_algorithm = ImplementationAlgorithm.SUCCESSIVE


# === Strategy Classes ===
# These classes combine target and implementation mixins to create concrete strategies.
# Each strategy is defined by exactly one target mixin and one implementation mixin.
# The 'pass' statements demonstrate the power of the mixin architecture -
# all functionality is inherited from the mixins.

class NCCBitArrayStrategy(NCCMixin, BitArrayMixin):
    """Strategy for NCC calculation using BitArray implementation.

    Combines:
    - NCCMixin: Provides NCC target-specific logic
    - BitArrayMixin: Provides BitArray implementation logic
    - CalculationStrategy: Abstract base class interface
    """
    pass


class MSCCBitArrayStrategy(MSCCMixin, BitArrayMixin):
    """Strategy for MSCC calculation using BitArray implementation.

    Combines:
    - MSCCMixin: Provides MSCC target-specific logic
    - BitArrayMixin: Provides BitArray implementation logic
    - CalculationStrategy: Abstract base class interface
    """
    pass


class BothBitArrayStrategy(BothMixin, BitArrayMixin):
    """Strategy for both NCC and MSCC calculation using BitArray implementation.

    Combines:
    - BothMixin: Provides BOTH target-specific logic
    - BitArrayMixin: Provides BitArray implementation logic
    - CalculationStrategy: Abstract base class interface
    """
    pass


class NCCSuccessiveStrategy(NCCMixin, SuccessiveMixin):
    """Strategy for NCC calculation using Successive implementation.

    Combines:
    - NCCMixin: Provides NCC target-specific logic
    - SuccessiveMixin: Provides Successive implementation logic
    - CalculationStrategy: Abstract base class interface
    """
    pass


class MSCCSuccessiveStrategy(MSCCMixin, SuccessiveMixin):
    """Strategy for MSCC calculation using Successive implementation.

    Combines:
    - MSCCMixin: Provides MSCC target-specific logic
    - SuccessiveMixin: Provides Successive implementation logic
    - CalculationStrategy: Abstract base class interface
    """
    pass


class BothSuccessiveStrategy(BothMixin, SuccessiveMixin):
    """Strategy for both NCC and MSCC calculation using Successive implementation.

    Combines:
    - BothMixin: Provides BOTH target-specific logic
    - SuccessiveMixin: Provides Successive implementation logic
    - CalculationStrategy: Abstract base class interface
    """
    pass


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

        # Explicit combination strategy mapping
        # This table provides O(1) lookup for all valid target-implementation combinations.
        # Each entry maps a (target, implementation) tuple to its specific strategy instance.
        # The table ensures that only supported combinations are available and provides
        # type-safe strategy selection at runtime.
        self._strategies = {
            # NCC strategies
            (CalculationTarget.NCC, ImplementationAlgorithm.BITARRAY): NCCBitArrayStrategy(),
            (CalculationTarget.NCC, ImplementationAlgorithm.SUCCESSIVE): NCCSuccessiveStrategy(),
            # MSCC strategies
            (CalculationTarget.MSCC, ImplementationAlgorithm.BITARRAY): MSCCBitArrayStrategy(),
            (CalculationTarget.MSCC, ImplementationAlgorithm.SUCCESSIVE): MSCCSuccessiveStrategy(),
            # BOTH strategies
            (CalculationTarget.BOTH, ImplementationAlgorithm.BITARRAY): BothBitArrayStrategy(),
            (CalculationTarget.BOTH, ImplementationAlgorithm.SUCCESSIVE): BothSuccessiveStrategy(),
        }

    def set_strategy(self, strategy: CalculationStrategy) -> None:
        """Set the calculation strategy.

        Args:
            strategy: New calculation strategy to use
        """
        self._strategy = strategy
        logger.debug(f"Strategy changed to {strategy.__class__.__name__}")

    def set_strategy_by_concepts(self, target: CalculationTarget, implementation: ImplementationAlgorithm) -> None:
        """Set strategy based on calculation target and implementation algorithm.

        Args:
            target: What to calculate
            implementation: How to calculate

        Raises:
            ValueError: If combination is not supported
        """
        combination = (target, implementation)
        if combination not in self._strategies:
            raise ValueError(f"Unsupported combination: target={target.value}, implementation={implementation.value}")

        self._strategy = self._strategies[combination]
        logger.debug(f"Strategy set to target={target.value}, implementation={implementation.value}")

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

        return self._strategy.requires_mappability

    def supports_skip_ncc(self) -> bool:
        """Check if current strategy supports skipping NCC.

        Returns:
            True if skip_ncc is supported

        Raises:
            RuntimeError: If no strategy is set
        """
        if self._strategy is None:
            raise RuntimeError("No calculation strategy set")

        return self._strategy.supports_skip_ncc

    @staticmethod
    def create_from_config(config: CalculationConfig) -> 'CalculationContext':
        """Create context with strategy based on configuration.

        Args:
            config: Calculation configuration

        Returns:
            CalculationContext with appropriate strategy set
        """
        context = CalculationContext()

        # Use concept-based approach
        if config.target is not None and config.implementation is not None:
            context.set_strategy_by_concepts(config.target, config.implementation)
        else:
            raise ValueError("Configuration must specify both target and implementation")

        return context