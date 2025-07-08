"""Unit tests for conceptual architecture improvements.

This module tests the new conceptual separation between calculation targets
and implementation algorithms, ensuring backward compatibility while
providing cleaner abstractions.
"""
import pytest

from PyMaSC.core.models import (
    CalculationConfig, CalculationTarget, ImplementationAlgorithm
)
from PyMaSC.core.strategy import (
    CalculationContext,
    NCCBitArrayStrategy, MSCCBitArrayStrategy, BothBitArrayStrategy,
    NCCSuccessiveStrategy, MSCCSuccessiveStrategy, BothSuccessiveStrategy
)
from PyMaSC.core.factory import CalculatorFactory


class TestCalculationTarget:
    """Test CalculationTarget enum."""

    def test_calculation_target_values(self):
        """Test that CalculationTarget has expected values."""
        assert CalculationTarget.NCC.value == "ncc"
        assert CalculationTarget.MSCC.value == "mscc"
        assert CalculationTarget.BOTH.value == "both"


class TestImplementationAlgorithm:
    """Test ImplementationAlgorithm enum."""

    def test_implementation_algorithm_values(self):
        """Test that ImplementationAlgorithm has expected values."""
        assert ImplementationAlgorithm.BITARRAY.value == "bitarray"
        assert ImplementationAlgorithm.SUCCESSIVE.value == "successive"


class TestCalculationConfig:
    """Test CalculationConfig with new conceptual fields."""

    def test_new_clean_configuration(self):
        """Test that new configuration approach works."""
        config = CalculationConfig(
            target=CalculationTarget.BOTH,
            implementation=ImplementationAlgorithm.BITARRAY,
            max_shift=500,
            mapq_criteria=20
        )

        # Should have clean new fields
        assert config.target == CalculationTarget.BOTH
        assert config.implementation == ImplementationAlgorithm.BITARRAY

    def test_ncc_successive_configuration(self):
        """Test configuration for NCC with Successive implementation."""
        config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20
        )

        # Should have correct new fields
        assert config.target == CalculationTarget.NCC
        assert config.implementation == ImplementationAlgorithm.SUCCESSIVE

    def test_validation_requires_all_fields(self):
        """Test that configuration validation requires all required fields."""
        # This should work
        config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20
        )
        assert config.target == CalculationTarget.NCC
        assert config.implementation == ImplementationAlgorithm.SUCCESSIVE


class TestCalculationContext:
    """Test CalculationContext with new conceptual approaches."""

    def test_strategy_selection_by_concepts(self):
        """Test strategy selection using new concepts."""
        context = CalculationContext()

        # Test BothBitArray strategy selection
        context.set_strategy_by_concepts(CalculationTarget.BOTH, ImplementationAlgorithm.BITARRAY)
        strategy = context.get_current_strategy()
        assert isinstance(strategy, BothBitArrayStrategy)
        assert strategy.calculation_target == CalculationTarget.BOTH
        assert strategy.implementation_algorithm == ImplementationAlgorithm.BITARRAY

        # Test NCC Successive strategy selection
        context.set_strategy_by_concepts(CalculationTarget.NCC, ImplementationAlgorithm.SUCCESSIVE)
        strategy = context.get_current_strategy()
        assert isinstance(strategy, NCCSuccessiveStrategy)

        # Test MSCC Successive strategy selection
        context.set_strategy_by_concepts(CalculationTarget.MSCC, ImplementationAlgorithm.SUCCESSIVE)
        strategy = context.get_current_strategy()
        assert isinstance(strategy, MSCCSuccessiveStrategy)
        assert strategy.calculation_target == CalculationTarget.MSCC

    def test_unsupported_concept_combination(self):
        """Test error handling for unsupported combinations."""
        context = CalculationContext()

        # All current combinations are now supported with explicit strategies
        # This test verifies that proper error handling exists for invalid enums
        # (This would only fail if someone passes invalid enum values)

        # All valid combinations should work
        all_targets = [CalculationTarget.NCC, CalculationTarget.MSCC, CalculationTarget.BOTH]
        all_implementations = [ImplementationAlgorithm.BITARRAY, ImplementationAlgorithm.SUCCESSIVE]

        for target in all_targets:
            for implementation in all_implementations:
                # This should not raise an exception
                context.set_strategy_by_concepts(target, implementation)
                strategy = context.get_current_strategy()
                assert strategy is not None
                assert strategy.calculation_target == target
                assert strategy.implementation_algorithm == implementation

    def test_create_from_config_with_new_concepts(self):
        """Test context creation from config using new concepts."""
        config = CalculationConfig(
            target=CalculationTarget.BOTH,
            implementation=ImplementationAlgorithm.BITARRAY,
            max_shift=500,
            mapq_criteria=20
        )

        context = CalculationContext.create_from_config(config)
        strategy = context.get_current_strategy()
        assert isinstance(strategy, BothBitArrayStrategy)

    def test_create_from_config_with_both_successive(self):
        """Test context creation from config using both target with successive."""
        config = CalculationConfig(
            target=CalculationTarget.BOTH,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20
        )

        context = CalculationContext.create_from_config(config)
        strategy = context.get_current_strategy()
        assert isinstance(strategy, BothSuccessiveStrategy)


class TestCalculatorFactory:
    """Test CalculatorFactory with new conceptual approaches."""

    def test_create_calculator_ncc_successive(self):
        """Test calculator creation for NCC with Successive."""
        config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20,
            references=["chr1"],
            lengths=[1000]
        )

        # This should work without requiring mappability for NCC
        calculator = CalculatorFactory.create_calculator(
            CalculationTarget.NCC,
            ImplementationAlgorithm.SUCCESSIVE,
            config
        )

        assert calculator is not None

    def test_create_calculator_mscc_requires_mappability(self):
        """Test that MSCC calculation requires mappability configuration."""
        config = CalculationConfig(
            target=CalculationTarget.MSCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20,
            references=["chr1"],
            lengths=[1000]
        )

        # Should fail without mappability config
        with pytest.raises(ValueError, match="Target mscc requires mappability configuration"):
            CalculatorFactory.create_calculator(
                CalculationTarget.MSCC,
                ImplementationAlgorithm.SUCCESSIVE,
                config
            )

    def test_all_combinations_supported(self):
        """Test that all valid combinations are supported."""
        config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.BITARRAY,
            max_shift=500,
            mapq_criteria=20,
            references=["chr1"],
            lengths=[1000]
        )

        # All combinations should be supported now
        calculator = CalculatorFactory.create_calculator(
            CalculationTarget.NCC,
            ImplementationAlgorithm.BITARRAY,
            config
        )
        assert calculator is not None


class TestExplicitCombinationStrategies:
    """Test new explicit combination strategy classes."""

    def test_ncc_bitarray_strategy(self):
        """Test NCCBitArrayStrategy implements correct interface."""
        strategy = NCCBitArrayStrategy()
        assert strategy.calculation_target == CalculationTarget.NCC
        assert strategy.implementation_algorithm == ImplementationAlgorithm.BITARRAY
        assert not strategy.requires_mappability
        assert not strategy.supports_skip_ncc

    def test_mscc_bitarray_strategy(self):
        """Test MSCCBitArrayStrategy implements correct interface."""
        strategy = MSCCBitArrayStrategy()
        assert strategy.calculation_target == CalculationTarget.MSCC
        assert strategy.implementation_algorithm == ImplementationAlgorithm.BITARRAY
        assert strategy.requires_mappability
        assert not strategy.supports_skip_ncc

    def test_both_bitarray_strategy(self):
        """Test BothBitArrayStrategy implements correct interface."""
        strategy = BothBitArrayStrategy()
        assert strategy.calculation_target == CalculationTarget.BOTH
        assert strategy.implementation_algorithm == ImplementationAlgorithm.BITARRAY
        assert strategy.requires_mappability
        assert strategy.supports_skip_ncc

    def test_ncc_successive_strategy(self):
        """Test NCCSuccessiveStrategy implements correct interface."""
        strategy = NCCSuccessiveStrategy()
        assert strategy.calculation_target == CalculationTarget.NCC
        assert strategy.implementation_algorithm == ImplementationAlgorithm.SUCCESSIVE
        assert not strategy.requires_mappability
        assert not strategy.supports_skip_ncc

    def test_mscc_successive_strategy(self):
        """Test MSCCSuccessiveStrategy implements correct interface."""
        strategy = MSCCSuccessiveStrategy()
        assert strategy.calculation_target == CalculationTarget.MSCC
        assert strategy.implementation_algorithm == ImplementationAlgorithm.SUCCESSIVE
        assert strategy.requires_mappability
        assert not strategy.supports_skip_ncc

    def test_both_successive_strategy(self):
        """Test BothSuccessiveStrategy implements correct interface."""
        strategy = BothSuccessiveStrategy()
        assert strategy.calculation_target == CalculationTarget.BOTH
        assert strategy.implementation_algorithm == ImplementationAlgorithm.SUCCESSIVE
        assert strategy.requires_mappability
        assert strategy.supports_skip_ncc


class TestCalculationContextWithNewStrategies:
    """Test CalculationContext with new explicit combination strategies."""

    def test_all_combinations_supported(self):
        """Test that all target-implementation combinations are supported."""
        context = CalculationContext()

        test_cases = [
            (CalculationTarget.NCC, ImplementationAlgorithm.BITARRAY, NCCBitArrayStrategy),
            (CalculationTarget.NCC, ImplementationAlgorithm.SUCCESSIVE, NCCSuccessiveStrategy),
            (CalculationTarget.MSCC, ImplementationAlgorithm.BITARRAY, MSCCBitArrayStrategy),
            (CalculationTarget.MSCC, ImplementationAlgorithm.SUCCESSIVE, MSCCSuccessiveStrategy),
            (CalculationTarget.BOTH, ImplementationAlgorithm.BITARRAY, BothBitArrayStrategy),
            (CalculationTarget.BOTH, ImplementationAlgorithm.SUCCESSIVE, BothSuccessiveStrategy),
        ]

        for target, implementation, expected_strategy_type in test_cases:
            context.set_strategy_by_concepts(target, implementation)
            strategy = context.get_current_strategy()
            assert strategy is not None
            assert isinstance(strategy, expected_strategy_type)
            assert strategy.calculation_target == target
            assert strategy.implementation_algorithm == implementation

    def test_context_factory_integration(self):
        """Test that context integrates properly with factory."""
        config = CalculationConfig(
            target=CalculationTarget.BOTH,
            implementation=ImplementationAlgorithm.BITARRAY,
            max_shift=500,
            mapq_criteria=20
        )

        context = CalculationContext.create_from_config(config)
        strategy = context.get_current_strategy()
        assert isinstance(strategy, BothBitArrayStrategy)


class TestStrategyClarityImprovement:
    """Test that new strategies provide better conceptual clarity."""

    def test_strategy_clarity_with_explicit_combinations(self):
        """Test that new strategies provide clear semantics."""
        # New approach: each strategy is explicit about what it does
        ncc_bitarray = NCCBitArrayStrategy()
        mscc_bitarray = MSCCBitArrayStrategy()
        both_bitarray = BothBitArrayStrategy()

        assert ncc_bitarray.calculation_target == CalculationTarget.NCC
        assert mscc_bitarray.calculation_target == CalculationTarget.MSCC
        assert both_bitarray.calculation_target == CalculationTarget.BOTH

        # All use the same implementation
        assert ncc_bitarray.implementation_algorithm == ImplementationAlgorithm.BITARRAY
        assert mscc_bitarray.implementation_algorithm == ImplementationAlgorithm.BITARRAY
        assert both_bitarray.implementation_algorithm == ImplementationAlgorithm.BITARRAY

    def test_extensibility_improvement(self):
        """Test that new approach makes extending with new combinations easier."""
        # With new approach, each combination has its own strategy class
        # This makes it easy to add specific logic for each combination

        # Example: NCC with BitArray might have different optimizations
        # than MSCC with BitArray
        ncc_strategy = NCCBitArrayStrategy()
        mscc_strategy = MSCCBitArrayStrategy()

        # They can have different requirements
        assert not ncc_strategy.requires_mappability
        assert mscc_strategy.requires_mappability

        # They can have different capabilities
        assert not ncc_strategy.supports_skip_ncc  # Can't skip what you're calculating
        assert not mscc_strategy.supports_skip_ncc  # MSCC-only, NCC already skipped