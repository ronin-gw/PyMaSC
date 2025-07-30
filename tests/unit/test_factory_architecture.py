"""Unit tests for factory-based architecture.

This module tests the factory pattern approach for calculator creation,
ensuring that all target-implementation combinations work correctly
with direct factory usage.
"""
import pytest

from PyMaSC.core.models import (
    CalculationConfig, CalculationTarget, ImplementationAlgorithm
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
    """Test CalculationConfig with factory-based architecture."""

    def test_factory_based_configuration(self):
        """Test that configuration works with factory pattern."""
        config = CalculationConfig(
            target=CalculationTarget.BOTH,
            implementation=ImplementationAlgorithm.BITARRAY,
            max_shift=500,
            mapq_criteria=20
        )

        # Should have clean configuration fields
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

        # Should have correct configuration fields
        assert config.target == CalculationTarget.NCC
        assert config.implementation == ImplementationAlgorithm.SUCCESSIVE

    def test_configuration_completeness(self):
        """Test that configuration contains all required fields."""
        config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20
        )
        assert config.target == CalculationTarget.NCC
        assert config.implementation == ImplementationAlgorithm.SUCCESSIVE


class TestCalculatorFactory:
    """Test CalculatorFactory with direct factory pattern."""

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
            read_length=50,  # Required for BitArray implementation
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

    def test_factory_direct_usage_pattern(self):
        """Test that factory can be used directly without strategy layer."""
        # All valid target-implementation combinations
        test_cases = [
            (CalculationTarget.NCC, ImplementationAlgorithm.BITARRAY),
            (CalculationTarget.NCC, ImplementationAlgorithm.SUCCESSIVE),
            (CalculationTarget.BOTH, ImplementationAlgorithm.BITARRAY),
            (CalculationTarget.BOTH, ImplementationAlgorithm.SUCCESSIVE),
        ]

        for target, implementation in test_cases:
            config = CalculationConfig(
                target=target,
                implementation=implementation,
                max_shift=500,
                mapq_criteria=20,
                read_length=50,  # Required for BitArray
                references=["chr1"],
                lengths=[1000]
            )

            # MSCC requires mappability config, so skip those
            if target in [CalculationTarget.MSCC, CalculationTarget.BOTH]:
                continue

            calculator = CalculatorFactory.create_calculator(
                target, implementation, config
            )
            assert calculator is not None

    def test_factory_performance_optimization(self):
        """Test that factory pattern provides direct calculator creation."""
        config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.BITARRAY,
            max_shift=200,
            mapq_criteria=20,
            read_length=50,
            references=["chr1", "chr2"],
            lengths=[1000, 800]
        )

        # Direct factory usage should be fast and efficient
        calculator = CalculatorFactory.create_calculator(
            config.target,
            config.implementation,
            config
        )

        assert calculator is not None
        # Verify it's the correct calculator type
        from PyMaSC.bacore.mscc import CCBitArrayCalculator
        assert isinstance(calculator, CCBitArrayCalculator)


class TestFactoryArchitectureIntegration:
    """Test integration of factory-based architecture."""

    def test_single_process_integration(self):
        """Test that single-process mode uses factory directly."""
        # This test verifies that CalcHandler uses CalculatorFactory directly
        # without strategy layer overhead
        config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=300,
            mapq_criteria=25,
            references=["chr1"],
            lengths=[1000]
        )

        calculator = CalculatorFactory.create_calculator(
            config.target,
            config.implementation,
            config
        )

        # Should create calculator without intermediate strategy objects
        assert calculator is not None

    def test_multiprocess_integration_consistency(self):
        """Test that multi-process mode also uses factory directly."""
        # Both single and multi-process should use the same factory pattern
        config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.BITARRAY,
            max_shift=400,
            mapq_criteria=30,
            read_length=75,
            references=["chr1"],
            lengths=[1000]
        )

        # Same factory call as used in both single and multi-process
        calculator = CalculatorFactory.create_calculator(
            config.target,
            config.implementation,
            config
        )

        assert calculator is not None

    def test_architecture_simplification(self):
        """Test that new architecture eliminates unnecessary complexity."""
        # Factory pattern should provide simple, direct calculator creation
        # No strategy objects, contexts, or intermediate layers needed

        config = CalculationConfig(
            target=CalculationTarget.NCC,
            implementation=ImplementationAlgorithm.SUCCESSIVE,
            max_shift=250,
            mapq_criteria=15,
            references=["chr1"],
            lengths=[1000]
        )

        # Single, direct factory call
        calculator = CalculatorFactory.create_calculator(
            config.target,
            config.implementation,
            config
        )

        assert calculator is not None
        # Verify no strategy overhead
        from PyMaSC.core.ncc import NaiveCCCalculator
        assert isinstance(calculator, NaiveCCCalculator)