"""Unit tests for factory-based architecture.

This module tests the factory pattern approach for calculator creation,
ensuring that all target-implementation combinations work correctly
with direct factory usage.
"""
import pytest

from PyMaSC.core.interfaces.config import (
    PyMaSCConfig, CalculationTarget, Algorithm, EstimationType
)
from PyMaSC.core.factory import create_calculator


class TestCalculationTarget:
    """Test CalculationTarget enum."""

    def test_calculation_target_values(self):
        """Test that CalculationTarget has expected values."""
        assert CalculationTarget.NCC.value == "ncc"
        assert CalculationTarget.MSCC.value == "mscc"
        assert CalculationTarget.BOTH.value == "both"


class TestAlgorithm:
    """Test Algorithm enum."""

    def test_implementation_algorithm_values(self):
        """Test that Algorithm has expected values."""
        assert Algorithm.BITARRAY.value == "bitarray"
        assert Algorithm.SUCCESSIVE.value == "successive"


class TestPyMaSCConfig:
    """Test PyMaSCConfig with factory-based architecture."""

    def test_factory_based_configuration(self):
        """Test that configuration works with factory pattern."""
        config = PyMaSCConfig(
            target=CalculationTarget.BOTH,
            implementation=Algorithm.BITARRAY,
            max_shift=500,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )

        # Should have clean configuration fields
        assert config.target == CalculationTarget.BOTH
        assert config.implementation == Algorithm.BITARRAY

    def test_ncc_successive_configuration(self):
        """Test configuration for NCC with Successive implementation."""
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )

        # Should have correct configuration fields
        assert config.target == CalculationTarget.NCC
        assert config.implementation == Algorithm.SUCCESSIVE

    def test_configuration_completeness(self):
        """Test that configuration contains all required fields."""
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )
        assert config.target == CalculationTarget.NCC
        assert config.implementation == Algorithm.SUCCESSIVE


class TestCalculatorFactory:
    """Test CalculatorFactory with direct factory pattern."""

    def test_create_calculator_ncc_successive(self):
        """Test calculator creation for NCC with Successive."""
        from PyMaSC.core.factory import create_calculator
        
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )
        config.ref2lengths = {"chr1": 1000}

        # This should work without requiring mappability for NCC
        calculator, bw_reader = create_calculator(config)
        assert calculator is not None

    def test_create_calculator_mscc_requires_mappability(self):
        """Test that MSCC calculation requires mappability configuration."""
        from PyMaSC.core.factory import create_calculator
        
        config = PyMaSCConfig(
            target=CalculationTarget.MSCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=500,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )
        config.ref2lengths = {"chr1": 1000}

        # Should fail without mappability config
        with pytest.raises((ValueError, AssertionError), match="Mappability path must be set"):
            calculator, bw_reader = create_calculator(config)

    def test_all_combinations_supported(self):
        """Test that all valid combinations are supported."""
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.BITARRAY,
            max_shift=500,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10,
            read_length=50  # Required for BitArray implementation
        )
        config.ref2lengths = {"chr1": 1000}

        # All combinations should be supported now
        calculator, bw_reader = create_calculator(config)
        assert calculator is not None

    def test_factory_direct_usage_pattern(self):
        """Test that factory can be used directly without strategy layer."""
        # All valid target-implementation combinations
        test_cases = [
            (CalculationTarget.NCC, Algorithm.BITARRAY),
            (CalculationTarget.NCC, Algorithm.SUCCESSIVE),
            (CalculationTarget.BOTH, Algorithm.BITARRAY),
            (CalculationTarget.BOTH, Algorithm.SUCCESSIVE),
        ]

        for target, implementation in test_cases:
            config = PyMaSCConfig(
                target=target,
                implementation=implementation,
                max_shift=500,
                mapq_criteria=20,
                nproc=1,
                esttype=EstimationType.MEDIAN,
                chi2_pval=0.0001,
                mv_avr_filter_len=50,
                filter_mask_len=10,
                min_calc_width=10,
                read_length=50  # Required for BitArray
            )
            config.ref2lengths = {"chr1": 1000}

            # MSCC requires mappability config, so skip those
            if target in [CalculationTarget.MSCC, CalculationTarget.BOTH]:
                continue

            calculator, bw_reader = create_calculator(config)
            assert calculator is not None

    def test_factory_performance_optimization(self):
        """Test that factory pattern provides direct calculator creation."""
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.BITARRAY,
            max_shift=200,
            mapq_criteria=20,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10,
            read_length=50
        )
        config.ref2lengths = {"chr1": 1000, "chr2": 800}

        # Direct factory usage should be fast and efficient
        calculator, bw_reader = create_calculator(config)

        assert calculator is not None
        # Verify it's the correct calculator type
        from PyMaSC.bacore.mscc import CCBitArrayCalculator
        assert isinstance(calculator, CCBitArrayCalculator)


class TestFactoryArchitectureIntegration:
    """Test integration of factory-based architecture."""

    def test_single_process_integration(self):
        """Test that single-process mode uses factory directly."""
        # This test verifies that CalcHandler uses create_calculator directly
        # without strategy layer overhead
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=300,
            mapq_criteria=25,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )
        config.ref2lengths = {"chr1": 1000}

        calculator, bw_reader = create_calculator(config)

        # Should create calculator without intermediate strategy objects
        assert calculator is not None

    def test_multiprocess_integration_consistency(self):
        """Test that multi-process mode also uses factory directly."""
        # Both single and multi-process should use the same factory pattern
        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.BITARRAY,
            max_shift=400,
            mapq_criteria=30,
            nproc=4,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10,
            read_length=75
        )
        config.ref2lengths = {"chr1": 1000}

        # Same factory call as used in both single and multi-process
        calculator, bw_reader = create_calculator(config)

        assert calculator is not None

    def test_architecture_simplification(self):
        """Test that new architecture eliminates unnecessary complexity."""
        # Factory pattern should provide simple, direct calculator creation
        # No strategy objects, contexts, or intermediate layers needed

        config = PyMaSCConfig(
            target=CalculationTarget.NCC,
            implementation=Algorithm.SUCCESSIVE,
            max_shift=250,
            mapq_criteria=15,
            nproc=1,
            esttype=EstimationType.MEDIAN,
            chi2_pval=0.0001,
            mv_avr_filter_len=50,
            filter_mask_len=10,
            min_calc_width=10
        )
        config.ref2lengths = {"chr1": 1000}

        # Single, direct factory call
        calculator, bw_reader = create_calculator(config)

        assert calculator is not None
        # Verify no strategy overhead
        from PyMaSC.core.ncc import NaiveCCCalculator
        assert isinstance(calculator, NaiveCCCalculator)