"""Tests for interface modules in PyMaSC.core.interfaces.

This module tests the abstract base classes and data models that define
the contracts for PyMaSC components.
"""
import unittest
from typing import Any
from dataclasses import dataclass
import numpy as np

from PyMaSC.interfaces.result import (
    AbstractDataclass, CorrelationResultModel, NCCResultModel, MSCCResultModel,
    BothChromResultModel, NCCGenomeWideResultModel, MSCCGenomeWideResultModel,
    BothGenomeWideResultModel
)


class TestAbstractDataclass(unittest.TestCase):
    """Test AbstractDataclass behavior."""

    def test_abstract_dataclass_cannot_be_instantiated(self):
        """Test that AbstractDataclass cannot be instantiated directly."""
        with self.assertRaises(TypeError) as cm:
            AbstractDataclass()
        self.assertIn("Cannot instantiate abstract class", str(cm.exception))

    def test_abstract_dataclass_subclass_cannot_be_instantiated(self):
        """Test that direct subclasses of AbstractDataclass cannot be instantiated."""
        @dataclass
        class DirectSubclass(AbstractDataclass):
            value: int = 1

        with self.assertRaises(TypeError) as cm:
            DirectSubclass()
        self.assertIn("Cannot instantiate abstract class", str(cm.exception))

    def test_abstract_dataclass_concrete_subclass_can_be_instantiated(self):
        """Test that concrete subclasses can be instantiated."""
        @dataclass
        class ConcreteBase(AbstractDataclass):
            pass

        @dataclass
        class ConcreteSubclass(ConcreteBase):
            value: int = 1

        # Should not raise
        instance = ConcreteSubclass()
        self.assertEqual(instance.value, 1)


class TestResultModels(unittest.TestCase):
    """Test result model interfaces."""

    def test_ncc_result_model_structure(self):
        """Test NCCResultModel structure."""
        # NCCResultModel should have specific integer sum types
        cc_data = np.zeros(5, dtype=np.float64)
        
        result = NCCResultModel(
            max_shift=4,
            read_len=36,
            genomelen=3000000000,
            forward_sum=50000,  # int
            reverse_sum=48000,  # int
            cc=cc_data
        )
        
        self.assertEqual(result.max_shift, 4)
        self.assertEqual(result.read_len, 36)
        self.assertEqual(result.genomelen, 3000000000)
        self.assertEqual(result.forward_sum, 50000)
        self.assertEqual(result.reverse_sum, 48000)
        self.assertIs(result.cc, cc_data)

    def test_mscc_result_model_structure(self):
        """Test MSCCResultModel structure."""
        # MSCCResultModel should have array sum types
        cc_data = np.zeros(5, dtype=np.float64)
        forward_array = np.array([5, 10, 15, 20, 25], dtype=np.int64)
        reverse_array = np.array([4, 8, 12, 16, 20], dtype=np.int64)
        mappable_len = (100, 200, 300, 400, 500)
        
        result = MSCCResultModel(
            max_shift=4,
            read_len=5,
            genomelen=3000000000,
            forward_sum=forward_array,  # array
            reverse_sum=reverse_array,  # array
            cc=cc_data,
            mappable_len=mappable_len
        )
        
        self.assertEqual(result.max_shift, 4)
        self.assertEqual(result.read_len, 5)
        self.assertEqual(result.genomelen, 3000000000)
        np.testing.assert_array_equal(result.forward_sum, forward_array)
        np.testing.assert_array_equal(result.reverse_sum, reverse_array)
        self.assertEqual(result.mappable_len, mappable_len)

    def test_both_chrom_result_model_structure(self):
        """Test BothChromResultModel structure."""
        # Create mock NCC and MSCC results
        ncc_result = NCCResultModel(
            max_shift=4,
            read_len=36,
            genomelen=1000000,
            forward_sum=500,
            reverse_sum=480,
            cc=np.zeros(5, dtype=np.float64)
        )
        
        mscc_result = MSCCResultModel(
            max_shift=4,
            read_len=5,
            genomelen=1000000,
            forward_sum=np.array([5, 10, 15, 20, 25], dtype=np.int64),
            reverse_sum=np.array([4, 8, 12, 16, 20], dtype=np.int64),
            cc=np.zeros(5, dtype=np.float64),
            mappable_len=(100, 200, 300, 400, 500)
        )
        
        both_result = BothChromResultModel(
            chrom=ncc_result,
            mappable_chrom=mscc_result
        )
        
        self.assertEqual(both_result.chrom, ncc_result)
        self.assertEqual(both_result.mappable_chrom, mscc_result)

    def test_genome_wide_result_models(self):
        """Test genome-wide result model structures."""
        # Test NCCGenomeWideResultModel
        ncc_genome = NCCGenomeWideResultModel(
            genomelen=3000000000,
            forward_sum=50000,
            reverse_sum=48000,
            chroms={}  # Empty mapping for testing
        )
        
        self.assertEqual(ncc_genome.genomelen, 3000000000)
        self.assertEqual(ncc_genome.forward_sum, 50000)
        self.assertEqual(ncc_genome.reverse_sum, 48000)
        self.assertEqual(ncc_genome.chroms, {})
        
        # Test MSCCGenomeWideResultModel
        mscc_genome = MSCCGenomeWideResultModel(
            genomelen=3000000000,
            chroms={}
        )
        
        self.assertEqual(mscc_genome.genomelen, 3000000000)
        self.assertEqual(mscc_genome.chroms, {})
        
        # Test BothGenomeWideResultModel
        both_genome = BothGenomeWideResultModel(
            genomelen=3000000000,
            forward_sum=50000,
            reverse_sum=48000,
            chroms={},
            mappable_chroms={}
        )
        
        self.assertEqual(both_genome.genomelen, 3000000000)
        self.assertEqual(both_genome.forward_sum, 50000)
        self.assertEqual(both_genome.reverse_sum, 48000)
        self.assertEqual(both_genome.chroms, {})
        self.assertEqual(both_genome.mappable_chroms, {})


class TestInterfaceCompatibility(unittest.TestCase):
    """Test interface compatibility and type annotations."""

    def test_correlation_result_model_type_annotations(self):
        """Test CorrelationResultModel type annotations through concrete subclasses."""
        # CorrelationResultModel is abstract, so test through concrete implementations
        
        # NCCResultModel accepts int values for forward_sum and reverse_sum
        ncc_result = NCCResultModel(
            max_shift=4,
            read_len=36,
            genomelen=1000000,
            forward_sum=500,  # int
            reverse_sum=480,  # int
            cc=np.zeros(5, dtype=np.float64)
        )
        self.assertEqual(ncc_result.forward_sum, 500)
        self.assertEqual(ncc_result.reverse_sum, 480)
        
        # MSCCResultModel accepts array values for forward_sum and reverse_sum
        mscc_result = MSCCResultModel(
            max_shift=4,
            read_len=36,
            genomelen=1000000,
            forward_sum=np.array([5, 10, 15], dtype=np.int64),  # array
            reverse_sum=np.array([4, 8, 12], dtype=np.int64),  # array
            cc=np.zeros(5, dtype=np.float64),
            mappable_len=None
        )
        np.testing.assert_array_equal(mscc_result.forward_sum, np.array([5, 10, 15]))
        np.testing.assert_array_equal(mscc_result.reverse_sum, np.array([4, 8, 12]))

    def test_optional_mappable_len(self):
        """Test that mappable_len is optional in MSCCResultModel."""
        # Should work with mappable_len
        result1 = MSCCResultModel(
            max_shift=4,
            read_len=5,
            genomelen=1000000,
            forward_sum=np.array([5, 10], dtype=np.int64),
            reverse_sum=np.array([4, 8], dtype=np.int64),
            cc=np.zeros(5, dtype=np.float64),
            mappable_len=(100, 200)
        )
        self.assertEqual(result1.mappable_len, (100, 200))
        
        # Should work without mappable_len
        result2 = MSCCResultModel(
            max_shift=4,
            read_len=5,
            genomelen=1000000,
            forward_sum=np.array([5, 10], dtype=np.int64),
            reverse_sum=np.array([4, 8], dtype=np.int64),
            cc=np.zeros(5, dtype=np.float64),
            mappable_len=None
        )
        self.assertIsNone(result2.mappable_len)


if __name__ == '__main__':
    unittest.main()