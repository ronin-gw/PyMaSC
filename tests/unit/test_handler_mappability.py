"""Tests for PyMaSC.handler.mappability module.

This module tests the mappability handler functionality including
error handling for BigWig file access, JSON statistics validation,
and worker process management.
"""
import unittest
from unittest.mock import Mock, patch, mock_open, MagicMock
import logging
import json
import tempfile
import os
from pathlib import Path
import numpy as np

from PyMaSC.handler.mappability import (
    MappabilityHandler, MappabilityCalcWorker, NumpyEncoder,
    BWIOError, JSONIOError, NeedUpdate
)


class TestBWIOError(unittest.TestCase):
    """Test BWIOError exception class."""

    def test_bwio_error_inheritance(self):
        """Test that BWIOError inherits from IOError."""
        self.assertTrue(issubclass(BWIOError, IOError))

    def test_bwio_error_instantiation(self):
        """Test BWIOError can be instantiated and raised."""
        error = BWIOError("Test error message")
        self.assertIsInstance(error, BWIOError)
        self.assertIsInstance(error, IOError)
        
        with self.assertRaises(BWIOError):
            raise error


class TestJSONIOError(unittest.TestCase):
    """Test JSONIOError exception class."""

    def test_json_io_error_inheritance(self):
        """Test that JSONIOError inherits from IOError."""
        self.assertTrue(issubclass(JSONIOError, IOError))

    def test_json_io_error_instantiation(self):
        """Test JSONIOError can be instantiated and raised."""
        error = JSONIOError("Test JSON error")
        self.assertIsInstance(error, JSONIOError)
        self.assertIsInstance(error, IOError)


class TestNeedUpdate(unittest.TestCase):
    """Test NeedUpdate exception class."""

    def test_need_update_inheritance(self):
        """Test that NeedUpdate inherits from Exception."""
        self.assertTrue(issubclass(NeedUpdate, Exception))

    def test_need_update_instantiation(self):
        """Test NeedUpdate can be instantiated and raised."""
        error = NeedUpdate("Update required")
        self.assertIsInstance(error, NeedUpdate)
        
        with self.assertRaises(NeedUpdate):
            raise error


class TestNumpyEncoder(unittest.TestCase):
    """Test NumpyEncoder JSON serialization."""

    def setUp(self):
        """Set up test fixtures."""
        self.encoder = NumpyEncoder()

    def test_numpy_encoder_handles_numpy_integer(self):
        """Test encoding of numpy integer types."""
        np_int = np.int32(42)
        result = self.encoder.default(np_int)
        self.assertEqual(result, 42)
        self.assertIsInstance(result, int)
        self.assertNotIsInstance(result, np.integer)

    def test_numpy_encoder_handles_numpy_float(self):
        """Test encoding of numpy float types."""
        np_float = np.float64(3.14159)
        result = self.encoder.default(np_float)
        self.assertAlmostEqual(result, 3.14159)
        self.assertIsInstance(result, float)
        self.assertNotIsInstance(result, np.floating)

    def test_numpy_encoder_handles_various_numpy_types(self):
        """Test encoding of various numpy numeric types."""
        test_values = [
            (np.int8(10), 10),
            (np.int16(100), 100),
            (np.int32(1000), 1000),
            (np.int64(10000), 10000),
            (np.float32(1.5), 1.5),
            (np.float64(2.5), 2.5),
        ]
        
        for np_value, expected in test_values:
            with self.subTest(np_type=type(np_value).__name__):
                result = self.encoder.default(np_value)
                self.assertEqual(result, expected)

    def test_numpy_encoder_raises_for_unsupported_types(self):
        """Test that unsupported types raise TypeError."""
        unsupported_obj = object()
        
        with self.assertRaises(TypeError):
            self.encoder.default(unsupported_obj)

    def test_numpy_encoder_json_serialization_integration(self):
        """Test full JSON serialization with numpy types."""
        data = {
            "int_value": np.int32(42),
            "float_value": np.float64(3.14159),
            "list_with_numpy": [np.int32(1), np.float64(2.5)],
            "regular_types": {"str": "test", "int": 100}
        }
        
        json_str = json.dumps(data, cls=NumpyEncoder)
        loaded_data = json.loads(json_str)
        
        # Verify the data was properly serialized and deserialized
        self.assertEqual(loaded_data["int_value"], 42)
        self.assertAlmostEqual(loaded_data["float_value"], 3.14159)
        self.assertEqual(loaded_data["list_with_numpy"], [1, 2.5])
        self.assertEqual(loaded_data["regular_types"], {"str": "test", "int": 100})




class TestMappabilityHandlerStaticMethods(unittest.TestCase):
    """Test static methods of MappabilityHandler."""

    def test_calc_mappable_len_required_shift_size_normal_case(self):
        """Test shift size calculation for normal case."""
        readlen = 36
        max_shift = 500
        result = MappabilityHandler.calc_mappable_len_required_shift_size(readlen, max_shift)
        expected = max_shift - readlen + 1  # 500 - 36 + 1 = 465
        self.assertEqual(result, expected)

    def test_calc_mappable_len_required_shift_size_small_shift(self):
        """Test shift size calculation when max_shift is small."""
        readlen = 36
        max_shift = 50  # Less than 2*readlen - 1 = 71
        result = MappabilityHandler.calc_mappable_len_required_shift_size(readlen, max_shift)
        expected = readlen  # Should return readlen
        self.assertEqual(result, expected)

    def test_calc_mappable_len_required_shift_size_edge_case(self):
        """Test shift size calculation at the boundary."""
        readlen = 36
        max_shift = 71  # Exactly 2*readlen - 1
        result = MappabilityHandler.calc_mappable_len_required_shift_size(readlen, max_shift)
        expected = readlen  # Should return readlen since max_shift is not > 2*readlen - 1
        self.assertEqual(result, expected)

    def test_from_config_method_exists(self):
        """Test that from_config class method exists."""
        # Just verify the method exists and is callable
        self.assertTrue(hasattr(MappabilityHandler, 'from_config'))
        self.assertTrue(callable(getattr(MappabilityHandler, 'from_config')))


if __name__ == '__main__':
    unittest.main()