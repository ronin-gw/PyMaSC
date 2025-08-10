"""Tests for PyMaSC.reader.stats module.

This module tests the statistics file reading and parsing functionality
including error handling for malformed data and missing fields.
"""
import unittest
from unittest.mock import Mock, patch, mock_open
import logging
import tempfile
import os
from pathlib import Path
from typing import Union, Literal
from dataclasses import dataclass

from PyMaSC.reader.stats import (
    _parse_stats_section, _parse_value, _get_field_type, _is_nanfield,
    _read_stats, load_stats
)
from PyMaSC.interfaces.output import (
    SUMMARY_LABELS, NCC_LABELS, MSCC_LABELS,
    SummaryItems, CorrelationItems, OutputStats,
    StatLabels, nanint, nanfloat
)


class TestParseStatsSection(unittest.TestCase):
    """Test _parse_stats_section function."""

    def setUp(self):
        """Set up test fixtures."""
        @dataclass
        class TestItems:
            required_field: str
            optional_nan_field: nanint
            regular_int_field: int
            
        @dataclass  
        class TestLabels:
            required_field: str = "Required Field"
            optional_nan_field: str = "Optional Field"
            regular_int_field: str = "Regular Field"

        self.test_dataclass = TestItems
        self.test_labels = TestLabels()

    def test_parse_stats_section_success_all_fields(self):
        """Test successful parsing when all fields are present."""
        raw_data = {
            "Required Field": "test_value",
            "Optional Field": "42",
            "Regular Field": "100"
        }
        
        result = _parse_stats_section(raw_data, self.test_labels, self.test_dataclass)
        
        self.assertEqual(result["required_field"], "test_value")
        self.assertEqual(result["optional_nan_field"], 42)
        self.assertEqual(result["regular_int_field"], 100)

    def test_parse_stats_section_missing_required_field_raises_error(self):
        """Test that missing required fields raise ValueError."""
        raw_data = {
            "Optional Field": "42",
            "Regular Field": "100"
            # Missing "Required Field"
        }
        
        with self.assertRaises(ValueError) as cm:
            _parse_stats_section(raw_data, self.test_labels, self.test_dataclass)
        
        self.assertIn("Required", str(cm.exception))
        self.assertIn("required_field", str(cm.exception))
        self.assertIn("missing", str(cm.exception))

    @patch('PyMaSC.reader.stats.logger')
    def test_parse_stats_section_missing_nan_field_logs_warning(self, mock_logger):
        """Test that missing nan-annotated fields log warning and set to 'nan'."""
        raw_data = {
            "Required Field": "test_value",
            # Missing "Optional Field" (nan field)
            "Regular Field": "100"
        }
        
        result = _parse_stats_section(raw_data, self.test_labels, self.test_dataclass)
        
        # Should set nan field to 'nan'
        self.assertEqual(result["optional_nan_field"], 'nan')
        # Should log warning
        mock_logger.warning.assert_called_once()
        warning_msg = mock_logger.warning.call_args[0][0]
        self.assertIn("Missing", warning_msg)
        self.assertIn("optional_nan_field", warning_msg)
        self.assertIn("nan", warning_msg)

    def test_parse_stats_section_with_section_name_in_error(self):
        """Test that section name appears in error messages."""
        raw_data = {}  # Empty data, missing required field
        
        with self.assertRaises(ValueError) as cm:
            _parse_stats_section(raw_data, self.test_labels, self.test_dataclass, "TEST_SECTION")
        
        self.assertIn("TEST_SECTION", str(cm.exception))

    @patch('PyMaSC.reader.stats.logger')
    def test_parse_stats_section_with_section_name_in_warning(self, mock_logger):
        """Test that section name appears in warning messages."""
        raw_data = {
            "Required Field": "test_value",
            "Regular Field": "100"
            # Missing optional nan field
        }
        
        _parse_stats_section(raw_data, self.test_labels, self.test_dataclass, "TEST_SECTION")
        
        mock_logger.warning.assert_called_once()
        warning_msg = mock_logger.warning.call_args[0][0]
        self.assertIn("TEST_SECTION", warning_msg)


class TestParseValue(unittest.TestCase):
    """Test _parse_value function."""

    def test_parse_value_nan_string_returns_nan(self):
        """Test that 'nan' string is returned as-is."""
        result = _parse_value("nan", "test_field", int)
        self.assertEqual(result, "nan")

    def test_parse_value_valid_int_conversion(self):
        """Test successful integer conversion."""
        result = _parse_value("42", "test_field", int)
        self.assertEqual(result, 42)
        self.assertIsInstance(result, int)

    def test_parse_value_valid_float_conversion(self):
        """Test successful float conversion."""
        result = _parse_value("3.14159", "test_field", float)
        self.assertAlmostEqual(result, 3.14159)
        self.assertIsInstance(result, float)

    def test_parse_value_valid_string_conversion(self):
        """Test string values pass through unchanged."""
        result = _parse_value("test_string", "test_field", str)
        self.assertEqual(result, "test_string")
        self.assertIsInstance(result, str)

    @patch('PyMaSC.reader.stats.logger')
    def test_parse_value_invalid_int_conversion_returns_nan(self, mock_logger):
        """Test that invalid int conversion returns 'nan' and logs warning."""
        result = _parse_value("not_a_number", "test_field", int)
        
        self.assertEqual(result, "nan")
        mock_logger.warning.assert_called_once()
        warning_msg = mock_logger.warning.call_args[0][0]
        self.assertIn("Failed to parse value", warning_msg)
        self.assertIn("not_a_number", warning_msg)
        self.assertIn("test_field", warning_msg)

    @patch('PyMaSC.reader.stats.logger')
    def test_parse_value_invalid_float_conversion_returns_nan(self, mock_logger):
        """Test that invalid float conversion returns 'nan' and logs warning."""
        result = _parse_value("invalid_float", "float_field", float)
        
        self.assertEqual(result, "nan")
        mock_logger.warning.assert_called_once()
        warning_msg = mock_logger.warning.call_args[0][0]
        self.assertIn("Failed to parse value", warning_msg)
        self.assertIn("invalid_float", warning_msg)
        self.assertIn("float_field", warning_msg)

    @patch('PyMaSC.reader.stats.logger')
    def test_parse_value_type_error_returns_nan(self, mock_logger):
        """Test that TypeError during conversion returns 'nan'."""
        # Simulate a type that raises TypeError
        def bad_type(value):
            raise TypeError("Cannot convert")
        
        result = _parse_value("test", "test_field", bad_type)
        
        self.assertEqual(result, "nan")
        mock_logger.warning.assert_called_once()


class TestGetFieldType(unittest.TestCase):
    """Test _get_field_type function."""

    def test_get_field_type_simple_annotation(self):
        """Test getting type from simple annotation."""
        @dataclass
        class SimpleClass:
            int_field: int
            str_field: str
            float_field: float
            
        self.assertEqual(_get_field_type("int_field", SimpleClass), int)
        self.assertEqual(_get_field_type("str_field", SimpleClass), str)
        self.assertEqual(_get_field_type("float_field", SimpleClass), float)

    def test_get_field_type_union_annotation_nanint(self):
        """Test getting type from Union annotation (nanint)."""
        @dataclass
        class UnionClass:
            nan_int_field: nanint  # Union[int, Literal["nan"]]
            
        result = _get_field_type("nan_int_field", UnionClass)
        self.assertEqual(result, int)

    def test_get_field_type_union_annotation_nanfloat(self):
        """Test getting type from Union annotation (nanfloat)."""
        @dataclass
        class UnionClass:
            nan_float_field: nanfloat  # Union[float, Literal["nan"]]
            
        result = _get_field_type("nan_float_field", UnionClass)
        self.assertEqual(result, float)

    def test_get_field_type_missing_field_returns_str(self):
        """Test that missing field returns str as default."""
        @dataclass
        class TestClass:
            existing_field: int
            
        result = _get_field_type("nonexistent_field", TestClass)
        self.assertEqual(result, str)

    def test_get_field_type_complex_union_with_none(self):
        """Test Union type handling with None."""
        @dataclass
        class ComplexClass:
            optional_int_field: Union[int, None]
            
        result = _get_field_type("optional_int_field", ComplexClass)
        self.assertEqual(result, int)


class TestIsNanfield(unittest.TestCase):
    """Test _is_nanfield function."""

    def test_is_nanfield_nanint_returns_true(self):
        """Test that nanint fields are identified correctly."""
        @dataclass
        class TestClass:
            nan_int_field: nanint
            
        result = _is_nanfield("nan_int_field", TestClass)
        self.assertTrue(result)

    def test_is_nanfield_nanfloat_returns_true(self):
        """Test that nanfloat fields are identified correctly."""
        @dataclass
        class TestClass:
            nan_float_field: nanfloat
            
        result = _is_nanfield("nan_float_field", TestClass)
        self.assertTrue(result)

    def test_is_nanfield_regular_int_returns_false(self):
        """Test that regular int fields return False."""
        @dataclass
        class TestClass:
            regular_int_field: int
            
        result = _is_nanfield("regular_int_field", TestClass)
        self.assertFalse(result)

    def test_is_nanfield_regular_float_returns_false(self):
        """Test that regular float fields return False."""
        @dataclass
        class TestClass:
            regular_float_field: float
            
        result = _is_nanfield("regular_float_field", TestClass)
        self.assertFalse(result)

    def test_is_nanfield_str_returns_false(self):
        """Test that string fields return False."""
        @dataclass
        class TestClass:
            str_field: str
            
        result = _is_nanfield("str_field", TestClass)
        self.assertFalse(result)

    def test_is_nanfield_missing_field_returns_false(self):
        """Test that missing field returns False."""
        @dataclass
        class TestClass:
            existing_field: int
            
        result = _is_nanfield("nonexistent_field", TestClass)
        self.assertFalse(result)


class TestReadStats(unittest.TestCase):
    """Test _read_stats function."""

    @patch('PyMaSC.reader.stats.logger')
    def test_read_stats_success(self, mock_logger):
        """Test successful stats file reading."""
        file_content = "Field1\tValue1\nField2\tValue2\nField3\tValue3\n"
        
        with patch("builtins.open", mock_open(read_data=file_content)):
            result = _read_stats("/test/path")
        
        expected = {
            "Field1": "Value1",
            "Field2": "Value2", 
            "Field3": "Value3"
        }
        self.assertEqual(result, expected)
        
        # Check that info message was logged
        mock_logger.info.assert_called_once()
        info_msg = mock_logger.info.call_args[0][0]
        self.assertIn("Load statistics", info_msg)
        self.assertIn("/test/path", info_msg)

    @patch('PyMaSC.reader.stats.logger')
    def test_read_stats_with_pathlib_path(self, mock_logger):
        """Test reading stats with pathlib.Path object."""
        file_content = "Name\tTest\n"
        path_obj = Path("/test/path.txt")
        
        with patch("builtins.open", mock_open(read_data=file_content)):
            result = _read_stats(path_obj)
        
        expected = {"Name": "Test"}
        self.assertEqual(result, expected)

    def test_read_stats_empty_file(self):
        """Test reading empty stats file."""
        with patch("builtins.open", mock_open(read_data="")):
            result = _read_stats("/test/empty")
        
        self.assertEqual(result, {})

    def test_read_stats_single_column_lines_cause_error(self):
        """Test that lines with single column cause ValueError."""
        file_content = "Field1\tValue1\nSingleColumnLine\nField2\tValue2\n"
        
        with patch("builtins.open", mock_open(read_data=file_content)):
            with self.assertRaises(ValueError) as cm:
                _read_stats("/test/path")
        
        # Should raise ValueError because dict() expects key-value pairs
        self.assertIn("dictionary update sequence", str(cm.exception))
        self.assertIn("length 1", str(cm.exception))


class TestLoadStats(unittest.TestCase):
    """Test load_stats function integration."""

    def test_load_stats_complete_file_success(self):
        """Test successful loading of complete stats file."""
        # Create a complete stats file content
        file_content = """Name\tTest Sample
Read length\t36
Expected library length\t150
Estimated library length\t65
Genome length\t3137454505
Forward reads\t622
Reverse reads\t670
Minimum NCC\t0.020135
NCC at read length\t0.117726
NCC at expected library length\tnan
NCC at estimated library length\t0.131668
NSC\tnan
RSC\tnan
Estimated NSC\t6.539169
Estimated RSC\t1.142857
FWHM\tnan
VSN\tnan
Estimated FWHM\t170
Estimated VSN\t0.034649
DMP length\t19906
Forward reads in DMP\t584
Reverse reads in DMP\t631
Minimum MSCC\t0.019834
MSCC at read length\t0.116031
MSCC at expected library length\tnan
MSCC at estimated library length\t0.129771
MSCC NSC\tnan
MSCC RSC\tnan
Estimated MSCC NSC\t6.544234
Estimated MSCC RSC\t1.145955
MSCC FWHM\tnan
MSCC VSN\tnan
Estimated MSCC FWHM\t168
Estimated MSCC VSN\t0.034122
"""

        with patch("builtins.open", mock_open(read_data=file_content)):
            result = load_stats("/test/complete_stats.txt")

        # Verify the structure
        self.assertIsInstance(result, OutputStats)
        
        # Check summary stats
        self.assertEqual(result.stats.name, "Test Sample")
        self.assertEqual(result.stats.read_len, 36)
        self.assertEqual(result.stats.expected_lib_len, 150)
        self.assertEqual(result.stats.est_lib_len, 65)
        
        # Check NCC stats
        self.assertEqual(result.ncc_stats.stats.genomelen, 3137454505)
        self.assertEqual(result.ncc_stats.stats.forward_reads, 622)
        self.assertEqual(result.ncc_stats.stats.reverse_reads, 670)
        self.assertAlmostEqual(result.ncc_stats.stats.min_cc, 0.020135, places=5)
        self.assertEqual(result.ncc_stats.stats.ccfl, "nan")  # Expected library length
        
        # Check MSCC stats  
        self.assertEqual(result.mscc_stats.stats.genomelen, 19906)  # DMP length
        self.assertEqual(result.mscc_stats.stats.forward_reads, 584)
        self.assertEqual(result.mscc_stats.stats.reverse_reads, 631)

    def test_load_stats_missing_required_summary_field(self):
        """Test error when required summary field is missing."""
        file_content = """Read length\t36
Expected library length\t150
Estimated library length\t65
"""  # Missing "Name" field

        with patch("builtins.open", mock_open(read_data=file_content)):
            with self.assertRaises(ValueError) as cm:
                load_stats("/test/incomplete.txt")
        
        self.assertIn("Required", str(cm.exception))
        self.assertIn("name", str(cm.exception))

    def test_load_stats_missing_required_summary_field_read_len(self):
        """Test error when required summary field 'Read length' is missing."""
        file_content = """Name\tTest
Expected library length\t150
Estimated library length\t65
Genome length\t3000000000
Forward reads\t622
Reverse reads\t670
Minimum NCC\t0.02
NCC at read length\t0.11
NCC at expected library length\tnan
NCC at estimated library length\t0.13
NSC\tnan
RSC\tnan
Estimated NSC\t6.5
Estimated RSC\t1.1
FWHM\tnan
VSN\tnan
Estimated FWHM\t170
Estimated VSN\t0.03
DMP length\t19906
Forward reads in DMP\t584
Reverse reads in DMP\t631
Minimum MSCC\t0.019
MSCC at read length\t0.116
MSCC at expected library length\tnan
MSCC at estimated library length\t0.129
MSCC NSC\tnan
MSCC RSC\tnan
Estimated MSCC NSC\t6.5
Estimated MSCC RSC\t1.1
MSCC FWHM\tnan
MSCC VSN\tnan
Estimated MSCC FWHM\t168
Estimated MSCC VSN\t0.034
"""  # Missing required "Read length" field

        with patch("builtins.open", mock_open(read_data=file_content)):
            with self.assertRaises(ValueError) as cm:
                load_stats("/test/incomplete_summary.txt")
        
        self.assertIn("Required", str(cm.exception))
        self.assertIn("read_len", str(cm.exception))

    @patch('PyMaSC.reader.stats.logger')
    def test_load_stats_missing_optional_ncc_fields_log_warnings(self, mock_logger):
        """Test that missing NCC fields log warnings but don't raise errors (all are nan-annotated)."""
        file_content = """Name\tTest
Read length\t36
Expected library length\t150
Estimated library length\t65
DMP length\t19906
Forward reads in DMP\t584
Reverse reads in DMP\t631
Minimum MSCC\t0.019
MSCC at read length\t0.116
MSCC at expected library length\tnan
MSCC at estimated library length\t0.129
MSCC NSC\tnan
MSCC RSC\tnan
Estimated MSCC NSC\t6.5
Estimated MSCC RSC\t1.1
MSCC FWHM\tnan
MSCC VSN\tnan
Estimated MSCC FWHM\t168
Estimated MSCC VSN\t0.034
"""  # Missing most NCC fields

        with patch("builtins.open", mock_open(read_data=file_content)):
            result = load_stats("/test/missing_ncc.txt")

        # Should succeed but log warnings for missing NCC fields
        self.assertIsInstance(result, OutputStats)
        
        # Check that warnings were logged for missing NCC fields
        warning_calls = [call for call in mock_logger.warning.call_args_list 
                        if "Missing NCC field" in str(call)]
        self.assertGreater(len(warning_calls), 0, "Should have logged warnings for missing NCC fields")
        
        # All NCC fields should be set to 'nan'
        self.assertEqual(result.ncc_stats.stats.genomelen, 'nan')
        self.assertEqual(result.ncc_stats.stats.forward_reads, 'nan')
        self.assertEqual(result.ncc_stats.stats.min_cc, 'nan')

    @patch('PyMaSC.reader.stats.logger')
    def test_load_stats_with_invalid_numeric_values(self, mock_logger):
        """Test handling of invalid numeric values in stats file."""
        file_content = """Name\tTest
Read length\tinvalid_number
Expected library length\t150
Estimated library length\t65
Genome length\t3000000000
Forward reads\t622
Reverse reads\t670
Minimum NCC\t0.02
NCC at read length\t0.11
NCC at expected library length\tnan
NCC at estimated library length\t0.13
NSC\tnan
RSC\tnan
Estimated NSC\t6.5
Estimated RSC\t1.1
FWHM\tnan
VSN\tnan
Estimated FWHM\t170
Estimated VSN\t0.03
DMP length\t19906
Forward reads in DMP\t584
Reverse reads in DMP\t631
Minimum MSCC\t0.019
MSCC at read length\t0.116
MSCC at expected library length\tnan
MSCC at estimated library length\t0.129
MSCC NSC\tnan
MSCC RSC\tnan
Estimated MSCC NSC\t6.5
Estimated MSCC RSC\t1.1
MSCC FWHM\tnan
MSCC VSN\tnan
Estimated MSCC FWHM\t168
Estimated MSCC VSN\t0.034
"""

        with patch("builtins.open", mock_open(read_data=file_content)):
            result = load_stats("/test/invalid_numbers.txt")

        # Should still succeed but with 'nan' for invalid values
        self.assertEqual(result.stats.read_len, "nan")
        mock_logger.warning.assert_called()

    def test_load_stats_file_io_error(self):
        """Test that IOError from file reading is properly handled."""
        # The @catch_IOError decorator should handle IOError and re-raise it
        # after logging
        with patch("builtins.open", side_effect=IOError("File not found")):
            with self.assertRaises(IOError):
                load_stats("/nonexistent/file.txt")

    def test_load_stats_malformed_file_format(self):
        """Test handling of malformed file format."""
        file_content = "This is not a tab-delimited file\nNo tabs here either\n"
        
        with patch("builtins.open", mock_open(read_data=file_content)):
            with self.assertRaises(ValueError):
                # Should fail when trying to parse required fields
                load_stats("/test/malformed.txt")


if __name__ == '__main__':
    unittest.main()