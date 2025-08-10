"""Tests for PyMaSC.utils.output utilities.

This module tests the file output utilities including error handling
decorators and directory preparation functions.
"""
import unittest
from unittest.mock import Mock, patch, MagicMock
import logging
import os
import tempfile
import shutil
from pathlib import Path

from PyMaSC.utils.output import catch_IOError, prepare_outdir


class TestCatchIOError(unittest.TestCase):
    """Test catch_IOError decorator functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.mock_logger = Mock(spec=logging.Logger)
        self.decorator = catch_IOError(self.mock_logger)

    def test_catch_ioerror_successful_function_call(self):
        """Test that successful function calls pass through unchanged."""
        @self.decorator
        def dummy_function(x, y):
            return x + y

        result = dummy_function(2, 3)
        self.assertEqual(result, 5)
        self.mock_logger.error.assert_not_called()

    def test_catch_ioerror_handles_ioerror_with_filename_and_errno(self):
        """Test IOError handling with filename and errno."""
        @self.decorator
        def failing_function():
            error = IOError("Permission denied")
            error.filename = "/test/path.txt"
            error.errno = 13
            raise error

        with self.assertRaises(IOError):
            failing_function()

        self.mock_logger.error.assert_called_once()
        call_args = self.mock_logger.error.call_args[0][0]
        self.assertIn("Faild to output '/test/path.txt'", call_args)
        self.assertIn("[Errno 13]", call_args)

    def test_catch_ioerror_handles_ioerror_with_message_attribute(self):
        """Test IOError handling when error has message attribute."""
        @self.decorator
        def failing_function():
            error = IOError("Custom error message")
            error.filename = "/test/file.txt"
            error.errno = 2
            error.message = "File not found"
            raise error

        with self.assertRaises(IOError):
            failing_function()

        self.mock_logger.error.assert_called_once()
        call_args = self.mock_logger.error.call_args[0][0]
        self.assertIn("File not found", call_args)

    def test_catch_ioerror_handles_ioerror_without_message_attribute(self):
        """Test IOError handling when error lacks message attribute."""
        @self.decorator
        def failing_function():
            error = IOError("Basic error")
            error.filename = "/test/file.txt"
            error.errno = 5
            # Deliberately not setting message attribute
            raise error

        with self.assertRaises(IOError):
            failing_function()

        self.mock_logger.error.assert_called_once()
        call_args = self.mock_logger.error.call_args[0][0]
        # Should handle missing message attribute gracefully
        self.assertIn("[Errno 5]", call_args)

    def test_catch_ioerror_handles_indexerror(self):
        """Test IndexError handling."""
        @self.decorator
        def failing_function():
            raise IndexError("Index out of range")

        with self.assertRaises(IndexError):
            failing_function()

        self.mock_logger.error.assert_called_once()
        call_args = self.mock_logger.error.call_args[0][0]
        self.assertIn("Invalid input file", call_args)
        self.assertIn("Index out of range", call_args)

    def test_catch_ioerror_handles_stopiteration(self):
        """Test StopIteration handling."""
        @self.decorator
        def failing_function():
            raise StopIteration("No more items")

        with self.assertRaises(StopIteration):
            failing_function()

        self.mock_logger.error.assert_called_once()
        call_args = self.mock_logger.error.call_args[0][0]
        self.assertIn("Invalid input file", call_args)
        self.assertIn("No more items", call_args)

    def test_catch_ioerror_preserves_function_metadata(self):
        """Test that the decorator preserves function metadata."""
        @self.decorator
        def test_function():
            """Test function docstring."""
            pass

        self.assertEqual(test_function.__name__, "test_function")
        self.assertEqual(test_function.__doc__, "Test function docstring.")

    def test_catch_ioerror_handles_functions_with_args_and_kwargs(self):
        """Test decorator works with functions that have args and kwargs."""
        @self.decorator
        def complex_function(a, b, c=None, *args, **kwargs):
            if c == "error":
                raise IOError("Test error")
            return (a, b, c, args, kwargs)

        # Test successful call with complex arguments
        result = complex_function(1, 2, c="success", extra=True)
        self.assertEqual(result, (1, 2, "success", (), {"extra": True}))

        # Test error case
        with self.assertRaises(IOError):
            complex_function(1, 2, c="error")


class TestPrepareOutdir(unittest.TestCase):
    """Test prepare_outdir function functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.mock_logger = Mock(spec=logging.Logger)
        self.temp_dir = None

    def tearDown(self):
        """Clean up test fixtures."""
        if self.temp_dir and os.path.exists(self.temp_dir):
            shutil.rmtree(self.temp_dir)

    def test_prepare_outdir_existing_directory_success(self):
        """Test successful preparation of existing directory."""
        with tempfile.TemporaryDirectory() as temp_dir:
            result = prepare_outdir(temp_dir, self.mock_logger)
            self.assertTrue(result)
            self.mock_logger.critical.assert_not_called()

    def test_prepare_outdir_existing_directory_with_pathlib_path(self):
        """Test successful preparation with pathlib.Path object."""
        with tempfile.TemporaryDirectory() as temp_dir:
            path_obj = Path(temp_dir)
            result = prepare_outdir(path_obj, self.mock_logger)
            self.assertTrue(result)
            self.mock_logger.critical.assert_not_called()

    def test_prepare_outdir_path_exists_but_not_directory(self):
        """Test failure when path exists but is not a directory."""
        with tempfile.NamedTemporaryFile(delete=False) as temp_file:
            temp_file_path = temp_file.name

        try:
            result = prepare_outdir(temp_file_path, self.mock_logger)
            self.assertFalse(result)
            
            # Check that appropriate error messages were logged
            self.assertEqual(self.mock_logger.critical.call_count, 2)
            calls = self.mock_logger.critical.call_args_list
            self.assertIn("Specified path as a output directory is not directory", calls[0][0][0])
            self.assertEqual(calls[1][0][0], temp_file_path)
        finally:
            os.unlink(temp_file_path)

    def test_prepare_outdir_creates_new_directory_successfully(self):
        """Test successful creation of new directory."""
        with tempfile.TemporaryDirectory() as parent_dir:
            new_dir = os.path.join(parent_dir, "new_output_dir")
            
            result = prepare_outdir(new_dir, self.mock_logger)
            self.assertTrue(result)
            self.assertTrue(os.path.isdir(new_dir))
            
            # Check info message was logged
            self.mock_logger.info.assert_called_once()
            info_call = self.mock_logger.info.call_args[0][0]
            self.assertIn("Make output directory", info_call)
            self.assertIn(new_dir, info_call)

    def test_prepare_outdir_creates_nested_directory_successfully(self):
        """Test successful creation of nested directories."""
        with tempfile.TemporaryDirectory() as parent_dir:
            nested_dir = os.path.join(parent_dir, "level1", "level2", "output")
            
            result = prepare_outdir(nested_dir, self.mock_logger)
            self.assertTrue(result)
            self.assertTrue(os.path.isdir(nested_dir))

    @patch('pathlib.Path.mkdir')
    def test_prepare_outdir_directory_creation_ioerror(self, mock_mkdir):
        """Test failure when directory creation raises IOError."""
        mock_mkdir.side_effect = IOError("Permission denied")
        # Add errno and message attributes to match the real IOError behavior
        mock_mkdir.side_effect.errno = 13
        mock_mkdir.side_effect.message = "Permission denied"

        with tempfile.TemporaryDirectory() as parent_dir:
            new_dir = os.path.join(parent_dir, "protected_dir")
            # Remove the parent to ensure it doesn't exist
            os.rmdir(parent_dir)
            
            result = prepare_outdir(new_dir, self.mock_logger)
            self.assertFalse(result)
            
            # Check error message was logged
            self.mock_logger.critical.assert_called_once()
            call_args = self.mock_logger.critical.call_args[0][0]
            self.assertIn("Faild to make output directory", call_args)
            self.assertIn("[Errno 13]", call_args)

    @patch('pathlib.Path.mkdir')
    def test_prepare_outdir_directory_creation_ioerror_without_message(self, mock_mkdir):
        """Test directory creation failure when IOError lacks message attribute."""
        error = IOError("Access denied")
        error.errno = 5
        # Deliberately not setting message attribute
        mock_mkdir.side_effect = error

        with tempfile.TemporaryDirectory() as parent_dir:
            new_dir = os.path.join(parent_dir, "restricted_dir")
            # Remove parent to trigger mkdir
            os.rmdir(parent_dir)
            
            result = prepare_outdir(new_dir, self.mock_logger)
            self.assertFalse(result)
            
            # Should handle missing message attribute gracefully
            self.mock_logger.critical.assert_called_once()
            call_args = self.mock_logger.critical.call_args[0][0]
            self.assertIn("[Errno 5]", call_args)

    @patch('os.access')
    def test_prepare_outdir_directory_not_writable(self, mock_access):
        """Test failure when directory is not writable."""
        mock_access.return_value = False

        with tempfile.TemporaryDirectory() as temp_dir:
            result = prepare_outdir(temp_dir, self.mock_logger)
            self.assertFalse(result)
            
            # Check that os.access was called with correct parameters
            mock_access.assert_called_once_with(temp_dir, os.W_OK)
            
            # Check error message was logged
            self.mock_logger.critical.assert_called_once()
            call_args = self.mock_logger.critical.call_args[0][0]
            self.assertIn("Output directory", call_args)
            self.assertIn("is not writable", call_args)
            self.assertIn(temp_dir, call_args)

    @patch('os.access')
    @patch('pathlib.Path.mkdir')
    def test_prepare_outdir_complex_failure_scenario(self, mock_mkdir, mock_access):
        """Test complex failure scenario with multiple potential issues."""
        # First call succeeds (directory creation), second call fails (write check)
        mock_access.return_value = False
        
        with tempfile.TemporaryDirectory() as parent_dir:
            new_dir = os.path.join(parent_dir, "complex_test_dir")
            os.rmdir(parent_dir)  # Remove to trigger mkdir
            
            result = prepare_outdir(new_dir, self.mock_logger)
            self.assertFalse(result)
            
            # Both mkdir and access check should have been attempted
            mock_mkdir.assert_called_once()
            mock_access.assert_called_once()


if __name__ == '__main__':
    unittest.main()