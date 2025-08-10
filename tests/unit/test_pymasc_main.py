"""Tests for PyMaSC.pymasc main CLI module.

This module tests the main CLI functionality including error handling
for file processing, output handling, and command-line argument validation.
"""
import unittest
from unittest.mock import Mock, patch, MagicMock, call
import argparse
import logging
import sys
from pathlib import Path
from typing import List, Optional

from PyMaSC.pymasc import (
    _parse_args, prepare_output, set_readlen, output_results
)
from PyMaSC.core.exceptions import ReadUnsortedError, ReadsTooFew, InputUnseekable, NothingToCalc
from PyMaSC.handler.mappability import BWIOError, JSONIOError


class TestParseArgs(unittest.TestCase):
    """Test argument parsing and validation."""

    def setUp(self):
        """Set up test fixtures."""
        self.original_argv = sys.argv

    def tearDown(self):
        """Clean up test fixtures."""
        sys.argv = self.original_argv

    @patch('PyMaSC.pymasc.get_pymasc_parser')
    def test_parse_args_skip_ncc_without_mappability_raises_error(self, mock_get_parser):
        """Test that --skip-ncc without --mappable raises parser error."""
        mock_parser = Mock()
        mock_get_parser.return_value = mock_parser
        
        # Mock args that would cause the error
        mock_args = Mock()
        mock_args.skip_ncc = True
        mock_args.mappability = None
        mock_parser.parse_args.return_value = mock_args
        
        # Should call parser.error
        mock_parser.error.side_effect = SystemExit("Parser error")
        
        with self.assertRaises(SystemExit):
            _parse_args()
        
        mock_parser.error.assert_called_once_with(
            "argument --skip-ncc: -m/--mappable must be specified."
        )

    @patch('PyMaSC.pymasc.get_pymasc_parser')
    @patch('PyMaSC.pymasc.set_rootlogger')
    @patch('PyMaSC.pymasc.logging_version')
    @patch('PyMaSC.pymasc.logger')
    def test_parse_args_library_length_exceeds_max_shift_warning(self, mock_logger, 
                                                               mock_logging_version, 
                                                               mock_set_rootlogger, 
                                                               mock_get_parser):
        """Test warning when library_length exceeds max_shift."""
        mock_parser = Mock()
        mock_get_parser.return_value = mock_parser
        
        mock_args = Mock()
        mock_args.skip_ncc = False
        mock_args.mappability = "test.bw"
        mock_args.mappability_stats = None
        mock_args.library_length = 1000
        mock_args.max_shift = 500  # Less than library_length
        mock_args.color = False
        mock_args.log_level = "INFO"
        mock_parser.parse_args.return_value = mock_args
        
        result = _parse_args()
        
        # Should log error and set library_length to None
        mock_logger.error.assert_called_once()
        error_msg = mock_logger.error.call_args[0][0]
        self.assertIn("expected library length > max shift", error_msg)
        self.assertIsNone(result.library_length)

    @patch('PyMaSC.pymasc.get_pymasc_parser')
    @patch('PyMaSC.pymasc.set_rootlogger')
    @patch('PyMaSC.pymasc.logging_version')
    def test_parse_args_mappability_stats_same_as_mappability_reset(self, mock_logging_version, 
                                                                   mock_set_rootlogger, 
                                                                   mock_get_parser):
        """Test that mappability_stats is reset when same as mappability path."""
        mock_parser = Mock()
        mock_get_parser.return_value = mock_parser
        
        mock_args = Mock()
        mock_args.skip_ncc = False
        mock_args.mappability = "test.bw"
        mock_args.mappability_stats = "test.bw"  # Same as mappability
        mock_args.library_length = None
        mock_args.max_shift = 500
        mock_args.color = False
        mock_args.log_level = "INFO"
        mock_parser.parse_args.return_value = mock_args
        
        result = _parse_args()
        
        # Should reset mappability_stats to None
        self.assertIsNone(result.mappability_stats)


class TestPrepareOutput(unittest.TestCase):
    """Test prepare_output function."""

    def setUp(self):
        """Set up test fixtures."""
        self.mock_logger = Mock(spec=logging.Logger)

    @patch('PyMaSC.pymasc.prepare_outdir')
    @patch('PyMaSC.pymasc.sys.exit')
    def test_prepare_output_directory_creation_failure_exits(self, mock_exit, mock_prepare_outdir):
        """Test that prepare_output exits when directory creation fails."""
        mock_prepare_outdir.return_value = False
        
        with patch('PyMaSC.pymasc.logger', self.mock_logger):
            prepare_output(["test.bam"], [None], "/unwritable/dir")
        
        mock_exit.assert_called_once_with(1)

    @patch('PyMaSC.pymasc.prepare_outdir')
    def test_prepare_output_warns_about_existing_files(self, mock_prepare_outdir):
        """Test that prepare_output warns about files that will be overwritten."""
        mock_prepare_outdir.return_value = True
        
        # Create a temporary file to simulate existing output file
        import tempfile
        import os
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create an existing file that would be overwritten
            existing_file = os.path.join(temp_dir, "test.pdf")
            with open(existing_file, 'w') as f:
                f.write("existing")
            
            with patch('PyMaSC.pymasc.logger', self.mock_logger):
                result = prepare_output(["test.bam"], [None], temp_dir, (".pdf",))
        
        # Should warn about existing file
        self.mock_logger.warning.assert_called()
        warning_msg = self.mock_logger.warning.call_args[0][0]
        self.assertIn("will be overwritten", warning_msg)

    @patch('PyMaSC.pymasc.prepare_outdir')
    def test_prepare_output_uses_custom_names(self, mock_prepare_outdir):
        """Test that prepare_output uses custom names when provided."""
        mock_prepare_outdir.return_value = True
        
        # Use a temporary directory for the test
        import tempfile
        with tempfile.TemporaryDirectory() as temp_dir:
            result = prepare_output(["test.bam"], ["custom_name"], temp_dir, ())
            
            # Should return path with custom name
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0].name, "custom_name")


class TestSetReadlen(unittest.TestCase):
    """Test set_readlen function."""

    def setUp(self):
        """Set up test fixtures."""
        self.mock_logger = Mock(spec=logging.Logger)

    def test_set_readlen_uses_specified_length(self):
        """Test that set_readlen uses user-specified read length."""
        mock_args = Mock()
        mock_args.read_length = 50
        
        mock_handler1 = Mock()
        mock_handler2 = Mock()
        handlers = [mock_handler1, mock_handler2]
        
        result = set_readlen(mock_args, handlers)
        
        self.assertEqual(result, 50)
        self.assertEqual(mock_handler1.read_len, 50)
        self.assertEqual(mock_handler2.read_len, 50)

    def test_set_readlen_estimates_from_data(self):
        """Test read length estimation from data."""
        mock_args = Mock()
        mock_args.read_length = None
        mock_args.readlen_estimator = "MODE"
        
        mock_handler1 = Mock()
        mock_handler1.estimate_readlen.return_value = 36
        mock_handler2 = Mock()
        mock_handler2.estimate_readlen.return_value = 50
        handlers = [mock_handler1, mock_handler2]
        
        with patch('PyMaSC.pymasc.logger', self.mock_logger):
            result = set_readlen(mock_args, handlers)
        
        # Should use maximum read length
        self.assertEqual(result, 50)
        self.assertEqual(mock_handler1.read_len, 50)
        self.assertEqual(mock_handler2.read_len, 50)

    def test_set_readlen_warns_about_multiple_lengths(self):
        """Test warning when multiple read lengths are detected."""
        mock_args = Mock()
        mock_args.read_length = None
        mock_args.readlen_estimator = "MODE"
        
        mock_handler1 = Mock()
        mock_handler1.estimate_readlen.return_value = 36
        mock_handler2 = Mock()
        mock_handler2.estimate_readlen.return_value = 50  # Different length
        handlers = [mock_handler1, mock_handler2]
        
        with patch('PyMaSC.pymasc.logger', self.mock_logger):
            result = set_readlen(mock_args, handlers)
        
        # Should warn about multiple lengths
        self.mock_logger.warning.assert_called_once()
        warning_msg = self.mock_logger.warning.call_args[0][0]
        self.assertIn("multiple read length candidates", warning_msg)
        self.assertIn("50", warning_msg)  # Should mention max length

    def test_set_readlen_handles_estimation_failures(self):
        """Test handling of read length estimation failures."""
        mock_args = Mock()
        mock_args.read_length = None
        mock_args.readlen_estimator = "MODE"
        
        mock_handler1 = Mock()
        mock_handler1.estimate_readlen.side_effect = ValueError("Estimation failed")
        mock_handler2 = Mock()
        mock_handler2.estimate_readlen.return_value = 50
        handlers = [mock_handler1, mock_handler2]
        
        with patch('PyMaSC.pymasc.logger', self.mock_logger):
            result = set_readlen(mock_args, handlers)
        
        # Should remove failed handler and use remaining ones
        self.assertEqual(len(handlers), 1)  # Handler with ValueError should be removed
        self.assertEqual(result, 50)


class TestOutputResults(unittest.TestCase):
    """Test output_results function."""

    def setUp(self):
        """Set up test fixtures."""
        self.mock_logger = Mock(spec=logging.Logger)

    @patch('PyMaSC.pymasc.output_stats')
    @patch('PyMaSC.pymasc.output_nreads_table')
    @patch('PyMaSC.pymasc.output_cc')
    @patch('PyMaSC.pymasc.output_mscc')
    def test_output_results_handles_none_result(self, mock_output_mscc, mock_output_cc, 
                                              mock_output_nreads, mock_output_stats):
        """Test that output_results handles None result gracefully."""
        mock_args = Mock()
        output_basename = Path("/test/output")
        
        output_results(mock_args, output_basename, None)
        
        # Should not call any output functions for None result
        mock_output_stats.assert_not_called()
        mock_output_nreads.assert_not_called()
        mock_output_cc.assert_not_called()
        mock_output_mscc.assert_not_called()

    @patch('PyMaSC.pymasc.output_stats')
    @patch('PyMaSC.pymasc.output_nreads_table')
    @patch('PyMaSC.pymasc.output_cc')
    @patch('PyMaSC.pymasc.output_mscc')
    def test_output_results_calls_appropriate_outputs(self, mock_output_mscc, mock_output_cc,
                                                     mock_output_nreads, mock_output_stats):
        """Test that output_results calls appropriate output functions."""
        mock_args = Mock()
        mock_args.skip_plots = True  # Skip plots to avoid import testing
        
        mock_result = Mock()
        mock_result.whole_ncc_stats = Mock()  # NCC stats present
        mock_result.whole_mscc_stats = Mock()  # MSCC stats present
        
        output_basename = Path("/test/output")
        
        output_results(mock_args, output_basename, mock_result)
        
        # Should call all output functions
        mock_output_stats.assert_called_once_with(output_basename, mock_result)
        mock_output_nreads.assert_called_once_with(output_basename, mock_result)
        mock_output_cc.assert_called_once_with(output_basename, mock_result)
        mock_output_mscc.assert_called_once_with(output_basename, mock_result)

    @patch('PyMaSC.pymasc.output_stats')
    @patch('PyMaSC.pymasc.output_nreads_table')
    @patch('PyMaSC.pymasc.output_cc')
    @patch('PyMaSC.pymasc.output_mscc')
    def test_output_results_skips_missing_stats(self, mock_output_mscc, mock_output_cc,
                                              mock_output_nreads, mock_output_stats):
        """Test that output_results skips output for missing stats."""
        mock_args = Mock()
        mock_args.skip_plots = True
        
        mock_result = Mock()
        mock_result.whole_ncc_stats = None  # No NCC stats
        mock_result.whole_mscc_stats = Mock()  # MSCC stats present
        
        output_basename = Path("/test/output")
        
        output_results(mock_args, output_basename, mock_result)
        
        # Should skip CC output but do MSCC
        mock_output_stats.assert_called_once()
        mock_output_nreads.assert_called_once()
        mock_output_cc.assert_not_called()  # Should be skipped
        mock_output_mscc.assert_called_once()

    def test_output_results_handles_plotting_import_error(self):
        """Test that output_results handles matplotlib import errors gracefully."""
        mock_args = Mock()
        mock_args.skip_plots = False  # Enable plots
        
        mock_result = Mock()
        mock_result.whole_ncc_stats = None
        mock_result.whole_mscc_stats = None
        
        output_basename = Path("/test/output")
        
        # Mock the dynamic import to raise ImportError
        def mock_import_from(module, names, level=0):
            if module == "PyMaSC.output.figure":
                raise ImportError("No module named matplotlib")
            # For other imports, return a mock object
            return __import__(module, fromlist=names, level=level)
        
        # Patch output functions
        with patch('PyMaSC.pymasc.output_stats'), \
             patch('PyMaSC.pymasc.output_nreads_table'), \
             patch('PyMaSC.pymasc.logger', self.mock_logger):
            
            # We need to trigger the import error inside the try block
            # Let's patch the module import directly in sys.modules
            import sys
            original_modules = sys.modules.copy()
            
            # Remove the module if it exists to trigger import
            if 'PyMaSC.output.figure' in sys.modules:
                del sys.modules['PyMaSC.output.figure']
            
            try:
                # Mock __import__ to raise ImportError for our specific module
                def failing_import(name, *args, **kwargs):
                    if name == "PyMaSC.output.figure":
                        raise ImportError("No module named matplotlib")
                    # Use original import for others
                    return original_import(name, *args, **kwargs)
                
                original_import = __builtins__['__import__']
                __builtins__['__import__'] = failing_import
                
                output_results(mock_args, output_basename, mock_result)
                
            finally:
                # Restore original state
                __builtins__['__import__'] = original_import
                sys.modules.update(original_modules)
        
        # Should log error about skipping plots
        self.mock_logger.error.assert_called()
        error_msg = self.mock_logger.error.call_args[0][0]
        self.assertIn("Skip output plots", error_msg)

    @patch('PyMaSC.pymasc.output_stats')
    @patch('PyMaSC.pymasc.output_nreads_table')
    @patch('PyMaSC.output.figure.plot_figures')  
    def test_output_results_creates_plots_when_available(self, mock_plot_figures,
                                                        mock_output_nreads, mock_output_stats):
        """Test that output_results creates plots when matplotlib is available."""
        mock_args = Mock()
        mock_args.skip_plots = False
        
        mock_result = Mock()
        mock_result.whole_ncc_stats = None
        mock_result.whole_mscc_stats = None
        
        output_basename = Path("/test/output")
        
        output_results(mock_args, output_basename, mock_result)
        
        # Should call plot_figures with the result
        mock_plot_figures.assert_called_once()


if __name__ == '__main__':
    unittest.main()