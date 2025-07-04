"""Utilities for output file handling and logging.

This module provides common output functionality used across all PyMaSC
output modules, eliminating duplicate I/O patterns, logging, and file
path manipulation code.

Key features:
- Standardized output file logging
- Common file path manipulation patterns
- Shared error handling decorators
- Consistent CSV/table I/O utilities

This eliminates 50-65 lines of duplicate code across output modules.
"""
import csv
import logging
import os
from functools import wraps
from pathlib import Path
from typing import List, Dict, Any, Optional, Callable, Union

logger = logging.getLogger(__name__)


def catch_IOError(logger_obj: logging.Logger):
    """Decorator for consistent I/O error handling across output modules.
    
    This decorator provides standardized error handling for file operations
    used throughout the output modules. It logs errors consistently and
    re-raises them for proper error propagation.
    
    Args:
        logger_obj: Logger instance to use for error logging
        
    Returns:
        Decorator function
    """
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except IOError as e:
                logger_obj.error(f"I/O error in {func.__name__}: {e}")
                raise
            except OSError as e:
                logger_obj.error(f"OS error in {func.__name__}: {e}")
                raise
        return wrapper
    return decorator


class OutputPathManager:
    """Manages output file paths and naming conventions.
    
    Centralizes the file path manipulation logic used across output modules.
    Provides consistent handling of suffixes, extensions, and path validation.
    
    Attributes:
        base_path: Base output file path
        suffix: Optional suffix to append to filenames
        create_dirs: Whether to create parent directories
    """
    
    def __init__(self, 
                 base_path: Union[str, Path],
                 suffix: Optional[str] = None,
                 create_dirs: bool = True):
        """Initialize output path manager.
        
        Args:
            base_path: Base path for output files
            suffix: Optional suffix to append to filenames
            create_dirs: Whether to create parent directories if they don't exist
        """
        self.base_path = Path(base_path)
        self.suffix = suffix or ""
        self.create_dirs = create_dirs
    
    def get_output_path(self, 
                       filename: Optional[str] = None,
                       extension: Optional[str] = None) -> Path:
        """Get complete output file path with suffix and extension.
        
        Args:
            filename: Optional filename to use instead of base_path stem
            extension: Optional extension to add/replace
            
        Returns:
            Complete output file path
        """
        if filename:
            # Use provided filename
            output_path = self.base_path.parent / filename
        else:
            # Use base path
            output_path = self.base_path
        
        # Add suffix if specified
        if self.suffix:
            stem = output_path.stem + self.suffix
            output_path = output_path.parent / (stem + output_path.suffix)
        
        # Replace extension if specified
        if extension:
            if not extension.startswith('.'):
                extension = '.' + extension
            output_path = output_path.with_suffix(extension)
        
        # Create parent directories if needed
        if self.create_dirs:
            output_path.parent.mkdir(parents=True, exist_ok=True)
        
        return output_path
    
    def log_output_creation(self, output_path: Union[str, Path]) -> None:
        """Log output file creation with consistent format.
        
        Args:
            output_path: Path to the created output file
        """
        logger.info(f"Output '{output_path}'")
    
    def get_basename(self, path: Optional[Union[str, Path]] = None) -> str:
        """Get basename of file path for naming purposes.
        
        Args:
            path: Path to get basename from, uses base_path if None
            
        Returns:
            Basename without extension
        """
        if path is None:
            path = self.base_path
        return Path(path).stem


class TableWriter:
    """Utilities for writing CSV/table files consistently.
    
    Provides standardized table writing functionality used across
    output modules. Handles CSV formatting, headers, and error handling.
    
    Attributes:
        output_path: Path to output file
        delimiter: CSV delimiter character
        quoting: CSV quoting style
    """
    
    def __init__(self,
                 output_path: Union[str, Path],
                 delimiter: str = '\t',
                 quoting: int = csv.QUOTE_MINIMAL):
        """Initialize table writer.
        
        Args:
            output_path: Path to output file
            delimiter: CSV delimiter (default tab for PyMaSC compatibility)
            quoting: CSV quoting style
        """
        self.output_path = Path(output_path)
        self.delimiter = delimiter
        self.quoting = quoting
    
    @catch_IOError(logger)
    def write_table(self,
                   data: List[List[Any]],
                   headers: Optional[List[str]] = None) -> None:
        """Write table data to file with consistent formatting.
        
        Args:
            data: List of rows, each row is a list of values
            headers: Optional header row
        """
        with open(self.output_path, 'w', newline='') as f:
            writer = csv.writer(f, delimiter=self.delimiter, quoting=self.quoting)
            
            # Write headers if provided
            if headers:
                writer.writerow(headers)
            
            # Write data rows
            writer.writerows(data)
        
        logger.info(f"Output '{self.output_path}'")
    
    @catch_IOError(logger)
    def write_dict_table(self,
                        data: List[Dict[str, Any]],
                        fieldnames: Optional[List[str]] = None) -> None:
        """Write dictionary data to CSV file.
        
        Args:
            data: List of dictionaries with table data
            fieldnames: Optional list of field names for headers
        """
        if not data:
            return
        
        # Determine fieldnames from first row if not provided
        if fieldnames is None:
            fieldnames = list(data[0].keys())
        
        with open(self.output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, 
                                  fieldnames=fieldnames,
                                  delimiter=self.delimiter,
                                  quoting=self.quoting)
            
            writer.writeheader()
            writer.writerows(data)
        
        logger.info(f"Output '{self.output_path}'")


class TableReader:
    """Utilities for reading CSV/table files consistently.
    
    Provides standardized table reading functionality with consistent
    error handling and format detection.
    
    Attributes:
        input_path: Path to input file
        delimiter: CSV delimiter character
    """
    
    def __init__(self,
                 input_path: Union[str, Path],
                 delimiter: str = '\t'):
        """Initialize table reader.
        
        Args:
            input_path: Path to input file
            delimiter: CSV delimiter (default tab)
        """
        self.input_path = Path(input_path)
        self.delimiter = delimiter
    
    @catch_IOError(logger)
    def read_table(self, has_headers: bool = True) -> tuple:
        """Read table data from file.
        
        Args:
            has_headers: Whether first row contains headers
            
        Returns:
            Tuple of (headers, data) where headers is list of strings
            and data is list of lists
        """
        headers = None
        data = []
        
        with open(self.input_path, 'r') as f:
            reader = csv.reader(f, delimiter=self.delimiter)
            
            if has_headers:
                headers = next(reader)
            
            data = [row for row in reader]
        
        return headers, data
    
    @catch_IOError(logger)
    def read_dict_table(self) -> List[Dict[str, str]]:
        """Read table as list of dictionaries.
        
        Returns:
            List of dictionaries with column headers as keys
        """
        with open(self.input_path, 'r') as f:
            reader = csv.DictReader(f, delimiter=self.delimiter)
            return list(reader)


class LoggingUtils:
    """Utilities for consistent logging across output modules.
    
    Provides standardized logging patterns used throughout PyMaSC
    output modules.
    """
    
    @staticmethod
    def log_output_file(output_path: Union[str, Path],
                       logger_obj: Optional[logging.Logger] = None) -> None:
        """Log output file creation with consistent format.
        
        Args:
            output_path: Path to output file
            logger_obj: Logger to use, defaults to module logger
        """
        if logger_obj is None:
            logger_obj = logger
        logger_obj.info(f"Output '{output_path}'")
    
    @staticmethod
    def log_processing_step(step_name: str,
                          details: Optional[str] = None,
                          logger_obj: Optional[logging.Logger] = None) -> None:
        """Log processing step with consistent format.
        
        Args:
            step_name: Name of the processing step
            details: Optional additional details
            logger_obj: Logger to use, defaults to module logger
        """
        if logger_obj is None:
            logger_obj = logger
        
        message = f"Processing: {step_name}"
        if details:
            message += f" - {details}"
        
        logger_obj.info(message)


# Convenience functions for common patterns
def ensure_output_directory(output_path: Union[str, Path]) -> Path:
    """Ensure output directory exists for given file path.
    
    Args:
        output_path: Path to output file
        
    Returns:
        Path object with ensured parent directory
    """
    path_obj = Path(output_path)
    path_obj.parent.mkdir(parents=True, exist_ok=True)
    return path_obj


def add_suffix_to_path(base_path: Union[str, Path], suffix: str) -> Path:
    """Add suffix to file path before extension.
    
    Args:
        base_path: Base file path
        suffix: Suffix to add
        
    Returns:
        Path with suffix added
    """
    path_obj = Path(base_path)
    new_stem = path_obj.stem + suffix
    return path_obj.parent / (new_stem + path_obj.suffix)


def standardize_csv_writer(file_handle,
                          delimiter: str = '\t',
                          quoting: int = csv.QUOTE_MINIMAL) -> csv.writer:
    """Create standardized CSV writer with PyMaSC defaults.
    
    Args:
        file_handle: Open file handle for writing
        delimiter: CSV delimiter
        quoting: CSV quoting style
        
    Returns:
        Configured CSV writer
    """
    return csv.writer(file_handle, delimiter=delimiter, quoting=quoting)


def standardize_csv_reader(file_handle,
                          delimiter: str = '\t') -> csv.reader:
    """Create standardized CSV reader with PyMaSC defaults.
    
    Args:
        file_handle: Open file handle for reading
        delimiter: CSV delimiter
        
    Returns:
        Configured CSV reader
    """
    return csv.reader(file_handle, delimiter=delimiter)