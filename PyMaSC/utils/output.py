"""File output utilities with error handling.

Provides utility functions for file and directory operations
with error handling and logging.

Key functionality:
- catch_IOError(): Decorator for I/O error handling
- prepare_outdir(): Output directory creation and validation
- Logging of file operation errors
- Permission checking and directory creation
"""
import os
from functools import wraps
from pathlib import Path
from typing import Any, Callable, TypeVar, Union
import logging


F = TypeVar('F', bound=Callable[..., Any])


def catch_IOError(logger: logging.Logger) -> Callable[[F], F]:
    """Decorator for I/O error handling.

    Creates a decorator that wraps functions to provide consistent
    error handling for I/O operations. Logs detailed error information
    and re-raises exceptions for proper error propagation.

    Args:
        logger: Logger instance for error reporting

    Returns:
        Decorator function that wraps I/O operations

    Example:
        @catch_IOError(logger)
        def write_file(path, data):
            with open(path, 'w') as f:
                f.write(data)
    """
    def _inner(func: F) -> F:
        @wraps(func)
        def _io_func(*args: Any, **kwargs: Any) -> Any:
            try:
                return func(*args, **kwargs)
            except IOError as e:
                logger.error("Faild to output '{}':\n[Errno {}] {}".format(
                    e.filename, e.errno, getattr(e, "message", ''))
                )
                raise e
            except (IndexError, StopIteration) as e:
                logger.error("Invalid input file: {}".format(e))
                raise e
        return _io_func  # type: ignore
    return _inner


def prepare_outdir(outdir: Union[str, Path], logger: logging.Logger) -> bool:
    """Create and validate output directory.

    Ensures the output directory exists, is writable, and is actually
    a directory. Creates the directory if it doesn't exist and validates
    permissions for file operations.

    Args:
        outdir: Path to output directory
        logger: Logger instance for status messages

    Returns:
        True if directory is ready for use, False if preparation failed

    Note:
        Logs detailed error messages for all failure scenarios
    """
    outdir_path = Path(outdir)
    if outdir_path.exists():
        if not outdir_path.is_dir():
            logger.critical("Specified path as a output directory is not directory.")
            logger.critical(str(outdir))
            return False
    else:
        logger.info("Make output directory: {}".format(outdir))
        try:
            outdir_path.mkdir(parents=True, exist_ok=True)
        except IOError as e:
            logger.critical("Faild to make output directory: [Errno {}] {}".format(e.errno, getattr(e, "message", "")))
            return False

    if not os.access(str(outdir), os.W_OK):
        logger.critical("Output directory '{}' is not writable.".format(outdir))
        return False

    return True
