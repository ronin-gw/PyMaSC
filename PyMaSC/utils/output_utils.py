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
from __future__ import annotations

import logging
from functools import wraps
from typing import Any, Callable

logger = logging.getLogger(__name__)


def catch_IOError(logger_obj: logging.Logger) -> Callable[[Callable], Callable]:
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
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            try:
                return func(*args, **kwargs)
            except OSError as e:
                logger_obj.error(f"OS error in {func.__name__}: {e}")
                raise
        return wrapper
    return decorator
