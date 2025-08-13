"""PyMaSC package initialization with version and common utilities.

Provides package constants, version information, and common utilities
for PyMaSC command-line tools.

The module provides:
- VERSION: Package version string
- WEBSITE_URL: Official PyMaSC website
- logging_version(): Version logging utility
- entrypoint(): Decorator for main entry point exception handling
"""
import logging
import sys
import traceback
import multiprocessing
from functools import wraps
from typing import Any, Callable

VERSION = "1.0.0"
WEBSITE_URL = "https://pymasc.sb.ecei.tohoku.ac.jp/"

logger = logging.getLogger(__name__)


def logging_version(logger: Any) -> None:
    """Log PyMaSC and Python version information.

    Outputs version information including PyMaSC version and Python details
    for debugging and support purposes.

    Args:
        logger: Logger instance to use for version output
    """
    logger.info("PyMaSC version {} with Python{}.{}.{}".format(
                *[VERSION] + list(sys.version_info[:3])))
    for line in sys.version.split('\n'):
        logger.debug(line)


def _ensure_spawn() -> None:
    """Ensure the multiprocessing start method is set to `spawn`."""
    try:
        current = multiprocessing.get_start_method(allow_none=True)
    except TypeError:
        current = None
    if current != "spawn":
        try:
            multiprocessing.set_start_method("spawn")
        except RuntimeError:
            logger.warning(
                "Failed to set multiprocessing start method to 'spawn'. "
                "This may cause issues with multiprocessing in PyMaSC."
            )


def entrypoint(logger: Any) -> Callable[[Callable[[], None]], Callable[[], None]]:
    """Decorator for main entry point exception handling.

    Provides consistent exception handling and cleanup for all PyMaSC
    command-line entry points. Handles KeyboardInterrupt gracefully
    and ensures proper logging of completion status.

    Args:
        logger: Logger instance for status messages

    Returns:
        Decorator function that wraps main entry point functions

    Example:
        @entrypoint(logger)
        def main():
            # Main application logic
            pass
    """
    def _entrypoint_wrapper_base(main_func: Callable[[], None]) -> Callable[[], None]:
        @wraps(main_func)
        def _inner() -> None:
            try:
                _ensure_spawn()
                main_func()
                logger.info("PyMASC finished.")
            except KeyboardInterrupt:
                sys.stderr.write("\r\033[K")
                sys.stderr.flush()
                logger.info("Got KeyboardInterrupt. bye")
                if 0 < logger.level <= logging.DEBUG:
                    traceback.print_exc()
        return _inner
    return _entrypoint_wrapper_base
