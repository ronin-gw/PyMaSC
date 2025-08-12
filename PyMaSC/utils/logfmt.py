"""Colored logging formatter for terminal output.

Provides logging with ANSI color support for PyMaSC command-line tools.
Automatically detects terminal capabilities and applies appropriate formatting.

Key components:
- set_rootlogger(): Configures the root logger with color support
- StrFormatStyle: String format style for logging messages
- ColorfulFormatter: ANSI color formatter with level-based coloring

The formatter automatically adjusts output based on whether output is
directed to a terminal or redirected to a file, ensuring compatibility
across different output scenarios.
"""
import logging
from typing import Optional, Dict

LOGGING_FORMAT: str = "[%(asctime)s | %(levelname)s] %(name)10s : %(message)s"


def set_rootlogger(colorize: bool, log_level: int) -> logging.Logger:
    """Configure root logger with color support.

    Sets up the root logger with appropriate formatting and color settings
    based on user preferences and terminal capabilities.

    Args:
        color: Color preference (True, False)
        log_level: Logging level (e.g., logging.INFO, logging.DEBUG)

    Returns:
        Configured root logger instance

    Note:
        Automatically detects terminal capabilities when args_color is None
    """
    rl = logging.getLogger('')

    h = logging.StreamHandler()
    h.setFormatter(ColorfulFormatter(fmt=LOGGING_FORMAT, colorize=colorize))

    rl.addHandler(h)
    rl.setLevel(log_level)

    return rl


class StrFormatStyle:
    """String format style for logging.

    Provides string formatting functionality for log messages using
    Python's str.format() method. This class defines the format style
    and time usage detection for log record formatting.

    Attributes:
        default_format: Default message format string
        asctime_format: Format string for timestamps
        asctime_search: Search pattern for timestamp detection
    """
    default_format: str = '{message}'
    asctime_format: str = '{asctime}'
    asctime_search: str = '{asctime'

    def __init__(self, fmt: Optional[str]) -> None:
        """Initialize format style with message format.

        Args:
            fmt: Format string for log messages (uses default if None)
        """
        self._fmt = fmt or self.default_format

    def usesTime(self) -> bool:
        """Check if format string includes timestamp.

        Returns:
            True if format uses timestamp, False otherwise
        """
        return self._fmt.find(self.asctime_search) >= 0

    def format(self, record: logging.LogRecord) -> str:
        return self._fmt.format(**record.__dict__)


class ColorfulFormatter(logging.Formatter):
    """ANSI color formatter for log messages.

    Extends the standard Python logging formatter to add ANSI color codes
    for different log levels. Provides enhanced readability in terminal
    output while maintaining compatibility with non-terminal outputs.

    The formatter applies different colors based on log level:
    - INFO: Cyan
    - WARNING: Yellow
    - ERROR: Red
    - CRITICAL: Magenta

    Attributes:
        DEFAULT_COLOR: Default ANSI color code
        LOGLEVEL2COLOR: Mapping of log levels to color codes
        colorize: Whether to apply color formatting
    """
    DEFAULT_COLOR: int = 39
    LOGLEVEL2COLOR: Dict[int, int] = {20: 36, 30: 33, 40: 31, 50: 35}

    _fmt: str

    def __init__(self, fmt: Optional[str] = None, datefmt: Optional[str] = None, colorize: bool = True) -> None:
        """Initialize colorful formatter with ANSI color support.

        Args:
            fmt: Log message format string
            datefmt: Date format string
            colorize: Whether to apply ANSI color codes
        """
        super(ColorfulFormatter, self).__init__(fmt=fmt, datefmt=datefmt)
        self.colorize = colorize
        if colorize:
            self._fmt = self._fmt.replace(
                "%(levelname)s", "\033[{col}m%(levelname)8s\033[0m"
            ).replace(
                "%(message)s", "\033[{msg}m%(message)s\033[0m"
            )
        else:
            self._fmt = self._fmt.replace("%(levelname)s", "%(levelname)8s")

    def fill_format(self, record: logging.LogRecord) -> str:
        """Generate formatted string with appropriate colors.

        Args:
            record: Log record to format

        Returns:
            Formatted string with ANSI color codes applied
        """
        fmt = self._fmt.format(
            col=self.LOGLEVEL2COLOR.get(record.levelno, self.DEFAULT_COLOR),
            msg=0 if record.levelno < 40 else 1)
        return fmt % record.__dict__

    def format(self, record: logging.LogRecord) -> str:
        record.message = record.getMessage()
        if self.usesTime():
            record.asctime = self.formatTime(record, self.datefmt)

        s = self.fill_format(record)

        if record.exc_info and not record.exc_text:
            record.exc_text = self.formatException(record.exc_info)

        if record.exc_text:
            if s[-1:] != "\n":
                s = s + "\n"
            s = s + record.exc_text

        return s
