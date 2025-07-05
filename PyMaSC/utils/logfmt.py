"""Colored logging formatter for terminal output.

This module provides enhanced logging functionality with ANSI color support
for PyMaSC command-line tools. It automatically detects terminal capabilities
and applies appropriate formatting for improved readability.

Key components:
- set_rootlogger(): Configures the root logger with color support
- StrFormatStyle: String format style for logging messages
- ColorfulFormatter: ANSI color formatter with level-based coloring

The formatter automatically adjusts output based on whether output is
directed to a terminal or redirected to a file, ensuring compatibility
across different output scenarios.
"""
import logging
import sys

LOGGING_FORMAT = "[%(asctime)s | %(levelname)s] %(name)10s : %(message)s"


def set_rootlogger(args_color, log_level):
    """Configure root logger with color support.

    Sets up the root logger with appropriate formatting and color settings
    based on user preferences and terminal capabilities.

    Args:
        args_color: Color preference ('TRUE', 'FALSE', or None for auto-detect)
        log_level: Logging level (e.g., logging.INFO, logging.DEBUG)

    Returns:
        Configured root logger instance

    Note:
        Automatically detects terminal capabilities when args_color is None
    """
    if args_color == "TRUE":
        colorize = True
    elif args_color == "FALSE":
        colorize = False
    else:
        colorize = sys.stderr.isatty()

    rl = logging.getLogger('')

    h = logging.StreamHandler()
    h.setFormatter(ColorfulFormatter(fmt=LOGGING_FORMAT, colorize=colorize))

    rl.addHandler(h)
    rl.setLevel(log_level)

    return rl


class StrFormatStyle(object):
    """String format style for logging.

    Provides string formatting functionality for log messages using
    Python's str.format() method. This class defines the format style
    and time usage detection for log record formatting.

    Attributes:
        default_format: Default message format string
        asctime_format: Format string for timestamps
        asctime_search: Search pattern for timestamp detection
    """
    default_format = '{message}'
    asctime_format = '{asctime}'
    asctime_search = '{asctime'

    def __init__(self, fmt):
        self._fmt = fmt or self.default_format

    def usesTime(self):
        return self._fmt.find(self.asctime_search) >= 0

    def format(self, record):
        print(self._fmt)
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
    DEFAULT_COLOR = 39
    LOGLEVEL2COLOR = {20: 36, 30: 33, 40: 31, 50: 35}

    @classmethod
    def bleach(cls):
        cls._fmt = cls._monofmt

    def __init__(self, fmt=None, datefmt=None, colorize=True):
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

    def fill_format(self, record):
        fmt = self._fmt.format(
            col=self.LOGLEVEL2COLOR.get(record.levelno, self.DEFAULT_COLOR),
            msg=0 if record.levelno < 40 else 1)
        return fmt % record.__dict__

    def format(self, record):
        record.message = record.getMessage()
        if self.usesTime():
            record.asctime = self.formatTime(record, self.datefmt)
        try:
            s = self.fill_format(record)
        except UnicodeDecodeError as e:
            try:
                record.name = record.name.decode('utf-8')
                s = self.fill_format(record)
            except UnicodeDecodeError:
                raise e
        if record.exc_info:
            if not record.exc_text:
                record.exc_text = self.formatException(record.exc_info)
        if record.exc_text:
            if s[-1:] != "\n":
                s = s + "\n"
            try:
                s = s + record.exc_text
            except UnicodeError:
                s = s + record.exc_text.decode(sys.getfilesystemencoding(), 'replace')
        return s
