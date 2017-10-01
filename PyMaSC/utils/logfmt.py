from __future__ import print_function

import logging
import sys

LOGGING_FORMAT = "[%(asctime)s | %(levelname)s] %(name)10s : %(message)s"


def set_rootlogger(colorize, log_level):
    rl = logging.getLogger('')

    h = logging.StreamHandler()
    h.setFormatter(ColorfulFormatter(fmt=LOGGING_FORMAT, colorize=colorize))

    rl.addHandler(h)
    rl.setLevel(log_level)

    return rl


class StrFormatStyle(object):
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
