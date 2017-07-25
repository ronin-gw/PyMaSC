import logging
import sys


class StrFormatStyle(object):
    default_format = '{message}'
    asctime_format = '{asctime}'
    asctime_search = '{asctime'

    def __init__(self, fmt):
        self._fmt = fmt or self.default_format

    def usesTime(self):
        return self._fmt.find(self.asctime_search) >= 0

    def format(self, record):
        print self._fmt
        return self._fmt.format(**record.__dict__)


class ColorfulFormatter(logging.Formatter):
    DEFAULT_COLOR = 39
    LOGLEVEL2COLOR = {20: 36, 30: 33, 40: 35, 50: 31}

    def __init__(self, fmt=None, datefmt=None):
        super(ColorfulFormatter, self).__init__(fmt=fmt, datefmt=datefmt)
        self._fmt = self._fmt.replace(
            "%(levelname)s", "\033[{col}m%(levelname)8s\033[0m"
        ).replace(
            "%(message)s", "\033[{msg}m%(message)s\033[0m"
        )

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
