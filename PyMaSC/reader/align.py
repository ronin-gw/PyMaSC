import logging
import sys

from reader.sam import SamReadGenerator, SamStreamReadGenerator
from reader.bam import BamReadGenerator


logger = logging.getLogger(__name__)
GZIP_MAGIC_NUMBER = '\x1f\x8b'


class InputUnseekable(Exception):
    pass


def get_read_generator_and_init_target(path, fmt=None):
    if path == '-':
        if not fmt:
            raise InputUnseekable
        path = None
        path_or_stream = sys.stdin
    else:
        if not fmt:
            fmt = _guess_format(path)
        path_or_stream = path

    if fmt == "BAM":
        return BamReadGenerator, path_or_stream
    elif fmt == "SAM":
        if path:
            return SamReadGenerator, path
        else:
            return SamStreamReadGenerator, path_or_stream


def _guess_format(path):
    with open(path, "rb") as f:
        if f.read(2) == GZIP_MAGIC_NUMBER:
            fmt = "BAM"
        else:
            fmt = "SAM"

        try:
            f.seek(0)
        except IOError as e:
            if e.errno == 29:
                logger.error("Faild to file seek back after format prediction.")
                raise InputUnseekable
            else:
                raise e

    return fmt
