import logging
import sys

from pysam import AlignmentFile

logger = logging.getLogger(__name__)
GZIP_MAGIC_NUMBER = '\x1f\x8b'


class InputUnseekable(Exception):
    pass


class BamReadGenerator(AlignmentFile):
    def next(self):
        while True:
            read = super(BamReadGenerator, self).next()
            # https://github.com/pysam-developers/pysam/issues/268
            try:
                return (read.is_reverse, read.is_duplicate, read.reference_name,
                        read.reference_start + 1, read.mapping_quality, read.query_length)
            except ValueError:
                pass


class SamStreamReadGenerator(object):
    def __iter__(self):
        return self

    def __init__(self, stream):
        self.stream = stream
        self._is_header_parsed = False
        self._buff = None
        self.references = []
        self.lengths = []

        self._parse_header()

    def _parse_header(self):
        for line in self:
            if not line.startswith('@'):
                self._buff = line
                break
            elif line.startswith("@SQ"):
                _, refname, length = line.split('\t')
                self.references.append(refname[3:])
                self.lengths.append(int(length[3:]))
        self.references = tuple(self.references)
        self.lengths = tuple(self.lengths)

        self._is_header_parsed = True

    def _next_source(self):
        return next(self.stream)

    def next(self):
        if not self._is_header_parsed:
            return self._next_source()

        if self._buff is None:
            raise StopIteration

        _, flag, chrom, pos, mapq, _, _, _, seq, _ = self._buff.split('\t', 9)
        try:
            self._buff = self._next_source()
        except StopIteration:
            self._buff = None

        flag = int(flag)
        return flag & 16, flag & 1024, chrom, int(pos), int(mapq), len(seq)


class SamReadGenerator(SamStreamReadGenerator, file):
    def __init__(self, *args, **kwargs):
        file.__init__(self, *args, **kwargs)
        SamStreamReadGenerator.__init__(self, self)

    def _next_source(self):
        return file.next(self)


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
