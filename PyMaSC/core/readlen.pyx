import logging

from pysam.libcalignmentfile cimport AlignmentFile
from pysam.libcalignedsegment cimport AlignedSegment

from PyMaSC.utils.progress import ReadCountProgressBar

logger = logging.getLogger(__name__)


def _mean(c):
    return int(round(
        sum(length * freq for length, freq in c.items()) / float(sum(c.values()))
    ))


def _median(c):
    num = sum(c.values())
    target = num / 2
    _sum = 0

    if num % 2:
        for l in sorted(c):
            _sum += c[l]
            if target <= _sum:
                return l
    else:
        length = sorted(c)
        for i, l in enumerate(length):
            _sum += c[l]
            if target < _sum:
                return l
            elif target == _sum:
                return int(round((l + float(length[i+1])) / 2))


def _mode(c):
    return [k for k, v in sorted(c.items(), key=lambda x: x[1])][-1]


ESTFUNCTIONS = dict(
    MEAN=_mean, MEDIAN=_median, MODE=_mode, MIN=min, MAX=max
)


def estimate_readlen(path, esttype, unsigned int mapq_criteria):
    estfunc = ESTFUNCTIONS[esttype]

    cdef bint is_duplicate
    cdef str chrom, _chrom
    cdef unsigned int mapq
    cdef unsigned long pos, readlen

    cdef dict counter = {}

    cdef dict chrom2length
    cdef unsigned long length
    cdef AlignedSegment read

    with AlignmentFile(path) as stream:
        chrom = None
        chrom2length = dict(zip(stream.references, stream.lengths))
        length = sum(chrom2length.values())

        progress = ReadCountProgressBar()
        progress.set_genome(length)

        for read in stream:
            try:
                _chrom = read.reference_name
            except ValueError:
                continue
            if chrom != _chrom:
                chrom = _chrom
                progress.set_chrom(chrom2length[chrom], chrom)
            progress.update(read.reference_start)

            if not read.is_duplicate and read.mapping_quality >= mapq_criteria:
                readlen = read.query_length + 1
                if readlen in counter:
                    counter[readlen] += 1
                else:
                    counter[readlen] = 1
    progress.finish()

    length = estfunc(counter)

    logger.info("Estimated read length = {}".format(length))
    return length
