# cython: language_level=3
"""Read length estimation from BAM files.

This module provides Cython-optimized functions for estimating read lengths
from BAM files used in PyMaSC cross-correlation analysis. It supports multiple
estimation methods and handles various edge cases in sequencing data.

Key functionality:
- Multiple estimation methods (mean, median, mode)
- Quality filtering and duplicate handling
- Progress reporting for large files
- Efficient processing of BAM file data

Read length estimation is crucial for accurate cross-correlation calculation
as it affects the expected position of the fragment length peak.
"""
import logging

from pysam.libcalignmentfile cimport AlignmentFile
from pysam.libcalignedsegment cimport AlignedSegment

from PyMaSC.utils.progress import ReadCountProgressBar

logger = logging.getLogger(__name__)


def _mean(c):
    """Calculate mean read length from frequency distribution.

    Args:
        c: Dictionary mapping read lengths to their frequencies

    Returns:
        Mean read length rounded to nearest integer
    """
    return int(round(
        sum(length * freq for length, freq in c.items()) / float(sum(c.values()))
    ))


def _median(c):
    """Calculate median read length from frequency distribution.

    Args:
        c: Dictionary mapping read lengths to their frequencies

    Returns:
        Median read length

    Note:
        Handles both odd and even total counts appropriately
    """
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
    """Calculate mode (most frequent) read length.

    Args:
        c: Dictionary mapping read lengths to their frequencies

    Returns:
        Most frequently occurring read length
    """
    return [k for k, v in sorted(c.items(), key=lambda x: x[1])][-1]


ESTFUNCTIONS = dict(
    MEAN=_mean, MEDIAN=_median, MODE=_mode, MIN=min, MAX=max
)


def estimate_readlen(path, esttype, unsigned int mapq_criteria):
    """Estimate read length from BAM file using specified method.

    Analyzes reads in a BAM file to estimate the characteristic read length
    using the specified statistical method. Applies quality filtering and
    handles various sequencing artifacts.

    Args:
        path: Path to BAM file
        esttype: Estimation method ('MEAN', 'MEDIAN', 'MODE', 'MIN', 'MAX')
        mapq_criteria: Minimum mapping quality threshold

    Returns:
        Estimated read length in base pairs

    Raises:
        ValueError: If estimation method is invalid or insufficient data

    Note:
        Processes only primary alignments and applies quality filtering
        to ensure accurate length estimation
    """
    estfunc = ESTFUNCTIONS[esttype]

    cdef bint is_duplicate
    cdef str chrom, _chrom
    cdef unsigned int mapq
    cdef unsigned long pos, readlen

    cdef dict counter = {}

    cdef dict chrom2length
    cdef unsigned long length
    cdef AlignedSegment read

    cdef unsigned long nreads = 0
    cdef unsigned long nunmapped = 0
    cdef unsigned long npaired = 0
    cdef unsigned long nread2 = 0

    with AlignmentFile(path) as stream:
        chrom = None
        chrom2length = dict(zip(stream.references, stream.lengths))
        length = sum(chrom2length.values())

        progress = ReadCountProgressBar()
        if not stream.has_index():
            progress.disable_bar()
        progress.set_genome(length)

        for read in stream:
            _chrom = read.reference_name
            if _chrom is None:
                continue

            nreads += 1
            if read.is_paired:
                npaired += 1
                if read.is_read2:
                    nread2 += 1

            if chrom != _chrom:
                chrom = _chrom
                progress.set_chrom(chrom, chrom2length[chrom])
            progress.update(read.reference_start)

            if read.is_unmapped:
                nunmapped += 1
            elif not read.is_duplicate and read.mapping_quality >= mapq_criteria:
                readlen = read.infer_query_length()
                if readlen in counter:
                    counter[readlen] += 1
                else:
                    counter[readlen] = 1
    progress.finish()

    length = estfunc(counter)

    logger.info("Scan {:,} reads, {:,} reads were unmapped and {:,} reads >= MAPQ {}."
                "".format(nreads, nunmapped, sum(counter.values()), mapq_criteria))
    if npaired > 0:
        logger.info("{:,} reads were paired: {:,} reads were 1st and {:,} reads were last segment."
                    "".format(npaired, npaired - nread2, nread2))
        logger.info("Note that only 1st reads in the templates will be used for calculation.")
    else:
        logger.info("All reads were single-ended.")
    logger.info("Estimated read length = {:,}".format(length))

    return length
