"""Single-process calculator for non-indexed BAM files.

This module provides a single-process implementation for cross-correlation
calculation when BAM files are not indexed or when multiprocessing is not
desired. It handles the complete workflow sequentially across all chromosomes.

Key functionality:
- Sequential processing of all chromosomes
- Support for both NCC and MSCC calculations
- Progress reporting for long-running operations
- Efficient memory usage for single-process execution

This implementation is used as a fallback when parallel processing
is not available or when working with streaming/unseekable inputs.
"""
import logging

from PyMaSC.core.mappability import BWFeederWithMappableRegionSum
from PyMaSC.core.ncc import NaiveCCCalculator, ReadUnsortedError
from PyMaSC.core.mscc import MSCCCalculator
from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


class SingleProcessCalculator(object):
    """Single-process cross-correlation calculator.
    
    Implements cross-correlation calculation in a single process,
    processing all chromosomes sequentially. This class is used when
    BAM files are not indexed or when multiprocessing is disabled.
    
    The calculator handles both naive cross-correlation and MSCC
    calculations depending on the availability of mappability data.
    It provides progress reporting and efficient memory management
    for single-threaded execution.
    
    Attributes:
        align_file: Opened BAM file object
        mapq_criteria: Minimum mapping quality threshold
        max_shift: Maximum shift distance for correlation
        references: List of chromosome names to process
        lengths: List of chromosome lengths
        read_len: Read length for calculations
        mappability_handler: Optional mappability correction handler
        skip_ncc: Whether to skip naive cross-correlation
        nccc: NaiveCCCalculator instance (if not skipped)
        mscc: MSCCCalculator instance (if mappability available)
    """
    def __init__(self, align_file, mapq_criteria, max_shift, references, lengths,
                 read_len, mappability_handler=None, skip_ncc=False):
        self.align_file = align_file
        self.mapq_criteria = mapq_criteria
        self.max_shift = max_shift
        self.references = references
        self.lengths = lengths
        self.read_len = read_len

        self.mappability_handler = mappability_handler
        self.skip_ncc = skip_ncc

        # for single process progress bar
        self._chr = None
        self._progress = ProgressBar()

        #
        if self.mappability_handler:
            self._bwfeeder = BWFeederWithMappableRegionSum(
                self.mappability_handler.path,
                self.read_len,
                self.mappability_handler.chrom2mappable_len
            )
        else:
            self._bwfeeder = None

        # set calculators
        if self.skip_ncc:
            self.nccc = None
            self.mscc = MSCCCalculator(self.max_shift, self.read_len,
                                       self.references, self.lengths, self._bwfeeder)
            self._stepping_calc = self._feed2mscc
        elif self._bwfeeder:
            self.nccc = NaiveCCCalculator(self.max_shift, self.references, self.lengths)
            self.mscc = MSCCCalculator(self.max_shift, self.read_len,
                                       self.references, self.lengths, self._bwfeeder)
            self._stepping_calc = self._feed2both
        else:
            self.nccc = NaiveCCCalculator(self.max_shift, self.references, self.lengths)
            self.mscc = None
            self._stepping_calc = self._feed2nccc

    def _feed2nccc(self, is_reverse, chrom, pos, readlen):
        if is_reverse:
            self.nccc.feed_reverse_read(chrom, pos, readlen)
        else:
            self.nccc.feed_forward_read(chrom, pos, readlen)

    def _feed2mscc(self, is_reverse, chrom, pos, readlen):
        if is_reverse:
            self.mscc.feed_reverse_read(chrom, pos, readlen)
        else:
            self.mscc.feed_forward_read(chrom, pos, readlen)

    def _feed2both(self, is_reverse, chrom, pos, readlen):
        if is_reverse:
            self.nccc.feed_reverse_read(chrom, pos, readlen)
            self.mscc.feed_reverse_read(chrom, pos, readlen)
        else:
            self.nccc.feed_forward_read(chrom, pos, readlen)
            self.mscc.feed_forward_read(chrom, pos, readlen)

    def run(self):
        _ref2genomelen = dict(zip(self.references, self.lengths))

        for read in self.align_file:
            chrom = read.reference_name
            if chrom is None or chrom not in self.references:
                continue

            if chrom != self._chr:
                self._progress.clean()
                self._progress.disable_bar()

                self._feed_read(read)

                self._chr = chrom
                self._progress.enable_bar()
                self._progress.set(chrom, _ref2genomelen[chrom])
                self._progress.update(read.reference_start)
            else:
                self._progress.update(read.reference_start)
                self._feed_read(read)

        self._progress.clean()
        self.align_file.close()
        self._deconstruct()

    def _feed_read(self, read):
        if (read.is_read2 or read.mapping_quality < self.mapq_criteria or
                read.is_unmapped or read.is_duplicate):
            return

        try:
            self._stepping_calc(read.is_reverse, read.reference_name,
                                read.reference_start + 1, read.infer_query_length())
        except ReadUnsortedError:
            logger.error("Input alignment file must be sorted.")
            self.align_file.close()
            self._bwfeeder.close()
            self._progress.clean()
            raise

    def _deconstruct(self):
        if not self.skip_ncc:
            self.nccc.finishup_calculation()

        if self.mscc:
            self.mscc.finishup_calculation()
            self._bwfeeder.close()
            self.ref2mappable_len = self._bwfeeder.chrom2mappable_len
            self.mappability_handler.chrom2mappable_len = self._bwfeeder.chrom2mappable_len
            self.mappability_handler.mappable_len = self._bwfeeder.mappable_len
