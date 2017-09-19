import logging

from PyMaSC.reader.align import get_read_generator_and_init_target
from PyMaSC.core.readlen import estimate_readlen
from PyMaSC.core.ncc import NaiveCCCalculator, ReadUnsortedError
from PyMaSC.core.mscc import MSCCCalculator
from PyMaSC.handler.result import CCResult, ReadsTooFew
from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


class CCCalcHandler(object):
    default_chi2_p_thresh = 0.05

    def __init__(self, path, fmt, esttype, max_shift, mapq_criteria, filter_len, readlen=None, bwfeeder=None, chi2_pval=None):
        self.path = path
        self.max_shift = max_shift
        self.mapq_criteria = mapq_criteria
        self.filter_len = filter_len
        self.bwfeeder = bwfeeder
        self.chi2_p_thresh = chi2_pval if chi2_pval else self.default_chi2_p_thresh

        # Raise `InputUnseekable`, possibly.
        self.align_parser, self.stream = get_read_generator_and_init_target(path, fmt)
        if readlen:
            self.read_len = readlen
        else:
            self.read_len = estimate_readlen(self.align_parser, self.stream, esttype, mapq_criteria)

        if self.read_len > max_shift:
            logger.warning("Read lengh ({}) seems to be longer than shift size ({}).".format(
                self.read_len, max_shift
            ))

        # to display progress bar
        self._chr = None
        self._progress = ProgressBar()

    def _set_calculator(self, references, lengths):
        self.nccc = NaiveCCCalculator(self.max_shift, references, lengths)

        if self.bwfeeder:
            self.mscc = MSCCCalculator(self.max_shift, self.read_len, references, lengths, self.bwfeeder)
            self._stepping_calc = self._feed2both
        else:
            self.mscc = None
            self._stepping_calc = self._feed2nccc

    def _feed2nccc(self, is_reverse, chrom, pos, readlen):
        if is_reverse:
            self.nccc.feed_reverse_read(chrom, pos, readlen)
        else:
            self.nccc.feed_forward_read(chrom, pos, readlen)

    def _feed2both(self, is_reverse, chrom, pos, readlen):
        if is_reverse:
            self.nccc.feed_reverse_read(chrom, pos, readlen)
            self.mscc.feed_reverse_read(chrom, pos, readlen)
        else:
            self.nccc.feed_forward_read(chrom, pos, readlen)
            self.mscc.feed_forward_read(chrom, pos, readlen)

    def _valid_read_generator(self):
        with self.align_parser(self.stream) as alignfile:
            self._set_calculator(alignfile.references, alignfile.lengths)
            _ref2genomelen = {r: l for r, l in zip(alignfile.references, alignfile.lengths)}

            for is_reverse, is_duplicate, chrom, pos, mapq, readlen in alignfile:
                if chrom != self._chr:
                    self._progress.clean()
                    self._progress.disable_bar()

                    if mapq < self.mapq_criteria or is_duplicate:
                        continue
                    yield is_reverse, chrom, pos, readlen

                    self._chr = chrom
                    self._progress.enable_bar()
                    self._progress.set(_ref2genomelen[chrom])
                    self._progress.update(pos)
                else:
                    self._progress.update(pos)
                    if mapq < self.mapq_criteria or is_duplicate:
                        continue
                    yield is_reverse, chrom, pos, readlen
            self._progress.clean()

    def run_calcuration(self):
        for is_reverse, chrom, pos, readlen in self._valid_read_generator():
            try:
                self._stepping_calc(is_reverse, chrom, pos, readlen)
            except ReadUnsortedError:
                logger.error("Input alignment file must be sorted.")
                return None

        self.nccc.finishup_calculation()
        if self.mscc:
            self.mscc.finishup_calculation()

        try:
            return CCResult(self.nccc, self.mscc, self.path, self.read_len, self.mapq_criteria, self.bwfeeder, self.filter_len)
        except ReadsTooFew:
            return None
