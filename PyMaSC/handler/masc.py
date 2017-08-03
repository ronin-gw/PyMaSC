import logging

from PyMaSC.reader.align import get_read_generator_and_init_target
from PyMaSC.core.masc import CCCalculator, ReadUnsortedError
from PyMaSC.handler.result import CCResult, ReadsTooFew

logger = logging.getLogger(__name__)


class CCCalcHandler(object):
    default_chi2_p_thresh = 0.05

    def __init__(self, path, fmt, max_shift, mapq_criteria, bwfeeder=None, chi2_pval=None):
        self.path = path
        self.max_shift = max_shift
        self.mapq_criteria = mapq_criteria
        self.bwfeeder = bwfeeder
        self.chi2_p_thresh = chi2_pval if chi2_pval else self.default_chi2_p_thresh

        self._chr = None
        self._wigbuff = []
        self._iter_stopped = False

        # Raise `InputUnseekable`, possibly.
        self.align_parser, self.stream = get_read_generator_and_init_target(path, fmt)

    def _valid_read_generator(self):
        with self.align_parser(self.stream) as alignfile:
            self.ccc = CCCalculator(self.max_shift, alignfile.references, alignfile.lengths,
                                    self.default_chi2_p_thresh, self.bwfeeder is not None)

            for is_reverse, is_duplicate, chrom, pos, mapq, readlen in alignfile:
                if mapq < self.mapq_criteria or is_duplicate:
                    continue

                if is_reverse:
                    read_pos = pos + readlen - 1
                else:
                    read_pos = pos

                # if pos > 100000 and chrom == "chr1":
                #     continue
                yield is_reverse, chrom, pos, readlen, self._get_alignability(chrom, read_pos, is_reverse)

        self.ccc._progress.clean()
        if self.bwfeeder and not self._iter_stopped:
            self.bwfeeder.enable_progress_bar()
            self._feeder.throw(GeneratorExit)

    def _get_alignability(self, chrom, pos, is_reverse):
        if self.bwfeeder is None:
            return False

        if self._chr != chrom:
            if self._chr is not None:
                self.ccc._progress.clean()
                if not self._iter_stopped:
                    self.bwfeeder.enable_progress_bar()
                    try:
                        self._feeder.throw(GeneratorExit)
                    except StopIteration:
                        pass
                    self.bwfeeder.disable_progress_bar()

            self._chr = chrom
            self._wigbuff = []
            self._iter_stopped = False
            try:
                self._feeder = self.bwfeeder.feed(chrom)
            except KeyError as e:
                logger.info("Mappability for '{}' not found. Skip calc mappability sensitive CC.".format(e.args[0]))
                self._iter_stopped = True

        if self._iter_stopped and not self._wigbuff:
            return False

        is_mappable = None
        i = 0
        for i, (start, end) in enumerate(self._wigbuff):
            if start < pos <= end:
                is_mappable = True
                break
            elif pos <= start:
                is_mappable = False
                break

        if not is_reverse:
            self._wigbuff = self._wigbuff[i:]
        if is_mappable is not None:
            return is_mappable

        while True:
            try:
                _chrom, start, end, _val = next(self._feeder)
            except StopIteration:
                self._iter_stopped = True
                return False

            self._wigbuff.append((start, end))

            if start < pos <= end:
                return True
            elif pos <= start:
                return False

    def run_calcuration(self):
        if self.bwfeeder is not None:
            self.bwfeeder.disable_progress_bar()

        for is_reverse, chrom, pos, readlen, is_mappable in self._valid_read_generator():
            try:
                if is_reverse:
                    self.ccc.feed_reverse_read(chrom, pos, readlen, is_mappable)
                else:
                    self.ccc.feed_forward_read(chrom, pos, readlen, is_mappable)
            except ReadUnsortedError:
                logger.error("Input read must be sorted.")
                return None

        if self.bwfeeder is not None:
            self.bwfeeder.enable_progress_bar()
            self.bwfeeder._calc_alignability()

        self.ccc.finishup_calculation()
        try:
            return CCResult(self.ccc, self.path, self.mapq_criteria, self.bwfeeder)
        except ReadsTooFew:
            return None
