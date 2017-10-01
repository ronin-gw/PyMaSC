import logging
from multiprocessing import Queue, Lock

import pysam

from PyMaSC.handler.masc_worker import NaiveCCCalcWorker, MSCCCalcWorker, NCCandMSCCCalcWorker
from PyMaSC.core.readlen import estimate_readlen
from PyMaSC.core.mappability import BWFeederWithMappableRegionSum
from PyMaSC.core.ncc import NaiveCCCalculator, ReadUnsortedError
from PyMaSC.core.mscc import MSCCCalculator
from PyMaSC.utils.progress import ProgressBar, MultiLineProgressManager

logger = logging.getLogger(__name__)


class InputUnseekable(Exception):
    pass


class SingleProcessCalculator(object):
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

    def _set_calculator(self):
        if self.skip_ncc:
            assert self.bwfeeder

            self.nccc = None
            self.mscc = MSCCCalculator(self.max_shift, self.read_len,
                                       self.references, self.lengths, self.bwfeeder)
            self._stepping_calc = self._feed2mscc
        elif self.bwfeeder:
            self.nccc = NaiveCCCalculator(self.max_shift, self.references, self.lengths)
            self.mscc = MSCCCalculator(self.max_shift, self.read_len,
                                       self.references, self.lengths, self.bwfeeder)
            self._stepping_calc = self._feed2both
        else:
            self.nccc = NaiveCCCalculator(self.max_shift, self.references, self.lengths)
            self.mscc = None
            self._stepping_calc = self._feed2nccc

    def _run_singleprocess_calculation(self):
        _ref2genomelen = dict(zip(self.references, self.lengths))
        if self.mappability_handler:
            self.bwfeeder = BWFeederWithMappableRegionSum(
                self.mappability_handler.path,
                self.read_len,
                self.mappability_handler.chrom2mappable_len
            )
        else:
            self.bwfeeder = None
        self._set_calculator()

        for read in self.align_file:
            # https://github.com/pysam-developers/pysam/issues/268
            try:
                chrom = read.reference_name
            except ValueError:
                continue
            pos = read.reference_start + 1

            if chrom != self._chr:
                self._progress.clean()
                self._progress.disable_bar()

                if (read.is_read2 or read.mapping_quality < self.mapq_criteria or
                        read.is_unmapped or read.is_duplicate):
                    continue

                self._feed_read(read.is_reverse, read.reference_name,
                                pos, read.query_length)

                self._chr = chrom
                self._progress.enable_bar()
                self._progress.set(chrom, _ref2genomelen[chrom])
                self._progress.update(pos)
            else:
                self._progress.update(pos)
                if read.mapping_quality < self.mapq_criteria or read.is_duplicate:
                    continue

                self._feed_read(read.is_reverse, read.reference_name,
                                pos, read.query_length)

        self._progress.clean()
        self.align_file.close()

        if not self.skip_ncc:
            self.nccc.finishup_calculation()
            self.ref2forward_sum = self.nccc.ref2forward_sum
            self.ref2reverse_sum = self.nccc.ref2reverse_sum
            self.ref2ccbins = self.nccc.ref2ccbins

        if self.mscc:
            self.mscc.finishup_calculation()
            self.ref2mappable_len = self.bwfeeder.chrom2mappable_len
            self.mappable_ref2forward_sum = self.mscc.ref2forward_sum
            self.mappable_ref2reverse_sum = self.mscc.ref2reverse_sum
            self.mappable_ref2ccbins = self.mscc.ref2ccbins
            self.mappability_handler.chrom2mappable_len = self.bwfeeder.chrom2mappable_len

    def _feed_read(self, is_reverse, chrom, pos, readlen):
        try:
            self._stepping_calc(is_reverse, chrom, pos, readlen)
        except ReadUnsortedError:
            self._progress.clean()
            logger.error("Input alignment file must be sorted.")
            raise


class CCCalcHandler(SingleProcessCalculator):
    def __init__(self, path, esttype, max_shift, mapq_criteria, nworker=1, skip_ncc=False):
        self.path = path
        self.esttype = esttype
        self.max_shift = max_shift
        self.mapq_criteria = mapq_criteria
        self.nworker = nworker
        self.skip_ncc = skip_ncc

        #
        try:
            self.align_file = pysam.AlignmentFile(path)
        except ValueError:
            logger.error("File has no sequences defined.")
            raise
        self.references = self.align_file.references
        self.lengths = self.align_file.lengths

        if not self.align_file.has_index() and self.nworker > 1:
            logger.error("Need indexed alignment file for multi-processng. "
                         "Calculation will be executed by single process.")
            self.nworker = 1
        elif self.nworker > 1:
            self.align_file.close()

        #
        self.ref2forward_sum = {}
        self.ref2reverse_sum = {}
        self.ref2ccbins = {}

        #
        self.mappable_ref2forward_sum = {}
        self.mappable_ref2reverse_sum = {}
        self.mappable_ref2ccbins = {}
        self.ref2mappable_len = {}

        #
        self.read_len = None
        self.mappability_handler = None

        # for single process progress bar
        self._chr = None
        self._progress = ProgressBar()

    def set_readlen(self, readlen=None):
        if readlen:
            self.read_len = readlen
        elif self.path == '-':
            logger.error("Cannot execute read length checking for unseekable input.")
            raise InputUnseekable
        else:
            logger.info("Check read length... : " + self.path)
            self.read_len = estimate_readlen(self.path, self.esttype, self.mapq_criteria)

        if self.read_len > self.max_shift:
            logger.warning("Read lengh ({}) seems to be longer than shift size ({}).".format(
                self.read_len, self.max_shift
            ))

    def set_mappability_handler(self, mappability_handler):
        self.mappability_handler = mappability_handler

    def run_calcuration(self):
        if self.nworker > 1:
            self._run_multiprocess_calcuration()
        else:
            self._run_singleprocess_calculation()

    def _run_multiprocess_calcuration(self):
        _chrom2finished = {c: False for c in self.references}

        order_queue = Queue()
        report_queue = Queue()
        logger_lock = Lock()
        progress = MultiLineProgressManager()

        if self.mappability_handler is None:
            workers = [
                NaiveCCCalcWorker(
                    order_queue, report_queue, logger_lock, self.path,
                    self.mapq_criteria, self.max_shift, self.references, self.lengths
                ) for _ in range(min(self.nworker, len(self.references)))
            ]
        elif self.skip_ncc:
            workers = [
                MSCCCalcWorker(
                    order_queue, report_queue, logger_lock, self.path,
                    self.mapq_criteria, self.max_shift, self.references, self.lengths,
                    self.mappability_handler.path, self.read_len,
                    self.mappability_handler.chrom2mappable_len
                ) for _ in range(min(self.nworker, len(self.references)))
            ]
        else:
            workers = [
                NCCandMSCCCalcWorker(
                    order_queue, report_queue, logger_lock, self.path,
                    self.mapq_criteria, self.max_shift, self.references, self.lengths,
                    self.mappability_handler.path, self.read_len,
                    self.mappability_handler.chrom2mappable_len
                ) for _ in range(min(self.nworker, len(self.references)))
            ]

        for w in workers:
            w.start()
        for chrom in self.references:
            order_queue.put(chrom)
        for _ in range(self.nworker):
            order_queue.put(None)

        try:
            while True:
                chrom, obj = report_queue.get()
                if chrom is None:  # update progress
                    chrom, body = obj
                    with logger_lock:
                        progress.update(chrom, body)
                else:
                    _m_len, (_f_sum, _r_sum, _ccbins), (_mf_sum, _mr_sum, _m_ccbins) = obj

                    if _m_len is not None:
                        self.mappability_handler.chrom2mappable_len[chrom] = _m_len
                        self.mappability_handler.chrom2is_called[chrom] = True

                    if None not in (_f_sum, _r_sum, _ccbins):
                        self.ref2forward_sum[chrom] = _f_sum
                        self.ref2reverse_sum[chrom] = _r_sum
                        self.ref2ccbins[chrom] = _ccbins

                    if None not in (_m_len, _mf_sum, _mr_sum, _m_ccbins):
                        self.mappable_ref2forward_sum[chrom] = _mf_sum
                        self.mappable_ref2reverse_sum[chrom] = _mr_sum
                        self.mappable_ref2ccbins[chrom] = _m_ccbins

                    _chrom2finished[chrom] = True
                    if all(_chrom2finished.values()):
                        break

                    with logger_lock:
                        progress.erase(chrom)
        except:
            for w in workers:
                if w.is_alive():
                    w.terminate()
            raise

        progress.clean()

        if self.mappability_handler is not None and not self.mappability_handler.is_called:
            self.mappability_handler.is_called = all(self.mappability_handler.chrom2is_called.values())
            self.mappability_handler.calc_mappability()
