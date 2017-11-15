import logging
from multiprocessing import Queue, Lock

import pysam

from PyMaSC.handler.masc_noindex_worker import SingleProcessCalculator
from PyMaSC.handler.masc_worker import NaiveCCCalcWorker, MSCCCalcWorker, NCCandMSCCCalcWorker
from PyMaSC.core.readlen import estimate_readlen
from PyMaSC.utils.progress import MultiLineProgressManager, ProgressBar, ProgressHook

logger = logging.getLogger(__name__)


class InputUnseekable(Exception):
    pass


class CCCalcHandler(object):
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
            # to avoid interfering mappability progress bar with multiline progress bar
            ProgressBar.global_switch = False
            ProgressHook.global_switch = True
            self._run_multiprocess_calcuration()
        else:
            self._run_singleprocess_calculation()

    def _run_singleprocess_calculation(self):
        worker = SingleProcessCalculator(
            self.align_file, self.mapq_criteria, self.max_shift, self.references,
            self.lengths, self.read_len, self.mappability_handler, self.skip_ncc
        )

        worker.run()

        if not self.skip_ncc:
            self.ref2forward_sum = worker.nccc.ref2forward_sum
            self.ref2reverse_sum = worker.nccc.ref2reverse_sum
            self.ref2ccbins = worker.nccc.ref2ccbins

        if worker.mscc:
            self.mappable_ref2forward_sum = worker.mscc.ref2forward_sum
            self.mappable_ref2reverse_sum = worker.mscc.ref2reverse_sum
            self.mappable_ref2ccbins = worker.mscc.ref2ccbins

        self._calc_unsolved_mappabilty()

    def _run_multiprocess_calcuration(self):
        _chrom2finished = {c: False for c in self.references}

        order_queue = Queue()
        report_queue = Queue()
        logger_lock = Lock()
        progress = MultiLineProgressManager()

        worker_args = [order_queue, report_queue, logger_lock, self.path,
                       self.mapq_criteria, self.max_shift, self.references, self.lengths]
        if self.mappability_handler:
            worker_args += [self.mappability_handler.path, self.read_len,
                            self.mappability_handler.chrom2mappable_len]
            if self.skip_ncc:
                worker_class = MSCCCalcWorker
            else:
                worker_class = NCCandMSCCCalcWorker
        else:
            worker_class = NaiveCCCalcWorker

        workers = [worker_class(*worker_args)
                   for _ in range(min(self.nworker, len(self.references)))]

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
                    mappable_len, cc_stats, masc_stats = obj
                    self._receive_results(chrom, mappable_len, cc_stats, masc_stats)

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
        self._calc_unsolved_mappabilty()

    def _receive_results(self, chrom, mappable_len, cc_stats, masc_stats):
        f_sum, r_sum, ccbins = cc_stats
        mf_sum, mr_sum, mccbins = masc_stats

        if mappable_len is not None:
            self.mappability_handler.chrom2mappable_len[chrom] = mappable_len
            self.mappability_handler.chrom2is_called[chrom] = True

        if None not in (f_sum, r_sum, ccbins):
            self.ref2forward_sum[chrom] = f_sum
            self.ref2reverse_sum[chrom] = r_sum
            self.ref2ccbins[chrom] = ccbins

        if None not in (mappable_len, mf_sum, mr_sum, mccbins):
            self.mappable_ref2forward_sum[chrom] = mf_sum
            self.mappable_ref2reverse_sum[chrom] = mr_sum
            self.mappable_ref2ccbins[chrom] = mccbins

    def _calc_unsolved_mappabilty(self):
        if self.mappability_handler is not None:
            if not self.mappability_handler.is_called:
                self.mappability_handler.is_called = all(
                    self.mappability_handler.chrom2is_called.values()
                )
                self.mappability_handler.calc_mappability()
            self.ref2mappable_len = self.mappability_handler.chrom2mappable_len
