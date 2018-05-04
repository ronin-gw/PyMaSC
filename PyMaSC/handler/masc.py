import logging
from multiprocessing import Queue, Lock

import pysam

from PyMaSC.handler.masc_noindex_worker import SingleProcessCalculator
from PyMaSC.handler.masc_worker import NaiveCCCalcWorker, MSCCCalcWorker, NCCandMSCCCalcWorker
from PyMaSC.core.readlen import estimate_readlen
from PyMaSC.utils.progress import MultiLineProgressManager, ProgressBar, ProgressHook
from PyMaSC.utils.calc import exec_worker_pool

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
        self.lengths = list(self.align_file.lengths)

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
            logger.error("Read lengh ({}) seems to be longer than shift size ({}).".format(
                self.read_len, self.max_shift
            ))
            raise ValueError

    def set_mappability_handler(self, mappability_handler):
        self.mappability_handler = mappability_handler
        bw_chromsizes = self.mappability_handler.chromsizes

        for i, reference in enumerate(self.references):
            if reference not in bw_chromsizes:
                logger.debug("mappability for '{}' not found".format(reference))
                continue
            self._compare_refsize(reference, bw_chromsizes[reference])

    def _compare_refsize(self, reference, bw_chr_size):
        i = self.references.index(reference)
        bam_chr_size = self.lengths[i]
        if bw_chr_size != bam_chr_size:
            logger.warning("'{}' reference length mismatch: SAM/BAM -> {:,}, "
                           "BigWig -> {:,}".format(reference, bam_chr_size, bw_chr_size))
            if bam_chr_size < bw_chr_size:
                logger.warning("Use longer length '{:d}' for '{}' anyway".format(
                               bw_chr_size, reference))
                self.lengths[i] = bw_chr_size

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
        self._order_queue = Queue()
        self._report_queue = Queue()
        self._logger_lock = Lock()

        worker_args = [self._order_queue, self._report_queue, self._logger_lock, self.path,
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

        self._run_calculation_with_workers(workers)
        self._calc_unsolved_mappabilty()

    def _run_calculation_with_workers(self, workers):
        _chrom2finished = {c: False for c in self.references}
        progress = MultiLineProgressManager()

        with exec_worker_pool(workers, self.references, self._order_queue):
            while True:
                chrom, obj = self._report_queue.get()
                if chrom is None:  # update progress
                    chrom, body = obj
                    with self._logger_lock:
                        progress.update(chrom, body)
                else:
                    mappable_len, cc_stats, masc_stats = obj
                    self._receive_results(chrom, mappable_len, cc_stats, masc_stats)

                    _chrom2finished[chrom] = True
                    if all(_chrom2finished.values()):
                        break

                    with self._logger_lock:
                        progress.erase(chrom)

        progress.clean()

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
