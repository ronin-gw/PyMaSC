import logging
from multiprocessing import Queue, Lock

from PyMaSC.handler.masc import CCCalcHandler
from PyMaSC.bacore.mscc import CCBitArrayCalculator
from PyMaSC.handler.masc_noindex_worker import SingleProcessCalculator
from PyMaSC.handler.masc_worker import CalcWorkerBase
from PyMaSC.reader.bigwig import BigWigReader
from PyMaSC.utils.progress import ProgressBar, ProgressHook

logger = logging.getLogger(__name__)


class BACalcHandler(CCCalcHandler):
    def _run_singleprocess_calculation(self):
        worker = BASingleProcessCalculator(
            self.align_file, self.mapq_criteria, self.max_shift, self.references,
            self.lengths, self.read_len, self.mappability_handler, self.skip_ncc
        )

        worker.run()

        if not self.skip_ncc:
            self.ref2forward_sum = worker.calculator.ref2forward_sum
            self.ref2reverse_sum = worker.calculator.ref2reverse_sum
            self.ref2ccbins = worker.calculator.ref2ccbins

        if self.mappability_handler:
            self.mappable_ref2forward_sum = worker.calculator.ref2mappable_forward_sum
            self.mappable_ref2reverse_sum = worker.calculator.ref2mappable_reverse_sum
            self.mappable_ref2ccbins = worker.calculator.ref2mascbins

            mappable_len = {k: v for k, v in worker.calculator.ref2mappable_len.items() if v is not None}
            self.ref2mappable_len = mappable_len
            self.mappability_handler.chrom2mappable_len = mappable_len
            self.mappability_handler.mappable_len = list(map(sum, zip(*mappable_len.values())))

    def _run_multiprocess_calcuration(self):
        self._order_queue = Queue()
        self._report_queue = Queue()
        self._logger_lock = Lock()

        workers = [
            BACalcWorker(
                self._order_queue, self._report_queue, self._logger_lock,
                self.path, self.mapq_criteria, self.max_shift, self.read_len,
                self.references, self.lengths, self.mappability_handler, self.skip_ncc
            ) for _ in range(min(self.nworker, len(self.references)))
        ]

        self._run_calculation_with_workers(workers)
        self._calc_unsolved_mappabilty()


class BASingleProcessCalculator(SingleProcessCalculator):
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
            self._bwfeeder = BigWigReader(mappability_handler.path)
        else:
            self._bwfeeder = None

        #
        self.calculator = CCBitArrayCalculator(
            self.max_shift, self.read_len, self.references, self.lengths,
            self._bwfeeder, self.skip_ncc
        )

    def _stepping_calc(self, is_reverse, chrom, pos, readlen):
        if is_reverse:
            self.calculator.feed_reverse_read(chrom, pos, readlen)
        else:
            self.calculator.feed_forward_read(chrom, pos, readlen)

    def _deconstruct(self):
        self.calculator.finishup_calculation()
        if self._bwfeeder:
            self._bwfeeder.close()


class BACalcWorker(CalcWorkerBase):
    def __init__(self, order_queue, report_queue, logger_lock,
                 path, mapq_criteria, max_shift, read_len, references, lengths,
                 mappability_handler=None, skip_ncc=False):
        super(BACalcWorker, self).__init__(
            order_queue, report_queue, logger_lock, path,
            mapq_criteria, max_shift, references, lengths
        )

        self.skip_ncc = skip_ncc
        self._bwfeeder = BigWigReader(mappability_handler.path) if mappability_handler else None
        self.calculator = CCBitArrayCalculator(
            max_shift, read_len, references, lengths, self._bwfeeder, self.skip_ncc,
            logger_lock, ProgressHook(report_queue)
        )

    def _put_result_to_report_queue(self, chrom):
        self.calculator.flush()

        self.report_queue.put((
            chrom, (
                self.calculator.ref2mappable_len.get(chrom), (
                    self.calculator.ref2forward_sum[chrom],
                    self.calculator.ref2reverse_sum[chrom],
                    self.calculator.ref2ccbins[chrom],
                ), (
                    self.calculator.ref2mappable_forward_sum[chrom],
                    self.calculator.ref2mappable_reverse_sum[chrom],
                    self.calculator.ref2mascbins[chrom],
                )
            )
        ))

    def _deconstruct(self):
        if self._bwfeeder:
            self._bwfeeder.close()
