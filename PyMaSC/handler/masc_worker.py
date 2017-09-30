import logging
import os
from multiprocessing import Process

from pysam import AlignmentFile

from PyMaSC.core.mappability import BWFeederWithMappableRegionSum
from PyMaSC.core.ncc import NaiveCCCalculator, ReadUnsortedError
from PyMaSC.core.mscc import MSCCCalculator
from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.progress import ProgressHook

logger = logging.getLogger(__name__)


class CalcWorkerBase(Process):
    def __init__(self, order_queue, report_queue, logger_lock, path,
                 mapq_criteria, max_shift, references, lengths):
        super(CalcWorkerBase, self).__init__()

        self.order_queue = order_queue
        self.report_queue = report_queue
        self.logger_lock = logger_lock

        self.path = path
        self.mapq_criteria = mapq_criteria

        self._ref2genomelen = dict(zip(references, lengths))

        self.calculator = None

        self._progress = ProgressHook(report_queue)
        self._chr = None

    def run(self):
        with self.logger_lock:
            logger.debug("{}: Hello. My pid is {}.".format(self.name, os.getpid()))

        with AlignmentFile(self.path) as alignfile:
            while True:
                chrom = self.order_queue.get()
                if chrom is None:
                    break
                elif chrom != self._chr:
                    self._chr = chrom
                    self._progress.set(chrom, self._ref2genomelen[chrom])

                with self.logger_lock:
                    logger.debug("{}: Process {}...".format(self.name, chrom))

                for read in alignfile.fetch(chrom):
                    self._progress.update(read.pos)

                    if (read.is_read2 or read.mapping_quality < self.mapq_criteria or
                            read.is_unmapped or read.is_duplicate):
                        continue

                    try:
                        self._feed_read(read)
                    except ReadUnsortedError:
                        logger.error("Input alignment file must be sorted.")
                        raise

                self._put_result_to_report_queue(chrom)

        with self.logger_lock:
            logger.debug("{}: Goodbye.".format(self.name))
        self._deconstructor()

    def _feed_read(self, read):
        if read.is_reverse:
            self.calculator.feed_reverse_read(
                read.reference_name,
                read.reference_start + 1,
                read.query_length
            )
        else:
            self.calculator.feed_forward_read(
                read.reference_name,
                read.reference_start + 1,
                read.query_length
            )

    def _put_result_to_report_queue(self, chrom):
        pass

    def _deconstructor(self):
        pass


class NaiveCCCalcWorker(CalcWorkerBase):
    def __init__(self, order_queue, report_queue, logger_lock, path,
                 mapq_criteria, max_shift, references, lengths):
        super(NaiveCCCalcWorker, self).__init__(
            order_queue, report_queue, logger_lock, path,
            mapq_criteria, max_shift, references, lengths
        )

        self.calculator = NaiveCCCalculator(max_shift, references, lengths, logger_lock)

    def _put_result_to_report_queue(self, chrom):
        self.calculator.flush()

        self.report_queue.put((
            chrom, (
                None, (
                    self.calculator.ref2forward_sum[chrom],
                    self.calculator.ref2reverse_sum[chrom],
                    self.calculator.ref2ccbins[chrom]
                ), (
                    None,
                    None,
                    None
                )
            )
        ))


class MSCCCalcWorker(NaiveCCCalcWorker):
    def __init__(self, order_queue, report_queue, logger_lock, path,
                 mapq_criteria, max_shift, references, lengths,
                 mappable_path, read_len, chrom2mappable_len):
        super(MSCCCalcWorker, self).__init__(
            order_queue, report_queue, logger_lock, path,
            mapq_criteria, max_shift, references, lengths
        )

        mappability_shift = MappabilityHandler.calc_mappable_len_required_shift_size(read_len, max_shift)
        self.bwfeeder = BWFeederWithMappableRegionSum(mappable_path, mappability_shift, chrom2mappable_len)
        self.calculator = MSCCCalculator(max_shift, read_len, references, lengths, self.bwfeeder, logger_lock)

        self.bwfeeder.disable_progress_bar()
        self._progress.enable_bar = ProgressHook._pass

    def _put_result_to_report_queue(self, chrom):
        self.calculator.flush()
        self.bwfeeder.flush()

        self.report_queue.put((
            chrom, (
                self.bwfeeder.chrom2mappable_len.get(chrom), (
                    None,
                    None,
                    None
                ), (
                    self.calculator.ref2forward_sum[chrom],
                    self.calculator.ref2reverse_sum[chrom],
                    self.calculator.ref2ccbins[chrom],
                )
            )
        ))

    def _deconstructor(self):
        self.bwfeeder.close()


class NCCandMSCCCalcWorker(MSCCCalcWorker):
    def __init__(self, order_queue, report_queue, logger_lock, path,
                 mapq_criteria, max_shift, references, lengths,
                 mappable_path, read_len, chrom2mappable_len):
        super(NCCandMSCCCalcWorker, self).__init__(
            order_queue, report_queue, logger_lock, path,
            mapq_criteria, max_shift, references, lengths,
            mappable_path, read_len, chrom2mappable_len
        )

        self.ncc_calculator = NaiveCCCalculator(max_shift, references, lengths, logger_lock)
        self.mscc_calculator = self.calculator

        self.bwfeeder.disable_progress_bar()
        self._progress.enable_bar = ProgressHook._pass

    def _feed_read(self, read):
        if read.is_reverse:
            self.ncc_calculator.feed_reverse_read(
                read.reference_name,
                read.reference_start + 1,
                read.query_length
            )
            self.mscc_calculator.feed_reverse_read(
                read.reference_name,
                read.reference_start + 1,
                read.query_length
            )
        else:
            self.ncc_calculator.feed_forward_read(
                read.reference_name,
                read.reference_start + 1,
                read.query_length
            )
            self.mscc_calculator.feed_forward_read(
                read.reference_name,
                read.reference_start + 1,
                read.query_length
            )

    def _put_result_to_report_queue(self, chrom):
        self.ncc_calculator.flush()
        self.mscc_calculator.flush()
        self.bwfeeder.flush()

        self.report_queue.put((
            chrom, (
                self.bwfeeder.chrom2mappable_len.get(chrom), (
                    self.ncc_calculator.ref2forward_sum[chrom],
                    self.ncc_calculator.ref2reverse_sum[chrom],
                    self.ncc_calculator.ref2ccbins[chrom],
                ), (
                    self.mscc_calculator.ref2forward_sum[chrom],
                    self.mscc_calculator.ref2reverse_sum[chrom],
                    self.mscc_calculator.ref2ccbins[chrom],
                )
            )
        ))
