import logging
import os
from multiprocessing import Process, Queue, Lock

from PyMaSC.reader.align import get_read_generator_and_init_target
from PyMaSC.core.readlen import estimate_readlen
from PyMaSC.core.mappability import BWFeederWithMappableRegionSum
from PyMaSC.core.ncc import NaiveCCCalculator, ReadUnsortedError
from PyMaSC.core.mscc import MSCCCalculator
from PyMaSC.handler.mappability import MappabilityHandler
from PyMaSC.utils.progress import ProgressHook

logger = logging.getLogger(__name__)


class CCCalcHandler(object):
    def __init__(self, path, fmt, esttype, max_shift, mapq_criteria, threads=1):
        self.path = path
        self.max_shift = max_shift
        self.mapq_criteria = mapq_criteria

        self.esttype = esttype
        self.mapq_criteria = mapq_criteria
        self.threads = threads

        # Raise `InputUnseekable`, possibly.
        self.align_parser, self.stream = get_read_generator_and_init_target(path, fmt)

        #
        self.references = self.lengths = None

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
        else:
            logger.info("Check read length... : " + self.path)
            self.read_len = estimate_readlen(self.align_parser, self.stream,
                                             self.esttype, self.mapq_criteria)

        if self.read_len > self.max_shift:
            logger.warning("Read lengh ({}) seems to be longer than shift size ({}).".format(
                self.read_len, self.max_shift
            ))

    def set_mappability_handler(self, mappability_handler):
        self.mappability_handler = mappability_handler

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
            _ref2genomelen = dict(zip(alignfile.references, alignfile.lengths))

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
        with self.align_parser(self.stream) as alignfile:
            self.references = alignfile.references
            self.lengths = alignfile.lengths
        _chrom2finished = {c: False for c in self.references}

        order_queue = Queue()
        report_queue = Queue()
        logger_lock = Lock()

        if self.mappability_handler is None:
            workers = [
                NaiveCCCalcWorker(
                    order_queue, report_queue, logger_lock, self.align_parser, self.stream,
                    self.mapq_criteria, self.max_shift, self.references, self.lengths
                ) for _ in range(min(self.threads, len(self.references)))
            ]
        else:
            workers = [
                NCCandMSCCCalcWorker(
                    order_queue, report_queue, logger_lock, self.align_parser, self.stream,
                    self.mapq_criteria, self.max_shift, self.references, self.lengths,
                    self.mappability_handler.path, self.read_len,
                    self.mappability_handler.chrom2mappable_len
                ) for _ in range(min(self.threads, len(self.references)))
            ]

        for w in workers:
            w.start()
        for chrom in self.references:
            order_queue.put(chrom)
        for _ in range(self.threads):
            order_queue.put(None)

        try:
            while True:
                chrom, obj = report_queue.get()
                if chrom is None:  # update progress
                    pass
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
        except:
            for w in workers:
                if w.is_alive():
                    w.terminate()
            raise

        if self.mappability_handler is not None and not self.mappability_handler.is_called:
            self.mappability_handler.is_called = all(self.mappability_handler.chrom2is_called.values())
            self.mappability_handler.calc_mappability()


class CalcWorkerBase(Process):
    def __init__(self, order_queue, report_queue, logger_lock, align_parser, stream,
                 mapq_criteria, max_shift, references, lengths):
        super(CalcWorkerBase, self).__init__()

        self.order_queue = order_queue
        self.report_queue = report_queue
        self.logger_lock = logger_lock

        self.align_parser = align_parser
        self.stream = stream
        self.mapq_criteria = mapq_criteria

        self._ref2genomelen = dict(zip(references, lengths))

        self.calculator = None

        self._progress = ProgressHook(report_queue)
        self._chr = None

    def run(self):
        with self.logger_lock:
            logger.debug("{}: Hello. My pid is {}.".format(self.name, os.getpid()))

        with self.align_parser(self.stream) as alignfile:
            while True:
                chrom = self.order_queue.get()
                if chrom is None:
                    break
                elif chrom != self._chr:
                    self._chr = chrom
                    self._progress.set(self._ref2genomelen[chrom])

                with self.logger_lock:
                    logger.debug("{}: Process {}...".format(self.name, chrom))

                for read in alignfile.fetch(chrom):
                    self._progress.update(read.pos)

                    if read.mapping_quality < self.mapq_criteria or read.is_duplicate:
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
    def __init__(self, order_queue, report_queue, logger_lock, align_parser, stream,
                 mapq_criteria, max_shift, references, lengths):
        super(NaiveCCCalcWorker, self).__init__(
            order_queue, report_queue, logger_lock, align_parser, stream,
            mapq_criteria, max_shift, references, lengths
        )

        self.calculator = NaiveCCCalculator(max_shift, references, lengths)

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
    def __init__(self, order_queue, report_queue, logger_lock, align_parser, stream,
                 mapq_criteria, max_shift, references, lengths,
                 mappable_path, read_len, chrom2mappable_len):
        super(MSCCCalcWorker, self).__init__(
            order_queue, report_queue, logger_lock, align_parser, stream,
            mapq_criteria, max_shift, references, lengths
        )

        mappability_shift = MappabilityHandler.calc_mappable_len_required_shift_size(read_len, max_shift)
        self.bwfeeder = BWFeederWithMappableRegionSum(mappable_path, mappability_shift, chrom2mappable_len)
        self.calculator = MSCCCalculator(max_shift, read_len, references, lengths, self.bwfeeder)

        self.bwfeeder._progress = self._progress

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
    def __init__(self, order_queue, report_queue, logger_lock, align_parser, stream,
                 mapq_criteria, max_shift, references, lengths,
                 mappable_path, read_len, chrom2mappable_len):
        super(NCCandMSCCCalcWorker, self).__init__(
            order_queue, report_queue, logger_lock, align_parser, stream,
            mapq_criteria, max_shift, references, lengths,
            mappable_path, read_len, chrom2mappable_len
        )

        self.ncc_calculator = NaiveCCCalculator(max_shift, references, lengths)
        self.mscc_calculator = self.calculator

        self.bwfeeder._progress = self._progress

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
