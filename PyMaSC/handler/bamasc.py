"""BitArray-based cross-correlation calculation handler.

This module provides BitArray-based implementations of cross-correlation
calculation handlers for memory-efficient processing. The BitArray approach
uses binary arrays to represent genomic positions, offering faster processing
and lower memory usage compared to the successive algorithm approach.

Key components:
- BACalcHandler: Main handler extending CCCalcHandler for BitArray algorithm
- BASingleProcessCalculator: Single-process BitArray calculator
- BACalcWorker: Multiprocess BitArray worker for parallel processing

The BitArray implementation is the default algorithm for PyMaSC due to its
performance advantages, though it requires more memory (~250MB per worker).

NOTE: This module is deprecated. Use PyMaSC.core.worker.UnifiedWorker instead.
Legacy workers are available for backward compatibility through compatibility wrappers.
"""
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
    """BitArray calculation handler extending base handler.

    This class extends CCCalcHandler to use the BitArray-based cross-correlation
    algorithm. It provides both single-process and multiprocess execution modes
    using the CCBitArrayCalculator for efficient computation.

    The BitArray approach represents genomic positions as binary arrays,
    enabling fast bitwise operations for cross-correlation calculation.
    Inherits all functionality from CCCalcHandler while providing BitArray-specific
    implementations of the calculation methods.
    """
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

    def _run_multiprocess_calculation(self):
        """Execute calculation using multiprocess mode with BitArray workers.

        Overrides parent method to use BitArray-specific workers while
        leveraging the parent's queue creation and workflow coordination.
        """
        # Use parent's queue creation
        self._create_worker_queues()

        # Create BitArray-specific workers
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
    """Single-process BitArray calculator.

    This class implements single-process cross-correlation calculation
    using the BitArray algorithm. It extends SingleProcessCalculator
    to use CCBitArrayCalculator for the core computation.

    The calculator processes reads sequentially across all chromosomes,
    feeding them to the BitArray calculator which maintains binary
    arrays for efficient position tracking and correlation computation.

    Attributes:
        calculator: CCBitArrayCalculator instance for core computation
        _bwfeeder: BigWigReader for mappability data (if available)
        _progress: Progress bar for user feedback
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
    """Multiprocess BitArray worker.

    This class implements a worker process for parallel cross-correlation
    calculation using the BitArray algorithm. Each worker processes one
    chromosome independently and reports results back to the main process.

    The worker uses inter-process communication queues for coordination
    and result reporting, enabling efficient parallel processing of
    multiple chromosomes simultaneously.

    Attributes:
        calculator: CCBitArrayCalculator instance for this worker
        _bwfeeder: BigWigReader for mappability data (if available)
    """
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


# Import compatibility wrapper for backward compatibility
from PyMaSC.core.worker_compat import BACalcWorker as _BACalcWorker

# Replace existing class with compatibility wrapper
BACalcWorker = _BACalcWorker
