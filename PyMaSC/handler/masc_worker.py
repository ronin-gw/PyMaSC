"""Multiprocessing worker classes for parallel cross-correlation calculation.

This module provides worker process implementations for parallel processing
of cross-correlation calculations across multiple chromosomes. Each worker
handles one chromosome independently, enabling efficient parallel computation
for large genomic datasets.

Key components:
- CalcWorkerBase: Base class for all calculation workers
- NaiveCCCalcWorker: Worker for naive cross-correlation calculation
- MSCCCalcWorker: Worker for mappability-sensitive cross-correlation
- NCCandMSCCCalcWorker: Worker for both NCC and MSCC calculations

The workers use inter-process communication to coordinate with the main
process and report progress and results efficiently.

NOTE: This module is deprecated. Use PyMaSC.core.worker.UnifiedWorker instead.
Legacy workers are available for backward compatibility through compatibility wrappers.
"""
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
    """Base class for multiprocessing calculation workers.
    
    Provides common functionality for all PyMaSC calculation workers,
    including inter-process communication, BAM file handling, and
    basic workflow coordination.
    
    This class extends multiprocessing.Process to enable parallel
    processing of different chromosomes simultaneously.
    
    Attributes:
        order_queue: Queue for receiving work orders from main process
        report_queue: Queue for sending results back to main process
        logger_lock: Lock for thread-safe logging operations
        path: Path to BAM file being processed
        mapq_criteria: Minimum mapping quality threshold
        calculator: Calculation engine instance
    """
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
                    self._progress.update(read.reference_start)

                    try:
                        self._feed_read(read)
                    except ReadUnsortedError:
                        logger.error("Input alignment file must be sorted.")
                        raise

                self._put_result_to_report_queue(chrom)

        with self.logger_lock:
            logger.debug("{}: Goodbye.".format(self.name))
        self._deconstruct()

    def _feed_read(self, read):
        if (read.is_read2 or read.mapping_quality < self.mapq_criteria or
                read.is_unmapped or read.is_duplicate):
            return

        if read.is_reverse:
            self.calculator.feed_reverse_read(
                read.reference_name,
                read.reference_start + 1,
                read.infer_query_length()
            )
        else:
            self.calculator.feed_forward_read(
                read.reference_name,
                read.reference_start + 1,
                read.infer_query_length()
            )

    def _put_result_to_report_queue(self, chrom):
        pass

    def _deconstruct(self):
        pass


class NaiveCCCalcWorker(CalcWorkerBase):
    """Worker process for naive cross-correlation calculation.
    
    Specializes CalcWorkerBase to perform naive cross-correlation (NCC)
    calculations for assigned chromosomes. Uses NaiveCCCalculator to
    compute cross-correlation without mappability correction.
    
    This worker is used when only basic cross-correlation analysis
    is required, without the complexity of mappability adjustment.
    """
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
    """Worker process for mappability-sensitive cross-correlation.
    
    Extends NaiveCCCalcWorker to perform MSCC calculations that
    incorporate mappability correction from BigWig files. This worker
    handles the additional complexity of mappability data integration.
    
    Used when mappability correction is required to address genomic
    regions with variable mappability that can confound standard
    cross-correlation analysis.
    """
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

    def _deconstruct(self):
        self.bwfeeder.close()


class NCCandMSCCCalcWorker(MSCCCalcWorker):
    """Worker process for both NCC and MSCC calculations.
    
    Performs both naive cross-correlation and mappability-sensitive
    cross-correlation calculations in a single worker process.
    This provides comprehensive analysis including both standard
    and mappability-corrected results.
    
    Used when complete cross-correlation analysis is required,
    enabling comparison between standard and mappability-corrected
    approaches within the same workflow.
    """
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
        if (read.is_read2 or read.mapping_quality < self.mapq_criteria or
                read.is_unmapped or read.is_duplicate):
            return

        elif read.is_reverse:
            self.ncc_calculator.feed_reverse_read(
                read.reference_name,
                read.reference_start + 1,
                read.infer_query_length()
            )
            self.mscc_calculator.feed_reverse_read(
                read.reference_name,
                read.reference_start + 1,
                read.infer_query_length()
            )
        else:
            self.ncc_calculator.feed_forward_read(
                read.reference_name,
                read.reference_start + 1,
                read.infer_query_length()
            )
            self.mscc_calculator.feed_forward_read(
                read.reference_name,
                read.reference_start + 1,
                read.infer_query_length()
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


# Import compatibility wrappers for backward compatibility
# These replace the original classes with unified worker implementations
from PyMaSC.core.worker_compat import (
    NaiveCCCalcWorker as _NaiveCCCalcWorker,
    MSCCCalcWorker as _MSCCCalcWorker,
    NCCandMSCCCalcWorker as _NCCandMSCCCalcWorker
)

# Replace existing classes with compatibility wrappers
# This ensures existing code continues to work without modification
NaiveCCCalcWorker = _NaiveCCCalcWorker
MSCCCalcWorker = _MSCCCalcWorker
NCCandMSCCCalcWorker = _NCCandMSCCCalcWorker
