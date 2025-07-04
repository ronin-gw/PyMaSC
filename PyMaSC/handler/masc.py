"""Main calculation handler implementing cross-correlation algorithms.

This module contains the core calculation handler that orchestrates cross-correlation
analysis with optional mappability correction. It supports both single-process and
multiprocess execution modes and manages the complete workflow from BAM file reading
to cross-correlation computation.

The handler coordinates:
- BAM file reading and validation
- Chromosome filtering and reference size validation
- Read length estimation or validation
- Mappability correction setup
- Cross-correlation calculation execution
- Result aggregation and storage
"""
import logging
from multiprocessing import Queue, Lock

import pysam

from PyMaSC.handler.masc_noindex_worker import SingleProcessCalculator
from PyMaSC.handler.masc_worker import NaiveCCCalcWorker, MSCCCalcWorker, NCCandMSCCCalcWorker
from PyMaSC.core.readlen import estimate_readlen
from PyMaSC.utils.progress import MultiLineProgressManager, ProgressBar, ProgressHook
from PyMaSC.utils.calc import exec_worker_pool, filter_chroms

logger = logging.getLogger(__name__)


class InputUnseekable(Exception):
    """Exception raised when input stream cannot be rewound for read length estimation.
    
    This exception is raised when trying to estimate read length from unseekable
    input streams (like stdin) that cannot be rewound for multiple passes.
    """
    pass


class NothingToCalc(Exception):
    """Exception raised when no chromosomes match the filtering criteria.
    
    This exception is raised when chromosome filtering results in no target
    chromosomes for analysis, typically due to overly restrictive include/exclude
    chromosome patterns.
    """
    pass


class CCCalcHandler(object):
    """Main cross-correlation calculation handler.
    
    This class orchestrates the complete cross-correlation analysis workflow,
    including file validation, chromosome filtering, read length estimation,
    mappability correction setup, and calculation execution.
    
    The handler supports both single-process and multiprocess execution modes,
    automatically falling back to single-process for unindexed BAM files.
    
    Attributes:
        path: Input BAM file path
        esttype: Read length estimation method
        max_shift: Maximum shift distance for correlation calculation
        mapq_criteria: Minimum mapping quality threshold
        nworker: Number of worker processes
        skip_ncc: Whether to skip naive cross-correlation calculation
        references: List of target chromosome names
        lengths: List of chromosome lengths
        read_len: Estimated or specified read length
        mappability_handler: Optional mappability correction handler
    """
    def __init__(self, path, esttype, max_shift, mapq_criteria, nworker=1, skip_ncc=False, chromfilter=None):
        """Initialize calculation handler with file and parameter validation.
        
        Args:
            path: Path to input BAM file
            esttype: Read length estimation method ('mean' or 'median')
            max_shift: Maximum shift distance for correlation calculation
            mapq_criteria: Minimum mapping quality threshold
            nworker: Number of worker processes (default: 1)
            skip_ncc: Whether to skip naive cross-correlation (default: False)
            chromfilter: Chromosome filtering patterns (default: None)
            
        Raises:
            ValueError: If BAM file has no sequences defined
            NothingToCalc: If no chromosomes match filtering criteria
        """
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

        #
        target_references = filter_chroms(self.align_file.references, chromfilter)
        if not target_references:
            logger.error("There is no targeted chromosomes.")
            raise NothingToCalc
        self.references = []
        self.lengths = []
        need_warning = False
        for reference, length in zip(self.align_file.references, self.align_file.lengths):
            if reference in target_references:
                self.references.append(reference)
                self.lengths.append(length)

        #
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
        """Set or estimate read length for cross-correlation calculation.
        
        Either uses the provided read length or estimates it from the BAM file
        using the specified estimation method.
        
        Args:
            readlen: Explicit read length to use, or None to estimate
            
        Raises:
            InputUnseekable: If read length estimation is needed but input is unseekable
            ValueError: If read length exceeds maximum shift distance
        """
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
        """Configure mappability correction handler.
        
        Sets up the mappability handler for MSCC calculation and validates
        chromosome size consistency between BAM and BigWig files.
        
        Args:
            mappability_handler: MappabilityHandler instance for correction
            
        Note:
            Automatically adjusts chromosome lengths if BigWig has longer sequences
        """
        self.mappability_handler = mappability_handler
        bw_chromsizes = self.mappability_handler.chromsizes

        for i, reference in enumerate(self.references):
            if reference not in bw_chromsizes:
                logger.debug("mappability for '{}' not found".format(reference))
                continue
            self._compare_refsize(reference, bw_chromsizes[reference])

    def _compare_refsize(self, reference, bw_chr_size):
        """Compare and validate chromosome sizes between BAM and BigWig files.
        
        Checks for size mismatches and adjusts BAM chromosome lengths if the
        BigWig file has longer sequences for the same chromosome.
        
        Args:
            reference: Chromosome name
            bw_chr_size: Chromosome size from BigWig file
            
        Note:
            Logs warnings for size mismatches and automatically uses longer length
        """
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
        """Execute cross-correlation calculation workflow.
        
        Chooses between single-process or multiprocess execution based on
        the number of workers and BAM file indexing status.
        
        Note:
            Automatically disables single-line progress bars for multiprocess mode
        """
        if self.nworker > 1:
            # to avoid interfering mappability progress bar with multiline progress bar
            ProgressBar.global_switch = False
            ProgressHook.global_switch = True
            self._run_multiprocess_calcuration()
        else:
            self._run_singleprocess_calculation()

    def _run_singleprocess_calculation(self):
        """Execute calculation using single-process mode.
        
        Creates a single process calculator worker and runs the complete
        cross-correlation analysis for all chromosomes sequentially.
        
        Note:
            Used when nworker=1 or when BAM file is not indexed
        """
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
        """Execute calculation using multiprocess mode.
        
        Sets up inter-process communication queues and launches worker processes
        for parallel chromosome processing. Manages worker coordination and
        result aggregation.
        
        Note:
            Requires indexed BAM files and automatically falls back to single-process
            if the file is not indexed
        """
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
        """Coordinate worker processes and collect results.
        
        Manages the worker process pool, handles progress reporting,
        and collects calculation results from all workers.
        
        Args:
            workers: List of worker process instances
            
        Note:
            Uses multiline progress display for concurrent chromosome processing
        """
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
        """Process and store results from worker processes.
        
        Unpacks worker results and stores them in the appropriate handler
        attributes for both naive cross-correlation and MSCC results.
        
        Args:
            chrom: Chromosome name
            mappable_len: Calculated mappable length for the chromosome
            cc_stats: Tuple of (forward_sum, reverse_sum, ccbins) for naive CC
            masc_stats: Tuple of (forward_sum, reverse_sum, ccbins) for MSCC
        """
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
        """Complete any remaining mappability calculations.
        
        Finalizes mappability calculations if they haven't been completed
        during the main calculation workflow and updates the mappable
        length dictionary.
        
        Note:
            Called at the end of both single-process and multiprocess workflows
        """
        if self.mappability_handler is not None:
            if not self.mappability_handler.is_called:
                self.mappability_handler.is_called = all(
                    self.mappability_handler.chrom2is_called.values()
                )
                self.mappability_handler.calc_mappability()
            self.ref2mappable_len = self.mappability_handler.chrom2mappable_len
