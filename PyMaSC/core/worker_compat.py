"""Compatibility wrappers for legacy worker classes.

This module provides backward compatibility for existing code that uses
the old worker classes. The wrappers delegate to the new UnifiedWorker
while maintaining the same interface as the original workers.

This ensures that existing code continues to work without modification
while new code can use the unified architecture.
"""
import warnings
from typing import Any, Dict, Optional
from multiprocessing import Queue, Lock

from PyMaSC.core.worker import UnifiedWorker
from PyMaSC.core.models import (
    CalculationConfig, MappabilityConfig, WorkerConfig, AlgorithmType
)


class LegacyWorkerBase:
    """Base class for legacy worker compatibility wrappers.

    Provides common functionality for creating UnifiedWorker instances
    with the same interface as the original worker classes.
    """

    def __init__(self, order_queue: Queue, report_queue: Queue, logger_lock: Lock,
                 bam_path: str, mapq_criteria: int, max_shift: int,
                 references: list, lengths: list, **kwargs):
        """Initialize legacy worker wrapper.

        Args:
            order_queue: Queue for receiving work orders
            report_queue: Queue for sending results
            logger_lock: Lock for thread-safe logging
            bam_path: Path to BAM file
            mapq_criteria: Minimum mapping quality threshold
            max_shift: Maximum shift distance
            references: List of reference chromosomes
            lengths: List of chromosome lengths
            **kwargs: Additional arguments specific to worker type
        """
        self._order_queue = order_queue
        self._report_queue = report_queue
        self._logger_lock = logger_lock
        self._bam_path = bam_path
        self._mapq_criteria = mapq_criteria
        self._max_shift = max_shift
        self._references = references
        self._lengths = lengths
        self._kwargs = kwargs

        # Create unified worker
        self._unified_worker = self._create_unified_worker()

    def _create_unified_worker(self) -> UnifiedWorker:
        """Create UnifiedWorker instance with appropriate configuration."""
        raise NotImplementedError("Subclasses must implement _create_unified_worker")

    def __getattr__(self, name: str) -> Any:
        """Delegate attribute access to unified worker."""
        return getattr(self._unified_worker, name)

    def run(self) -> None:
        """Run the worker process."""
        self._unified_worker.run()

    @property
    def name(self) -> str:
        """Get worker name."""
        return self._unified_worker.name


class NaiveCCCalcWorker(LegacyWorkerBase):
    """Compatibility wrapper for NaiveCCCalcWorker.

    Provides the same interface as the original NaiveCCCalcWorker
    while delegating to UnifiedWorker internally.
    """

    def __init__(self, order_queue: Queue, report_queue: Queue, logger_lock: Lock,
                 bam_path: str, mapq_criteria: int, max_shift: int,
                 references: list, lengths: list):
        """Initialize NaiveCCCalcWorker compatibility wrapper."""
        warnings.warn(
            "NaiveCCCalcWorker is deprecated. Use UnifiedWorker instead.",
            DeprecationWarning,
            stacklevel=2
        )
        super().__init__(order_queue, report_queue, logger_lock, bam_path,
                        mapq_criteria, max_shift, references, lengths)

    def _create_unified_worker(self) -> UnifiedWorker:
        """Create UnifiedWorker configured for NCC calculation."""
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.NAIVE_CC,
            max_shift=self._max_shift,
            mapq_criteria=self._mapq_criteria,
            references=self._references,
            lengths=self._lengths,
            skip_ncc=False
        )

        worker_config = WorkerConfig(
            calculation_config=calc_config,
            mappability_config=None,
            progress_enabled=True,
            logger_lock=self._logger_lock
        )

        return UnifiedWorker(
            worker_config,
            self._order_queue,
            self._report_queue,
            self._logger_lock,
            self._bam_path
        )


class MSCCCalcWorker(LegacyWorkerBase):
    """Compatibility wrapper for MSCCCalcWorker.

    Provides the same interface as the original MSCCCalcWorker
    while delegating to UnifiedWorker internally.
    """

    def __init__(self, order_queue: Queue, report_queue: Queue, logger_lock: Lock,
                 bam_path: str, mapq_criteria: int, max_shift: int,
                 references: list, lengths: list, mappable_path: str,
                 read_len: int, chrom2mappable_len: Dict[str, int]):
        """Initialize MSCCCalcWorker compatibility wrapper."""
        warnings.warn(
            "MSCCCalcWorker is deprecated. Use UnifiedWorker instead.",
            DeprecationWarning,
            stacklevel=2
        )
        super().__init__(order_queue, report_queue, logger_lock, bam_path,
                        mapq_criteria, max_shift, references, lengths,
                        mappable_path=mappable_path, read_len=read_len,
                        chrom2mappable_len=chrom2mappable_len)

    def _create_unified_worker(self) -> UnifiedWorker:
        """Create UnifiedWorker configured for MSCC calculation."""
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.MSCC,
            max_shift=self._max_shift,
            mapq_criteria=self._mapq_criteria,
            references=self._references,
            lengths=self._lengths,
            read_length=self._kwargs['read_len'],
            skip_ncc=True  # MSCC only
        )

        map_config = MappabilityConfig(
            mappability_path=self._kwargs['mappable_path'],
            read_len=self._kwargs['read_len']
        )

        worker_config = WorkerConfig(
            calculation_config=calc_config,
            mappability_config=map_config,
            progress_enabled=True,
            logger_lock=self._logger_lock
        )

        return UnifiedWorker(
            worker_config,
            self._order_queue,
            self._report_queue,
            self._logger_lock,
            self._bam_path
        )


class NCCandMSCCCalcWorker(LegacyWorkerBase):
    """Compatibility wrapper for NCCandMSCCCalcWorker.

    Provides the same interface as the original NCCandMSCCCalcWorker
    while delegating to UnifiedWorker internally.
    """

    def __init__(self, order_queue: Queue, report_queue: Queue, logger_lock: Lock,
                 bam_path: str, mapq_criteria: int, max_shift: int,
                 references: list, lengths: list, mappable_path: str,
                 read_len: int, chrom2mappable_len: Dict[str, int]):
        """Initialize NCCandMSCCCalcWorker compatibility wrapper."""
        warnings.warn(
            "NCCandMSCCCalcWorker is deprecated. Use UnifiedWorker instead.",
            DeprecationWarning,
            stacklevel=2
        )
        super().__init__(order_queue, report_queue, logger_lock, bam_path,
                        mapq_criteria, max_shift, references, lengths,
                        mappable_path=mappable_path, read_len=read_len,
                        chrom2mappable_len=chrom2mappable_len)

    def _create_unified_worker(self) -> UnifiedWorker:
        """Create UnifiedWorker configured for both NCC and MSCC calculation."""
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.MSCC,
            max_shift=self._max_shift,
            mapq_criteria=self._mapq_criteria,
            references=self._references,
            lengths=self._lengths,
            read_length=self._kwargs['read_len'],
            skip_ncc=False  # Include both NCC and MSCC
        )

        map_config = MappabilityConfig(
            mappability_path=self._kwargs['mappable_path'],
            read_len=self._kwargs['read_len']
        )

        worker_config = WorkerConfig(
            calculation_config=calc_config,
            mappability_config=map_config,
            progress_enabled=True,
            logger_lock=self._logger_lock
        )

        return UnifiedWorker(
            worker_config,
            self._order_queue,
            self._report_queue,
            self._logger_lock,
            self._bam_path
        )


class BACalcWorker(LegacyWorkerBase):
    """Compatibility wrapper for BACalcWorker.

    Provides the same interface as the original BACalcWorker
    while delegating to UnifiedWorker internally.
    """

    def __init__(self, order_queue: Queue, report_queue: Queue, logger_lock: Lock,
                 bam_path: str, mapq_criteria: int, max_shift: int, read_len: int,
                 references: list, lengths: list, mappability_handler: Optional[Any] = None,
                 skip_ncc: bool = False):
        """Initialize BACalcWorker compatibility wrapper."""
        warnings.warn(
            "BACalcWorker is deprecated. Use UnifiedWorker instead.",
            DeprecationWarning,
            stacklevel=2
        )
        super().__init__(order_queue, report_queue, logger_lock, bam_path,
                        mapq_criteria, max_shift, references, lengths,
                        read_len=read_len, mappability_handler=mappability_handler,
                        skip_ncc=skip_ncc)

    def _create_unified_worker(self) -> UnifiedWorker:
        """Create UnifiedWorker configured for BitArray calculation."""
        calc_config = CalculationConfig(
            algorithm=AlgorithmType.BITARRAY,
            max_shift=self._max_shift,
            mapq_criteria=self._mapq_criteria,
            references=self._references,
            lengths=self._lengths,
            read_length=self._kwargs['read_len'],
            skip_ncc=self._kwargs['skip_ncc']
        )

        map_config = None
        mappability_handler = self._kwargs.get('mappability_handler')
        if mappability_handler:
            map_config = MappabilityConfig(
                mappability_path=mappability_handler.path,
                read_len=self._kwargs['read_len']
            )

        worker_config = WorkerConfig(
            calculation_config=calc_config,
            mappability_config=map_config,
            progress_enabled=True,
            logger_lock=self._logger_lock
        )

        return UnifiedWorker(
            worker_config,
            self._order_queue,
            self._report_queue,
            self._logger_lock,
            self._bam_path
        )


# Convenience function for backward compatibility
def create_legacy_worker(worker_type: str, *args, **kwargs) -> Any:
    """Create a legacy worker using the old interface.

    This function provides backward compatibility for code that
    creates workers using string identifiers.

    Args:
        worker_type: Type of worker ('ncc', 'mscc', 'ncc_mscc', 'bitarray')
        *args, **kwargs: Arguments for worker creation

    Returns:
        Compatible worker instance
    """
    if worker_type == 'ncc':
        return NaiveCCCalcWorker(*args, **kwargs)
    elif worker_type == 'mscc':
        return MSCCCalcWorker(*args, **kwargs)
    elif worker_type == 'ncc_mscc':
        return NCCandMSCCCalcWorker(*args, **kwargs)
    elif worker_type == 'bitarray':
        return BACalcWorker(*args, **kwargs)
    else:
        raise ValueError(f"Unknown worker type: {worker_type}")