"""Abstract interfaces for cross-correlation calculators.

Defines the calculator protocol hierarchy for NCC, MSCC, and composite calculations.
"""
from abc import ABC, abstractmethod
from typing import Optional
from multiprocessing.synchronize import Lock

from ..result import (
    ChromResult, NCCResult, MSCCResult,
    GenomeWideResult, NCCGenomeWideResult, MSCCGenomeWideResult
)


class CrossCorrelationCalculator(ABC):
    """Base interface for cross-correlation calculators.

    All implementations must provide methods for:
    - Processing forward and reverse strand reads
    - Finalizing calculations
    - Flushing intermediate results
    - Accessing result data
    """

    max_shift: int
    logger_lock: Optional[Lock]

    @abstractmethod
    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a forward strand read for cross-correlation calculation.

        Args:
            chrom: Chromosome name
            pos: Read position (1-based)
            readlen: Read length
        """
        pass

    @abstractmethod
    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a reverse strand read for cross-correlation calculation.

        Args:
            chrom: Chromosome name
            pos: Read position (1-based)
            readlen: Read length
        """
        pass

    @abstractmethod
    def finishup_calculation(self) -> None:
        """Complete cross-correlation calculation for all chromosomes.

        This method should be called after all reads have been processed
        to finalize the correlation calculations.
        """
        pass

    @abstractmethod
    def flush(self, chrom: Optional[str] = None) -> None:
        """Flush intermediate results.

        Used in multiprocessing to ensure results are available
        for collection before worker process termination.

        Args:
            chrom: Specific chromosome to flush (None for all).
                  In multiprocessing, ensures result availability
                  for chromosomes with no reads processed.
        """
        pass

    @abstractmethod
    def get_result(self, chrom: str) -> ChromResult:
        pass

    @abstractmethod
    def get_whole_result(self) -> GenomeWideResult:
        pass


class NCCCalculatorModel(CrossCorrelationCalculator):
    """Calculator interface for Naive Cross-Correlation (NCC) analysis."""

    @abstractmethod
    def get_result(self, chrom: str) -> NCCResult:
        """Get NCC results for a specific chromosome."""
        pass

    @abstractmethod
    def get_whole_result(self) -> NCCGenomeWideResult:
        """Get genome-wide NCC results."""
        pass


class MSCCCalculatorModel(CrossCorrelationCalculator):
    """Calculator interface for Mappability-Sensitive Cross-Correlation (MSCC) analysis."""

    @abstractmethod
    def get_result(self, chrom: str) -> MSCCResult:
        """Get MSCC results for a specific chromosome."""
        pass

    @abstractmethod
    def get_whole_result(self) -> MSCCGenomeWideResult:
        """Get genome-wide MSCC results."""
        pass


class BothCalculatorModel(CrossCorrelationCalculator):
    """Calculator interface for combined NCC and MSCC analysis."""

    read_len: int
    skip_ncc: bool

    @abstractmethod
    def get_result(self, chrom: str) -> ChromResult:
        """Get combined results for a specific chromosome."""
        pass

    @abstractmethod
    def get_whole_result(self) -> GenomeWideResult:
        """Get genome-wide combined results."""
        pass
