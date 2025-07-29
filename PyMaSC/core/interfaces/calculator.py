from abc import ABC, abstractmethod
from typing import Optional
from multiprocessing.synchronize import Lock

# Calculators use result models that support cc calculation from ccbins.
from ..result import (
    ChromResult, NCCResult, MSCCResult, BothChromResult,
    GenomeWideResult, NCCGenomeWideResult, MSCCGenomeWideResult, BothGenomeWideResult
)


class CrossCorrelationCalculator(ABC):
    """Abstract base class for all cross-correlation calculators.

    This interface defines the common contract that all correlation calculation
    implementations must follow, enabling polymorphic usage across different
    algorithms (NCC, MSCC, BitArray, etc.).

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
        """Flush any intermediate calculation results.

        Used in multiprocessing scenarios to ensure results are available
        for collection before worker process termination.

        If `chrom` is provided, flush results for that specific chromosome;
        especially in multiprocessing, this argument is used to ensure that
        result is available for the chromosome that has no reads processed.
        """
        pass

    @abstractmethod
    def get_result(self, chrom: str) -> ChromResult:
        pass

    @abstractmethod
    def get_whole_result(self) -> GenomeWideResult:
        pass


class NCCCalculatorModel(CrossCorrelationCalculator):
    @abstractmethod
    def get_result(self, chrom: str) -> NCCResult:
        pass

    @abstractmethod
    def get_whole_result(self) -> NCCGenomeWideResult:
        pass


class MSCCCalculatorModel(CrossCorrelationCalculator):
    @abstractmethod
    def get_result(self, chrom: str) -> MSCCResult:
        pass

    @abstractmethod
    def get_whole_result(self) -> MSCCGenomeWideResult:
        pass


class BothCalculatorModel(CrossCorrelationCalculator):
    read_len: int
    skip_ncc: bool

    @abstractmethod
    def get_result(self, chrom: str) -> BothChromResult:
        pass

    @abstractmethod
    def get_whole_result(self) -> BothGenomeWideResult:
        pass
