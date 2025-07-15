from abc import ABC, abstractmethod

from .result import GenomeWideResult


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
    def flush(self) -> None:
        """Flush any intermediate calculation results.

        Used in multiprocessing scenarios to ensure results are available
        for collection before worker process termination.
        """
        pass

    @abstractmethod
    def get_result(self) -> GenomeWideResult:
        pass