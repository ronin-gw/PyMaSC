"""Type stubs for ncc.pyx - Naive Cross-Correlation calculation module."""

from typing import Optional, List, Any
from .interfaces.calculator import CrossCorrelationCalculator
from .interfaces.result import NCCGenomeWideResult


class ReadUnsortedError(IndexError):
    """Exception raised when input reads are not sorted by position."""
    pass


class NaiveCCCalculator(CrossCorrelationCalculator):
    """Cython implementation of naive cross-correlation calculation."""

    # Public readonly attributes
    max_shift: int
    logger_lock: Optional[Any]

    def __init__(
        self,
        max_shift: int,
        references: List[str],
        lengths: List[int],
        logger_lock: Optional[Any] = None
    ) -> None:
        """Initialize NCC calculator.

        Args:
            max_shift: Maximum shift distance for correlation calculation
            references: List of chromosome names
            lengths: List of chromosome lengths
            logger_lock: Optional lock for thread-safe logging
        """
        ...

    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a forward strand read for cross-correlation calculation.

        Args:
            chrom: Chromosome name
            pos: Read start position (1-based)
            readlen: Read length in base pairs
        """
        ...

    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a reverse strand read for cross-correlation calculation.

        Args:
            chrom: Chromosome name
            pos: Read start position (1-based)
            readlen: Read length in base pairs
        """
        ...

    def flush(self) -> None:
        """Finalize cross-correlation calculation for current chromosome."""
        ...

    def finishup_calculation(self) -> None:
        """Complete cross-correlation calculation for all chromosomes."""
        ...

    def get_result(self) -> NCCGenomeWideResult:
        """Get the genome-wide NCC calculation results.

        Returns:
            NCCGenomeWideResult: Complete results for all chromosomes including
                per-chromosome NCCResult objects and genome-wide statistics.
        """
        ...