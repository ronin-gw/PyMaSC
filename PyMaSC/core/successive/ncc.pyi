"""Type stubs for ncc.pyx - Naive Cross-Correlation calculation module."""

from typing import Optional, Sequence, Any
from multiprocessing.synchronize import Lock

from PyMaSC.interfaces.calculator import NCCCalculatorModel
from PyMaSC.result import NCCResult, NCCGenomeWideResult


class ReadUnsortedError(IndexError):
    """Exception raised when input reads are not sorted by position."""
    pass


class NaiveCCCalculator(NCCCalculatorModel):
    """Naive cross-correlation calculator."""

    max_shift: int
    logger_lock: Optional[Lock]

    def __init__(
        self,
        max_shift: int,
        references: Sequence[str],
        lengths: Sequence[int],
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

    def flush(self, chrom: Optional[str] = None) -> None:
        """Finalize cross-correlation calculation for current chromosome."""
        ...

    def finishup_calculation(self) -> None:
        """Complete cross-correlation calculation for all chromosomes."""
        ...

    def get_result(self, chrom: str) -> NCCResult:
        ...

    def get_whole_result(self) -> NCCGenomeWideResult:
        """Get genome-wide NCC results."""
        ...