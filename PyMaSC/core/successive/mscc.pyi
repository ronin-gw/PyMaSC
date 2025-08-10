"""Type stubs for mscc.pyx - Mappability-Sensitive Cross-Correlation calculation module."""

from typing import Optional, List, Sequence
from multiprocessing.synchronize import Lock

from PyMaSC.interfaces.calculator import MSCCCalculatorModel
from PyMaSC.result import MSCCResult, MSCCGenomeWideResult

from PyMaSC.reader.bigwig import BigWigReader


class MSCCCalculator(MSCCCalculatorModel):
    """Mappability-Sensitive Cross-Correlation calculator."""

    max_shift: int
    read_len: int
    logger_lock: Optional[Lock]

    def __init__(
        self,
        max_shift: int,
        read_len: int,
        references: List[str],
        lengths: Sequence[int],
        bwfeeder: BigWigReader,
        logger_lock: Optional[Lock] = None
    ) -> None:
        """Initialize MSCC calculator with mappability support.

        Args:
            max_shift: Maximum shift distance for correlation calculation
            read_len: Read length for mappability window calculation
            references: List of chromosome names
            lengths: List of chromosome lengths
            bwfeeder: BigWig reader for mappability data
            logger_lock: Optional lock for thread-safe logging
        """
        ...

    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a forward strand read with mappability weighting.

        Args:
            chrom: Chromosome name
            pos: Read start position (1-based)
            readlen: Read length in base pairs
        """
        ...

    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a reverse strand read with mappability weighting.

        Args:
            chrom: Chromosome name
            pos: Read start position (1-based)
            readlen: Read length in base pairs
        """
        ...

    def flush(self, chrom: Optional[str] = None) -> None:
        """Finalize MSCC calculation for current chromosome."""
        ...

    def finishup_calculation(self) -> None:
        """Complete MSCC calculation for all chromosomes."""
        ...

    def get_result(self, chrom: str) -> MSCCResult:
        ...

    def get_whole_result(self) -> MSCCGenomeWideResult:
        """Get the genome-wide MSCC calculation results.

        Returns:
            MSCCGenomeWideResult: Complete results for all chromosomes including
                per-chromosome MSCCResult objects and genome-wide statistics.
        """
        ...