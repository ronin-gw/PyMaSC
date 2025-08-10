"""Type stubs for bacore/mscc.pyx - BitArray-based MSCC calculation module."""

from typing import Optional, List, Sequence, Any
from multiprocessing.synchronize import Lock

from PyMaSC.interfaces.calculator import BothCalculatorModel
from PyMaSC.result import BothChromResult, GenomeWideResult
from PyMaSC.reader.bigwig import BigWigReader


class CCBitArrayCalculator(BothCalculatorModel):
    """MSCC calculator using bit array operations.

    Computes both NCC and MSCC statistics simultaneously.
    """

    MAPPABILITY_THRESHOLD: float
    EXTRA_ALLOCATE_SIZE: int

    max_shift: int
    read_len: int
    skip_ncc: bool
    logger_lock: Optional[Lock]

    def __init__(
        self,
        max_shift: int,
        read_len: int,
        references: List[str],
        lengths: Sequence[int],
        bwfeeder: Optional[BigWigReader] = None,
        skip_ncc: bool = False,
        logger_lock: Optional[Lock] = None,
        progress_bar: Optional[Any] = None
    ) -> None:
        """Initialize MSCC calculator.

        Args:
            max_shift: Maximum shift distance for cross-correlation
            read_len: Expected read length for calculations
            references: List of chromosome names
            lengths: List of chromosome lengths
            bwfeeder: BigWig feeder for mappability data (optional)
            skip_ncc: Whether to skip naive cross-correlation calculation
            logger_lock: Lock for thread-safe logging (optional)
            progress_bar: Progress reporting object (optional)
        """
        ...

    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a forward strand read at specified position.

        Marks the read's 5' end position in the forward strand bit array
        for cross-correlation calculation. Handles chromosome switching
        and position validation automatically.

        Args:
            chrom: Chromosome name
            pos: 1-based genomic position of read's 5' end
            readlen: Length of the read in base pairs

        Note:
            Automatically triggers correlation calculation when switching
            between chromosomes. Ignores duplicate reads at same position.
        """
        ...

    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Process a reverse strand read at specified position.

        Marks the read's 3' end position (corresponding to 5' end of
        the reverse complement) in the reverse strand bit array for
        cross-correlation calculation.

        Args:
            chrom: Chromosome name
            pos: 1-based genomic position of read's 5' end
            readlen: Length of the read in base pairs

        Note:
            For reverse reads, marks position at pos + readlen - 1 to
            represent the 5' end of the reverse complement strand.
            Prevents duplicate counting at the same effective position.
        """
        ...

    def flush(self, chrom: Optional[str] = None) -> None:
        """Flush any intermediate calculation results.

        Processes buffered data and calculates cross-correlation for
        the current chromosome. Safe to call multiple times.
        """
        ...

    def finishup_calculation(self) -> None:
        """Complete cross-correlation calculation for all chromosomes.

        Finalizes the calculation by processing any remaining buffered
        data and calculating mappability statistics for chromosomes
        with no aligned reads. This ensures complete genome-wide
        statistics are available.

        Note:
            Must be called after all reads have been processed to
            ensure complete and accurate cross-correlation results.
        """
        ...

    def get_result(self, chrom: str) -> BothChromResult:
        ...

    def get_whole_result(self) -> GenomeWideResult:
        """Get the genome-wide NCC and MSCC calculation results.

        Returns:
            BothGenomeWideResult: Complete results for all chromosomes including
                per-chromosome NCCResult and MSCCResult objects and genome-wide statistics.
        """
        ...
