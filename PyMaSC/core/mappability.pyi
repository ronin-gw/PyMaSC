"""Type stubs for mappability module.

This module provides Cython-optimized implementations for calculating
mappability statistics from BigWig files. Mappability data is essential
for MSCC calculations to correct for genomic regions with variable
mappability that can confound cross-correlation analysis.

Key components:
- MappableLengthCalculator: Core mappability calculation engine
- BWFeederWithMappableRegionSum: Specialized feeder for on-demand calculation
- ContinueCalculation: Exception for flow control in generators

The module efficiently processes BigWig mappability data and calculates
shift-dependent mappability statistics required for MSCC analysis.
"""
from typing import Dict, List, Optional, Tuple, Generator
from multiprocessing.synchronize import Lock

from PyMaSC.reader.bigwig import BigWigReader
from PyMaSC.utils.progress import ProgressHook


class ContinueCalculation(Exception):
    """Exception for flow control in mappability calculation generators.

    Used to signal that a generator should stop yielding data points
    but continue internal calculation processes. This allows for
    flexible control flow in mappability data processing.
    """
    pass


class MappableLengthCalculator(BigWigReader):
    """Core mappability calculation engine for BigWig data.

    Calculates mappability statistics from BigWig files for use in
    MSCC cross-correlation analysis. Processes mappability data
    efficiently using Cython optimization and maintains statistics
    for different shift distances.

    The calculator processes BigWig intervals and computes cumulative
    mappability statistics required for correcting cross-correlation
    analysis in regions with variable mappability.

    Attributes:
        MAPPABILITY_THRESHOLD: Minimum mappability value to consider
        max_shift: Maximum shift distance for calculations
        chrom2is_called: Per-chromosome calculation status
        chrom2mappable_len: Per-chromosome mappability arrays
        mappable_len: Genome-wide mappability statistics
        is_called: Whether all chromosomes have been processed
    """
    MAPPABILITY_THRESHOLD: float
    max_shift: int
    chrom2is_called: Dict[str, bool]
    chrom2mappable_len: Dict[str, Tuple[int, ...]]
    mappable_len: List[int]
    is_called: bool
    _progress: ProgressHook

    def __init__(self, path: str, max_shift: int = 0, logger_lock: Optional[Lock] = None) -> None:
        """Initialize mappability calculator.

        Args:
            path: Path to BigWig mappability file
            max_shift: Maximum shift distance for calculations
            logger_lock: Optional lock for thread-safe logging
        """
        ...

    def enable_progress_bar(self) -> None: ...

    def disable_progress_bar(self) -> None: ...

    def flush(self) -> None: ...

    def calc_mappability(self, chrom: Optional[str] = None) -> None:
        """Calculate mappability statistics for specified chromosome(s).

        Processes BigWig mappability data and computes shift-dependent
        mappability statistics. Can process a single chromosome or
        all unprocessed chromosomes.

        Args:
            chrom: Chromosome name to process (None for all unprocessed)

        Note:
            Updates chrom2mappable_len and mappable_len attributes
            with calculated statistics
        """
        ...

    def get_mappable_len(
        self,
        chrom: Optional[str] = None,
        shift_from: Optional[int] = None,
        shift_to: Optional[int] = None,
        force: bool = False
    ) -> Optional[List[int]]:
        """Retrieve mappability statistics for specified parameters.

        Returns mappability statistics for a chromosome or genome-wide,
        optionally filtered by shift distance range. Can trigger
        calculation if not already computed.

        Args:
            chrom: Chromosome name (None for genome-wide)
            shift_from: Start of shift range (None for beginning)
            shift_to: End of shift range (None for end)
            force: Whether to force calculation if not computed

        Returns:
            List of mappability values for specified range

        Raises:
            KeyError: If chromosome data not calculated and force=False
        """
        ...


class BWFeederWithMappableRegionSum(MappableLengthCalculator):
    """Specialized BigWig feeder with on-demand mappability calculation.

    Extends MappableLengthCalculator to provide streaming access to
    BigWig data while calculating mappability statistics on-demand.
    This class is used in MSCC calculations where mappability data
    is needed incrementally during cross-correlation computation.

    The feeder can operate in two modes:
    1. Using pre-calculated mappability statistics (if available)
    2. Calculating mappability on-demand while streaming data

    Attributes:
        Inherits all attributes from MappableLengthCalculator
    """

    def __init__(
        self,
        path: str,
        max_shift: int,
        chrom2mappable_len: Optional[Dict[str, Tuple[int, ...]]] = None
    ) -> None:
        """Initialize BigWig feeder with mappability calculation.

        Args:
            path: Path to BigWig mappability file
            max_shift: Maximum shift distance for calculations
            chrom2mappable_len: Pre-calculated mappability statistics (optional)
        """
        ...

    def feed(self, chrom: str) -> Generator[Tuple[int, int, float], None, None]:
        """Feed BigWig data for specified chromosome.

        Returns a generator that yields BigWig intervals for the
        specified chromosome. If mappability statistics haven't
        been calculated for this chromosome, performs on-demand
        calculation while yielding data.

        Args:
            chrom: Chromosome name to process

        Returns:
            Generator yielding (start, end, value) tuples for BigWig intervals

        Note:
            Automatically chooses between pre-calculated and on-demand
            calculation modes based on available data
        """
        ...
