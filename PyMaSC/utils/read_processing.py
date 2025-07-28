"""Utilities for read processing and filtering.

This module provides common read processing functionality used across
all PyMaSC calculation workers and handlers, eliminating duplicate
read filtering and processing logic.

Key features:
- Standardized read filtering criteria
- Common read validation patterns
- Shared read feeding logic
- Consistent error handling for read processing

This eliminates 60-75 lines of duplicate code across workers.
"""
import logging
from typing import Protocol, runtime_checkable

from pysam import AlignedSegment

logger = logging.getLogger(__name__)


@runtime_checkable
class ReadCalculator(Protocol):
    """Protocol for objects that can receive reads."""

    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Feed a forward strand read to the calculator."""
        ...

    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        """Feed a reverse strand read to the calculator."""
        ...


class ReadFilter:
    """Standardized read filtering based on PyMaSC criteria.

    Provides consistent read filtering logic used across all workers
    and calculation modules. Centralizes the quality and type criteria
    to ensure consistent behavior.

    Attributes:
        mapq_criteria: Minimum mapping quality threshold
        skip_read2: Whether to skip second reads in pairs
        skip_unmapped: Whether to skip unmapped reads
        skip_duplicates: Whether to skip duplicate reads
    """

    def __init__(self,
                 mapq_criteria: int = 0,
                 skip_read2: bool = True,
                 skip_unmapped: bool = True,
                 skip_duplicates: bool = True):
        """Initialize read filter with criteria.

        Args:
            mapq_criteria: Minimum mapping quality threshold
            skip_read2: Whether to skip second reads in pairs
            skip_unmapped: Whether to skip unmapped reads
            skip_duplicates: Whether to skip duplicate reads
        """
        self.mapq_criteria = mapq_criteria
        self.skip_read2 = skip_read2
        self.skip_unmapped = skip_unmapped
        self.skip_duplicates = skip_duplicates

    def should_skip_read(self, read: AlignedSegment) -> bool:
        """Check if a read should be skipped based on quality criteria.

        This is the standardized read filtering logic used throughout PyMaSC.
        It replaces the duplicate filtering code found in 7+ locations.

        Args:
            read: Read to check

        Returns:
            True if read should be skipped, False if it should be processed
        """
        # Skip second reads in pairs if configured
        if self.skip_read2 and read.is_read2:
            return True

        # Skip reads below mapping quality threshold
        if read.mapping_quality < self.mapq_criteria:
            return True

        # Skip unmapped reads if configured
        if self.skip_unmapped and read.is_unmapped:
            return True

        # Skip duplicate reads if configured
        if self.skip_duplicates and read.is_duplicate:
            return True

        return False


class ReadProcessor:
    """
    """

    calculator: ReadCalculator
    read_filter: ReadFilter

    def __init__(self, calculator: ReadCalculator, read_filter: ReadFilter) -> None:
        """
        """
        self.calculator = calculator
        self.read_filter = read_filter

    def process_read(self, read: AlignedSegment) -> bool:
        """Process a single read through filtering and feeding.

        This standardizes the read processing pattern used across all workers.
        Replaces duplicate logic found in multiple worker implementations.

        Args:
            read: Read to process

        Returns:
            True if read was processed, False if skipped

        Raises:
            ReadUnsortedError: If reads appear to be unsorted (implementation specific)
        """
        # Apply quality filtering
        if self.read_filter.should_skip_read(read):
            return False

        # Extract read information
        chrom = read.reference_name
        if chrom is None:
            return False  # Skip reads without reference

        pos = read.reference_start + 1  # Convert to 1-based coordinate
        readlen = read.infer_query_length()
        if readlen is None:
            return False  # Skip reads without queryable length

        # Feed to appropriate calculator method based on strand
        if read.is_reverse:
            self.feed_reverse_read(chrom, pos, readlen)
        else:
            self.feed_forward_read(chrom, pos, readlen)

        return True

    def feed_forward_read(self, chrom: str, pos: int, readlen: int) -> None:
        return self.calculator.feed_forward_read(chrom, pos, readlen)

    def feed_reverse_read(self, chrom: str, pos: int, readlen: int) -> None:
        return self.calculator.feed_reverse_read(chrom, pos, readlen)


def create_read_processor(calculator: ReadCalculator, mapq_criteria: int) -> ReadProcessor:
    """Create standard read processor with PyMaSC filter settings.

    Args:
        calculator: Calculator to feed reads to
        mapq_criteria: Minimum mapping quality threshold

    Returns:
        ReadProcessor with standard configuration
    """
    filter_obj = ReadFilter(mapq_criteria)
    return ReadProcessor(calculator, filter_obj)
