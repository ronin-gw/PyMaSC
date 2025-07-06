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


class ReadUnsortedError(Exception):
    """Exception raised when reads are not sorted as expected."""
    pass


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

    def get_filter_stats(self) -> dict:
        """Get current filter configuration as dictionary.

        Returns:
            Dictionary with current filter settings
        """
        return {
            'mapq_criteria': self.mapq_criteria,
            'skip_read2': self.skip_read2,
            'skip_unmapped': self.skip_unmapped,
            'skip_duplicates': self.skip_duplicates
        }


class ReadProcessor:
    """Standardized read processing for feeding to calculators.

    Handles the common pattern of processing reads and feeding them
    to calculation objects. Centralizes the read coordinate conversion
    and strand-specific feeding logic.

    Attributes:
        calculator: Calculator object to feed reads to
        filter: ReadFilter for quality filtering
        position_offset: Offset to add to read positions (default 1 for 1-based)
    """

    def __init__(self,
                 calculator: ReadCalculator,
                 read_filter: ReadFilter,
                 position_offset: int = 1):
        """Initialize read processor.

        Args:
            calculator: Calculator to feed reads to
            read_filter: Filter for read quality checking
            position_offset: Offset for coordinate conversion (1 for 1-based coordinates)
        """
        self.calculator = calculator
        self.filter = read_filter
        self.position_offset = position_offset

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
        if self.filter.should_skip_read(read):
            return False

        # Extract read information
        chrom = read.reference_name
        if chrom is None:
            return False  # Skip reads without reference

        pos = read.reference_start + self.position_offset
        readlen = read.infer_query_length()
        if readlen is None:
            return False  # Skip reads without queryable length

        # Feed to appropriate calculator method based on strand
        try:
            if read.is_reverse:
                self.calculator.feed_reverse_read(chrom, pos, readlen)
            else:
                self.calculator.feed_forward_read(chrom, pos, readlen)

            return True

        except Exception as e:
            # Convert certain exceptions to ReadUnsortedError for consistency
            if "sorted" in str(e).lower() or "order" in str(e).lower():
                raise ReadUnsortedError("Input alignment file must be sorted") from e
            else:
                # Re-raise other exceptions as-is
                raise


class DualReadProcessor:
    """Read processor that feeds reads to two calculators simultaneously.

    Used for calculations that need both NCC and MSCC results.
    Eliminates duplicate read feeding logic in combined workers.

    Attributes:
        primary_calculator: Primary calculator (typically MSCC)
        secondary_calculator: Secondary calculator (typically NCC)
        filter: ReadFilter for quality checking
        position_offset: Offset for coordinate conversion
    """

    def __init__(self,
                 primary_calculator: ReadCalculator,
                 secondary_calculator: ReadCalculator,
                 read_filter: ReadFilter,
                 position_offset: int = 1):
        """Initialize dual read processor.

        Args:
            primary_calculator: Primary calculator to feed reads to
            secondary_calculator: Secondary calculator to feed reads to
            read_filter: Filter for read quality checking
            position_offset: Offset for coordinate conversion
        """
        self.primary_calculator = primary_calculator
        self.secondary_calculator = secondary_calculator
        self.filter = read_filter
        self.position_offset = position_offset

    def process_read(self, read: AlignedSegment) -> bool:
        """Process read through both calculators.

        Args:
            read: Read to process

        Returns:
            True if read was processed, False if skipped
        """
        # Apply quality filtering
        if self.filter.should_skip_read(read):
            return False

        # Extract read information
        chrom = read.reference_name
        if chrom is None:
            return False  # Skip reads without reference

        pos = read.reference_start + self.position_offset
        readlen = read.infer_query_length()
        if readlen is None:
            return False  # Skip reads without queryable length

        # Feed to both calculators
        try:
            if read.is_reverse:
                self.primary_calculator.feed_reverse_read(chrom, pos, readlen)
                self.secondary_calculator.feed_reverse_read(chrom, pos, readlen)
            else:
                self.primary_calculator.feed_forward_read(chrom, pos, readlen)
                self.secondary_calculator.feed_forward_read(chrom, pos, readlen)

            return True

        except Exception as e:
            if "sorted" in str(e).lower() or "order" in str(e).lower():
                raise ReadUnsortedError("Input alignment file must be sorted") from e
            else:
                raise


# Convenience functions for backward compatibility and easy adoption
def create_standard_filter(mapq_criteria: int) -> ReadFilter:
    """Create standard PyMaSC read filter.

    Args:
        mapq_criteria: Minimum mapping quality threshold

    Returns:
        ReadFilter with standard PyMaSC settings
    """
    return ReadFilter(
        mapq_criteria=mapq_criteria,
        skip_read2=True,
        skip_unmapped=True,
        skip_duplicates=True
    )


def create_read_processor(calculator: ReadCalculator, mapq_criteria: int) -> ReadProcessor:
    """Create standard read processor with PyMaSC filter settings.

    Args:
        calculator: Calculator to feed reads to
        mapq_criteria: Minimum mapping quality threshold

    Returns:
        ReadProcessor with standard configuration
    """
    filter_obj = create_standard_filter(mapq_criteria)
    return ReadProcessor(calculator, filter_obj)


def create_dual_processor(primary_calc: ReadCalculator, 
                         secondary_calc: ReadCalculator,
                         mapq_criteria: int) -> DualReadProcessor:
    """Create dual read processor with standard PyMaSC filter settings.

    Args:
        primary_calc: Primary calculator
        secondary_calc: Secondary calculator  
        mapq_criteria: Minimum mapping quality threshold

    Returns:
        DualReadProcessor with standard configuration
    """
    filter_obj = create_standard_filter(mapq_criteria)
    return DualReadProcessor(primary_calc, secondary_calc, filter_obj)


# Legacy compatibility - reproduce exact filtering logic from old code
def legacy_should_skip_read(read: AlignedSegment, mapq_criteria: int) -> bool:
    """Legacy read filtering function for exact backward compatibility.

    This function reproduces the exact filtering logic from the original
    worker implementations to ensure identical behavior during migration.

    Args:
        read: Read to check
        mapq_criteria: Minimum mapping quality threshold

    Returns:
        True if read should be skipped
    """
    return (read.is_read2 or 
            read.mapping_quality < mapq_criteria or
            read.is_unmapped or 
            read.is_duplicate)