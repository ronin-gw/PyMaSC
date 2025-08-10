"""Type stubs for bigwig.pyx - BigWig file reader interface for mappability data."""

from typing import Optional, List, Tuple, Union, Iterator, Dict, Any
from os import PathLike


class BigWigFileIterator:
    """Iterator for BigWig intervals.

    Yields interval tuples containing genomic coordinates and values.
    """

    def __init__(self, intervals: Optional[List[Tuple[int, int, float]]]) -> None:
        """Initialize iterator with interval data.

        Args:
            intervals: List of interval tuples (start, end, value) or None
        """
        ...

    def __iter__(self) -> Iterator[Tuple[int, int, float]]:
        """Return iterator object for Python iteration protocol.

        Returns:
            Self as the iterator object
        """
        ...

    def __next__(self) -> Tuple[int, int, float]:
        """Get next interval for Python iteration protocol.

        Returns:
            Tuple of (start, end, value) for the next interval

        Raises:
            StopIteration: When no more intervals are available
        """
        ...


class BigWigReader:
    """BigWig file reader for mappability data.

    Reads BigWig files containing mappability data for MSCC calculations.

    Attributes:
        path: Path to the BigWig file
        closed: Whether the file has been closed
        chromsizes: Dictionary mapping chromosome names to lengths
    """

    path: str
    closed: bool
    file: Any  # pyBigWig file object
    chromsizes: Dict[str, int]

    def __init__(self, path: Union[str, PathLike[str]]) -> None:
        """Initialize BigWig reader.

        Args:
            path: Path to BigWig file to read

        Raises:
            IOError: If the specified file does not exist
            RuntimeError: If the BigWig file cannot be opened
        """
        ...

    def fetch(self, valfilter: float, chrom: str) -> BigWigFileIterator:
        """Fetch intervals from specified chromosome.

        Args:
            valfilter: Minimum value threshold for intervals
            chrom: Chromosome name to fetch data from

        Returns:
            BigWigFileIterator containing filtered intervals

        Note:
            Only intervals with values >= valfilter are included

        Raises:
            KeyError: If chromosome is not found in the BigWig file
        """
        ...

    def disable_progress_bar(self) -> None:
        """Disable progress bar (compatibility method)."""
        ...

    def close(self) -> None:
        """Close the BigWig file."""
        ...
