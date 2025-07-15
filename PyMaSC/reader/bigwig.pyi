"""Type stubs for bigwig.pyx - BigWig file reader interface for mappability data."""

from typing import Optional, List, Tuple, Union, Iterator, Dict, Any
from os import PathLike


class BigWigFileIterator:
    """Iterator for BigWig intervals with original API compatibility.
    
    Provides iteration over BigWig intervals while maintaining compatibility
    with the original BigWigFile API. The iterator yields interval structures
    containing genomic coordinates and associated values.
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
    """pyBigWig-based BigWig file reader with original API compatibility.
    
    Provides efficient access to BigWig files containing mappability data
    for MSCC calculations. The class wraps pyBigWig functionality while
    maintaining compatibility with the original bx-python-based API.
    
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
        
        Opens the BigWig file and extracts chromosome size information
        for efficient data access during MSCC calculations.
        
        Args:
            path: Path to BigWig file to read
            
        Raises:
            IOError: If the specified file does not exist
            RuntimeError: If the BigWig file cannot be opened
        """
        ...
    
    def fetch(self, valfilter: float, chrom: str) -> BigWigFileIterator:
        """Fetch intervals from BigWig file for specified chromosome.
        
        Retrieves all intervals from the specified chromosome that meet
        the value filter criteria. Returns an iterator for efficient
        processing of large datasets.
        
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
        """Disable progress bar for compatibility with MSCC calculations.
        
        This method provides compatibility with the original bx-python API
        but has no actual effect since pyBigWig doesn't have progress bars.
        """
        ...
    
    def close(self) -> None:
        """Close the BigWig file and release resources.
        
        Closes the underlying pyBigWig file object and marks the
        reader as closed to prevent further operations.
        """
        ...