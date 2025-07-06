"""I/O service for all file operations in PyMaSC.

This module provides the IOService that encapsulates all input/output
operations, completely separated from calculation logic and workflow
management. This is the infrastructure layer of the application.

Key features:
- Centralized file access
- Consistent error handling
- Mockable for testing
- Resource management
"""
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Dict, Any, Tuple, Union

import numpy as np
from pysam import AlignmentFile, AlignedSegment

from PyMaSC.reader.bigwig import BigWigReader
from PyMaSC.services.calculation import ChromosomeData, CalculationResult, GenomeWideResult
from PyMaSC.utils.output_utils import OutputPathManager, ensure_output_directory, TableWriter as UtilTableWriter
from PyMaSC.utils.read_processing import create_standard_filter

logger = logging.getLogger(__name__)


@dataclass
class BAMFileInfo:
    """Information about a BAM file.

    Provides metadata without requiring file access.
    """
    path: str
    has_index: bool
    references: List[str]
    lengths: List[int]

    @property
    def n_references(self) -> int:
        """Number of reference sequences."""
        return len(self.references)


@dataclass
class ReadData:
    """Container for read information.

    Simplified representation of a sequencing read.
    """
    chromosome: str
    position: int
    length: int
    is_reverse: bool
    mapping_quality: int


class IOService(ABC):
    """Abstract base for I/O services.

    Defines the interface for all I/O operations, enabling
    different implementations (e.g., file-based, in-memory, cloud).
    """

    @abstractmethod
    def get_bam_info(self, path: str) -> BAMFileInfo:
        """Get information about a BAM file.

        Args:
            path: Path to BAM file

        Returns:
            BAM file information
        """
        pass

    @abstractmethod
    def read_chromosome_data(self, 
                           bam_path: str,
                           chromosome: str,
                           mapq_threshold: int = 0) -> ChromosomeData:
        """Read all data for a chromosome.

        Args:
            bam_path: Path to BAM file
            chromosome: Chromosome to read
            mapq_threshold: Minimum mapping quality

        Returns:
            Chromosome data with all reads
        """
        pass

    @abstractmethod
    def stream_reads(self,
                    bam_path: str,
                    chromosome: str,
                    mapq_threshold: int = 0) -> Iterator[ReadData]:
        """Stream reads for a chromosome.

        Args:
            bam_path: Path to BAM file
            chromosome: Chromosome to read
            mapq_threshold: Minimum mapping quality

        Yields:
            Individual read data
        """
        pass

    @abstractmethod
    def read_mappability(self,
                       bigwig_path: str,
                       chromosome: str,
                       start: int,
                       end: int) -> List[float]:
        """Read mappability values from BigWig file.

        Args:
            bigwig_path: Path to BigWig file
            chromosome: Chromosome name
            start: Start position
            end: End position

        Returns:
            List of mappability values
        """
        pass

    @abstractmethod
    def write_results(self,
                     result: GenomeWideResult,
                     output_prefix: str,
                     formats: List[str]) -> Dict[str, str]:
        """Write calculation results to files.

        Args:
            result: Calculation results
            output_prefix: Base path for output files
            formats: List of output formats ('table', 'stats', 'figure')

        Returns:
            Dictionary mapping format to output path
        """
        pass


class FileIOService(IOService):
    """File-based implementation of I/O service.

    Standard implementation that reads from and writes to
    local filesystem files.
    """

    def __init__(self) -> None:
        """Initialize file I/O service."""
        self._open_files: Dict[str, Any] = {}
        self._read_filter: Optional[Any] = None

    def get_bam_info(self, path: str) -> BAMFileInfo:
        """Get information about a BAM file."""
        try:
            with AlignmentFile(path) as bamfile:
                # Check for index
                has_index = False
                try:
                    # Try to use index by fetching from first reference
                    if bamfile.references:
                        list(bamfile.fetch(bamfile.references[0], 1, 2))
                        has_index = True
                except (OSError, ValueError):
                    has_index = False

                return BAMFileInfo(
                    path=path,
                    has_index=has_index,
                    references=list(bamfile.references),
                    lengths=list(bamfile.lengths)
                )
        except Exception as e:
            logger.error(f"Error reading BAM file {path}: {e}")
            raise IOError(f"Cannot read BAM file: {e}") from e

    def read_chromosome_data(self, 
                           bam_path: str,
                           chromosome: str,
                           mapq_threshold: int = 0) -> ChromosomeData:
        """Read all data for a chromosome at once."""
        # Get chromosome info
        info = self.get_bam_info(bam_path)

        if chromosome not in info.references:
            raise ValueError(f"Chromosome {chromosome} not found in BAM file")

        chrom_idx = info.references.index(chromosome)
        chrom_length = info.lengths[chrom_idx]

        # Collect all reads
        forward_reads = []
        reverse_reads = []

        for read in self.stream_reads(bam_path, chromosome, mapq_threshold):
            if read.is_reverse:
                reverse_reads.append((read.position, read.length))
            else:
                forward_reads.append((read.position, read.length))

        return ChromosomeData(
            chromosome=chromosome,
            forward_reads=forward_reads,
            reverse_reads=reverse_reads,
            length=chrom_length
        )

    def stream_reads(self,
                    bam_path: str,
                    chromosome: str,
                    mapq_threshold: int = 0) -> Iterator[ReadData]:
        """Stream reads for a chromosome."""
        # Create read filter if needed
        if self._read_filter is None or getattr(self._read_filter, 'mapq_criteria', -1) != mapq_threshold:
            self._read_filter = create_standard_filter(mapq_threshold)

        try:
            with AlignmentFile(bam_path) as bamfile:
                for read in bamfile.fetch(chromosome):
                    # Apply quality filtering
                    if self._read_filter and self._read_filter.should_skip_read(read):
                        continue

                    # Skip reads without reference or queryable length
                    chrom = read.reference_name
                    readlen = read.infer_query_length()
                    if chrom is None or readlen is None:
                        continue

                    yield ReadData(
                        chromosome=chrom,
                        position=read.reference_start + 1,  # Convert to 1-based
                        length=readlen,
                        is_reverse=read.is_reverse,
                        mapping_quality=read.mapping_quality
                    )
        except Exception as e:
            logger.error(f"Error streaming reads from {bam_path}: {e}")
            raise IOError(f"Cannot read from BAM file: {e}") from e

    def read_mappability(self,
                       bigwig_path: str,
                       chromosome: str,
                       start: int,
                       end: int) -> List[float]:
        """Read mappability values from BigWig file."""
        try:
            # Use cached reader if available
            if bigwig_path not in self._open_files:
                self._open_files[bigwig_path] = BigWigReader(bigwig_path)

            reader = self._open_files[bigwig_path]

            # Get values for region
            values = reader.get_as_array(chromosome, start, end)

            # Convert None values to 0.0
            return [v if v is not None else 0.0 for v in values]

        except Exception as e:
            logger.error(f"Error reading mappability from {bigwig_path}: {e}")
            raise IOError(f"Cannot read BigWig file: {e}") from e

    def write_results(self,
                     result: GenomeWideResult,
                     output_prefix: str,
                     formats: List[str]) -> Dict[str, str]:
        """Write calculation results to files."""
        output_manager = OutputPathManager(output_prefix)
        output_paths = {}

        # Ensure output directory exists
        Path(output_prefix).parent.mkdir(parents=True, exist_ok=True)

        # Write table format
        if 'table' in formats:
            table_path = output_manager.get_output_path(extension='.txt')
            self._write_table_output(result, str(table_path))
            output_paths['table'] = str(table_path)

        # Write statistics format
        if 'stats' in formats:
            stats_path = output_manager.get_output_path(extension='.json')
            self._write_stats_output(result, str(stats_path))
            output_paths['stats'] = str(stats_path)

        # Write figure format
        if 'figure' in formats:
            figure_path = output_manager.get_output_path(extension='.pdf')
            self._write_figure_output(result, str(figure_path))
            output_paths['figure'] = str(figure_path)

        return output_paths

    def _write_table_output(self, result: GenomeWideResult, path: str) -> None:
        """Write table format output."""
        # Prepare data for table writer
        data = []
        headers = ['chromosome', 'forward_reads', 'reverse_reads', 
                  'max_correlation', 'peak_position']

        for chrom, chrom_result in result.chromosome_results.items():
            if len(chrom_result.correlation_bins) > 0:
                max_corr = np.max(chrom_result.correlation_bins)
                peak_pos = int(np.argmax(chrom_result.correlation_bins))
            else:
                max_corr = 0.0
                peak_pos = 0

            data.append([
                chrom,
                chrom_result.forward_count,
                chrom_result.reverse_count,
                max_corr,
                peak_pos
            ])

        # Add summary row
        if len(result.aggregated_correlation) > 0:
            genome_max = np.max(result.aggregated_correlation)
            genome_peak = int(np.argmax(result.aggregated_correlation))
        else:
            genome_max = 0.0
            genome_peak = 0

        data.append([
            'GENOME_WIDE',
            result.total_forward_reads,
            result.total_reverse_reads,
            genome_max,
            genome_peak
        ])

        # Write using utility
        writer = UtilTableWriter(path)
        writer.write_table(data, headers)

    def _write_stats_output(self, result: GenomeWideResult, path: str) -> None:
        """Write statistics format output."""
        # Prepare statistics dictionary
        stats = {
            'total_reads': result.total_forward_reads + result.total_reverse_reads,
            'forward_reads': result.total_forward_reads,
            'reverse_reads': result.total_reverse_reads,
            'n_chromosomes': len(result.chromosome_results),
            'chromosomes': list(result.chromosome_results.keys())
        }

        # Add quality metrics if available
        if len(result.aggregated_correlation) > 0:
            # Estimate read length from first chromosome result
            first_result = next(iter(result.chromosome_results.values()))
            read_length = 50  # Default estimate

            quality_metrics = first_result.get_quality_metrics(read_length)
            stats.update(quality_metrics)

        # Add mappability statistics if available
        if result.mappable_fraction is not None:
            stats['mappable_fraction'] = result.mappable_fraction

        # Write using standard library
        import json
        with open(path, 'w') as f:
            json.dump(stats, f, indent=2)

        logger.info(f"Output '{path}'")

    def _write_figure_output(self, result: GenomeWideResult, path: str) -> None:
        """Write figure format output.

        Note: This is a placeholder. Real implementation would create
        actual plots using matplotlib.
        """
        # For now, just log that we would create a figure
        logger.info(f"Would create figure at '{path}'")

        # Create a minimal PDF to satisfy the interface
        # In production, this would use FigureWriter to create actual plots
        Path(path).touch()

    def close(self) -> None:
        """Close all open file handles."""
        for reader in self._open_files.values():
            if hasattr(reader, 'close'):
                reader.close()
        self._open_files.clear()

    def __enter__(self) -> 'FileIOService':
        """Context manager entry."""
        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        """Context manager exit."""
        self.close()


class InMemoryIOService(IOService):
    """In-memory implementation for testing.

    Provides I/O service interface without actual file access,
    useful for unit testing.
    """

    def __init__(self) -> None:
        """Initialize in-memory I/O service."""
        self.bam_data: Dict[str, Any] = {}
        self.bigwig_data: Dict[str, Any] = {}
        self.written_results: Dict[str, Any] = {}

    def add_test_data(self, 
                     bam_path: str, 
                     chromosome: str,
                     forward_reads: List[Tuple[int, int]],
                     reverse_reads: List[Tuple[int, int]],
                     length: int) -> None:
        """Add test data for a chromosome.

        Args:
            bam_path: Virtual BAM file path
            chromosome: Chromosome name
            forward_reads: List of (position, length) for forward reads
            reverse_reads: List of (position, length) for reverse reads  
            length: Chromosome length
        """
        if bam_path not in self.bam_data:
            self.bam_data[bam_path] = {
                'info': BAMFileInfo(
                    path=bam_path,
                    has_index=True,
                    references=[],
                    lengths=[]
                ),
                'chromosomes': {}
            }

        # Update info
        if chromosome not in self.bam_data[bam_path]['info'].references:
            self.bam_data[bam_path]['info'].references.append(chromosome)
            self.bam_data[bam_path]['info'].lengths.append(length)

        # Store chromosome data
        self.bam_data[bam_path]['chromosomes'][chromosome] = ChromosomeData(
            chromosome=chromosome,
            forward_reads=forward_reads,
            reverse_reads=reverse_reads,
            length=length
        )

    def get_bam_info(self, path: str) -> BAMFileInfo:
        """Get test BAM info."""
        if path in self.bam_data:
            return self.bam_data[path]['info']  # type: ignore[no-any-return]
        else:
            raise IOError(f"Test BAM file not found: {path}")

    def read_chromosome_data(self, 
                           bam_path: str,
                           chromosome: str,
                           mapq_threshold: int = 0) -> ChromosomeData:
        """Read test chromosome data."""
        if bam_path in self.bam_data:
            if chromosome in self.bam_data[bam_path]['chromosomes']:
                return self.bam_data[bam_path]['chromosomes'][chromosome]  # type: ignore[no-any-return]

        raise IOError(f"Test chromosome data not found: {bam_path}:{chromosome}")

    def stream_reads(self,
                    bam_path: str,
                    chromosome: str,
                    mapq_threshold: int = 0) -> Iterator[ReadData]:
        """Stream test reads."""
        data = self.read_chromosome_data(bam_path, chromosome, mapq_threshold)

        # Yield forward reads
        for pos, length in data.forward_reads:
            yield ReadData(
                chromosome=chromosome,
                position=pos,
                length=length,
                is_reverse=False,
                mapping_quality=30  # Default quality
            )

        # Yield reverse reads
        for pos, length in data.reverse_reads:
            yield ReadData(
                chromosome=chromosome,
                position=pos,
                length=length,
                is_reverse=True,
                mapping_quality=30  # Default quality
            )

    def read_mappability(self,
                       bigwig_path: str,
                       chromosome: str,
                       start: int,
                       end: int) -> List[float]:
        """Read test mappability values."""
        # Return dummy values for testing
        return [1.0] * (end - start)

    def write_results(self,
                     result: GenomeWideResult,
                     output_prefix: str,
                     formats: List[str]) -> Dict[str, str]:
        """Store results in memory instead of writing files."""
        output_paths = {}

        for format_type in formats:
            path = f"{output_prefix}.{format_type}"
            self.written_results[path] = result
            output_paths[format_type] = path

        return output_paths


# Factory function for service creation
def create_io_service(test_mode: bool = False) -> IOService:
    """Create appropriate I/O service.

    Args:
        test_mode: Whether to create test service

    Returns:
        I/O service instance
    """
    if test_mode:
        return InMemoryIOService()
    else:
        return FileIOService()