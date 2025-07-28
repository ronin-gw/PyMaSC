"""BAM file processing and validation.

This module provides specialized functionality for BAM file handling,
including validation, chromosome filtering, and metadata extraction.
It separates BAM file concerns from calculation logic.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import List, Optional, Dict, Tuple, Union
import os

import pysam

from PyMaSC.utils.calc import filter_chroms

logger = logging.getLogger(__name__)


@dataclass
class BAMMetadata:
    """Container for BAM file metadata.

    Attributes:
        path: Path to BAM file
        has_index: Whether BAM file has an index
        references: List of reference sequence names
        lengths: List of reference sequence lengths
        filtered_references: References after filtering
        filtered_lengths: Lengths after filtering
    """
    path: str
    has_index: bool
    references: List[str]
    lengths: List[int]
    filtered_references: List[str]
    filtered_lengths: List[int]


class BAMValidationError(ValueError):
    """Raised when BAM file validation fails."""
    pass


class BAMFileProcessor:
    """Handles BAM file processing and validation.

    This class encapsulates all BAM file-specific operations,
    providing a clean interface for file validation, chromosome
    filtering, and metadata extraction.
    """

    def __init__(self, path: Union[str, os.PathLike[str]]):
        """Initialize BAM processor.

        Args:
            path: Path to BAM file
        """
        self.path = str(path)
        self._align_file: Optional[pysam.AlignmentFile] = None
        self._metadata: Optional[BAMMetadata] = None

    def validate_and_open(self, chromfilter: Optional[List[Tuple[bool, List[str]]]] = None) -> BAMMetadata:
        """Validate BAM file and extract metadata.

        Args:
            chromfilter: Optional chromosome filter in format [(include_flag, patterns)]

        Returns:
            BAM file metadata

        Raises:
            BAMValidationError: If validation fails
        """
        # Open BAM file
        try:
            self._align_file = pysam.AlignmentFile(self.path)
        except ValueError as e:
            logger.error("File has no sequences defined.")
            raise BAMValidationError("File has no sequences defined.") from e

        # Extract basic metadata
        all_references = list(self._align_file.references)
        all_lengths = list(self._align_file.lengths)

        if not all_references:
            raise BAMValidationError("BAM file has no reference sequences")

        # Apply chromosome filtering
        target_references = filter_chroms(all_references, chromfilter)
        if not target_references:
            logger.error("There is no targeted chromosomes.")
            raise BAMValidationError("No chromosomes match filter criteria")

        # Build filtered lists
        filtered_refs = []
        filtered_lens = []
        for ref, length in zip(all_references, all_lengths):
            if ref in target_references:
                filtered_refs.append(ref)
                filtered_lens.append(length)

        # Check for index
        has_index = self._align_file.has_index()

        # Create metadata
        self._metadata = BAMMetadata(
            path=self.path,
            has_index=has_index,
            references=all_references,
            lengths=all_lengths,
            filtered_references=filtered_refs,
            filtered_lengths=filtered_lens
        )

        return self._metadata

    def validate_chromosome_sizes(self, external_sizes: Dict[str, int]) -> Dict[str, int]:
        """Validate and reconcile chromosome sizes with external source.

        Args:
            external_sizes: External chromosome sizes (e.g., from BigWig)

        Returns:
            Updated chromosome sizes
        """
        if not self._metadata:
            raise RuntimeError("Must call validate_and_open first")

        updated_sizes = {}

        for i, ref in enumerate(self._metadata.filtered_references):
            bam_size = self._metadata.filtered_lengths[i]

            if ref not in external_sizes:
                logger.debug(f"External size for '{ref}' not found")
                updated_sizes[ref] = bam_size
                continue

            external_size = external_sizes[ref]

            if external_size != bam_size:
                logger.warning(
                    f"'{ref}' reference length mismatch: "
                    f"SAM/BAM -> {bam_size:,}, External -> {external_size:,}"
                )
                if bam_size < external_size:
                    logger.warning(
                        f"Use longer length '{external_size:d}' for '{ref}' anyway"
                    )
                    updated_sizes[ref] = external_size
                else:
                    updated_sizes[ref] = bam_size
            else:
                updated_sizes[ref] = bam_size

        return updated_sizes

    def check_multiprocess_compatibility(self, require_index: bool = True) -> bool:
        """Check if BAM file is compatible with multiprocessing.

        Args:
            require_index: Whether index is required

        Returns:
            True if compatible, False otherwise
        """
        if not self._metadata:
            raise RuntimeError("Must call validate_and_open first")

        if require_index and not self._metadata.has_index:
            logger.error(
                "Need indexed alignment file for multi-processing. "
                "Calculation will be executed by single process."
            )
            return False

        return True

    def close(self) -> None:
        """Close BAM file if open."""
        if self._align_file and not self._align_file.closed:
            self._align_file.close()
            self._align_file = None

    def reopen(self) -> pysam.AlignmentFile:
        """Reopen BAM file for reading.

        Returns:
            Opened alignment file
        """
        if self._align_file and not self._align_file.closed:
            return self._align_file

        self._align_file = pysam.AlignmentFile(self.path)
        return self._align_file

    @property
    def metadata(self) -> Optional[BAMMetadata]:
        """Get BAM metadata if available."""
        return self._metadata

    @property
    def is_stdin(self) -> bool:
        """Check if reading from stdin."""
        return self.path == '-'

    def __enter__(self) -> 'BAMFileProcessor':
        """Context manager entry."""
        return self

    def __exit__(self, _exc_type: Optional[type], _exc_val: Optional[BaseException], _exc_tb: Optional[object]) -> None:
        """Context manager exit."""
        self.close()
