"""BAM file processing and validation.

This module provides specialized functionality for BAM file handling,
including validation, chromosome filtering, and metadata extraction.
It separates BAM file concerns from calculation logic.
"""
from __future__ import annotations

import logging
import weakref
from typing import (
    runtime_checkable, Protocol, Literal, Optional,
    Any, Iterable, Iterator, List, Tuple, Dict
)

import sys
if sys.version_info >= (3, 9):
    from typing import Self
else:
    from typing_extensions import Self

import pysam

from PyMaSC.utils.calc import filter_chroms

logger = logging.getLogger(__name__)


class BAMValidationError(ValueError):
    """Raised when BAM file validation fails."""
    pass


class BAMNoReadsError(BAMValidationError):
    """Raised when no reads are found in the BAM file."""
    pass


class BAMNoTargetChroms(BAMValidationError):
    """Raised when no target chromosomes are found in the BAM file."""
    pass


@runtime_checkable
class AlignmentLike(Protocol):
    def __enter__(self) -> Self: ...

    def __exit__(self, exc_type: Optional[type], exc: Optional[BaseException], tb: Any) -> bool: ...

    def close(self) -> None: ...

    def has_index(self) -> bool: ...

    @property
    def closed(self) -> bool: ...

    @property
    def references(self) -> Iterable[str]: ...

    @property
    def lengths(self) -> Iterable[int]: ...

    def __iter__(self) -> Iterator: ...

    def fetch(self, *args: Any, **kwargs: Any) -> pysam.IteratorRow: ...


class BAMFileProcessor(AlignmentLike):
    _filtered_references: Optional[Tuple[str, ...]]
    _filtered_lengths: Optional[Tuple[int, ...]]

    def __init__(self, path: str, *args: Any, **kwargs: Any):
        self._path = path
        self._af = pysam.AlignmentFile(path, *args, **kwargs)
        self._finalizer = weakref.finalize(self, self._safe_close, self._af)
        self._filtered_references = None
        self._filtered_lengths = None

    def close(self) -> None:
        if self._finalizer and self._finalizer.alive:
            self._finalizer()
        if self._af and not self._af.closed:
            self._af.close()

    @staticmethod
    def _safe_close(af: pysam.AlignmentFile) -> None:
        try:
            if af and not af.closed:
                af.close()
        except Exception:
            pass

    def __enter__(self) -> Self:
        return self

    def __exit__(self, exc_type: Optional[type], exc: Optional[BaseException], tb: Optional[Any]) -> Literal[False]:
        self.close()
        return False

    def has_index(self) -> bool:
        return self._af.has_index()

    @property
    def closed(self) -> bool:
        return self._af.closed

    @property
    def references(self) -> Tuple[str, ...]:
        return self._af.references

    @property
    def lengths(self) -> Tuple[int, ...]:
        return self._af.lengths

    def __iter__(self) -> Iterator[pysam.AlignedSegment]:
        return iter(self._af)

    def fetch(self, *args: Any, **kwargs: Any) -> pysam.IteratorRow:
        return self._af.fetch(*args, **kwargs)

    def apply_chromfilter(self, chromfilter: Optional[List[Tuple[bool, List[str]]]] = None) -> Tuple[Tuple[str, ...], Tuple[int, ...]]:
        """
        """
        all_references = list(self._af.references)
        all_lengths = list(self._af.lengths)

        if not all_references:
            logger.error("BAM file has no reference sequences")
            raise BAMNoReadsError

        target_references = filter_chroms(all_references, chromfilter)
        if not target_references:
            logger.error("There is no targeted chromosomes.")
            raise BAMNoTargetChroms

        filtered_ref: Tuple[str, ...]
        filtered_length: Tuple[int, ...]
        filtered_ref, filtered_length = zip(
            *((ref, length) for ref, length in zip(all_references, all_lengths)
              if ref in target_references)
        )

        self._filtered_references = filtered_ref
        self._filtered_lengths = filtered_length

        return filtered_ref, filtered_length

    def validate_chromosome_sizes(self, external_sizes: Dict[str, int]) -> Dict[str, int]:
        """Validate and reconcile chromosome sizes with external source.

        Args:
            external_sizes: External chromosome sizes (e.g., from BigWig)

        Returns:
            Updated chromosome sizes
        """
        if self._filtered_references is None or self._filtered_lengths is None:
            raise RuntimeError("Must call validate_and_open first")

        filtered_sizes = dict(zip(self._filtered_references, self._filtered_lengths))
        updated_sizes = {}

        for ref, bam_size in filtered_sizes.items():
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
        if require_index and not self.has_index():
            logger.error(
                "Need indexed alignment file for multi-processing. "
                "Calculation will be executed by single process."
            )
            return False

        return True
