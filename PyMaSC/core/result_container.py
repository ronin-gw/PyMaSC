"""Type-safe containers for cross-correlation results.

This module provides strongly-typed container classes for organizing
cross-correlation calculation results. It addresses the type confusion
in the original implementation where values could be scalars or arrays.

Key improvements:
- Clear separation between NCC (scalar) and MSCC (array) data
- Type-safe access to results
- Validation of data consistency
- Clear data flow between components
"""
from __future__ import annotations

from abc import ABC
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Any

import numpy as np


@dataclass
class NCCData:
    """Container for Naive Cross-Correlation (NCC) data.

    NCC data uses scalar values for read counts since it doesn't
    account for mappability variations across shift positions.

    Attributes:
        forward_sum: Total forward read count (scalar)
        reverse_sum: Total reverse read count (scalar)
        ccbins: Cross-correlation values array
        genomelen: Total genome length
    """
    forward_sum: int
    reverse_sum: int
    ccbins: np.ndarray
    genomelen: int

    def __post_init__(self):
        """Validate data consistency."""
        if self.forward_sum < 0 or self.reverse_sum < 0:
            raise ValueError("Read counts must be non-negative")
        if self.genomelen <= 0:
            raise ValueError("Genome length must be positive")
        if self.ccbins is None:
            raise ValueError("Cross-correlation bins cannot be None")
        if len(self.ccbins) == 0:
            raise ValueError("Cross-correlation bins cannot be empty")


@dataclass
class MSCCData:
    """Container for Mappability-Sensitive Cross-Correlation (MSCC) data.

    MSCC data uses arrays for read counts as they vary by shift position
    based on mappability considerations.

    Attributes:
        forward_sum: Forward read counts by shift position (array)
        reverse_sum: Reverse read counts by shift position (array)
        ccbins: Mappability-adjusted cross-correlation values
        mappable_len: Mappable genome length (scalar or array)
    """
    forward_sum: np.ndarray
    reverse_sum: np.ndarray
    ccbins: np.ndarray
    mappable_len: Union[int, np.ndarray, None]

    def __post_init__(self):
        """Validate data consistency."""
        # Check for None values first
        if self.ccbins is None:
            raise ValueError("MSCC cross-correlation bins cannot be None")
        if self.forward_sum is None:
            raise ValueError("MSCC forward sum cannot be None")
        if self.reverse_sum is None:
            raise ValueError("MSCC reverse sum cannot be None")

        # Ensure all arrays have the same length
        lengths = [len(self.forward_sum), len(self.reverse_sum), len(self.ccbins)]
        if len(set(lengths)) > 1:
            raise ValueError(f"Array lengths must match: {lengths}")

        # Validate mappable_len type and value with more flexible checking
        if isinstance(self.mappable_len, np.ndarray):
            if len(self.mappable_len) == 0:
                raise ValueError("Mappable length array cannot be empty")
        elif isinstance(self.mappable_len, (list, tuple)):
            # Convert list/tuple to array
            self.mappable_len = np.array(self.mappable_len)
            if len(self.mappable_len) == 0:
                raise ValueError("Mappable length array cannot be empty")
        elif isinstance(self.mappable_len, (int, np.integer, float, np.floating)):
            # Accept numeric types and convert to int
            self.mappable_len = int(self.mappable_len)
            if self.mappable_len <= 0:
                raise ValueError("Mappable length must be positive")
        elif self.mappable_len is None:
            # Allow None for missing data
            pass
        else:
            raise ValueError(f"Mappable length must be int, array, or None, got {type(self.mappable_len)}: {self.mappable_len}")


@dataclass
class ChromosomeResult:
    """Container for per-chromosome calculation results.

    Stores both NCC and MSCC results for a single chromosome,
    handling cases where one or both analyses are performed.

    Attributes:
        chromosome: Chromosome name
        length: Chromosome length
        ncc_data: Optional NCC results
        mscc_data: Optional MSCC results
    """
    chromosome: str
    length: int
    ncc_data: Optional[NCCData] = None
    mscc_data: Optional[MSCCData] = None

    def __post_init__(self):
        """Validate that at least one analysis type is present."""
        if self.ncc_data is None and self.mscc_data is None:
            raise ValueError(f"Chromosome {self.chromosome} has no analysis data")
        if self.length <= 0:
            raise ValueError(f"Chromosome {self.chromosome} length must be positive")

    @property
    def has_ncc(self) -> bool:
        """Check if NCC data is available."""
        return self.ncc_data is not None

    @property
    def has_mscc(self) -> bool:
        """Check if MSCC data is available."""
        return self.mscc_data is not None

    @property
    def max_shift(self) -> int:
        """Get maximum shift distance from available data."""
        if self.ncc_data:
            return len(self.ncc_data.ccbins) - 1
        elif self.mscc_data:
            return len(self.mscc_data.ccbins) - 1
        return 0


@dataclass
class ResultContainer:
    """Main container for all cross-correlation results.

    Provides organized storage and access to calculation results
    across all chromosomes, with support for both NCC and MSCC analyses.

    Attributes:
        read_length: Read length in base pairs
        references: List of analyzed chromosome names
        chromosome_results: Per-chromosome results
        skip_ncc: Whether NCC calculation was skipped
        has_mappability: Whether mappability data was used
    """
    read_length: int
    references: List[str]
    chromosome_results: Dict[str, ChromosomeResult] = field(default_factory=dict)
    skip_ncc: bool = False
    has_mappability: bool = False

    def __post_init__(self):
        """Validate container consistency."""
        if self.read_length <= 0:
            raise ValueError("Read length must be positive")
        if not self.references:
            raise ValueError("References list cannot be empty")

    def add_chromosome_result(self, result: ChromosomeResult) -> None:
        """Add a chromosome result to the container.

        Args:
            result: ChromosomeResult to add

        Raises:
            ValueError: If chromosome is not in references list
        """
        if result.chromosome not in self.references:
            raise ValueError(f"Chromosome {result.chromosome} not in references")
        self.chromosome_results[result.chromosome] = result

        # Update flags based on available data
        if result.has_mscc:
            self.has_mappability = True

    def get_ncc_data(self) -> Dict[str, NCCData]:
        """Get all available NCC data by chromosome.

        Returns:
            Dictionary mapping chromosome names to NCC data
        """
        return {
            chrom: result.ncc_data
            for chrom, result in self.chromosome_results.items()
            if result.ncc_data is not None
        }

    def get_mscc_data(self) -> Dict[str, MSCCData]:
        """Get all available MSCC data by chromosome.

        Returns:
            Dictionary mapping chromosome names to MSCC data
        """
        return {
            chrom: result.mscc_data
            for chrom, result in self.chromosome_results.items()
            if result.mscc_data is not None
        }

    def get_total_reads(self) -> Dict[str, int]:
        """Calculate total read counts by chromosome.

        Returns:
            Dictionary with forward, reverse, and total counts
        """
        forward_total = 0
        reverse_total = 0

        for result in self.chromosome_results.values():
            if result.ncc_data:
                forward_total += result.ncc_data.forward_sum
                reverse_total += result.ncc_data.reverse_sum

        return {
            'forward': forward_total,
            'reverse': reverse_total,
            'total': forward_total + reverse_total
        }

    def get_genome_length(self) -> int:
        """Calculate total genome length across all chromosomes.

        Returns:
            Total genome length in base pairs
        """
        return sum(result.length for result in self.chromosome_results.values())

    def validate_completeness(self) -> None:
        """Validate that all expected chromosomes have results.

        Raises:
            ValueError: If any chromosome is missing results
        """
        missing = set(self.references) - set(self.chromosome_results.keys())
        if missing:
            raise ValueError(f"Missing results for chromosomes: {missing}")


def create_from_handler_data(
    handler_data: Dict[str, Any],
    read_length: int,
    references: List[str],
    lengths: List[int]
) -> ResultContainer:
    """Create ResultContainer from handler output data.

    Factory function to create a properly structured ResultContainer
    from the dictionary-based output of calculation handlers.

    Args:
        handler_data: Dictionary with handler calculation results
        read_length: Read length in base pairs
        references: List of chromosome names
        lengths: List of chromosome lengths

    Returns:
        Populated ResultContainer instance
    """
    container = ResultContainer(
        read_length=read_length,
        references=references,
        skip_ncc=handler_data.get('skip_ncc', False)
    )

    # Create chromosome results
    ref_to_length = dict(zip(references, lengths))

    for ref in references:
        ncc_data = None
        mscc_data = None

        # Extract NCC data if available
        if (not container.skip_ncc and
            ref in handler_data.get('ref2ccbins', {}) and
            handler_data['ref2ccbins'][ref] is not None):
            ncc_data = NCCData(
                forward_sum=handler_data['ref2forward_sum'].get(ref, 0),
                reverse_sum=handler_data['ref2reverse_sum'].get(ref, 0),
                ccbins=handler_data['ref2ccbins'][ref],
                genomelen=ref_to_length[ref]
            )

        # Extract MSCC data if available
        mappable_forward = handler_data.get('mappable_ref2forward_sum', {}).get(ref)
        mappable_reverse = handler_data.get('mappable_ref2reverse_sum', {}).get(ref)
        mappable_ccbins = handler_data.get('mappable_ref2ccbins', {}).get(ref)

        if (mappable_ccbins is not None and
            mappable_forward is not None and
            mappable_reverse is not None):
            mscc_data = MSCCData(
                forward_sum=mappable_forward,
                reverse_sum=mappable_reverse,
                ccbins=mappable_ccbins,
                mappable_len=handler_data.get('ref2mappable_len', {}).get(ref, None)
            )

        if ncc_data or mscc_data:
            result = ChromosomeResult(
                chromosome=ref,
                length=ref_to_length[ref],
                ncc_data=ncc_data,
                mscc_data=mscc_data
            )
            container.add_chromosome_result(result)

    return container
