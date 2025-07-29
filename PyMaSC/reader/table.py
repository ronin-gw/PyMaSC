"""Table reading utilities for PyMaSC.

This module provides functions to read various PyMaSC table formats including
NCC (Normalized Cross-Correlation) tables, MSCC (Mappability-Sensitive Cross-Correlation)
tables, and read count tables.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict
from os import PathLike

import numpy as np
import numpy.typing as npt

from PyMaSC.output.table import TableIO, NReadsIO

logger = logging.getLogger(__name__)


@dataclass
class CCData:
    """Intermediate representation for table data.

    Holds raw data loaded from various table formats before
    conversion to final statistics objects.
    """
    whole: npt.NDArray[np.float64]
    chrom: Dict[str, npt.NDArray[np.float64]]


@dataclass
class NreadData:
    forward_sum: Dict[str, int]
    reverse_sum: Dict[str, int]
    mappable_forward: Dict[str, npt.NDArray[np.int64]]
    mappable_reverse: Dict[str, npt.NDArray[np.int64]]


def load_cc_table(path: PathLike[str]) -> CCData:
    """Load CC table from a tab-delimited file.

    Reads a CC table in the format produced by PyMaSC, with columns for
    'shift', 'whole' (genome-wide), and individual chromosomes.

    Args:
        path: Path to the CC table file

    Returns:
        Tuple of (whole_genome_cc, chromosome_cc_dict) where:
        - whole_genome_cc: numpy array of genome-wide cross-correlation values
        - chromosome_cc_dict: dictionary mapping chromosome names to their CC arrays

    Raises:
        IOError: If file cannot be read
        ValueError: If required 'whole' column is missing

    Example:
        >>> whole_cc, chrom_cc = load_ncc_table('results_cc.tab')
        >>> print(f"Max correlation: {whole_cc.max()}")
        >>> print(f"Chromosomes: {list(chrom_cc.keys())}")
    """
    with TableIO(path, 'r') as tab:
        table = tab.read(float)

    # Extract and validate whole genome column
    if "whole" not in table:
        raise ValueError("Required column 'whole' not found in NCC table")

    whole_cc = np.array(table.pop("whole"), dtype=np.float64)

    # Convert remaining columns to numpy arrays
    chrom_cc = {chrom: np.array(values, dtype=np.float64) for chrom, values in table.items()}

    return CCData(whole=whole_cc, chrom=chrom_cc)


def load_nreads_table(path: PathLike[str]) -> NreadData:
    """Load read count table from a tab-delimited file.

    Reads an nreads table containing forward and reverse read counts
    for each chromosome and shift position. The table may contain:
    - A 'raw' row with total read counts
    - Multiple rows with mappable read counts by shift distance

    Args:
        path: Path to the nreads table file

    Returns:
        Tuple of (forward_sum, reverse_sum, mappable_forward, mappable_reverse):
        - forward_sum: Total forward reads by chromosome
        - reverse_sum: Total reverse reads by chromosome
        - mappable_forward: Forward reads by shift and chromosome
        - mappable_reverse: Reverse reads by shift and chromosome

    Example:
        >>> fw_sum, rv_sum, fw_map, rv_map = load_nreads_table('results_nreads.tab')
        >>> print(f"Total forward reads: {fw_sum.get('whole', 0)}")
        >>> print(f"Chr1 mappable forward at shift 0: {fw_map['chr1'][0]}")
    """
    logger.info(f"Loading nreads table from '{path}'")

    with NReadsIO(path, 'r') as tab:
        forward_sum, reverse_sum, mappable_forward, mappable_reverse = tab.read()

    # Convert lists to numpy arrays for mappable counts
    mappable_forward_arrays = {
        chrom: np.array(counts, dtype=np.int64)
        for chrom, counts in mappable_forward.items()
    }
    mappable_reverse_arrays = {
        chrom: np.array(counts, dtype=np.int64)
        for chrom, counts in mappable_reverse.items()
    }

    return NreadData(
        forward_sum=forward_sum,
        reverse_sum=reverse_sum,
        mappable_forward=mappable_forward_arrays,
        mappable_reverse=mappable_reverse_arrays
    )
