"""Tab-delimited table I/O for cross-correlation and read count data.

This module provides comprehensive I/O functionality for PyMaSC data tables,
including cross-correlation tables, MSCC tables, and read count tables.
All tables use tab-delimited format with standardized column structures
for easy integration with external analysis tools.

Key components:
- TableIO: Generic tab-delimited table I/O handler
- NReadsIO: Specialized I/O for read count tables with forward/reverse pairs
- Table loading functions for cross-correlation and MSCC data
- Table output functions for all result types

The module handles both reading and writing of tables with proper error
handling and logging for debugging and validation.
"""
from __future__ import annotations, print_function

import logging
from functools import partial
import csv
from collections import defaultdict
from copy import copy
import os
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, TextIO, Tuple, Union

import numpy as np

from PyMaSC.utils.output import catch_IOError

CCOUTPUT_SUFFIX = "_cc.tab"
MSCCOUTPUT_SUFFIX = "_mscc.tab"
NREADOUTPUT_SUFFIX = "_nreads.tab"

logger = logging.getLogger(__name__)


class TableIO(object):
    """Generic tab-delimited table I/O handler.

    Provides basic functionality for reading and writing tab-delimited tables
    with standardized format and error handling. Used as base class for
    specialized table handlers.

    Attributes:
        DIALECT: CSV dialect for tab-delimited format
        path: File path for I/O operations
        mode: File access mode ('r' or 'w')
        fp: File pointer for operations
    """
    DIALECT = "excel-tab"

    def __init__(self, path: os.PathLike[str], mode: str = 'r') -> None:
        """Initialize table I/O handler.

        Args:
            path: File path for table operations
            mode: File access mode ('r' for read, 'w' for write)

        Raises:
            NotImplemented: If mode is not 'r' or 'w'
        """
        self.path = path
        self.mode = mode
        self.fp: Optional[TextIO] = None

        if mode not in ('r', 'w'):
            raise NotImplementedError(f"Unsupported mode: {mode}")

    def open(self) -> None:
        """Open file for table operations."""
        self.fp = open(self.path, self.mode, newline='')

    def __enter__(self) -> 'TableIO':
        """Context manager entry.

        Returns:
            Self for use in with statements
        """
        self.open()
        return self

    def __exit__(self, ex_type: Optional[type], ex_value: Optional[Exception], trace: Optional[Any]) -> None:
        """Context manager exit with file cleanup.

        Args:
            ex_type: Exception type
            ex_value: Exception value
            trace: Exception traceback
        """
        if self.fp is not None:
            self.fp.close()

    close = __exit__

    def read(self, typefunc: Callable[[str], Any]) -> Dict[str, List[Any]]:
        """Read table data with type conversion.

        Args:
            typefunc: Function to convert string values to desired type

        Returns:
            Dictionary mapping column names to converted data lists
        """
        assert self.fp is not None, "File not opened"
        tab = csv.reader(self.fp, dialect=TableIO.DIALECT)
        header = next(tab)[1:]
        table_tuples = dict(zip(header, zip(*(map(typefunc, row[1:]) for row in tab))))
        # Convert tuples to lists for consistency
        table = {k: list(v) for k, v in table_tuples.items()}
        return table

    def write(self, header: List[str], body: Any) -> None:
        """Write table data to file.

        Args:
            header: Column headers (excluding 'shift' column)
            body: Iterable of data rows
        """
        assert self.fp is not None, "File not opened"
        tab = csv.writer(self.fp, dialect=TableIO.DIALECT)
        tab.writerow(("shift", ) + tuple(header))
        tab.writerows((i, ) + tuple(row) for i, row in enumerate(body))


@catch_IOError(logger)
def _load_table(path: os.PathLike[str], logfmt: str) -> Dict[str, List[float]]:
    """Generic table loading with logging.

    Args:
        path: File path to load
        logfmt: Format string for logging message

    Returns:
        Dictionary mapping column names to data arrays

    Note:
        Removes 'whole' column from results and warns if missing
    """
    logger.info(logfmt.format(path))

    with TableIO(path) as tab:
        table = tab.read(float)

    if "whole" in table:
        table.pop("whole")
    else:
        logger.warning("Mandatory column 'whole' not found")

    return table


load_cc = partial(_load_table, logfmt="Load CC table from '{}'")
load_masc = partial(_load_table, logfmt="Load MSCC table from '{}'")


@catch_IOError(logger)
def _output_cctable(outfile: os.PathLike[str], stats_result: Any, suffix: str, target_attr: str) -> None:
    """Output cross-correlation tables using new StatisticsResult interface.

    Args:
        outfile: Base output file path
        stats_result: StatisticsResult object containing data
        suffix: File suffix to append
        target_attr: Attribute name ('cc' or 'masc') to extract
    """
    outfile_path = Path(outfile)
    # Create: /path/to/file.ext -> /path/to/file_cc.tab or /path/to/file_mscc.tab
    stem_with_suffix = outfile_path.stem + suffix.replace('.tab', '')
    outfile_with_suffix = outfile_path.parent / f"{stem_with_suffix}.tab"
    logger.info("Output '{}'".format(outfile_with_suffix))

    def _extract_correlation_data(stats_result, target_attr: str):
        """Extract correlation data from StatisticsResult.

        Args:
            stats_result: StatisticsResult object
            target_attr: 'cc' for NCC or 'masc' for MSCC

        Returns:
            Tuple of (whole_genome_cc, ref2cc_dict)
        """
        # Extract whole-genome correlation
        genome_stats = stats_result.genome_wide_stats
        if target_attr == "cc" and genome_stats.cc is not None:
            whole_cc = genome_stats.cc.cc
        elif target_attr == "masc" and genome_stats.masc is not None:
            whole_cc = genome_stats.masc.cc
        else:
            # Return empty array if no data available
            whole_cc = np.array([])

        # Extract per-chromosome correlations
        ref2cc = {}
        for chrom in stats_result.references:
            chrom_stats = stats_result.chromosome_stats.get(chrom)
            if chrom_stats is None:
                continue

            if target_attr == "cc" and chrom_stats.cc is not None:
                cc_data = chrom_stats.cc.cc
                if cc_data is not None and not np.isnan(cc_data).all():
                    ref2cc[chrom] = cc_data
            elif target_attr == "masc" and chrom_stats.masc is not None:
                cc_data = chrom_stats.masc.cc
                if cc_data is not None and not np.isnan(cc_data).all():
                    ref2cc[chrom] = cc_data

        return whole_cc, ref2cc

    # Extract correlation data using new interface
    cc, ref2cc = _extract_correlation_data(stats_result, target_attr)
    keys = sorted(ref2cc.keys())

    with TableIO(outfile_with_suffix, 'w') as tab:
        tab.write(["whole", ] + keys, ([cc] + [ref2cc[k][i] for k in keys] for i, cc in enumerate(cc)))


output_cc = partial(_output_cctable, suffix=CCOUTPUT_SUFFIX, target_attr="cc")
output_mscc = partial(_output_cctable, suffix=MSCCOUTPUT_SUFFIX, target_attr="masc")


class NReadsIO(TableIO):
    """Specialized I/O for read count tables.

    Handles reading and writing of read count tables that contain
    forward-reverse read pair counts in 'forward-reverse' format.
    Supports both raw read counts and mappable read counts.

    The class processes tables with the following structure:
    - Optional 'raw' row with total read counts
    - Multiple rows with mappable read counts by shift distance
    - Columns represent different chromosomes or genome-wide ('whole')
    """
    @staticmethod
    def make_forward_reverse_tables(header: List[str], body: List[List[str]]) -> Tuple[Dict[str, List[int]], Dict[str, List[int]]]:
        """Parse forward-reverse read count pairs from table rows.

        Args:
            header: Column names
            body: Table rows containing 'forward-reverse' formatted data

        Returns:
            Tuple of (forward_dict, reverse_dict) with counts by column
        """
        forward = defaultdict(list)
        reverse = defaultdict(list)

        for row in body:
            for key, pair in zip(header, row[1:]):
                if isinstance(pair, str) and '-' in pair:
                    f, r = map(int, pair.split('-'))
                    forward[key].append(f)
                    reverse[key].append(r)

        return forward, reverse

    def read(self) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, List[int]], Dict[str, List[int]]]:
        """Read read count table with forward/reverse separation.

        Returns:
            Tuple of (forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum)
            Each is a dictionary mapping chromosome names to read counts
        """
        assert self.fp is not None, "File not opened"
        tab = csv.reader(self.fp, dialect=TableIO.DIALECT)

        header = next(tab)[1:]

        first = next(tab)
        if first[0] == "raw":
            forward_sum_lists, reverse_sum_lists = self.make_forward_reverse_tables(header, [first])
            forward_sum = {k: v[0] for k, v in forward_sum_lists.items()}
            reverse_sum = {k: v[0] for k, v in reverse_sum_lists.items()}
            mscc_iter = list(tab)
        else:
            forward_sum = reverse_sum = {}
            mscc_iter = [first] + list(tab)

        mappable_forward_sum, mappable_reverse_sum = self.make_forward_reverse_tables(header, mscc_iter)
        if not mappable_forward_sum:
            mappable_forward_sum = mappable_reverse_sum = {}

        return forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum

    @staticmethod
    def _make_table(rowname: Union[str, int], forward: Union[List[int], Tuple[Any, ...]], reverse: Union[List[int], Tuple[Any, ...]]) -> List[Union[str, int]]:
        """Create table row with forward-reverse formatted data.

        Args:
            rowname: Name for the first column
            forward: Forward read counts tuple
            reverse: Reverse read counts tuple

        Returns:
            List of row data with 'forward-reverse' formatted pairs
        """
        return [rowname] + ["{}-{}".format(f, r) for f, r in zip(forward, reverse)]

    def write(self, header: List[str], forward_sum: Optional[Dict[str, int]], reverse_sum: Optional[Dict[str, int]], mappable_forward_sum: Optional[Dict[str, List[int]]], mappable_reverse_sum: Optional[Dict[str, List[int]]]) -> None:
        """Write read count table with forward/reverse data.

        Args:
            header: Column names
            forward_sum: Total forward read counts by chromosome
            reverse_sum: Total reverse read counts by chromosome
            mappable_forward_sum: Mappable forward read counts by shift and chromosome
            mappable_reverse_sum: Mappable reverse read counts by shift and chromosome
        """
        # Handle None or empty mappable dictionaries gracefully
        if mappable_forward_sum:
            mappable_forward_sum = {k: v for k, v in mappable_forward_sum.items() if v is not None}
        else:
            mappable_forward_sum = {}

        if mappable_reverse_sum:
            mappable_reverse_sum = {k: v for k, v in mappable_reverse_sum.items() if v is not None}
        else:
            mappable_reverse_sum = {}

        assert self.fp is not None, "File not opened"
        tab = csv.writer(self.fp, dialect=TableIO.DIALECT)
        tab.writerow(("shift", ) + tuple(header))

        if forward_sum and reverse_sum:
            tab.writerow(self._make_table("raw", [forward_sum.get(col, 0) for col in header],
                                                 [reverse_sum.get(col, 0) for col in header]))

        if mappable_forward_sum and mappable_reverse_sum:
            shiftsize = len(mappable_forward_sum["whole"])
            for i, (forward, reverse) in enumerate(zip(zip(*[mappable_forward_sum.get(col, [0] * shiftsize) for col in header]),
                                                       zip(*[mappable_reverse_sum.get(col, [0] * shiftsize) for col in header]))):
                tab.writerow(self._make_table(i, forward, reverse))


@catch_IOError(logger)
def load_nreads_table(path: os.PathLike[str]) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, List[int]], Dict[str, List[int]]]:
    """Load read count table from file.

    Args:
        path: File path to read count table

    Returns:
        Tuple of read count dictionaries (forward, reverse, mappable_forward, mappable_reverse)

    Raises:
        KeyError: If no valid data could be loaded
    """
    logger.info("Load Nreads table from '{}'".format(path))

    with NReadsIO(path) as tab:
        result = tab.read()
        forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum = result

    for d in [forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum]:
        if isinstance(d, dict) and "whole" in d:
            d.pop("whole")
        elif d:
            logger.warning("Mandatory column 'whole' not found")

    if all(not d for d in [forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum]):
        logger.critical("Nothing to load.")
        raise KeyError

    return forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum


@catch_IOError(logger)
def output_nreads_table(outfile: os.PathLike[str], stats_result: Any) -> None:
    """Output read count table from analysis results using new StatisticsResult interface.

    Args:
        outfile: Base output file path
        stats_result: StatisticsResult object containing read count data
    """
    outfile_path = Path(outfile)
    # Create: /path/to/file.ext -> /path/to/file_nreads.tab
    stem_with_suffix = outfile_path.stem + NREADOUTPUT_SUFFIX.replace('.tab', '')
    outfile_with_suffix = outfile_path.parent / f"{stem_with_suffix}.tab"
    logger.info("Output '{}'".format(outfile_with_suffix))

    def _extract_read_counts(stats_result, count_type: str, is_mappable: bool = False):
        """Extract read counts from StatisticsResult.

        Args:
            stats_result: StatisticsResult object
            count_type: 'forward_sum' or 'reverse_sum'
            is_mappable: Whether to extract mappable read counts

        Returns:
            Dictionary with chromosome counts and 'whole' genome count.
            For non-mappable: values are integers
            For mappable: values are lists (arrays) indexed by shift distance
        """
        result_dict = {}

        if is_mappable:
            # For mappable counts, we need arrays indexed by shift distance
            # First determine the max shift to create consistent array sizes
            max_shift = 300  # Default fallback
            genome_stats = stats_result.genome_wide_stats
            if genome_stats and hasattr(genome_stats, 'max_shift'):
                max_shift = genome_stats.max_shift
            elif genome_stats and genome_stats.masc and hasattr(genome_stats.masc, 'cc') and genome_stats.masc.cc is not None:
                max_shift = len(genome_stats.masc.cc) - 1

            array_size = max_shift + 1

            # Extract per-chromosome mappable counts as arrays
            for chrom in stats_result.references:
                chrom_stats = stats_result.chromosome_stats.get(chrom)
                if chrom_stats is None or chrom_stats.masc is None:
                    result_dict[chrom] = [0] * array_size
                    continue

                count_value = getattr(chrom_stats.masc, count_type, None)
                if count_value is not None and hasattr(count_value, '__iter__') and not isinstance(count_value, str):
                    # Convert to list and ensure proper length
                    count_list = list(count_value)
                    if len(count_list) == array_size:
                        result_dict[chrom] = count_list
                    elif len(count_list) < array_size:
                        # Pad with zeros if too short
                        result_dict[chrom] = count_list + [0] * (array_size - len(count_list))
                    else:
                        # Truncate if too long
                        result_dict[chrom] = count_list[:array_size]
                else:
                    # Single value or None - create array with that value at position 0
                    value = int(count_value) if count_value is not None else 0
                    result_dict[chrom] = [value] + [0] * (array_size - 1)

            # Extract whole-genome mappable count as array
            if genome_stats and genome_stats.masc is not None:
                whole_count = getattr(genome_stats.masc, count_type, None)
                if whole_count is not None and hasattr(whole_count, '__iter__') and not isinstance(whole_count, str):
                    count_list = list(whole_count)
                    if len(count_list) == array_size:
                        result_dict["whole"] = count_list
                    elif len(count_list) < array_size:
                        result_dict["whole"] = count_list + [0] * (array_size - len(count_list))
                    else:
                        result_dict["whole"] = count_list[:array_size]
                else:
                    value = int(whole_count) if whole_count is not None else 0
                    result_dict["whole"] = [value] + [0] * (array_size - 1)
            else:
                result_dict["whole"] = [0] * array_size

        else:
            # For non-mappable counts, return single integers
            # Extract per-chromosome counts
            for chrom in stats_result.references:
                chrom_stats = stats_result.chromosome_stats.get(chrom)
                if chrom_stats is None or chrom_stats.cc is None:
                    result_dict[chrom] = 0
                    continue

                count_value = getattr(chrom_stats.cc, count_type, 0)
                result_dict[chrom] = int(count_value) if count_value is not None else 0

            # Extract whole-genome count
            genome_stats = stats_result.genome_wide_stats
            if genome_stats and genome_stats.cc is not None:
                whole_count = getattr(genome_stats.cc, count_type, 0)
                result_dict["whole"] = int(whole_count) if whole_count is not None else 0
            else:
                result_dict["whole"] = 0

        return result_dict

    with NReadsIO(outfile_with_suffix, 'w') as tab:
        header = ["whole"] + sorted(stats_result.references)

        # Extract read counts using new interface
        forward_sum = _extract_read_counts(stats_result, "forward_sum", is_mappable=False)
        reverse_sum = _extract_read_counts(stats_result, "reverse_sum", is_mappable=False)

        # Extract mappable read counts if available
        if stats_result.has_mappability:
            mappable_forward = _extract_read_counts(stats_result, "forward_sum", is_mappable=True)
            mappable_reverse = _extract_read_counts(stats_result, "reverse_sum", is_mappable=True)
        else:
            mappable_forward = None
            mappable_reverse = None

        # NReadsIO has its own write() method with different signature
        tab.write(header, forward_sum, reverse_sum, mappable_forward, mappable_reverse)
