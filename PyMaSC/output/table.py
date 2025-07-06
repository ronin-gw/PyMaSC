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
from __future__ import print_function

import logging
from functools import partial
import csv
from itertools import chain
from collections import defaultdict
from copy import copy
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

    def __init__(self, path: str, mode: str = 'r') -> None:
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
        self.fp = open(self.path, self.mode, newline='')  # type: ignore

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
def _load_table(path: str, logfmt: str) -> Dict[str, List[float]]:
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
def _output_cctable(outfile: str, ccr: Any, suffix: str, target_attr: str) -> None:
    """Output cross-correlation tables.

    Args:
        outfile: Base output file path
        ccr: CCResult object containing data
        suffix: File suffix to append
        target_attr: Attribute name ('cc' or 'masc') to extract
    """
    outfile += suffix
    logger.info("Output '{}'".format(outfile))

    cc = getattr(ccr.whole, target_attr).cc
    ref2cc = {k: getattr(ccr.ref2stats[k], target_attr) for k in ccr.references}
    ref2cc = {k: v.cc for k, v in ref2cc.items() if v is not None and not np.isnan(v.cc).all()}
    keys = sorted(ref2cc.keys())

    with TableIO(outfile, 'w') as tab:
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

    def read(self) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, List[int]], Dict[str, List[int]]]:  # type: ignore[override]
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

    def write(self, header: List[str], forward_sum: Optional[Dict[str, int]], reverse_sum: Optional[Dict[str, int]], mappable_forward_sum: Optional[Dict[str, List[int]]], mappable_reverse_sum: Optional[Dict[str, List[int]]]) -> None:  # type: ignore[override]
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
def load_nreads_table(path: str) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, List[int]], Dict[str, List[int]]]:
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
        result = tab.read()  # type: ignore[call-arg]
        forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum = result

    for d in [forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum]:
        if isinstance(d, dict) and "whole" in d:
            d.pop("whole")
        elif d:
            logger.warning("Mandatory column 'whole' not found")

    if all(not d for d in [forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum]):
        logger.critical("Nothing to load.")
        raise KeyError

    return forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum  # type: ignore[return-value]


@catch_IOError(logger)
def output_nreads_table(outfile: str, ccr: Any) -> None:
    """Output read count table from analysis results.

    Args:
        outfile: Base output file path
        ccr: CCResult object containing read count data
    """
    outfile += NREADOUTPUT_SUFFIX
    logger.info("Output '{}'".format(outfile))

    def _make_dict_with_whole(add_target, target1, target2):
        if not add_target:
            return None
        d = copy(add_target)
        target_obj = getattr(ccr.whole, target1)
        if target_obj is None:
            return None
        d["whole"] = getattr(target_obj, target2)
        return d

    with NReadsIO(outfile, 'w') as tab:
        header = ["whole"] + sorted(ccr.references)
        forward_sum = _make_dict_with_whole(ccr.ref2forward_sum, "cc", "forward_sum")
        reverse_sum = _make_dict_with_whole(ccr.ref2reverse_sum, "cc", "reverse_sum")
        mappable_forward = _make_dict_with_whole(ccr.mappable_ref2forward_sum, "masc", "forward_sum")
        mappable_reverse = _make_dict_with_whole(ccr.mappable_ref2reverse_sum, "masc", "reverse_sum")
        # NReadsIO has its own write() method with different signature
        tab.write(header, forward_sum, reverse_sum, mappable_forward, mappable_reverse)  # type: ignore[call-arg]
