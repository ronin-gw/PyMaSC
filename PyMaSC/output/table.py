"""Tab-delimited table I/O for cross-correlation and read count data.

Provides I/O functionality for PyMaSC data tables including cross-correlation,
MSCC, and read count tables. Uses tab-delimited format with standardized
column structures for integration with external analysis tools.

Key components:
- TableIO: Generic tab-delimited table I/O handler
- NReadsIO: Specialized I/O for read count tables with forward/reverse pairs
- Table loading functions for cross-correlation and MSCC data
- Table output functions for all result types

Handles both reading and writing of tables with error handling
and logging for debugging and validation.
"""
from __future__ import annotations

import logging
import os
import csv

from collections import defaultdict
from functools import partial
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, TextIO, Tuple, Union, Literal, cast

import sys
if sys.version_info >= (3, 11):
    from typing import Self
else:
    from typing_extensions import Self

import numpy as np

from PyMaSC.interfaces.stats import GenomeWideStats, WholeGenomeStats, ChromosomeStats
from PyMaSC.utils.output import catch_IOError
from typing import TypeVar


CCOUTPUT_SUFFIX = "_cc.tab"
MSCCOUTPUT_SUFFIX = "_mscc.tab"
NREADOUTPUT_SUFFIX = "_nreads.tab"

logger = logging.getLogger(__name__)


class TableIOBase:
    path: os.PathLike[str]
    mode: str
    fp: Optional[TextIO]

    def open(self) -> None:
        """Open file for table operations."""
        self.fp = cast(TextIO, open(self.path, self.mode, newline=''))

    def __enter__(self) -> Self:
        """Context manager entry.

        Returns:
            Self for use in with statements
        """
        self.open()
        return self

    def __exit__(self, _ex_type: Optional[type], _ex_value: Optional[Exception], _trace: Optional[Any]) -> None:
        """Context manager exit with file cleanup.

        Args:
            ex_type: Exception type
            ex_value: Exception value
            trace: Exception traceback
        """
        if self.fp is not None:
            self.fp.close()

    close = __exit__

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


class TableIO(TableIOBase):
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

    T = TypeVar('T')

    def read(self, typefunc: Callable[[str], T]) -> Dict[str, List[T]]:
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


load_cc: Callable[[os.PathLike[str]], Dict[str, List[float]]] = partial(_load_table, logfmt="Load CC table from '{}'")
load_masc: Callable[[os.PathLike[str]], Dict[str, List[float]]] = partial(_load_table, logfmt="Load MSCC table from '{}'")


@catch_IOError(logger)
def _output_cctable(
    outfile: os.PathLike[str],
    stats: GenomeWideStats,
    suffix: str,
    target_attr: Union[Literal['ncc'], Literal['mscc']]
) -> None:
    """Output cross-correlation tables using new StatisticsResult interface.

    Args:
        outfile: Base output file path
        stats_result: StatisticsResult object containing data
        suffix: File suffix to append
        target_attr: Attribute name ('cc' or 'mscc') to extract
    """
    outfile_path = Path(outfile)
    # Create: /path/to/file.ext -> /path/to/file_cc.tab or /path/to/file_mscc.tab
    stem_with_suffix = outfile_path.stem + suffix.replace('.tab', '')
    outfile_with_suffix = outfile_path.parent / f"{stem_with_suffix}.tab"
    logger.info("Output '{}'".format(outfile_with_suffix))

    # Extract whole-genome correlation
    whole_cc: WholeGenomeStats = getattr(stats, f'whole_{target_attr}_stats')
    chrom_stats: Dict[str, ChromosomeStats] = getattr(stats, f'{target_attr}_stats')

    assert whole_cc is not None
    assert chrom_stats is not None

    # Extract correlations
    cc = whole_cc.cc
    ref2cc = {chrom: stats.cc for chrom, stats in chrom_stats.items()
              if not np.isnan(stats.cc).all()}
    keys = sorted(ref2cc.keys())

    with TableIO(outfile_with_suffix, 'w') as tab:
        tab.write(["whole", ] + keys, ([cc] + [ref2cc[k][i] for k in keys] for i, cc in enumerate(cc)))


output_cc: Callable[[os.PathLike[str], GenomeWideStats], None] = partial(_output_cctable, suffix=CCOUTPUT_SUFFIX, target_attr="ncc")
output_mscc: Callable[[os.PathLike[str], GenomeWideStats], None] = partial(_output_cctable, suffix=MSCCOUTPUT_SUFFIX, target_attr="mscc")


class NReadsIO(TableIOBase):
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

    def write(
        self,
        header: List[str],
        forward_sum: Optional[Dict[str, int]],
        reverse_sum: Optional[Dict[str, int]],
        mappable_forward_sum: Optional[Dict[str, List[int]]],
        mappable_reverse_sum: Optional[Dict[str, List[int]]]
    ) -> None:
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
def output_nreads_table(outfile: os.PathLike[str], stats: GenomeWideStats) -> None:
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

    #
    def _extract_read_counts(whole_stats: WholeGenomeStats, chromstats: Dict[str, ChromosomeStats]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
        forward_sum = {}
        forward_sum['whole'] = whole_stats.stats.forward_reads
        forward_sum.update(**{chrom: stat.stats.forward_reads for chrom, stat in chromstats.items()})

        reverse_sum = {}
        reverse_sum['whole'] = whole_stats.stats.reverse_reads
        reverse_sum.update(**{chrom: stat.stats.reverse_reads for chrom, stat in chromstats.items()})

        return forward_sum, reverse_sum

    forward_sum = reverse_sum = None
    if stats.whole_ncc_stats is not None:
        assert stats.ncc_stats is not None
        forward_sum, reverse_sum = _extract_read_counts(stats.whole_ncc_stats, stats.ncc_stats)

    mappable_forward = mappable_reverse = None
    if stats.whole_mscc_stats is not None:
        assert stats.mscc_stats is not None
        mappable_forward, mappable_reverse = _extract_read_counts(stats.whole_mscc_stats, stats.mscc_stats)

    with NReadsIO(outfile_with_suffix, 'w') as tab:
        header = ["whole"] + sorted(stats.references)
        tab.write(header, forward_sum, reverse_sum, mappable_forward, mappable_reverse)
