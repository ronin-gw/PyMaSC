"""Command-line argument parsing utilities for PyMaSC entry points.

Provides argument parsing for PyMaSC's command-line utilities:
pymasc, pymasc-plot, and pymasc-precalc.

Key functionality:
- Common argument groups shared across utilities
- Custom argparse actions for validation and type conversion
- Parser factories for each command-line tool
- Input validation and error handling
"""
import argparse
import logging
import sys
from pathlib import Path
from typing import Any, Optional, Sequence, Type, Union

import PyMaSC
from PyMaSC.interfaces.config import EstimationType

READLEN_ESTIMATION_TYPES = tuple(e.value for e in EstimationType)
EPILOG = (" \nVisit PyMaSC web site for more information and to get human genome "
          "mappability tracks\n" + PyMaSC.WEBSITE_URL + '\n ')

NEAR_READLEN_ERR_CRITERION = 5


def _make_upper(s: str) -> str:
    """Convert string to uppercase.

    Args:
        s: String to convert

    Returns:
        Uppercase version of the input string
    """
    return s.upper()


class StoreLoggingLevel(argparse.Action):
    """Custom argparse action for logging level arguments.

    Converts string logging level names (e.g., 'INFO', 'DEBUG') to their
    corresponding logging module constants for use in log level configuration.
    """
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: Union[Any, Sequence[Any], None],
        option_string: Optional[str] = None
    ) -> None:
        assert isinstance(values, str), "Logging level must be a string"
        setattr(namespace, self.dest, getattr(logging, values))


class ForceNaturalNumber(argparse.Action):
    """Custom argparse action for positive integer validation.

    Ensures that integer arguments are positive (> 0) and provides
    clear error messages for invalid values.
    """
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: Union[Any, Sequence[Any], None],
        option_string: Optional[str] = None
    ) -> None:
        assert isinstance(values, int), "Argument must be an integer"
        if values < 1:
            parser.error("argument {} must be > 0.".format('/'.join(self.option_strings)))
        setattr(namespace, self.dest, values)


class ToColorizeOption(argparse.Action):
    """Custom argparse action for colorization options.

    Converts 'TRUE'/'FALSE' strings to boolean values for colorization
    preferences in logging output.
    """
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: Union[Any, Sequence[Any], None],
        option_string: Optional[str] = None
    ) -> None:
        assert isinstance(values, str), "Colorization option must be a string"
        if values == "TRUE":
            colorize = True
        elif values == "FALSE":
            colorize = False
        else:
            colorize = sys.stderr.isatty()
        setattr(namespace, self.dest, colorize)


def make_multistate_append_action(key: bool) -> Type[argparse.Action]:
    """Create multi-state argument action factory.

    Creates a custom argparse action that appends (key, value) tuples
    to a list, allowing multiple related arguments to be collected.

    Args:
        key: Key to associate with values in the tuple list

    Returns:
        Custom argparse Action class for multi-state arguments
    """
    class _MultistateAppendAction(argparse.Action):
        def __call__(self, parser: argparse.ArgumentParser, namespace: argparse.Namespace, values: Any, option_string: Optional[str] = None) -> None:
            args = getattr(namespace, self.dest)
            args = [] if args is None else args
            args.append((key, values))
            setattr(namespace, self.dest, args)

    return _MultistateAppendAction


def add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add common arguments shared across all PyMaSC utilities.

    Adds logging, progress control, and color output arguments that
    are used by all PyMaSC command-line tools.

    Args:
        parser: ArgumentParser instance to add arguments to
    """
    parser.add_argument(
        "-v", "--log-level", type=_make_upper, default=logging.INFO,
        action=StoreLoggingLevel, choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        help="Set verbosity. (Default: INFO)"
    )
    parser.add_argument(
        "--disable-progress", action="store_true",
        help="Disable progress bar"
    )
    parser.add_argument(
        "--color", type=_make_upper, default=True, action=ToColorizeOption, choices=("TRUE", "FALSE"),
        help="Coloring log. (Default: auto)"
    )
    parser.add_argument(
        "--version", action="version", version="PyMaSC " + PyMaSC.VERSION
    )


def add_multiprocess_args(group: argparse._ArgumentGroup) -> None:
    """Add multiprocessing-related arguments.

    Args:
        group: Argument group to add multiprocessing options to
    """
    group.add_argument(
        "-p", "--process", type=int, default=1, action=ForceNaturalNumber,
        help="Number of worker process. (Default: 1)"
    )


def add_mappability_args(group: argparse._ArgumentGroup) -> None:
    """Add mappability file and statistics arguments.

    Adds arguments for specifying BigWig mappability files and
    pre-calculated mappability statistics for MSCC analysis.

    Args:
        group: Argument group to add mappability options to
    """
    group.add_argument(
        "-m", "--mappability", metavar="REGION_FILE", type=Path,
        help="BigWig format mappable region file."
    )
    group.add_argument(
        "--mappability-stats", type=Path,
        help="Read/Save path for mappability stats. "
             "(Default: [REGION_FILE]_mappability.json)"
    )


def add_shift_arg(group: argparse._ArgumentGroup) -> None:
    """Add maximum shift distance argument.

    Args:
        group: Argument group to add shift argument to
    """
    group.add_argument(
        "-d", "--max-shift", type=int, action=ForceNaturalNumber, default=1000,
        help="PyMaSC calculate CC with reverse strand shift from 1 to [MAX_SHIFT] bases. (Default: 1000)"
    )


def add_liblen_arg(group: argparse._ArgumentGroup) -> None:
    """Add library length argument for expected fragment size.

    Args:
        group: Argument group to add library length option to
    """
    group.add_argument(
        "-l", "--library-length", type=int, action=ForceNaturalNumber,
        help="Your expected library length for input sample(s)."
    )


def add_chrom_filter_args(group: argparse._ArgumentGroup) -> None:
    """Add chromosome filtering arguments.

    Adds include and exclude chromosome pattern arguments for
    selective analysis of chromosome subsets.

    Args:
        group: Argument group to add chromosome filtering options to
    """
    group.add_argument(
        "-i", "--include-chrom", nargs='+', dest="chromfilter", metavar="CHROM",
        action=make_multistate_append_action(True),
        help="Include chromosomes to calculate. You can use Unix shell-style "
             "wildcards ('.', '*', '[]' and '[!]'). This option can be declared "
             "multiple times to include chromosomes specified in a just before "
             "-e/--exclude-chrom option. Note that this option is case-sensitive."
    )
    group.add_argument(
        "-e", "--exclude-chrom", nargs='+', dest="chromfilter", metavar="CHROM",
        action=make_multistate_append_action(False),
        help="Exclude chromosomes from calculation. You can use Unix shell-style "
             "wildcards ('.', '*', '[]' and '[!]'). This option can be declared "
             "multiple times to exclude chromosomes specified in a just before "
             "-i/--include-chrom option. Note that this option is case-sensitive."
    )


def add_result_proc_args(group: argparse._ArgumentGroup) -> None:
    """Add result processing arguments.

    Adds arguments for controlling result processing including smoothing,
    statistical testing, masking, and fragment length estimation.

    Args:
        group: Argument group to add result processing options to
    """
    group.add_argument(
        "--chi2-pval", type=float, default=0.05,
        help="p-value threshold for Chi-squared test to check strand specificity "
             "using number of reads mapped to each strands. (Default: 0.05)"
    )
    group.add_argument(
        "-w", "--smooth-window", type=int, default=15, action=ForceNaturalNumber,
        help="Moving average window size for smoothing MSCC "
             "to estimate library length. (Default: 15)"
    )
    group.add_argument(
        "--mask-size", type=int, default=NEAR_READLEN_ERR_CRITERION,
        help="If difference between a read length and the estimated library length "
             "is equal or less than the length specified by this option, "
             "PyMaSC masks correlation coefficients in the read length +/- specified length "
             "and try to estimate mean library length again. (Default: {}, Specify < 1 to disable)".format(NEAR_READLEN_ERR_CRITERION)
    )
    group.add_argument(
        "--bg-avr-width", type=int, action=ForceNaturalNumber, default=50,
        help="The minimum of coefficient will be calcurated as the median of the end of specified bases. (Default: 50bp)"
    )


def get_pymasc_parser() -> argparse.ArgumentParser:
    """Create main PyMaSC argument parser.

    Constructs the complete argument parser for the main pymasc command,
    including all argument groups and options for cross-correlation analysis.

    Returns:
        ArgumentParser configured for main PyMaSC analysis
    """
    parser = argparse.ArgumentParser(
        description="Estimation and visualization tool for library length, "
                    "NSC and RSC metrics with\nmappability sensitive cross-correlation calculation.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    add_common_args(parser)

    proc_args = parser.add_argument_group("Processing behaviors")
    add_multiprocess_args(proc_args)
    proc_args.add_argument(
        "--successive", action="store_true",
        help="Calc with successive algorithm instead of bit array implementation"
    )
    proc_args.add_argument(
        "--skip-ncc", action="store_true",
        help="Skip naive cross-correlation calculation. Mappability file must be specified."
    )
    proc_args.add_argument(
        "--skip-plots", action="store_true",
        help="Skip output figures."
    )

    input_args = parser.add_argument_group("Input alignment file arguments")
    input_args.add_argument(
        "reads", nargs="+", type=Path,
        help="SAM/BAM format mapped reads. Input must be sorted by positions."
    )
    input_args.add_argument(
        "-r", "--read-length", type=int, action=ForceNaturalNumber,
        help="Set read length manually and disable read length estimation."
    )
    input_args.add_argument(
        "--readlen-estimator", type=_make_upper,
        default="MEDIAN", choices=READLEN_ESTIMATION_TYPES,
        help="Select a representative value used to estimate a read length from "
             "observed read lengths. "
             "Choices: mean, median, mode, min, max (Default: median)"
    )
    add_liblen_arg(input_args)

    map_args = parser.add_argument_group("Input mappability file arguments")
    add_mappability_args(map_args)

    filter = parser.add_argument_group("Input file filtering arguments")
    filter.add_argument(
        "-q", "--mapq", type=int, default=1,
        help="Filter out reads which have less than specified "
             "SAM mapping quality score. (Default: 1)"
    )
    add_chrom_filter_args(filter)

    proc_params = parser.add_argument_group("PyMaSC parameters")
    add_shift_arg(proc_params)
    add_result_proc_args(proc_params)

    output = parser.add_argument_group("Output file arguments")
    output.add_argument(
        "-n", "--name", nargs='*', default=[],
        help="Output file base name(s). (Default: input file name without extension)"
    )
    output.add_argument(
        "-o", "--outdir", default='.', type=Path,
        help="Output directory. (Default: current directory)"
    )

    return parser


def get_precalc_parser() -> argparse.ArgumentParser:
    """Create pre-calculation argument parser.

    Constructs the argument parser for the pymasc-precalc command,
    which pre-calculates mappability statistics from BigWig files.

    Returns:
        ArgumentParser configured for mappability pre-calculation
    """
    parser = argparse.ArgumentParser(
        description="Pre-calculate mappability region statistics for "
                    "PyMaSC successive algorithm.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    add_common_args(parser)

    proc_args = parser.add_argument_group("Processing behaviors")
    add_multiprocess_args(proc_args)

    map_args = parser.add_argument_group("Input mappability file arguments")
    add_mappability_args(map_args)

    proc_params = parser.add_argument_group("PyMaSC parameters")
    add_shift_arg(proc_params)
    proc_params.add_argument(
        "-r", "--max-readlen", type=int, action=ForceNaturalNumber, default=1000,
        help="Set max read length to calculate mappable region length."
    )

    return parser


def get_plot_parser() -> argparse.ArgumentParser:
    """Create plotting argument parser.

    Constructs the argument parser for the pymasc-plot command,
    which generates plots from pre-calculated PyMaSC results.

    Returns:
        ArgumentParser configured for plot generation
    """
    parser = argparse.ArgumentParser(
        description="Plot figures from PyMaSC statistic outputs.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    add_common_args(parser)

    input_args = parser.add_argument_group("Input alignment file arguments")
    input_args.add_argument(
        "statfile", nargs='?', type=Path,
        help="A base path to the statistic files (*_stats.tab, *_cc.tab and *_masc.tab) "
             "to plot figures."
    )
    input_args.add_argument(
        "--stats", type=Path,
        help="For specifying path to a statistic file (*_stats.tab) separately."
    )
    input_args.add_argument(
        "--cc", type=Path,
        help="For specifying path to a cross-correlation table file (*_cc.tab) separately."
    )
    input_args.add_argument(
        "--masc", type=Path,
        help="For specifying path to a mappability sensitive cross-correlation file (*_stats.tab) separately."
    )
    input_args.add_argument(
        "--nreads", type=Path,
        help="For specifying path to a # of reads file (*_nreads.tab) separately."
    )
    input_args.add_argument(
        "-s", "--sizes", type=Path,
        help="A file to obtain length of chromosomes. Tab delimited files,"
             "like *.chrom.sizes and *.fai files, or SAM/BAM format files are acceptable."
    )
    input_args.add_argument(
        "-m", "--mappability-stats", type=Path,
        help="A JSON file to obtain mappable length of chromosomes generated by PyMaSC"
             "for a BigWig file."
    )

    filter = parser.add_argument_group("Chromosome filtering arguments")
    add_chrom_filter_args(filter)

    proc_params = parser.add_argument_group("PyMaSC parameters")
    add_result_proc_args(proc_params)
    add_liblen_arg(proc_params)

    output = parser.add_argument_group("Output file arguments")
    output.add_argument(
        "-n", "--name",
        help="Change output file base name. (Default: same as name field in input)"
    )
    output.add_argument(
        "-o", "--outdir", default='.', type=Path,
        help="Output directory. (Default: current directory)"
    )
    output.add_argument(
        "-f", "--force-overwrite", nargs="*", type=str.lower, choices=('all', 'stats', 'cc', 'mscc'), default=[],
        help="Overwrite specified files even if input and output path are same. "
             "(choices: 'all', 'stats', 'cc', 'mscc' (multiple choices are acceptable))"
    )

    return parser
