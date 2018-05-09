import argparse
import logging

import PyMaSC

READLEN_ESTIMATION_TYPES = ("MEAN", "MEDIAN", "MODE", "MIN", "MAX")


def _make_upper(s):
    return s.upper()


class StoreLoggingLevel(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, getattr(logging, values))


class ForceNaturalNumber(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values < 1:
            parser.error("argument {} must be > 0.".format('/'.join(self.option_strings)))
        setattr(namespace, self.dest, values)


def add_common_args(parser):
    parser.add_argument(
        "-v", "--log-level", nargs='?', type=_make_upper, default=logging.INFO,
        action=StoreLoggingLevel, choices=("DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"),
        help="Set verbosity. (Default: INFO)"
    )
    parser.add_argument(
        "--disable-progress", action="store_true",
        help="Disable progress bar"
    )
    parser.add_argument(
        "--color", nargs='?', type=_make_upper, default=None, choices=("TRUE", "FALSE"),
        help="Coloring log. (Default: auto)"
    )
    parser.add_argument(
        "--version", action="version", version="PyMaSC " + PyMaSC.VERSION
    )


def add_multiprocess_args(group):
    group.add_argument(
        "-p", "--process", nargs='?', type=int, default=1, action=ForceNaturalNumber,
        help="Number of worker process. (Default: 1)"
    )


def add_mappability_args(group):
    group.add_argument(
        "-m", "--mappability", nargs=1, metavar="REGION_FILE",
        help="BigWig format mappable region file."
    )
    group.add_argument(
        "--mappability-stats", nargs='?',
        help="Read/Save path for mappability stats. "
             "(Default: [REGION_FILE]_mappability.json)"
    )


def add_shift_arg(group):
    group.add_argument(
        "-d", "--max-shift", nargs='?', type=int, action=ForceNaturalNumber, default=1000,
        help="PyMaSC calculate CC with reverse strand shift from 1 to [MAX_SHIFT] bases. (Default: 1000)"
    )


def add_liblen_arg(group):
    group.add_argument(
        "-l", "--library-length", nargs='?', type=int, action=ForceNaturalNumber,
        help="Expected library length for input sample(s)."
    )


def add_result_proc_args(group):
    group.add_argument(
        "--chi2-pval", nargs='?', type=float, default=0.05,
        help="p-value threshold for Chi-squared test to check strand specificity "
             "using number of reads mapped to each strands. (Default: 0.05)"
    )
    group.add_argument(
        "-w", "--smooth-window", nargs='?', type=int, default=15, action=ForceNaturalNumber,
        help="Moving average window size for smoothing MSCC "
             "to estimate library length. (Default: 15)"
    )


def get_pymasc_parser():
    parser = argparse.ArgumentParser(
        description="Estimation and visualization tool for library length, "
                    "NSC and RSC metrics with mappability sensitive cross-correlation calculation."
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
        "reads", nargs="+",
        help="SAM/BAM format mapped reads. Input must be sorted by positions."
    )
    input_args.add_argument(
        "-r", "--read-length", nargs='?', type=int, action=ForceNaturalNumber,
        help="Set read length manually and disable read length estimation."
    )
    input_args.add_argument(
        "--readlen-estimator", nargs='?', type=_make_upper,
        default="MEDIAN", choices=READLEN_ESTIMATION_TYPES,
        help="Define the statistic used to estimate a read length from observed "
             "read lengths. Choices: mean, median, mode, min, max (Default: median)"
    )
    input_args.add_argument(
        "-q", "--mapq", nargs='?', type=int, default=1,
        help="Filter out reads which have less than specified "
             "SAM mapping quality score. (Default: 1)"
    )
    add_liblen_arg(input_args)

    map_args = parser.add_argument_group("Input mappability file arguments")
    add_mappability_args(map_args)

    proc_params = parser.add_argument_group("PyMaSC parameters")
    add_shift_arg(proc_params)
    add_result_proc_args(proc_params)

    output = parser.add_argument_group("Output file arguments")
    output.add_argument(
        "-n", "--name", nargs='*', default=[],
        help="Output file base name(s). (Default: input file name without extension)"
    )
    output.add_argument(
        "-o", "--outdir", nargs='?', default='.',
        help="Output directory. (Default: current directory)"
    )

    return parser


def get_precalc_parser():
    parser = argparse.ArgumentParser(
        description="Pre-calculate mappability region statistics for "
                    "PyMaSC successive algorithm."
    )

    add_common_args(parser)

    proc_args = parser.add_argument_group("Processing behaviors")
    add_multiprocess_args(proc_args)

    map_args = parser.add_argument_group("Input mappability file arguments")
    add_mappability_args(map_args)

    proc_params = parser.add_argument_group("PyMaSC parameters")
    add_shift_arg(proc_params)
    proc_params.add_argument(
        "-r", "--max-readlen", nargs='?', type=int, action=ForceNaturalNumber, default=1000,
        help="Set max read length to calculate mappable region length."
    )

    return parser


def get_plot_parser():
    parser = argparse.ArgumentParser(
        description="Plot figures from PyMaSC statistic outputs."
    )

    add_common_args(parser)

    input_args = parser.add_argument_group("Input alignment file arguments")
    input_args.add_argument(
        "statfile", nargs='?',
        help="A base path to the statistic files (*_stats.tab, *_cc.tab and *_masc.tab) "
             "to plot figures."
    )
    input_args.add_argument(
        "-s", "--stats", nargs='?',
        help="Specify path to a statistic file (*_stats.tab) separately."
    )
    input_args.add_argument(
        "-c", "--cc", nargs='?',
        help="Specify path to a cross-correlation table file (*_cc.tab) separately."
    )
    input_args.add_argument(
        "-m", "--masc", nargs='?',
        help="Specify path to a mappability sensitive cross-correlation file (*_stats.tab) separately."
    )

    proc_params = parser.add_argument_group("PyMaSC parameters")
    add_result_proc_args(proc_params)
    add_liblen_arg(proc_params)

    output = parser.add_argument_group("Output file arguments")
    output.add_argument(
        "-n", "--name", nargs='?',
        help="Change output file base name. (Default: same as name field in input)"
    )
    output.add_argument(
        "-o", "--outdir", nargs='?', default='.',
        help="Output directory. (Default: current directory)"
    )

    return parser
