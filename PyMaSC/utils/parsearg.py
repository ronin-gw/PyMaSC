import argparse
import logging


def _make_upper(s):
    return s.upper()


class StoreLoggingLevel(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, getattr(logging, values))


class ForceNaturalNumber(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values < 1:
            parser.error("argument -d/--max-shift: shift length must be > 0.")
        setattr(namespace, self.dest, values)


def add_common_args(parser):
    parser.add_argument("-v", "--log-level", nargs='?', type=_make_upper, default=logging.INFO, action=StoreLoggingLevel,
                        choices=("INFO", "WARNING", "ERROR", "CRITICAL"),
                        help="Set verbosity. (Default: INFO)")
    parser.add_argument("--color", nargs='?', type=_make_upper, default=None, choices=("TRUE", "FALSE"),
                        help="Coloring log. (Default: auto)")


def add_mappability_args(group):
    group.add_argument("-m", "--mappable", nargs=1,
                       help="BigWig/BigBed format mappable region file.")
    group.add_argument("--map-path", nargs=1,
                       help="Read/Save path for mappability stats. (Default: <region file>.json)")


def add_shift_arg(group):
    group.add_argument("-d", "--max-shift", nargs='?', type=int, action=ForceNaturalNumber, default=1000)


def get_parser():
    parser = argparse.ArgumentParser(description="Mappability-sensitive cross-correlation calculator for NGS data")

    add_common_args(parser)

    input_args = parser.add_argument_group("Input file arguments")
    input_args.add_argument("reads", nargs="+",
                            help="SAM/BAM format mapped reads. Input must be sorted.")
    input_args.add_argument("-f", "--format", nargs='?', type=_make_upper, default=None, choices=("BAM", "SAM"),
                            help="Specify input file type. (Default: auto)")
    add_mappability_args(input_args)

    params = parser.add_argument_group("Parameters")
    add_shift_arg(params)
    params.add_argument("-q", "--mapq", nargs='?', type=int, default=1,
                        help="Filter out reads which have less than specified MAPQ. (Default: 1)")
    params.add_argument("-p", "--chi2-pval", nargs='?', type=float, default=0.05,
                        help="Chi-squared test p-value threshold to check strand specificity. (Default: 0.05)")

    output = parser.add_argument_group("Output file arguments")
    output.add_argument("--outdir", nargs='?', default='.',
                        help="Output directory. (Default: current directory)")

    return parser