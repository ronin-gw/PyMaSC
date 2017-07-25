import argparse
import logging


def _make_upper(s):
    return s.upper()


class StoreLoggingLevel(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, getattr(logging, values))


def get_parser():
    parser = argparse.ArgumentParser(description="Mappability-sensitive cross-correlation calculator for NGS data")

    parser.add_argument("-v", "--log-level", nargs='?', type=_make_upper, default=logging.INFO, action=StoreLoggingLevel,
                        choices=("INFO", "WARNING", "ERROR", "CRITICAL"),
                        help="Set verbosity. (Default: INFO)")

    input_args = parser.add_argument_group("Input file arguments")
    input_args.add_argument("reads", nargs="+",
                            help="SAM/BAM format mapped reads.")
    input_args.add_argument("-f", "--format", nargs='?', type=_make_upper, default=None, choices=("BAM", "SAM"),
                            help="Specify input file type. (Default: auto)")
    input_args.add_argument("-m", "--mappable", nargs="+",
                            help="BED format mappable regions. Use auto detection if not specified.")

    params = parser.add_argument_group("Parameters")
    params.add_argument("-d", "--max-shift", nargs='?', type=int, default=1000)
    params.add_argument("-q", "--mapq", nargs='?', type=int, default=1,
                        help="Filter out reads which have less than specified MAPQ. (Default: 1)")

    output = parser.add_argument_group("Output file arguments")
    output.add_argument("-n", "--name", nargs=1,
                        help="Basename for output files. (Default: for each read file name)")
    output.add_argument("--outdir", nargs='?', default='.',
                        help="Output directory. (Default: the current directory)")

    return parser
