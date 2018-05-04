from __future__ import print_function

import os.path
import logging
from functools import partial

import numpy as np

from PyMaSC.utils.output import catch_IOError

CCOUTPUT_SUFFIX = "_cc.tab"
MSCCOUTPUT_SUFFIX = "_mscc.tab"
STATSFILE_SUFFIX = "_stats.tab"

logger = logging.getLogger(__name__)


@catch_IOError(logger)
def _output_cctable(outfile, ccr, suffix, target_attr):
    outfile += suffix
    logger.info("Output '{}'".format(outfile))

    with open(outfile, 'w') as f:
        keys = sorted(ccr.references)
        cc = getattr(ccr.whole, target_attr)
        ref2cc = {k: getattr(ccr.ref2stats[k], target_attr) for k in keys}
        keys = [k for k, v in ref2cc.items() if v is not None and not np.isnan(v).all()]

        print('\t'.join(["shift", "whole"] + keys), file=f)
        fmt = '{}\t' * (len(keys) + 1) + '{}'
        for i, cc in enumerate(cc):
            print(fmt.format(i, cc, *[ref2cc[k][i] for k in keys]), file=f)


output_cc = partial(_output_cctable, suffix=CCOUTPUT_SUFFIX, target_attr="cc")
output_mscc = partial(_output_cctable, suffix=MSCCOUTPUT_SUFFIX, target_attr="masc")


@catch_IOError(logger)
def _load_table(path, logfmt):
    logger.info(logfmt.format(path))

    with open(path) as f:
        header = f.readline().rstrip().split('\t')[1:]
        table = dict(zip(header, zip(*(map(float, l.split('\t')[1:]) for l in f))))
    if "whole" not in table:
        logger.critical()
    whole = table.pop("whole")
    return whole, table


load_cc = partial(_load_table, logfmt="Load CC table from '{}'")
load_masc = partial(_load_table, logfmt="Load MSCC table from '{}'")


STAT_ATTR = (
    ("Genome length", "genomelen"),
    ("Forward reads", "forward_sum"),
    ("Reverse reads", "reverse_sum"),
    ("Read length", "read_len"),
    ("Minimum CC", "cc_min"),
    ("CC at read length", "ccrl"),
    ("Expected library length", "library_len"),
    ("CC at library length", "ccfl"),
    ("NSC", "nsc"),
    ("RSC", "rsc"),
    ("Estimated library length", "est_lib_len"),
    ("Estimated CC at library length", "est_ccfl"),
    ("Estimated NSC", "est_nsc"),
    ("Estimated RSC", "est_rsc")
)


@catch_IOError(logger)
def output_stats(outfile, ccr):
    basename = os.path.basename(outfile)
    outfile += STATSFILE_SUFFIX
    logger.info("Output '{}'".format(outfile))

    with open(outfile, 'w') as f:
        print("{}\t{}".format("Name", basename), file=f)
        for row, attr in STAT_ATTR:
            print("{}\t{}".format(row, getattr(ccr.whole, attr, "nan") or "nan"), file=f)


@catch_IOError(logger)
def load_stats(path):
    logger.info("Load statistics from '{}'.".format(path))
    stat2attr = {k: v for k, v in STAT_ATTR if v in
                 ("genomelen", "forward_sum", "reverse_sum", "read_len", "library_len")}
    attrs = {}
    with open(path) as f:
        for l in f:
            row, val = l.split('\t', 1)
            if row == "Name":
                attrs["name"] = val.rstrip()
            if row in stat2attr:
                val = int(val) if val.strip().isdigit() else None
                attrs[stat2attr[row]] = val
    return attrs
