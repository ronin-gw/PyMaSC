from __future__ import print_function

import os.path
import logging

from PyMaSC.utils.output import catch_IOError

STATSFILE_SUFFIX = "_stats.tab"
STAT_ATTR = (
    ("Read length", "read_len"),
    ("Expected library length", "library_len"),
    ("Genome length", "genomelen"),
    ("Forward reads", "forward_sum"),
    ("Reverse reads", "reverse_sum"),
    ("Minimum CC", "cc_min"),
    ("CC at read length", "ccrl"),
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
