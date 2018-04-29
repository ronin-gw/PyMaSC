from __future__ import print_function

import os.path
import logging

import numpy as np

from PyMaSC.utils.output import catch_IOError

CCOUTPUT_SUFFIX = "_cc.tab"
MSCCOUTPUT_SUFFIX = "_mscc.tab"
STATSFILE_SUFFIX = "_stats.tab"

logger = logging.getLogger(__name__)


def _output_cctable(outfile, ccr, target_attr):
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


@catch_IOError(logger, "CC table")
def output_cc(outfile, ccr):
    outfile += CCOUTPUT_SUFFIX
    if ccr.skip_ncc:
        return None
    _output_cctable(outfile, ccr, "cc")


@catch_IOError(logger, "MSCC table")
def output_mscc(outfile, ccr):
    outfile += MSCCOUTPUT_SUFFIX
    if not ccr.calc_masc:
        return None
    _output_cctable(outfile, ccr, "masc")


@catch_IOError(logger, "stats")
def output_stats(outfile, ccr):
    basename = os.path.basename(outfile)
    outfile += STATSFILE_SUFFIX
    logger.info("Output '{}'".format(outfile))

    with open(outfile, 'w') as f:
        print("{}\t{}".format("name", basename), file=f)
        for row, attr in (
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
        ):
            print("{}\t{}".format(row, getattr(ccr.whole, attr, "nan") or "nan"), file=f)
