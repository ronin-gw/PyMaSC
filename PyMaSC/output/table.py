from __future__ import print_function

import logging
from functools import partial
import csv
from itertools import chain
from collections import defaultdict
from copy import copy

import numpy as np

from PyMaSC.utils.output import catch_IOError

CCOUTPUT_SUFFIX = "_cc.tab"
MSCCOUTPUT_SUFFIX = "_mscc.tab"
NREADOUTPUT_SUFFIX = "_nreads.tab"

logger = logging.getLogger(__name__)


class TableIO(object):
    DIALECT = "excel-tab"

    def __init__(self, path, mode='r'):
        self.path = path
        self.mode = mode
        self.fp = None

        if mode not in ('r', 'w'):
            raise NotImplemented

    def open(self):
        self.fp = open(self.path, self.mode, newline='')

    def __enter__(self):
        self.open()
        return self

    def __exit__(self, ex_type, ex_value, trace):
        self.fp.close()

    close = __exit__

    def read(self, typefunc):
        tab = csv.reader(self.fp, dialect=TableIO.DIALECT)
        header = next(tab)[1:]
        table = dict(zip(header, zip(*(map(typefunc, row[1:]) for row in tab))))
        return header, table

    def write(self, header, body):
        tab = csv.writer(self.fp, dialect=TableIO.DIALECT)
        tab.writerow(("shift", ) + tuple(header))
        tab.writerows((i, ) + tuple(row) for i, row in enumerate(body))


@catch_IOError(logger)
def _load_table(path, logfmt):
    logger.info(logfmt.format(path))

    with TableIO(path) as tab:
        _, table = tab.read(float)

    if "whole" in table:
        table.pop("whole")
    else:
        logger.warning("Mandatory column 'whole' not found")

    return table


load_cc = partial(_load_table, logfmt="Load CC table from '{}'")
load_masc = partial(_load_table, logfmt="Load MSCC table from '{}'")


@catch_IOError(logger)
def _output_cctable(outfile, ccr, suffix, target_attr):
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
    @staticmethod
    def make_forward_reverse_tables(header, body):
        forward = defaultdict(list)
        reverse = defaultdict(list)

        for row in body:
            for key, pair in zip(header, row[1:]):
                f, r = map(int, pair.split('-'))
                forward[key].append(f)
                reverse[key].append(r)

        return forward, reverse

    def read(self):
        tab = csv.reader(self.fp, dialect=TableIO.DIALECT)

        header = next(tab)[1:]

        first = next(tab)
        if first[0] == "raw":
            forward_sum, reverse_sum = self.make_forward_reverse_tables(header, [first])
            forward_sum = {k: v[0] for k, v in forward_sum.items()}
            reverse_sum = {k: v[0] for k, v in reverse_sum.items()}
            mscc_iter = tab
        else:
            forward_sum = reverse_sum = {}
            mscc_iter = chain([first], tab)

        mappable_forward_sum, mappable_reverse_sum = self.make_forward_reverse_tables(header, mscc_iter)
        if not mappable_forward_sum:
            mappable_forward_sum = mappable_reverse_sum = {}

        return forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum

    @staticmethod
    def _make_table(rowname, forward, reverse):
        return chain([rowname], ("{}-{}".format(f, r) for f, r in zip(forward, reverse)))

    def write(self, header, forward_sum, reverse_sum, mappable_forward_sum, mappable_reverse_sum):
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
def load_nreads_table(path):
    logger.info("Load Nreads table from '{}'".format(path))

    with NReadsIO(path) as tab:
        dicts = tab.read()

    for d in dicts:
        if "whole" in d:
            d.pop("whole")
        elif d:
            logger.warning("Mandatory column 'whole' not found")

    if all(not d for d in dicts):
        logger.critical("Nothing to load.")
        raise KeyError

    return dicts


@catch_IOError(logger)
def output_nreads_table(outfile, ccr):
    outfile += NREADOUTPUT_SUFFIX
    logger.info("Output '{}'".format(outfile))

    def _make_dict_with_whole(add_target, target1, target2):
        if not add_target:
            return None
        d = copy(add_target)
        d["whole"] = getattr(getattr(ccr.whole, target1), target2)
        return d

    with NReadsIO(outfile, 'w') as tab:
        tab.write(["whole"] + sorted(ccr.references),
                  _make_dict_with_whole(ccr.ref2forward_sum, "cc", "forward_sum"),
                  _make_dict_with_whole(ccr.ref2reverse_sum, "cc", "reverse_sum"),
                  _make_dict_with_whole(ccr.mappable_ref2forward_sum, "masc", "forward_sum"),
                  _make_dict_with_whole(ccr.mappable_ref2reverse_sum, "masc", "reverse_sum")
                  )
