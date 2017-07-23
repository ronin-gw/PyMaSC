#!/usr/bin/env python
import logging
import argparse

from pysam import AlignmentFile
import numpy

logger = logging.getLogger(__name__)


class CCCalculator(object):
    def __init__(self, max_shift, genomelen):
        """
        self._ccbins: summation bins to calc cc
        self._forward_buff: forward strand read mappability buff.
                            fixed length: `self.max_shift`
                            The last base points `self._deq_tail_pos.`
        self._reverse_buff: reverse strand read mappability buff.
                            The first base points `self._deq_tail_pos + 1`
        """
        self.max_shift = max_shift
        self.genomelen = genomelen

        self._chr = None
        # self._ccbins = [0] * max_shift
        self._ccbins = numpy.array([0] * max_shift)
        self._forward_buff = numpy.array([0] * self.max_shift)
        self._reverse_buff = numpy.array([], dtype=numpy.int64)
        self._init_pos_buff()

        self.forward_sum = self.reverse_sum = 0
        self.cc = []

    def _init_pos_buff(self):
        # self._forward_buff = [0] * self.max_shift
        self._forward_buff.fill(0)
        self._reverse_buff.fill(0)

        self._deq_tail_pos = self.max_shift
        self._last_pos = 0
        self._last_forward_pos = self._last_reverse_pos = 0

    def _flush(self):
        self._shift_with_update(len(self._reverse_buff) - 1)
        self._init_pos_buff()

    def _check_pos(self, chrom, pos):
        if chrom != self._chr:
            if self._chr is not None:
                self._flush()
            self._chr = chrom
            logger.info("Process {}...".format(chrom))

        if pos < self._last_pos:
            raise IndexError("Input read must be sorted.")
        self._last_pos = pos

    def feed_forward_read(self, chrom, pos):
        logger.debug("Got forward read {}".format(pos))
        self._check_pos(chrom, pos)
        if self._last_forward_pos == pos:
            logger.debug("Read is duplicated")
            return None
        self._last_forward_pos = pos

        logger.debug(self._reverse_buff)
        logger.debug(self._forward_buff)

        self.forward_sum += 1
        offset = pos - self._deq_tail_pos - 1
        if offset < 0:
            logger.debug("Update buff {}".format(offset))
            self._forward_buff[offset] = 1
        else:
            logger.debug("Shift buff {}".format(offset))
            self._shift_with_update(offset, update_forward=True)

        logger.debug(self._reverse_buff)
        logger.debug(self._forward_buff)

    def feed_reverse_read(self, chrom, pos, readlen):
        logger.debug("Got reverse read {}".format(pos))
        self._check_pos(chrom, pos)
        if self._last_reverse_pos == pos:
            logger.debug("Read is duplicated")
            return None
        self._last_reverse_pos = pos

        logger.debug(self._reverse_buff)
        logger.debug(self._forward_buff)

        self.reverse_sum += 1
        offset = pos - self._deq_tail_pos - 1
        rev_pos = pos + readlen - 1
        rev_offset = rev_pos - self._deq_tail_pos
        if rev_offset < 1:
            logger.debug("Update immediately {}".format(offset))
            self._update_cc(rev_offset)
        else:
            if offset > 0:
                logger.debug("Shift buff {}".format(offset))
                self._shift_with_update(offset - 1)
                rev_offset -= offset
            try:
                self._reverse_buff[rev_offset - 1] = 1
            except IndexError:
                self._reverse_buff = numpy.append(self._reverse_buff, [0] * (rev_offset - len(self._reverse_buff) - 1) + [1])

        logger.debug(self._reverse_buff)
        logger.debug(self._forward_buff)

    def _shift_with_update(self, offset, update_forward=False):
        remain_buff_len = len(self._reverse_buff) - offset - 1

        if remain_buff_len > 0:
            for rev_state in self._reverse_buff[:offset + (1 if update_forward else 0)]:
                if rev_state:
                    self._update_cc()
                # self._forward_buff = self._forward_buff[1:] + [0]
                self._forward_buff = numpy.concatenate((self._forward_buff[1:], [0]))

            if update_forward:
                if self._reverse_buff[offset]:
                    self._update_cc()
                self._forward_buff = numpy.concatenate((self._forward_buff[1:], [1]))

            self._reverse_buff = self._reverse_buff[offset+1:]
        else:
            for rev_state in self._reverse_buff:
                if rev_state:
                    self._update_cc()
                # self._forward_buff = self._forward_buff[1:] + [0]
                self._forward_buff = numpy.concatenate((self._forward_buff[1:], [0]))

            self._forward_buff.fill(0)
            if update_forward:
                self._forward_buff[0] = 1
            self._reverse_buff = numpy.array([], dtype=numpy.int64)

        self._deq_tail_pos += offset + 1

    def _update_cc(self, offset=0):
        if offset:
            for i, state in enumerate(self._forward_buff[offset-1::-1]):
                self._ccbins[i] += state
        else:
            self._ccbins += self._forward_buff[::-1]

    def calc_cc(self):
        self._flush()

        self.forward_mean = float(self.forward_sum) / self.genomelen
        self.reverse_mean = float(self.reverse_sum) / self.genomelen

        self.forward_var = self.forward_mean * (1 - self.forward_mean)
        self.reverse_var = self.reverse_mean * (1 - self.reverse_mean)

        sum_prod = self.forward_mean * self.reverse_mean
        var_geomean = (self.forward_var * self.reverse_var) ** 0.5

        # for i in xrange(self.max_shift):
        #     self.cc.append(
        #         ((self._ccbins[i] / float(self.genomelen - i - 1)) - sum_prod) / var_geomean
        #     )
        self.cc = (self._ccbins / (self.genomelen - numpy.array(range(self.max_shift)) - 1.) - sum_prod) / var_geomean

        return self.cc


def _main():
    parser = argparse.ArgumentParser(description="Mappability-sensitive cross-correlation calculator for NGS data")
    parser.add_argument("reads", nargs="+", help="SAM/BAM format mapped reads.")
    parser.add_argument("-d", "--max-shift", nargs='?', type=int, default=1000)
    parser.add_argument("-m", "--mappable", nargs="+", help="BED format mappable regions. Use auto detection if not specified.")
    parser.add_argument("-q", "--mapq", nargs='?', type=int, default=1, help="Filter out reads which have less than specified MAPQ. (default: 1)")
    parser.add_argument("-n", "--name", nargs=1, help="Basename for output files. (default: for each read file name)")
    parser.add_argument("--outdir", nargs='?', default='.', help="Output directory. (default: the current directory)")

    args = parser.parse_args()
    chrom = None

    with AlignmentFile(args.reads[0]) as f:
        ccc = CCCalculator(args.max_shift, sum(f.lengths))
        # ccc = CCCalculator(args.max_shift, f.lengths[0])

        for i, read in enumerate(f):
            if read.mapping_quality < args.mapq or read.is_duplicate:
                continue

            if read.is_reverse:
                ccc.feed_reverse_read(read.reference_name, read.pos+1, read.query_length)
            else:
                ccc.feed_forward_read(read.reference_name, read.pos+1)

            if chrom != read.reference_name:
                if chrom is not None:
                    break
                chrom = read.reference_name
            # elif not (i % 10000):
            #     sys.stderr.write("\r{}".format(i))

    # ccc = CCCalculator(args.max_shift, 12789634)
    # for i, read in enumerate(sys.stdin):
    #     _, flag, chrom, pos, mapq, _, _, _, seq, _ = read.split('\t', 9)
    #     flag = int(flag)
    #
    #     if int(mapq) < args.mapq or flag & 1024:
    #         continue
    #     if flag & 16:
    #         ccc.feed_reverse_read(chrom, int(pos), len(seq))
    #     else:
    #         ccc.feed_forward_read(chrom, int(pos))f

    # print ccc.calc_cc()
    # for i in tuple(ccc.calc_cc()):
    #     print i


if __name__ == "__main__":
    # logger.setLevel(logging.DEBUG)
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler())
    _main()
