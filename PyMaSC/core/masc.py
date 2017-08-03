import logging
import types

import numpy as np

from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


class ReadUnsortedError(IndexError):
    pass


class CCCalculator(object):
    # def next(self):
    #     pass

    @staticmethod
    def _get_shift_update_fun(obj, forward_buff_attr_name, reverse_buff_attr_name, ccbins_attr_name):
        def update_func(self, offset, update_forward=False):
            forward_buff = getattr(self, forward_buff_attr_name)
            reverse_buff = getattr(self, reverse_buff_attr_name)
            ccbins = getattr(self, ccbins_attr_name)

            remain_buff_len = forward_buff.size - offset - 1
            if remain_buff_len > 0:
                rotate = offset + (0 if update_forward else 1)
                i = -1
                for i, rev_state in enumerate(reverse_buff[:rotate]):
                    if rev_state:
                        ccbins[i:] += forward_buff[i:][::-1]
                forward_buff = np.concatenate((forward_buff[rotate:], np.repeat(0, rotate)))

                if update_forward:
                    if reverse_buff[offset:offset+1]:
                        ccbins += forward_buff[::-1]
                    forward_buff = np.concatenate((forward_buff[1:], [1]))
            else:
                for i, rev_state in enumerate(reverse_buff[:forward_buff.size]):
                    if rev_state:
                        ccbins[i:] += forward_buff[i:][::-1]
                forward_buff.fill(0)

                if update_forward:
                    forward_buff[-1] = 1

            setattr(self, forward_buff_attr_name, forward_buff)
            setattr(self, reverse_buff_attr_name, reverse_buff[offset+1:])
            setattr(self, ccbins_attr_name, ccbins)

        return types.MethodType(update_func, obj)

    def _get_pass_func(self):
        def pass_func(self, *args, **kwargs):
            pass

        return types.MethodType(pass_func, self)

    def __init__(self, max_shift, references, lengths, chi2_pval, calc_masc=False):
        """
        self._ccbins: summation bins to calc cc
        self._forward_buff: forward strand read mappability buff.
                            fixed length: `self.max_shift`
                            The last base points `self._deq_tail_pos.`
        self._reverse_buff: reverse strand read mappability buff.
                            The first base points `self._deq_tail_pos + 1`
        """
        self.max_shift = max_shift
        self.ref2genomelen = dict(zip(references, lengths))
        self.genomelen = sum(lengths)
        # self.genomelen = self.ref2genomelen["chr1"]
        self.chi2_p_thresh = chi2_pval
        self.calc_masc = calc_masc

        # internal buff
        self._chr = None

        # Naive CC
        # stats
        self.forward_sum = self.reverse_sum = 0
        self.reference2forward_sum = {ref: 0 for ref in references}
        self.reference2reverse_sum = {ref: 0 for ref in references}
        self.forward_read_len_sum = self.reverse_read_len_sum = 0
        self.ccbins = []
        self.reference2ccbins = {ref: None for ref in references}
        # internal buff
        self._forward_sum = self._reverse_sum = 0
        self._ccbins = np.repeat(0, max_shift)
        self._forward_buff = np.repeat(0, self.max_shift)
        self._reverse_buff = np.array([], dtype=np.int64)
        # calc func
        self._shift_with_update_ncc = self._get_shift_update_fun(
            self, "_forward_buff", "_reverse_buff", "_ccbins"
        )

        # MSCC
        if calc_masc:
            # stats
            self.mappable_forward_sum = self.mappable_reverse_sum = 0
            self.mappable_reference2forward_sum = {ref: 0 for ref in references}
            self.mappable_reference2reverse_sum = {ref: 0 for ref in references}
            self.mappable_forward_len_sum = self.mappable_reverse_len_sum = 0
            self.mappable_ccbins = []
            self.mappable_reference2ccbins = {ref: None for ref in references}
            # internal buff
            self._m_forward_sum = self._m_reverse_sum = 0
            self._m_ccbins = np.repeat(0, max_shift)
            self._m_forward_buff = np.repeat(0, self.max_shift)
            self._m_reverse_buff = np.array([], dtype=np.int64)
            # calc func
            self._shift_with_update_mscc = self._get_shift_update_fun(
                self, "_m_forward_buff", "_m_reverse_buff", "_m_ccbins"
            )
        else:
            # stats
            self.mappable_forward_sum = self.mappable_reverse_sum = None
            self.mappable_reference2forward_sum = None
            self.mappable_reference2reverse_sum = None
            self.mappable_forward_len_sum = self.mappable_reverse_len_sum = None
            self.mappable_ccbins = None
            self.mappable_reference2ccbins = None
            # calc func
            self._shift_with_update_mscc = self._get_pass_func()

        self._init_pos_buff()
        self._progress = ProgressBar()

    def _init_pos_buff(self):
        self._forward_sum = self._reverse_sum = 0
        self._ccbins.fill(0)
        self._forward_buff.fill(0)
        self._reverse_buff.fill(0)
        if self.calc_masc:
            self._m_forward_sum = self._m_reverse_sum = 0
            self._m_ccbins.fill(0)
            self._m_forward_buff.fill(0)
            self._m_reverse_buff.fill(0)

        self._deq_tail_pos = self.max_shift
        self._last_pos = 0
        self._last_forward_pos = self._last_reverse_pos = 0

    def _flush(self):
        self._shift_with_update_ncc(self._reverse_buff.size - 1)

        self.forward_sum += self._forward_sum
        self.reverse_sum += self._reverse_sum
        self.reference2forward_sum[self._chr] = self._forward_sum
        self.reference2reverse_sum[self._chr] = self._reverse_sum
        self.reference2ccbins[self._chr] = tuple(self._ccbins)

        if self.calc_masc:
            self.mappable_forward_sum += self._m_forward_sum
            self.mappable_reverse_sum += self._m_reverse_sum
            self.mappable_reference2forward_sum[self._chr] = self._m_forward_sum
            self.mappable_reference2reverse_sum[self._chr] = self._m_reverse_sum
            self._shift_with_update_mscc(self._m_reverse_buff.size - 1)
            self.mappable_reference2ccbins[self._chr] = tuple(self._m_ccbins)

        self._init_pos_buff()
        self._progress.clean()

    def _check_pos(self, chrom, pos):
        if chrom != self._chr:
            if self._chr is not None:
                self._flush()
            self._chr = chrom

            logger.info("Calc {}...".format(chrom))

            self._progress.set(self.ref2genomelen[chrom])

        if pos < self._last_pos:
            raise ReadUnsortedError

        self._progress.update(pos)

        self._last_pos = pos

    def feed_forward_read(self, chrom, pos, readlen, is_mappable=False):
        logger.debug("Got forward read {}".format(pos))
        self._check_pos(chrom, pos)
        if self._last_forward_pos == pos:
            logger.debug("Read is duplicated")
            return None
        self._last_forward_pos = pos

        logger.debug(self._reverse_buff)
        # logger.debug(self._m_reverse_buff)
        # logger.debug(self._m_forward_buff)
        logger.debug(self._forward_buff)

        self._forward_sum += 1
        self.forward_read_len_sum += readlen

        if is_mappable:
            self._m_forward_sum += 1

        offset = pos - self._deq_tail_pos - 1
        if offset < 0:
            logger.debug("Update buff {}".format(offset))
            self._forward_buff[offset] = 1
            if is_mappable:
                self._m_forward_buff[offset] = 1
        else:
            logger.debug("Shift buff {}".format(offset))
            self._shift_with_update(offset, update_forward=True, is_mappabile=is_mappable)

        logger.debug(self._reverse_buff)
        # logger.debug(self._m_reverse_buff)
        # logger.debug(self._m_forward_buff)
        logger.debug(self._forward_buff)

    def feed_reverse_read(self, chrom, pos, readlen, is_mappable=False):
        logger.debug("Got reverse read {}".format(pos))
        self._check_pos(chrom, pos)
        if self._last_reverse_pos == pos:
            logger.debug("Read is duplicated")
            return None
        self._last_reverse_pos = pos

        logger.debug(self._reverse_buff)
        # logger.debug(self._m_reverse_buff)
        # logger.debug(self._m_forward_buff)
        logger.debug(self._forward_buff)

        self._reverse_sum += 1
        self.reverse_read_len_sum += readlen
        offset = pos - self._deq_tail_pos - 1
        rev_pos = pos + readlen - 1

        if is_mappable:
            self._m_reverse_sum += 1

        rev_offset = rev_pos - self._deq_tail_pos
        if rev_offset < 1:
            logger.debug("Update immediately {}".format(offset))
            self._ccbins[:rev_offset-1] += self._forward_buff[rev_offset-2::-1]
            if is_mappable:
                self._m_ccbins[:rev_offset-1] += self._m_forward_buff[rev_offset-2::-1]
        else:
            if offset > 0:
                logger.debug("Shift buff {}".format(offset))
                self._shift_with_update(offset - 1, is_mappabile=is_mappable)
                rev_offset -= offset
            try:
                self._reverse_buff[rev_offset - 1] = 1
                if is_mappable:
                    self._m_reverse_buff[rev_offset - 1] = 1
            except IndexError:
                self._reverse_buff = np.append(
                    self._reverse_buff, [0] * (rev_offset - self._reverse_buff.size - 1) + [1]
                )
                if is_mappable:
                    self._m_reverse_buff = np.append(
                        self._m_reverse_buff, [0] * (rev_offset - self._m_reverse_buff.size - 1) + [1]
                    )

        logger.debug(self._reverse_buff)
        # logger.debug(self._m_reverse_buff)
        # logger.debug(self._m_forward_buff)
        logger.debug(self._forward_buff)

    def _shift_with_update(self, offset, update_forward=False, is_mappabile=False):
        self._shift_with_update_ncc(offset, update_forward)
        if is_mappabile:
            self._shift_with_update_mscc(offset, update_forward)

        self._deq_tail_pos += offset + 1

    def _update_cc(self, rev_offset, is_mappabile=False):
        self._update_ncc(rev_offset)
        if is_mappabile:
            self._update_mscc(rev_offset)

    def finishup_calculation(self):
        self._flush()

        for bins in zip(*[v for v in self.reference2ccbins.values() if v]):
            self.ccbins.append(sum(bins))

        if self.calc_masc:
            for bins in zip(*[v for v in self.mappable_reference2ccbins.values() if v]):
                self.mappable_ccbins.append(sum(bins))
