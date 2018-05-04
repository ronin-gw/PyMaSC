import logging

from libc.stdint cimport int_fast64_t as int64
from libc.string cimport strcmp, strcpy
from cython cimport boundscheck, wraparound

from .bitarray cimport bit_index_t, bitarray
from PyMaSC.reader.bx.bigwig_file cimport BigWigFile

from PyMaSC.core.ncc import ReadUnsortedError
from PyMaSC.utils.progress import ProgressBar

logger = logging.getLogger(__name__)


cdef class CCBitArrayCalculator(object):
    cdef public:
        double MAPPABILITY_THRESHOLD
        int64 EXTRA_ALLOCATE_SIZE

    cdef readonly:
        int64 max_shift
        int64 read_len
        dict ref2genomelen
        int64 genomelen
        int64 forward_sum, reverse_sum
        dict ref2forward_sum, ref2reverse_sum
        dict ref2mappable_forward_sum,  ref2mappable_reverse_sum
        int64 forward_read_len_sum, reverse_read_len_sum
        dict ref2ccbins, ref2mascbins, ref2mappable_len
        bint skip_ncc
    cdef:
        str _chr
        list _solved_chr
        bitarray _forward_array, _reverse_array
        object _bwfeeder, logger_lock, _progress

        int64 _last_pos, _last_forward_pos, _last_reverse_pos, _array_extend_size
        bint _buff_flashed, _allow_bitarray_progress

    def __init__(self, int64 max_shift, int64 read_len, references, lengths,
                 bwfeeder=None, skip_ncc=False, logger_lock=None, progress_bar=None):
        self.MAPPABILITY_THRESHOLD = 1.
        # bitarary extra allocate size for reads longer than self.read_len
        self.EXTRA_ALLOCATE_SIZE = 100

        self.max_shift = max_shift
        self.read_len = read_len
        self.ref2genomelen = dict(zip(references, lengths))
        self.genomelen = sum(lengths)
        # stats
        self.forward_sum = self.reverse_sum = 0
        self.ref2forward_sum = {ref: 0 for ref in references}
        self.ref2reverse_sum = {ref: 0 for ref in references}
        # {forward,reverse}_sum[d[0...i]] = Number of reads in Md
        self.ref2mappable_forward_sum = {ref: None for ref in references}
        self.ref2mappable_reverse_sum = {ref: None for ref in references}
        self.forward_read_len_sum = self.reverse_read_len_sum = 0
        self.ref2ccbins = {ref: None for ref in references}
        self.ref2mascbins = {ref: None for ref in references}
        self.ref2mappable_len = {ref: None for ref in references}

        # internal buff
        self._chr = ''
        self._solved_chr = []
        self._buff_flashed = False
        self._array_extend_size = self.read_len + self.max_shift + self.EXTRA_ALLOCATE_SIZE

        #
        self._bwfeeder = bwfeeder
        self.skip_ncc = skip_ncc
        self.logger_lock = logger_lock
        self._progress = progress_bar if progress_bar else ProgressBar()
        # Issue: bitarray progress stacks in multiprocessing
        self._allow_bitarray_progress = progress_bar is None

    def _logging_info(self, msg):
        if self.logger_lock:
            self.logger_lock.acquire()
        logger.info(msg)
        if self.logger_lock:
            self.logger_lock.release()

    def _init_progress(self, body, prefix, suffix, maxval):
        self._progress.body = body
        self._progress.fmt = "\r" + prefix + "{:<" + str(len(body)) + "}" + suffix
        self._progress.set(self._chr, maxval)

    def _init_mappability_progress(self):
        self._progress.enable_bar()
        self._init_progress("@@@@@^" * 12, '^', ' ', self.ref2genomelen[self._chr])

    def _init_bitarray_progress(self):
        if self._allow_bitarray_progress:
            self._init_progress("xxxxxX" * 12, 'X', ' ', self.max_shift)
        else:
            self._progress.disable_bar()

    def _init_buff(self):
        # To access by 1-based index, allocate arrays with 1 more size.
        # Additionally, extra spaces same or more than (read_length - 1)
        # is needed to contain reverse strand reads.
        chromsize = self.ref2genomelen[self._chr] + self._array_extend_size
        self._forward_array = bitarray(chromsize)
        self._reverse_array = bitarray(chromsize)

        self._last_pos = 0
        self._last_forward_pos = self._last_reverse_pos = 0

    def flush(self):
        if self._buff_flashed:
            return None

        # keep this method safe to call; for chroms which have no reads.
        if self._chr != '':
            self._calc_correlation()

        self._buff_flashed = True

    cdef inline int _calc_correlation(self) except -1:
        chromsize = self.ref2genomelen[self._chr] + 1

        cdef bit_index_t forward_sum, reverse_sum
        cdef list ccbin
        cdef list mappable_len, mappable_forward_sum, mappable_reverse_sum, mascbin
        cdef bitarray mappability, forward_mappability, reverse_mappability
        cdef bitarray mappable_forward = bitarray(chromsize)
        cdef bitarray mappable_reverse = bitarray(chromsize)
        cdef bitarray doubly_mappable_region = bitarray(chromsize)
        cdef list reverse_mappability_buff
        cdef int i

        # naive cross-correlation
        if not self.skip_ncc:
            forward_sum = self._forward_array.count()
            reverse_sum = self._reverse_array.count()
            self.forward_sum += forward_sum
            self.reverse_sum += reverse_sum
            self.ref2forward_sum[self._chr] = forward_sum
            self.ref2reverse_sum[self._chr] = reverse_sum
            ccbin = self.ref2ccbins[self._chr] = []

        # mappability-sensitive cross-correlation
        try:
            mappability = self._load_mappability()
        except KeyError as e:
            self._logging_info("Mappability for '{}' not found. "
                               "Skip calc mappability sensitive CC.".format(e.args[0]))
            mappability = None
        if mappability:
            mappable_len = self.ref2mappable_len[self._chr] = [None] * self.read_len
            mappable_forward_sum = self.ref2mappable_forward_sum[self._chr] = []
            mappable_reverse_sum = self.ref2mappable_reverse_sum[self._chr] = []
            mascbin = self.ref2mascbins[self._chr] = []

            forward_mappability = mappability
            reverse_mappability = mappability.clone()
            reverse_mappability_buff = reverse_mappability[:self.read_len - 1]
            reverse_mappability.rshift(self.read_len - 1, 0)

        # calc cross-correlation
        self._logging_info("Calculate cross-correlation for {}...".format(self._chr))
        self._init_bitarray_progress()

        for i in xrange(self.max_shift + 1):
            # calc mappability-sensitive cross-correlation
            if mappability:
                doubly_mappable_region.alloc_and(forward_mappability, reverse_mappability)
                if i < self.read_len:
                    mappable_len[self.read_len - i - 1] = doubly_mappable_region.count()
                elif i < self.read_len * 2 - 1:
                    # assert mappable_len[i - self.read_len + 1] == doubly_mappable_region.count()
                    pass
                else:
                    mappable_len.append(doubly_mappable_region.count())

                mappable_forward.alloc_and(self._forward_array, doubly_mappable_region)
                mappable_reverse.alloc_and(self._reverse_array, doubly_mappable_region)

                mappable_forward_sum.append(mappable_forward.count())
                mappable_reverse_sum.append(mappable_reverse.count())
                mascbin.append(mappable_forward.acount(mappable_reverse))

                reverse_mappability.lshift(
                    1,
                    reverse_mappability_buff.pop() if reverse_mappability_buff else 0
                )

            # calc naive cross-correlation
            if not self.skip_ncc:
                ccbin.append(self._forward_array.acount(self._reverse_array))

            self._reverse_array.rshift(1, 0)
            self._progress.update(i)

        self._progress.clean()

    cdef inline bitarray _load_mappability(self):
        cdef bitarray mappability
        cdef BigWigFile feeder
        cdef float _val

        if not self._bwfeeder:
            return None

        # KeyError raises if there is no mappability tracks for _chr
        feeder = self._bwfeeder.fetch(self.MAPPABILITY_THRESHOLD, self._chr)
        self._logging_info("Loading {} mappability to bit array...".format(self._chr))

        # Allocate same size as forward/reverse read array
        chromsize = self.ref2genomelen[self._chr] + self._array_extend_size
        mappability = bitarray(chromsize)
        self._init_mappability_progress()
        for begin, end, _val in feeder:
            mappability.set(begin + 1, end)
            self._progress.update(end)
        self._progress.update(self.ref2genomelen[self._chr])
        self._progress.clean()

        return mappability

    cdef inline _check_pos(self, str chrom, int64 pos):
        if chrom != self._chr:
            if self._chr != '':
                if chrom in self._solved_chr:
                    raise ReadUnsortedError
                self._solved_chr.append(self._chr)
                self.flush()
                self._buff_flashed = False
            self._chr = chrom
            self._init_buff()
            self._logging_info("Loading {} reads to bit array...".format(chrom))

        if pos < self._last_pos:
            raise ReadUnsortedError

        self._last_pos = pos

    @wraparound(False)
    @boundscheck(False)
    def feed_forward_read(self, str chrom, int64 pos, int64 readlen):
        self._check_pos(chrom, pos)

        if self._last_forward_pos == pos:
            return None  # duplicated read
        self._last_forward_pos = pos

        self.forward_read_len_sum += readlen
        self._forward_array[pos] = 1

    @wraparound(False)
    @boundscheck(False)
    def feed_reverse_read(self, str chrom, int64 pos, int64 readlen):
        self._check_pos(chrom, pos)

        if self._reverse_array[pos + readlen - 1] == 0:
            self._reverse_array[pos + readlen - 1] = 1
            self.reverse_read_len_sum += readlen

    def finishup_calculation(self):
        cdef bitarray mappability, forward_mappability, reverse_mappability
        cdef int i

        self.flush()

        # calc mappability for references that there is no reads aligned
        if not self._bwfeeder:
            return None

        for chrom, dmlength in self.ref2mappable_len.items():
            if dmlength is not None:
                continue

            self._chr = chrom
            try:
                mappability = self._load_mappability()
            except KeyError:
                continue

            mappable_len = self.ref2mappable_len[self._chr] = []
            forward_mappability = mappability
            reverse_mappability = mappability.clone()

            self._logging_info("Calc {} mappable length...".format(self._chr))
            for i in xrange(self.max_shift + 1):
                mappable_len.append(forward_mappability.acount(reverse_mappability))
                reverse_mappability.rshift(1, 0)
