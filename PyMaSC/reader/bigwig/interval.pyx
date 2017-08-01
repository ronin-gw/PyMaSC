#
# KentLib bigWigIntervalQuery iterable implementation by Cython
#
from libc.stdlib cimport free

from PyMaSC.reader.bigwig.kentlib.types cimport bits32, bits64
from PyMaSC.reader.bigwig.kentlib.files cimport bbiFile, cirTreeFile


cdef extern from "sig.h":
    int bigWigSig


cdef extern from "common.h":
    void *needLargeMem(size_t size)
    void freeMem(void *pt)
    void slFreeList(void *listPt)
    bits32 memReadBits32(
        char **pPt,
        boolean isSwapped
    )
    float memReadFloat(
        char **pPt,
        boolean isSwapped
    )
    void fileOffsetSizeFindGap(
        fileOffsetSize *list,
    	fileOffsetSize **pBeforeGap,
        fileOffsetSize **pAfterGap
    )


cdef extern from "udc.h":
    void udcSeek(
        udcFile *file,
        bits64 offset
    )
    void udcMustRead(
        udcFile *file,
        void *buf,
        bits64 size
    )


cdef extern from "bwgInternal.h":
    enum bwgSectionType:
        bwgTypeBedGraph = 1
        bwgTypeVariableStep = 2
        bwgTypeFixedStep = 3

    struct bwgSectionHead:
        pass

    void bwgSectionHeadFromMem(
        char **pPt,
        bwgSectionHead *head,
        boolean isSwapped
    )


cdef extern from "bbiFile.h":
    fileOffsetSize *bbiOverlappingBlocks(
        bbiFile *bbi,
        cirTreeFile *ctf,
    	char *chrom,
        bits32 start,
        bits32 end,
        bits32 *retChromId
    )
    void bbiAttachUnzoomedCir(bbiFile *bbi)


cdef extern from "zlibFace.h":
    size_t zUncompress(
        void *compressed,
    	size_t compressedSize,
    	void *uncompBuf,
    	size_t uncompBufSize
    )


cdef class BWIntervalGenerator:
    cdef init(self, bbiFile *bigwig, char *chrom, bits32 start, bits32 end):

        self.file = bigwig
        if (self.file.typeSig != bigWigSig):
            raise IOError("Trying to do bigWigIntervalQuery on a non big-wig file.")

        self.start = start
        self.end = end

        bbiAttachUnzoomedCir(self.file)
        self.blockList = bbiOverlappingBlocks(
            self.file, self.file.unzoomedCir, chrom, start, end, NULL
        )

        self.udc = self.file.udc
        self.isSwapped = self.file.isSwapped


        if (self.file.uncompressBufSize > 0):
            self.uncomp_buff = <char *>needLargeMem(self.file.uncompressBufSize)

        self.block = self.blockList
        self._init_contigblocks()
        self._init_block()
        self.yield_count = 0


    cdef _init_contigblocks(self):
        if self.block == NULL:
            raise StopIteration

        fileOffsetSizeFindGap(self.block, &self.beforeGap, &self.afterGap)
        cdef bits64 mergedOffset = self.block.offset
        cdef bits64 mergedSize = self.beforeGap.offset + self.beforeGap.size - mergedOffset
        udcSeek(self.udc, mergedOffset)

        self.mergedBuf = <char *>needLargeMem(mergedSize)
        udcMustRead(self.udc, self.mergedBuf, mergedSize)
        self.blockBuf = self.mergedBuf

    cdef int _init_block(self):
        if self.block == self.afterGap:
            return 1

        cdef int uncomp_size
        if self.uncomp_buff:
            self.block_point = self.uncomp_buff
            uncomp_size = zUncompress(self.blockBuf, self.block.size, self.uncomp_buff, self.file.uncompressBufSize)
            self.blockEnd = self.block_point + uncomp_size
        else:
            self.block_point = self.blockBuf
            self.blockEnd = self.block_point + self.block.size

        bwgSectionHeadFromMem(&self.block_point, &self.head, self.isSwapped)
        return 0

    cdef Interval _iter_block(self):
        cdef Interval el
        el.start = 0
        el.end = 0

        if self.head.type == bwgTypeBedGraph:
            while self.yield_count < self.head.itemCount:
                self.yield_count += 1
                self.s = memReadBits32(&self.block_point, self.isSwapped)
                self.e = memReadBits32(&self.block_point, self.isSwapped)
                val = memReadFloat(&self.block_point, self.isSwapped)

                if self.s < self.start:
                    self.s = self.start
                if self.e > self.end:
                    self.e = self.end
                if self.s < self.e:
                    el.start = self.s
                    el.end = self.e
                    el.val = val
                    return el

        elif self.head.type == bwgTypeVariableStep:
            while self.yield_count < self.head.itemCount:
                self.yield_count += 1
                self.s = memReadBits32(&self.block_point, self.isSwapped)
                self.e = self.s + self.head.itemSpan
                val = memReadFloat(&self.block_point, self.isSwapped)

                if self.s < self.start:
                    self.s = self.start
                if self.e > self.end:
                    self.e = self.end
                if self.s < self.e:
                    el.start = self.s
                    el.end = self.e
                    el.val = val
                    return el

        elif self.head.type == bwgTypeFixedStep:
            self.s = self.head.start
            self.e = self.s + self.head.itemSpan

            while self.yield_count < self.head.itemCount:
                self.yield_count += 1
                val = memReadFloat(&self.block_point, self.isSwapped)
                self.clippedS = self.s
                self.clippedE = self.e

                if self.clippedS < self.start:
                    self.clippedS = self.start
                if self.clippedE > self.end:
                    self.clippedE = self.end
                if self.clippedS < self.clippedE:
                    el.start = self.clippedS
                    el.end = self.clippedE
                    el.val = val
                    return el

                self.s += self.head.itemStep
                self.e += self.head.itemStep
        else:
            raise RuntimeError("Unknown bigWig section head type.")

        self.yield_count = 0
        return el

    cdef _update_block(self):
        assert self.block_point == self.blockEnd
        self.blockBuf += self.block.size
        self.block = self.block.next

    cdef _clean_contig(self):
        if self.mergedBuf != NULL:
            free(self.mergedBuf)
            self.mergedBuf = NULL

    cdef _clean_up(self):
        cdef fileOffsetSize *next

        if self.uncomp_buff != NULL:
            freeMem(self.uncomp_buff)
        if self.blockList != NULL:
            slFreeList(&self.blockList)

    def __iter__(self):
        return self

    def __next__(self):
        cdef Interval result

        while True:
            result = self._iter_block()
            if not (result.start == 0 and result.end == 0):
                return result.start, result.end, result.val
            self._update_block()
            while True:
                if self._init_block() == 0:
                    break
                self._clean_contig()
                self._init_contigblocks()

    def __dealloc__(self):
        self._clean_contig()
        self._clean_up()
