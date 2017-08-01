import os
import itertools

from libc.stdio cimport FILE
from libc.string cimport strcpy

from .kentlib.types cimport bits32, bits64, boolean
from .kentlib.structs cimport fileOffsetSize
from .kentlib.files cimport udcFile, cirTreeFile, bbiFile
from .interval cimport BWIntervalGenerator

cdef extern from "common.h":
    int FALSE

    FILE *mustOpen(char *fileName, char *mode)
    void carefulClose(FILE **pFile)


cdef extern from "bbiFile.h":
    struct bbiChromInfo:
        bbiChromInfo *next
        char *name
        bits32 id
        bits32 size

    bbiChromInfo *bbiChromList(bbiFile *bbi)
    void bbiFileClose(bbiFile **pBwf)
    void bbiChromInfoFreeList(bbiChromInfo **pList)
    fileOffsetSize *bbiOverlappingBlocks(
        bbiFile *bbi,
        cirTreeFile *ctf,
    	char *chrom,
        bits32 start,
        bits32 end,
        bits32 *retChromId
    )

cdef extern from "bigWig.h":
    bbiFile *bigWigFileOpen(char *fileName)
    int bigWigIntervalDump(
        bbiFile *bwf,
        char *chrom,
        bits32 start,
        bits32 end,
        int maxCount,
        FILE *out
    )
    void bigWigFileCreate(
        char *inName,
        char *chromSizes,
        int blockSize,
        int itemsPerSlot,
        boolean clipDontDie,
        boolean compress,
        char *outName
    )


def _assert_readable_file(path):
    if not (os.path.isfile(path) and os.access(path, os.R_OK)):
        raise IOError("Can't read file '{}'".format(path))


def _assert_writable_file(path):
    if not (os.path.isfile(path) and os.access(path, os.W_OK)):
        raise IOError("Can't write file '{}'".format(path))


def wig2bigwig(path, sizefile, output):
    if not path:
        raise IOError("Failed to convet wig to bigwig. Input file is not specified.")
    elif not sizefile:
        raise IOError("Failed to convet wig to bigwig. 'chrom_size' file required.")
    elif not output:
        raise IOError("Failed to convet wig to bigwig. Output file is not specified.")

    for f in (path, sizefile):
        _assert_readable_file(f)
    _assert_writable_file(output)

    bigWigFileCreate(path, sizefile, 256, 1024, FALSE, FALSE, output)


def bigwig2wig(path, output):
    if not path:
        raise IOError("Failed to convet bigwig to wig. Input file is not specified.")
    elif not output:
        raise IOError("Failed to convet bigwig to wig. Output file is not specified.")

    _assert_readable_file(path)
    _assert_writable_file(output)

    cdef bbiFile *bigwig = bigWigFileOpen(path)
    cdef FILE *wig = mustOpen(output, "w")

    cdef bbiChromInfo *chrom
    cdef bbiChromInfo *chromList = bbiChromList(bigwig)
    chrom = chromList
    while chrom != NULL:
        bigWigIntervalDump(bigwig, chrom.name, 0, chrom.size, 0, wig)
        chrom = chrom.next

    bbiChromInfoFreeList(&chromList)
    carefulClose(&wig)
    bbiFileClose(&bigwig)


cdef class BigWigFile:
    cdef:
        readonly char path[1024]
        readonly bint closed
        readonly dict chromsizes

        bbiFile *file
        bbiChromInfo *chrominfo

    def __init__(self, path, chrom_size=None):
        if not os.path.exists(path) and path != '-':
            raise IOError("input file '{0}' dose not exist.".format(path))
        elif path == '-':
            path = "stdin"

        cdef char *path_tmp
        prefix, ext = os.path.splitext(path)
        if ext == ".wig":
            bwfile = prefix + ".bw"
            if os.path.exists(bwfile):
                path_tmp = bwfile
            else:
                wig2bigwig(path, chrom_size, bwfile)
            path_tmp = bwfile
        else:
            path_tmp = path

        strcpy(self.path, path_tmp)

        self._open()
        self.closed = False

        self._set_chromsizes()

    def _open(self):
        self.file = bigWigFileOpen(self.path)
        self.chrominfo = bbiChromList(self.file)

    def _set_chromsizes(self):
        self.chromsizes = {}
        cdef bbiChromInfo *chrom = self.chrominfo

        while chrom != NULL:
            self.chromsizes[chrom.name] = chrom.size
            chrom = chrom.next

    def fetch(self, chrom=None, begin=None, end=None):
        if chrom is None:
            chroms = self.chromsizes.keys()
        else:
            if chrom not in self.chromsizes:
                raise AttributeError("Reference name '{}' not found.".format(chrom))
            chroms = [chrom]
        if begin is None or begin < 0:
            begin = 0

        _begin = begin
        iters = []
        for chrom in chroms:
            if end is None:
                _end = self.chromsizes[chrom]

            i = BWIntervalGenerator()
            i.init(self.file, chrom, _begin, _end)
            iters.append(i)

        return itertools.chain.from_iterable(iters)


    def close(self):
        if not self.closed:
            if self.file != NULL:
                bbiFileClose(&self.file)
            if self.chrominfo != NULL:
                bbiChromInfoFreeList(&self.chrominfo)
            self.closed = True

    def __dealloc__(self):
        self.close()
