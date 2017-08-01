from libc.time cimport time_t

from .types cimport boolean, bits16, bits32, bits64


# from udc.c
cdef struct connInfo:
    int socket
    bits64 offset
    int ctrlSocket

cdef struct udcRemoteFileInfo:
    bits64 updateTime
    bits64 size
    connInfo ci

ctypedef int (*UdcDataCallback)(char *url, bits64 offset, int size, void *buffer, connInfo *ci)
ctypedef boolean (*UdcInfoCallback)(char *url, udcRemoteFileInfo *retInfo)

cdef struct udcProtocol:
    udcProtocol *next
    UdcDataCallback fetchData
    UdcInfoCallback fetchInfo

cdef struct udcBitmap:
    udcBitmap *next
    bits32 blockSize
    bits64 remoteUpdate
    bits64 fileSize
    bits32 version
    bits64 localUpdate
    bits64 localAccess
    boolean isSwapped
    int fd


cdef extern from "udc.h":
    struct udcFile:
        udcFile *next
        char *url
        char *protocol
        udcProtocol *prot
        time_t updateTime
        bits64 size
        bits64 offset
        char *cacheDir
        char *bitmapFileName
        char *sparseFileName
        int fdSparse
        boolean sparseReadAhead
        char *sparseReadAheadBuf
        bits64 sparseRAOffset
        udcBitmap *bits
        bits64 startData
        bits64 endData
        bits32 bitmapVersion
        connInfo connInfo


cdef extern from "bPlusTree.h":
    struct bptFile:
        pass


cdef extern from "cirTree.h":
    struct cirTreeFile:
        pass


cdef extern from "bbiFile.h":
    struct bbiZoomLevel:
        pass

    struct bbiFile:
        bbiFile *next
        char *fileName
        udcFile *udc
        bits32 typeSig
        boolean isSwapped
        bptFile *chromBpt
        bits16 version
        bits16 zoomLevels
        bits64 chromTreeOffset
        bits64 unzoomedDataOffset
        bits64 unzoomedIndexOffset
        bits16 fieldCount
        bits16 definedFieldCount
        bits64 asOffset
        bits64 totalSummaryOffset
        bits32 uncompressBufSize
        cirTreeFile *unzoomedCir
        bbiZoomLevel *levelList
