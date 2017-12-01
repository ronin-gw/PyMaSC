import sys


if sys.version_info > (3, ):
    def tostr(s):
        return s
    from itertools import zip_longest
    from io import BytesIO as StringIO
    xrange = range
else:
    def tostr(s):
        return s.encode("utf-8")
    from itertools import izip_longest as zip_longest
    from cStringIO import StringIO as StringIO
    xrange = xrange
