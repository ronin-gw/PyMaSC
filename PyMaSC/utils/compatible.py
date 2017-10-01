import sys


if sys.version_info > (3, ):
    # def tostr(s):
    #     return bytes(s, "utf-8")
    from itertools import zip_longest
    from io import BytesIO as StringIO
else:
    # tostr = str
    from itertools import izip_longest as zip_longest
    from cStringIO import StringIO as StringIO
