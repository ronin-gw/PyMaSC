class SamStreamReadGenerator(object):
    def __iter__(self):
        return self

    def __init__(self, stream):
        self.stream = stream
        self._is_header_parsed = False
        self._buff = None
        self.references = []
        self.lengths = []

        self._parse_header()

    def _parse_header(self):
        for line in self:
            if not line.startswith('@'):
                self._buff = line
                break
            elif line.startswith("@SQ"):
                _, refname, length = line.split('\t')
                self.references.append(refname[3:])
                self.lengths.append(int(length[3:]))
        self.references = tuple(self.references)
        self.lengths = tuple(self.lengths)

        self._is_header_parsed = True

    def _next_source(self):
        return next(self.stream)

    def next(self):
        if not self._is_header_parsed:
            return self._next_source()

        if self._buff is None:
            raise StopIteration

        _, flag, chrom, pos, mapq, _, _, _, seq, _ = self._buff.split('\t', 9)
        try:
            self._buff = self._next_source()
        except StopIteration:
            self._buff = None

        flag = int(flag)
        return flag & 16, flag & 1024, chrom, int(pos), int(mapq), len(seq)


class SamReadGenerator(SamStreamReadGenerator, file):
    def __init__(self, *args, **kwargs):
        file.__init__(self, *args, **kwargs)
        SamStreamReadGenerator.__init__(self, self)

    def _next_source(self):
        return file.next(self)
