import sys


class ProgressBar(object):
    enable = False

    @classmethod
    def _pass(self, *args, **kwargs):
        pass

    def __init__(self, body="<1II1>" * 12, prefix='>', suffix='<', output=sys.stderr):
        self.body = body
        self.fmt = "\r" + prefix + "{:<" + str(len(body)) + "}" + suffix
        self.output = output

        if not self.enable:
            self.format = self.clean = self.set = self.update = self._pass

    def format(self, s):
        self.output.write(self.fmt.format(s))

    def clean(self):
        self.output.write("\r\033[K")
        self.output.flush()

    def set(self, maxval):
        self._unit = maxval / len(self.body)
        self.pos = 0
        self._next_update = self._unit
        self.format('')

    def update(self, val):
        if val > self._next_update:
            while val > self._next_update:
                self.pos += 1
                self._next_update += self._unit
            self.format(self.body[:self.pos])
