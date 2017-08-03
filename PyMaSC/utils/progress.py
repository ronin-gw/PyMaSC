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

        if self.enable:
            self.enable_bar()
        else:
            self.disable_bar()

    def enable_bar(self):
        self.format = self._format
        self.clean = self._clean
        self.update = self._update

    def disable_bar(self):
        self.format = self.clean = self.update = self._pass

    def _format(self, s):
        self.output.write(self.fmt.format(s))

    def _clean(self):
        self.output.write("\r\033[K")
        self.output.flush()

    def set(self, maxval):
        self._unit = maxval / len(self.body)
        self.pos = 0
        self._next_update = self._unit
        self.format('')

    def _update(self, val):
        if val > self._next_update:
            while val > self._next_update:
                self.pos += 1
                self._next_update += self._unit
            self.format(self.body[:self.pos])
