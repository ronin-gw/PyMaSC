import sys
import array
import fcntl
from termios import TIOCGWINSZ
from functools import partialmethod


class ProgressBase(object):
    global_switch = False

    @classmethod
    def _pass(self, *args, **kwargs):
        pass


class ProgressBar(ProgressBase):
    def __init__(self, body="<1II1>" * 12, prefix='>', suffix='<', output=sys.stderr):
        self.body = body
        self.fmt = "\r" + prefix + "{:<" + str(len(body)) + "}" + suffix
        self.output = output

        if self.global_switch:
            self.enable_bar()
        else:
            self.disable_bar()

    def enable_bar(self):
        if self.global_switch:
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

    def set(self, name, maxval):
        self._unit = float(maxval) / len(self.body)
        self.pos = 0
        self._next_update = self._unit
        self.format('')

    def _update(self, val):
        if val > self._next_update:
            while val > self._next_update:
                self.pos += 1
                self._next_update += self._unit
            self.format(self.body[:self.pos])


class ProgressHook(ProgressBar):
    def __init__(self, queue, body="<1II1>" * 12, prefix='>', suffix='<'):
        super(ProgressHook, self).__init__(body, prefix, suffix, queue)
        self.name = None

    def _clean(self):
        pass

    def set(self, name, maxval):
        self.name = name
        super(ProgressHook, self).set(name, maxval)

    def _format(self, s):
        self.output.put((
            None, (self.name, self.fmt.format(s))
        ))


class MultiLineProgressManager(ProgressBase):
    def __init__(self, output=sys.stderr):
        #
        self.output = output
        if not self.output.isatty():
            self.global_switch = False

        #
        if not self.global_switch:
            self.erase = self.clean = self.update = self._pass
            return None

        # get terminal size
        buf = array.array('H', ([0] * 4))
        stdfileno = self.output.fileno()
        fcntl.ioctl(stdfileno, TIOCGWINSZ, buf, 1)
        self.max_height, self.max_width = buf[:2]
        self.max_width -= 1

        #
        self.key2lineno = {}
        self.lineno2key = {}
        self.key2body = {}
        self.nlines = 0

    def _cursor_n(self, n, ctrl):
        if n < 1:
            return None
        self._write("\033[{}{}".format(n, ctrl))

    _down = partialmethod(_cursor_n, ctrl="E")
    _up = partialmethod(_cursor_n, ctrl="F")

    def _reset_line(self):
        self._write("\033[K")

    def _write(self, l):
        self.output.write(l[:self.max_width])
        self.output.flush()

    def _refresh_lines(self, from_, to):
        for i in range(from_, to):
            k = self.lineno2key[i]
            self._reset_line()
            self._write("{} {}".format(self.key2body[k], k))
            if i < self.nlines:
                self._write('\n')

    def update(self, key, body):
        if key not in self.key2lineno:
            self.nlines += 1
            lineno = self.key2lineno[key] = self.nlines
            self.lineno2key[lineno] = key
        else:
            lineno = self.key2lineno[key]
        self.key2body[key] = body

        self._refresh_lines(1, self.nlines + 1)

        self._up(self.nlines - 1)
        if self.nlines == 1:
            self._write("\033[G")

    def erase(self, key):
        try:
            lineno = self.key2lineno[key]
        except KeyError:
            return None

        self._refresh_lines(1, lineno)

        if lineno == self.nlines:
            self._reset_line()

        for i in range(lineno + 1, self.nlines + 1):
            k = self.lineno2key[i - 1] = self.lineno2key[i]
            self.key2lineno[k] -= 1

            self._write("{} {}".format(self.key2body[k], k))
            self._write('\n')

        self.nlines -= 1
        self._reset_line()
        self._up(self.nlines)

        if self.nlines == 1:
            self._write("\033[G")

        del self.key2lineno[key], self.key2body[key]

    def clean(self):
        self._reset_line()
        for i in range(self.nlines - 1):
            self._down(1)
            self._reset_line()
        for i in range(self.nlines - 1):
            self._up(1)


class ReadCountProgressBar(ProgressBar, MultiLineProgressManager):
    def __init__(self, g_body="^@@@@@@@@@" * 10, g_prefix='', g_suffix='^',
                 c_body="<1II1>" * 12, c_prefix='>', c_suffix='< {}', output=sys.stderr):
        MultiLineProgressManager.__init__(self, output)
        if not self.global_switch:
            self.set_chrom = self.set_genome = self.update = self.finish = self._pass
            return None

        self.genome_body = g_body
        self.genome_fmt = g_prefix + "{:<" + str(len(g_body)) + "}" + g_suffix
        self.chrom_body = c_body
        self.chrom_fmt = c_prefix + "{:<" + str(len(c_body)) + "}" + c_suffix

        if self.global_switch:
            self.enable_bar()
        else:
            self.disable_bar()

        self.output = output
        self._genome_offset = None

    def enable_bar(self):
        if self.global_switch:
            self.set_chrom = self._set_chrom
            self.set_genome = self._set_genome
            self.finish = self._finish
            self.update = self._update

    def disable_bar(self):
        self.set_chrom = self.set_genome = self.finish = self.update = self._pass

    def _set_chrom(self, maxval, name):
        if self._genome_offset is None:
            self._genome_offset = 0
        else:
            self._genome_offset += self._chrom_maxval
        self._chrom = name
        self._chrom_maxval = maxval
        self._chrom_unit = float(maxval) / len(self.chrom_body)
        self.chrom_pos = 0
        self._chrom_next_update = self._chrom_unit

        self._reset_line()
        self._write(self.chrom_fmt.format('', self._chrom))
        self._write('\n')
        self._reset_line()
        self._write(self.genome_fmt.format(self.genome_body[:self.genome_pos]))
        self._up(1)

    def _set_genome(self, maxval):
        self._genome_unit = float(maxval) / len(self.genome_body)
        self.genome_pos = 0
        self._genome_next_update = self._genome_unit

    def _update(self, val):
        if val > self._chrom_next_update:
            while val > self._chrom_next_update:
                self.chrom_pos += 1
                self._chrom_next_update += self._chrom_unit
            self._write(self.chrom_fmt.format(self.chrom_body[:self.chrom_pos], self._chrom))
            self._write('\n')

            if val + self._genome_offset > self._genome_next_update:
                while val + self._genome_offset > self._genome_next_update:
                    self.genome_pos += 1
                    self._genome_next_update += self._genome_unit
                self._write(self.genome_fmt.format(self.genome_body[:self.genome_pos]))

            self._up(1)
            self._write("\033[G")

    def _finish(self):
        self._down(1)
        self._reset_line()
        self._up(1)
        self._reset_line()
