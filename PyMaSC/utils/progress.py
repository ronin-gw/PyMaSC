"""Terminal progress bar and multiline progress display utilities.

Provides progress reporting for PyMaSC analysis operations. Supports both
single-line progress bars and multiline displays for parallel processing.

Key components:
- ProgressBase: Base class with global enable/disable control
- ProgressBar: Single-line terminal progress bar implementation
- ProgressHook: Progress reporting for multiprocessing environments
- MultiLineProgressManager: Concurrent progress display for multiple operations
- ReadCountProgressBar: Specialized progress bar for read counting operations

The progress system automatically detects terminal capabilities and can be
globally disabled for non-interactive environments or when output is redirected.
"""
from __future__ import annotations

import sys
import array
import fcntl
from abc import ABC, abstractmethod
from multiprocessing.queues import Queue
from termios import TIOCGWINSZ
from typing import Any, Dict, TextIO, Union, TypeVar, Generic, Callable


class ProgressBase(ABC):
    """Base class for progress indicators.

    Provides global enable/disable control for all progress indicators
    and defines the interface for progress reporting classes.

    Attributes:
        global_switch: Class-level flag to enable/disable all progress indicators
    """

    global_switch: bool = sys.stderr.isatty()

    @classmethod
    def _pass(cls, *args: Any, **kwargs: Any) -> None:
        pass

    @abstractmethod
    def enable_bar(self) -> None:
        pass

    @abstractmethod
    def disable_bar(self) -> None:
        pass


OutputT = TypeVar("OutputT", TextIO, Queue)


class SingleProgressBar(ProgressBase, Generic[OutputT], ABC):
    """Abstract base class for single-line progress bars.

    Provides common functionality for progress bars that display on a single line,
    with support for different output targets (TextIO or Queue).
    """
    format: Callable[[str], None]
    clean: Callable[[], None]
    update: Callable[[Union[int, float]], None]

    def _init_config(self, body: str, prefix: str, suffix: str) -> None:
        """Initialize progress bar configuration and enable/disable based on global switch.

        Args:
            body: Progress bar body string
            prefix: String prefix for the progress bar
            suffix: String suffix for the progress bar
        """
        self._init_format(body, prefix, suffix)
        if self.global_switch:
            self.enable_bar()
        else:
            self.disable_bar()

    def _init_format(self, body: str, prefix: str, suffix: str) -> None:
        """Initialize the format for the progress bar."""
        self.body = body
        self.fmt = "\r" + prefix + "{:<" + str(len(body)) + "}" + suffix

    def reset_progress(self, body: str, prefix: str, suffix: str, name: str, maxval: Union[int, float]) -> None:
        """Reset progress bar with new configuration and maximum value.

        Args:
            body: New progress bar body string
            prefix: New string prefix
            suffix: New string suffix
            name: Progress operation name
            maxval: Maximum value for progress completion
        """
        self._init_format(body, prefix, suffix)
        self.set(name, maxval)

    def enable_bar(self) -> None:
        """Enable progress bar display if global switch allows."""
        if self.global_switch:
            self.format = self._format
            self.clean = self._clean
            self.update = self._update

    def disable_bar(self) -> None:
        """Disable progress bar display by setting methods to no-ops."""
        self.format = self.clean = self.update = self._pass

    @abstractmethod
    def _format(self, s: str) -> None:
        pass

    @abstractmethod
    def _clean(self) -> None:
        pass

    def set(self, name: str, maxval: Union[int, float]) -> None:
        """Set progress parameters and initialize progress tracking.

        Args:
            name: Name of the progress operation
            maxval: Maximum value representing 100% completion
        """
        self._unit = float(maxval) / len(self.body)
        self.pos = 0
        self._next_update = self._unit

    def _update(self, val: Union[int, float]) -> None:
        if val > self._next_update:
            while val > self._next_update:
                self.pos += 1
                self._next_update += self._unit
            self.format(self.body[:self.pos])


class ProgressBar(SingleProgressBar[TextIO]):
    """Single-line progress bar implementation.

    Displays a visual progress bar in the terminal with customizable appearance.
    The bar automatically updates as progress is made and can be cleanly
    removed when operations complete.

    Attributes:
        body: String defining the progress bar appearance
        fmt: Format string for displaying the progress bar
        output: Output stream for progress display (default: stderr)
    """
    def __init__(
        self,
        output: TextIO = sys.stderr,
        body: str = "<1II1>" * 12,
        prefix: str = ">",
        suffix: str = "<",
    ) -> None:
        """Initialize progress bar.

        Args:
            output: Output stream for progress display
            body: Visual pattern for the progress bar
            prefix: Text displayed before the progress bar
            suffix: Text displayed after the progress bar
        """
        self._init_config(body, prefix, suffix)
        self.output = output

    def _format(self, s: str) -> None:
        self.output.write(self.fmt.format(s))

    def _clean(self) -> None:
        self.output.write("\r\033[K")
        self.output.flush()


class ProgressHook(SingleProgressBar[Queue]):
    """Progress reporting for multiprocessing.

    Extends ProgressBar to support progress reporting from worker processes
    back to the main process through inter-process communication queues.
    """
    def __init__(
        self,
        output: Queue,
        body: str = "<1II1>" * 12,
        prefix: str = ">",
        suffix: str = "<",
    ) -> None:
        self._init_config(body, prefix, suffix)
        self.output = output

    def _clean(self) -> None:
        pass

    def set(self, name: str, maxval: Union[int, float]) -> None:
        self.name = name
        super(ProgressHook, self).set(name, maxval)

    def _format(self, s: str) -> None:
        self.output.put((
            None, (self.name, self.fmt.format(s))
        ))


class MultiLineProgressManager(ProgressBase):
    """Multi-line progress display manager.

    Manages concurrent progress displays for multiple operations, allowing
    each operation to have its own progress line that can be independently
    updated without interfering with others.

    This class is essential for multiprocessing scenarios where multiple
    chromosomes are being processed simultaneously and each needs its own
    progress indicator.

    Attributes:
        _terminal_width: Detected terminal width for proper formatting
        _status: Dictionary tracking status of each progress line
    """
    erase: Callable[[str], None]
    clean: Callable[[], None]
    update: Callable[[str, str], None]

    def __init__(self, output: TextIO = sys.stderr) -> None:
        """Initialize multiline progress manager.

        Args:
            output: Output stream for progress display
        """
        #
        self.output = output
        if not self.output.isatty():
            self.global_switch = False

        #
        if not self.global_switch:
            self.disable_bar()
            return None
        else:
            self.enable_bar()

        # get terminal size
        buf = array.array('H', ([0] * 4))
        stdfileno = self.output.fileno()
        fcntl.ioctl(stdfileno, TIOCGWINSZ, buf, True)
        self.max_height, self.max_width = buf[:2]
        self.max_width -= 1
        self.key2lineno: Dict[str, int] = {}
        self.lineno2key: Dict[int, str] = {}
        self.key2body: Dict[str, str] = {}
        self.nlines = 0

    def enable_bar(self) -> None:
        if self.global_switch:
            self.erase = self._erase
            self.clean = self._clean
            self.update = self._update

    def disable_bar(self) -> None:
        self.erase = self.clean = self.update = self._pass

    def _cursor_n(self, n: int, ctrl: str) -> None:
        if n < 1:
            return None
        self._write(f"\033[{n}{ctrl}")

    def _down(self, n: int) -> None:
        self._cursor_n(n, "E")

    def _up(self, n: int) -> None:
        self._cursor_n(n, "F")

    def _reset_line(self) -> None:
        self._write("\033[K")

    def _write(self, line: str) -> None:
        self.output.write(line[: self.max_width])
        self.output.flush()

    def _refresh_lines(self, from_: int, to: int) -> None:
        for i in range(from_, to):
            k = self.lineno2key[i]
            self._reset_line()
            self._write(f"{self.key2body[k]} {k}")
            if i < self.nlines:
                self._write("\n")

    def _update(self, key: str, body: str) -> None:
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

    def _erase(self, key: str) -> None:
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

    def _clean(self) -> None:
        self._reset_line()
        for i in range(self.nlines - 1):
            self._down(1)
            self._reset_line()
        for i in range(self.nlines - 1):
            self._up(1)


class ReadCountProgressBar(MultiLineProgressManager):
    """Specialized progress bar for genome wide read counting operations.

    Combines single-line progress bar functionality with multiline management
    for read counting operations that may involve multiple chromosomes or
    processing stages.

    This class provides read-specific progress formatting and handles the
    coordination between read counting stages and chromosome processing.
    """
    set_genome: Callable[[int], None]
    set_chrom: Callable[[str, int], None]
    set: Callable[[str, int], None]

    def __init__(
        self,
        g_body: str = "^@@@@@@@@@" * 10,
        g_prefix: str = "",
        g_suffix: str = "^",
        c_body: str = "<1II1>" * 12,
        c_prefix: str = ">",
        c_suffix: str = "< {}",
        output: TextIO = sys.stderr,
    ) -> None:
        MultiLineProgressManager.__init__(self, output)

        self.genome_body = g_body
        self.genome_fmt = g_prefix + "{:<" + str(len(g_body)) + "}" + g_suffix
        self.chrom_body = c_body
        self.chrom_fmt = c_prefix + "{:<" + str(len(c_body)) + "}" + c_suffix

        if self.global_switch:
            self.enable_bar()
        else:
            self.disable_bar()

        self.output = output
        self._genome_offset = 0

    def _update(self, *args: Any, **kwargs: Any) -> None:
        """Compatible update method that handles both parent class signatures.

        Can be called as:
        - update(val) for ProgressBar compatibility
        - update(key, body) for MultiLineProgressManager compatibility
        """
        if len(args) == 1:
            self._update_wholegenome(args[0])
        elif len(args) == 2:
            # MultiLineProgressManager style: update(key, body)
            super()._update(args[0], args[1])
        else:
            self._pass(*args, **kwargs)

    def enable_bar(self) -> None:
        if self.global_switch:
            self.set = self.set_chrom = self._set_chrom
            self.set_genome = self._set_genome
            self.update = self._update
            self.finish = self._finish

    def disable_bar(self) -> None:
        self.set = self.set_chrom = self.set_genome = self.update = self.finish = self._pass

    def _set_chrom(self, name: str, maxval: int) -> None:
        self._chrom = name
        self._chrom_maxval = maxval
        if self._genome_offset != 0:
            self._genome_offset += self._chrom_maxval
        self._chrom_unit = float(maxval) / len(self.chrom_body)
        self.chrom_pos = 0
        self._chrom_next_update = self._chrom_unit
        self._reset_line()
        self._write(self.chrom_fmt.format("", self._chrom))
        self._write("\n")
        self._reset_line()
        self._write(self.genome_fmt.format(self.genome_body[: self.genome_pos]))
        self._up(1)

    def _set_genome(self, maxval: int) -> None:
        self._genome_unit = float(maxval) / len(self.genome_body)
        self.genome_pos = 0
        self._genome_next_update = self._genome_unit

    def _update_wholegenome(self, val: Union[int, float]) -> None:
        if val > self._chrom_next_update:
            while val > self._chrom_next_update:
                self.chrom_pos += 1
                self._chrom_next_update += self._chrom_unit
            self._write(self.chrom_fmt.format(self.chrom_body[: self.chrom_pos], self._chrom))
            self._write("\n")
            if val + self._genome_offset > self._genome_next_update:
                while val + self._genome_offset > self._genome_next_update:
                    self.genome_pos += 1
                    self._genome_next_update += self._genome_unit
                self._write(self.genome_fmt.format(self.genome_body[: self.genome_pos]))
            self._up(1)
            self._write("\033[G")

    def _finish(self) -> None:
        self._down(1)
        self._reset_line()
        self._up(1)
        self._reset_line()
