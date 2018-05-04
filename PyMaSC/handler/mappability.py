import logging
import os
import json
from multiprocessing import Process, Queue, Lock

import numpy as np

from PyMaSC.core.mappability import MappableLengthCalculator
from PyMaSC.utils.progress import ProgressHook, MultiLineProgressManager
from PyMaSC.utils.compatible import tostr, xrange
from PyMaSC.utils.output import prepare_outdir
from PyMaSC.utils.calc import exec_worker_pool

logger = logging.getLogger(__name__)


class BWIOError(IOError):
    pass


class JSONIOError(IOError):
    pass


class NeedUpdate(Exception):
    pass


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.long, np.float, np.float_)):
            return float(obj)
        elif isinstance(obj, (np.uint, np.int32, np.int64)):
            return int(obj)
        else:
            return super(self, NumpyEncoder).default(obj)


class MappabilityHandler(MappableLengthCalculator):
    @staticmethod
    def calc_mappable_len_required_shift_size(readlen, max_shift):
        return max_shift - readlen + 1 if max_shift > 2*readlen - 1 else readlen

    def __init__(self, path, max_shift=0, readlen=0, map_path=None, nworker=1):
        max_shift = self.calc_mappable_len_required_shift_size(readlen, max_shift)
        self.nworker = nworker

        if not os.access(path, os.R_OK):
            reason = "no such file." if os.path.isfile(path) else "file is unreadable."
            logger.critical("Failed to open '{}': {}".format(path, reason))
            raise BWIOError

        super(MappabilityHandler, self).__init__(path, max_shift)
        self.close()
        self._progress.disable_bar()
        self.need_save_stats = True

        if map_path:
            self.map_path = map_path
        else:
            self.map_path = os.path.splitext(path)[0] + "_mappability.json"

        if not os.path.exists(self.map_path):
            self._check_saving_directory_is_writable()
            logger.info("Calcurate mappable length with max shift size {}.".format(max_shift))
        elif not os.path.isfile(self.map_path):
            logger.critical("Specified path is not file: '{}'".format(self.map_path))
            raise JSONIOError
        elif not os.access(self.map_path, os.R_OK):
            logger.error("Failed to read '{}'".format(self.map_path))
        else:
            self._try_load_mappability_stats()
            if self.need_save_stats:
                self._check_stats_is_overwritable()
                logger.info("Calcurate mappable length with max shift size {}.".format(max_shift))
            else:
                logger.info("Use mappability stats read from '{}'".format(self.map_path))

    def _check_saving_directory_is_writable(self):
        dirname = os.path.dirname(self.map_path)
        dirname = dirname if dirname else '.'
        if not prepare_outdir(dirname, logger):
            raise JSONIOError

    def _try_load_mappability_stats(self):
        try:
            stats = self._read_mappability_stats()
        except IOError as e:
            logger.error("Failed to read '{}'".format(self.map_path))
            logger.error("[Errno {}] {}".format(e.errno, e.message))
        except (TypeError, OverflowError, ValueError, KeyError, IndexError) as e:
            logger.error("Failed to load json file: '{}'".format(self.map_path))
        except NeedUpdate:
            pass

        self._load_mappability_stats(stats)

    def _read_mappability_stats(self):
        with open(self.map_path) as f:
            stats = json.load(f)

        for k in ("max_shift", "__whole__", "references"):
            if k not in stats:
                logger.error("Mandatory key '{}' not found.".format(k))
                raise KeyError(k)

        if stats["max_shift"] < self.max_shift:
            logger.info("Specified shift length longer than former analysis. The stats will be updated.")
            raise NeedUpdate

        if stats["max_shift"] != len(stats["__whole__"]) - 1:
            logger.error("Max shift length for whole genome unmatched.")
            raise IndexError

        for ref in self.chromsizes:
            if ref not in stats["references"]:
                logger.error("Reference '{}' not found.".format(ref))
                raise KeyError(ref)

            if stats["max_shift"] != len(stats["references"][ref]) - 1:
                logger.error("Max shift length for 'ref' unmatched.".format(ref))
                raise IndexError

        return stats

    def _load_mappability_stats(self, stats):
        self.mappable_len = stats["__whole__"][:self.max_shift + 1]
        self.chrom2mappable_len = {ref: b[:self.max_shift + 1] for ref, b in stats["references"].items()}
        self.chrom2is_called = {ref: True for ref in self.chromsizes}
        self.is_called = True
        self.need_save_stats = False

    def _check_stats_is_overwritable(self):
        if not os.access(self.map_path, os.W_OK):
            logger.critical("Failed to overwrite '{}'".format(self.map_path))
            raise JSONIOError
        else:
            logger.warning("Existing file '{}' will be overwritten.".format(self.map_path))

    def save_mappability_stats(self):
        if not self.need_save_stats:
            return logger.info("Mappability stats updating is not required.")

        logger.info("Save mappable length to '{}'".format(self.map_path))
        try:
            with open(self.map_path, 'w') as f:
                json.dump({
                    "max_shift": self.max_shift,
                    "__whole__": self.mappable_len,
                    "references": self.chrom2mappable_len
                }, f, indent=4, sort_keys=True, cls=NumpyEncoder)
        except IOError as e:
            logger.error("Faild to output: {}\n[Errno {}] {}".format(
                         e.filename, e.errno, e.message))

        self.need_save_stats = False

    def calc_mappability(self):
        target_chroms = [tostr(c) for c, b in self.chrom2is_called.items() if b is False]
        if not target_chroms:
            return self._sumup_mappability()

        order_queue = Queue()
        report_queue = Queue()
        logger_lock = Lock()
        progress = MultiLineProgressManager()

        workers = [MappabilityCalcWorker(self.path, self.max_shift, order_queue, report_queue, logger_lock)
                   for _ in range(min(self.nworker, len(target_chroms)))]

        with exec_worker_pool(workers, target_chroms, order_queue):
            while not self.is_called:
                chrom, obj = report_queue.get()
                if chrom is None:  # update progress
                    chrom, body = obj
                    with logger_lock:
                        progress.update(chrom, body)
                else:
                    length = obj
                    self.chrom2mappable_len[chrom] = tuple(length)
                    self.chrom2is_called[chrom] = True
                    if all(self.chrom2is_called.values()):
                        self.is_called = True
                    with logger_lock:
                        progress.erase(chrom)

        progress.clean()
        self._sumup_mappability()

    def _sumup_mappability(self):
        for length in self.chrom2mappable_len.values():
            for i in xrange(self.max_shift + 1):
                self.mappable_len[i] += length[i]


class MappabilityCalcWorker(Process):
    def __init__(self, path, max_shift, order_queue, report_queue, logger_lock):
        super(MappabilityCalcWorker, self).__init__()

        self.calculator = MappableLengthCalculator(path, max_shift, logger_lock)
        self.calculator._progress.disable_bar()
        self.order_queue = order_queue
        self.report_queue = report_queue
        self.logger_lock = logger_lock
        self.calculator._progress = ProgressHook(report_queue)

    def run(self):
        with self.logger_lock:
            logger.debug("{}: Hello. My pid is {}.".format(self.name, os.getpid()))

        while True:
            chrom = self.order_queue.get()
            if chrom is None:
                break

            with self.logger_lock:
                logger.debug("{}: Process {}...".format(self.name, chrom))
            self.calculator.calc_mappability(chrom)
            self.report_queue.put((chrom, self.calculator.chrom2mappable_len[chrom]))

        with self.logger_lock:
            logger.debug("{}: Goodbye.".format(self.name))
        self.calculator.close()
