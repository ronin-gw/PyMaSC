import fnmatch
from itertools import groupby, chain
import logging

import numpy as np

logger = logging.getLogger(__name__)


def moving_avr_filter(arr, window):
    f = np.repeat(1, window) / float(window)
    avr = np.correlate(arr, f, mode="same")
    h_w = window // 2
    for i in range(h_w):
        avr[i] = np.average(arr[0:(h_w + i)])
        avr[-(i + 1)] = np.average(arr[-(h_w + i):])
    return avr


def filter_chroms(chroms, filters):
    if filters is None:
        logger.debug("There is no chromosome filters.")
        return chroms

    logger.debug("Filtering chromosome lists: " + repr(chroms))
    chroms = set(chroms)
    include_chroms = set()

    for i, (to_include, values) in enumerate(groupby(filters, key=lambda f: f[0])):
        patterns = set(chain(*(f[1] for f in values)))
        logger.debug("{} filters: {}".format("Include" if to_include else "Exclude", repr(patterns)))

        filtered_chroms = set.union(*(set(fnmatch.filter(chroms, p)) for p in patterns))
        logger.debug("Matches: " + repr(filtered_chroms))

        if not to_include:
            include_chroms |= chroms - filtered_chroms
        chroms = filtered_chroms

        logger.debug("current include chroms: " + repr(include_chroms))

    if to_include:
        include_chroms |= chroms

    logger.debug("Include chroms: " + repr(include_chroms))
    return include_chroms


class exec_worker_pool(object):
    def __init__(self, workers, tasks, task_queue):
        self.workers = workers
        self.tasks = tasks
        self.task_queue = task_queue

    def __enter__(self):
        for w in self.workers:
            w.start()
        for t in self.tasks:
            self.task_queue.put(t)
        for _ in range(len(self.workers)):
            self.task_queue.put(None)

    def __exit__(self, type, value, traceback):
        for w in self.workers:
            if w.is_alive():
                w.terminate()
