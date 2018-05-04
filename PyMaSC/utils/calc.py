import numpy as np


def moving_avr_filter(arr, window):
    f = np.repeat(1, window) / float(window)
    return np.correlate(arr, f, mode="same")


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
