import numpy as np


def moving_avr_filter(arr, window):
    f = np.repeat(1, window) / float(window)
    return np.correlate(arr, f, mode="same")
