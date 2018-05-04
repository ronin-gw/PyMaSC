import os
import sys
from functools import wraps


def catch_IOError(logger):
    def _inner(func):
        @wraps(func)
        def _io_func(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except IOError as e:
                logger.error("Faild to output '{}':\n[Errno {}] {}".format(
                    e.filename, e.errno, e.message)
                )
        return _io_func
    return _inner


def prepare_outdir(outdir, logger):
    if os.path.exists(outdir):
        if not os.path.isdir(outdir):
            logger.critical("Specified path as a output directory is not directory.")
            logger.critical(outdir)
            return False
    else:
        logger.info("Make output directory: {}".format(outdir))
        try:
            os.mkdir(outdir)
        except IOError as e:
            logger.critical("Faild to make output directory: [Errno {}] {}".format(e.errno, e.message))
            return False

    if not os.access(outdir, os.W_OK):
        logger.critical("Output directory '{}' is not writable.".format(outdir))
        return False

    return True
