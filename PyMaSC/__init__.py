import sys
from functools import wraps

VERSION = "0.2.4"


def logging_version(logger):
    logger.info("PyMaSC version {} with Python{}.{}.{}".format(
                *[VERSION] + list(sys.version_info[:3])))
    for line in sys.version.split('\n'):
        logger.debug(line)


def entrypoint(logger):
    def _entrypoint_wrapper_base(main_func):
        @wraps(main_func)
        def _inner():
            try:
                main_func()
                logger.info("PyMASC plot finished.")
            except KeyboardInterrupt:
                sys.stderr.write("\r\033[K")
                sys.stderr.flush()
                logger.info("Got KeyboardInterrupt. bye")
        return _inner
    return _entrypoint_wrapper_base
