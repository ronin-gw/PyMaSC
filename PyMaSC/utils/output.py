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
