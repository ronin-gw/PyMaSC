def catch_IOError(_logger, label=''):
    errmsg = "Faild to output" + ((' ' + label) if label else '') + ": {}:\n[Errno {}] {}"

    def _inner(func):
        def _main(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except IOError as e:
                _logger.error(errmsg.format(e.filename, e.errno, e.message))
        return _main

    return _inner
