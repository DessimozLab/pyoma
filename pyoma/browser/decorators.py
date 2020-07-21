from functools import wraps
import logging
from time import time


def timethis(level, name=None):
    """
    Decorator that reports the execution time. level is the logging
    level, name is the logger name. If name isn't specified,
    it default to the function's module

    :Example:
    @timethis(logging.INFO)
    def add(x,y):
        return x+y

    This would log the execution time with level info to the
    module's logger
    """

    def decorate(func):
        logname = name if name else func.__module__
        log = logging.getLogger(logname)

        @wraps(func)
        def wrapper(*args, **kwargs):
            start = time()
            result = func(*args, **kwargs)
            end = time()
            log.log(level, "time for {}: {}sec".format(func.__name__, end - start))
            return result

        return wrapper

    return decorate
