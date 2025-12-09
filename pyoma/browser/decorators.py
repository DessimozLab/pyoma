from functools import wraps
import logging
from time import time, perf_counter, process_time


def timethis(level=logging.INFO, name=None):
    """
    Decorator that reports the execution time. level is the logging
    level, name is the logger name. If name isn't specified,
    it default to the function's module

    :Example:
    @timethis(level=logging.INFO)
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
            start_wall = perf_counter()
            start_cpu = process_time()
            result = func(*args, **kwargs)
            end_wall = perf_counter()
            end_cpu = process_time()
            wall = end_wall - start_wall
            cpu = end_cpu - start_cpu
            efficiency = cpu / wall if wall > 0 else 0
            log.log(level, f"time for {func.__name__}: CPU={cpu:.6f}s, Wall={wall:.6f}s, Efficiency={efficiency:.2%}")
            return result

        return wrapper

    return decorate


def outdated_database_warning():
    def decorate(func):
        log = logging.getLogger(func.__module__)

        @wraps(func)
        def wrapper(*args, **kwargs):
            log.warning("outdated database warning in {}".format(func.__name__))
            result = func(*args, **kwargs)
            return result

        return wrapper

    return decorate
