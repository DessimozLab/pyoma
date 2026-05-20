from importlib.metadata import version as lib_version

_author__ = "Adrian Altenhoff"
__version__ = lib_version("pyoma")


def version():
    """returns the current library version of pyoma"""
    return __version__
