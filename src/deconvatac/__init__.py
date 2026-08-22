from importlib.metadata import PackageNotFoundError, version


from . import data, methods, pp, tl

__all__ = ["data", "methods", "pp", "tl"]

try:
    __version__ = version("deconvATAC")
except PackageNotFoundError:
    __version__ = "0+unknown"


#from . import pl, pp, tl

#__all__ = ["pl", "pp", "tl"]

#__version__ = version("deconvATAC")
