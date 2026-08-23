from __future__ import annotations

import importlib
from typing import Type

from .base import BaseDeconvolver


METHODS: dict[str, Type[BaseDeconvolver]] = {}

_BUILTIN_METHODS = {
    "cell2location": "deconvatac.methods.cell2location:Cell2LocationDeconvolver",
    "destvi": "deconvatac.methods.destvi:DestVIDeconvolver",
    "rctd": "deconvatac.methods.rctd:RCTDDeconvolver",
    "spatialdwls": "deconvatac.methods.spatialdwls:SpatialDWLSDeconvolver",
    "tangram": "deconvatac.methods.tangram:TangramDeconvolver",
    "nnls": "deconvatac.methods.nnls:NNLSDeconvolver",
    "shapemix": "deconvatac.methods.shapemix:ShapeMixDeconvolver",
}


def register_method(name: str, cls: Type[BaseDeconvolver]) -> None:
    """Register a deconvolution method adapter."""
    METHODS[name] = cls


def get_method(name: str) -> Type[BaseDeconvolver]:
    """Return a method adapter class by name."""
    normalized = name.lower()
    if normalized not in METHODS:
        if normalized not in _BUILTIN_METHODS:
            available = ", ".join(list_methods())
            raise KeyError(f"Unknown method '{name}'. Available methods: {available}")
        module_name, class_name = _BUILTIN_METHODS[normalized].split(":")
        module = importlib.import_module(module_name)
        METHODS[normalized] = getattr(module, class_name)
    return METHODS[normalized]


def list_methods() -> list[str]:
    """List registered and built-in method names."""
    return sorted(set(METHODS) | set(_BUILTIN_METHODS))
