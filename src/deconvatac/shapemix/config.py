"""Strict configuration for the frozen ShapeMix model-version-1 contract."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, fields
from numbers import Integral, Real
from typing import Any, Mapping


def _positive_float(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real number.")
    normalized = float(value)
    if not math.isfinite(normalized) or normalized <= 0:
        raise ValueError(f"{name} must be finite and strictly positive.")
    return normalized


def _positive_integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer.")
    normalized = int(value)
    if normalized <= 0:
        raise ValueError(f"{name} must be strictly positive.")
    return normalized


def _nonnegative_integer(value: Any, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer.")
    normalized = int(value)
    if normalized < 0:
        raise ValueError(f"{name} must be nonnegative.")
    return normalized


@dataclass(frozen=True)
class ShapeMixConfig:
    """Validated ShapeMix configuration.

    Statistical choices are fixed by model version 1.  Positive optimizer-size
    overrides are accepted for toy tests and explicitly labeled sensitivity
    runs; :meth:`validate_protocol_v1` distinguishes those configurations from
    the frozen primary settings.
    """

    use_shape: bool = True
    total_likelihood: str = "negative_binomial"
    conditional_shape_likelihood: str = "multinomial"
    signature_rate_pseudocount: float = 0.5
    signature_shape_concentration: float = 1.0
    exposure_mode: str = "absorbed_in_abundance"
    dispersion_mode: str = "reference_crossfit_global_scaled_by_abundance"
    dispersion_crossfit_folds: int = 2
    dispersion_alpha_floor: float = 1.0e-8
    background_mode: str = "none"
    abundance_prior: str = "gamma"
    abundance_prior_shape: float = 2.0
    abundance_prior_rate: float = 1.0
    optimizer: str = "adam"
    learning_rate: float = 0.03
    max_steps: int = 2_000
    patience: int = 100
    tolerance: float = 1.0e-5
    restarts: int = 3
    spot_batch_size: int = 64
    peak_chunk_size: int = 512
    seed: int = 0
    device: str = "cpu"
    dtype: str = "float32"

    def __post_init__(self) -> None:
        if not isinstance(self.use_shape, bool):
            raise TypeError("use_shape must be a boolean.")

        fixed_strings = {
            "total_likelihood": "negative_binomial",
            "conditional_shape_likelihood": "multinomial",
            "exposure_mode": "absorbed_in_abundance",
            "dispersion_mode": "reference_crossfit_global_scaled_by_abundance",
            "background_mode": "none",
            "abundance_prior": "gamma",
            "optimizer": "adam",
            "device": "cpu",
            "dtype": "float32",
        }
        for name, expected in fixed_strings.items():
            observed = getattr(self, name)
            if not isinstance(observed, str):
                raise TypeError(f"{name} must be a string.")
            if observed != expected:
                raise ValueError(
                    f"{name}={observed!r} is not supported by ShapeMix model version 1; "
                    f"expected {expected!r}."
                )

        fixed_floats = {
            "signature_rate_pseudocount": 0.5,
            "signature_shape_concentration": 1.0,
            "dispersion_alpha_floor": 1.0e-8,
            "abundance_prior_shape": 2.0,
            "abundance_prior_rate": 1.0,
        }
        for name, expected in fixed_floats.items():
            normalized = _positive_float(getattr(self, name), name)
            object.__setattr__(self, name, normalized)
            if normalized != expected:
                raise ValueError(
                    f"{name}={normalized!r} changes the frozen model-version-1 contract; "
                    f"expected {expected!r}."
                )

        folds = _positive_integer(
            self.dispersion_crossfit_folds, "dispersion_crossfit_folds"
        )
        object.__setattr__(self, "dispersion_crossfit_folds", folds)
        if folds != 2:
            raise ValueError("ShapeMix model version 1 requires exactly two dispersion folds.")

        for name in (
            "max_steps",
            "patience",
            "restarts",
            "spot_batch_size",
            "peak_chunk_size",
        ):
            object.__setattr__(self, name, _positive_integer(getattr(self, name), name))
        for name in ("learning_rate", "tolerance"):
            object.__setattr__(self, name, _positive_float(getattr(self, name), name))

        seed = _nonnegative_integer(self.seed, "seed")
        object.__setattr__(self, "seed", seed)
        if seed != 0:
            raise ValueError("ShapeMix protocol version 1 fixes the method seed at 0.")

    @classmethod
    def from_mapping(cls, value: Mapping[str, Any]) -> "ShapeMixConfig":
        """Parse direct parameters or a strict ``method/params`` configuration."""
        if not isinstance(value, Mapping):
            raise TypeError("ShapeMix configuration must be a mapping.")
        if any(not isinstance(key, str) for key in value):
            raise TypeError("ShapeMix configuration keys must be strings.")

        if "method" in value or "params" in value:
            unknown_outer = sorted(set(value).difference({"method", "params"}))
            if unknown_outer:
                raise ValueError(
                    f"Unknown ShapeMix configuration keys: {', '.join(unknown_outer)}."
                )
            if value.get("method") != "shapemix":
                raise ValueError("ShapeMix method configuration must declare method: shapemix.")
            params = value.get("params")
            if not isinstance(params, Mapping):
                raise TypeError("ShapeMix method configuration requires a params mapping.")
        else:
            params = value

        if any(not isinstance(key, str) for key in params):
            raise TypeError("ShapeMix parameter keys must be strings.")

        allowed = {field.name for field in fields(cls)}
        unknown = sorted(set(params).difference(allowed))
        if unknown:
            raise ValueError(f"Unknown ShapeMix parameter keys: {', '.join(unknown)}.")
        return cls(**dict(params))

    def to_dict(self) -> dict[str, Any]:
        """Return parameters in stable dataclass field order."""
        return {field.name: getattr(self, field.name) for field in fields(self)}

    @property
    def content_sha256(self) -> str:
        """Hash the complete normalized configuration using canonical JSON."""
        encoded = json.dumps(
            self.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    @property
    def is_protocol_v1(self) -> bool:
        """Return whether all settings equal the frozen primary configuration."""
        defaults = type(self)()
        return all(
            field.name == "use_shape"
            or getattr(self, field.name) == getattr(defaults, field.name)
            for field in fields(self)
        )

    def validate_protocol_v1(self) -> None:
        """Fail if an accepted optimizer override makes this a sensitivity run."""
        defaults = type(self)()
        changed = {
            field.name: (getattr(defaults, field.name), getattr(self, field.name))
            for field in fields(self)
            if field.name != "use_shape"
            and getattr(self, field.name) != getattr(defaults, field.name)
        }
        if changed:
            details = ", ".join(
                f"{name}: expected {expected!r}, observed {observed!r}"
                for name, (expected, observed) in changed.items()
            )
            raise ValueError(f"Configuration is not protocol-version-1 conforming: {details}.")


def validate_nested_ablation_configs(
    shape_aware: ShapeMixConfig,
    count_only: ShapeMixConfig,
) -> None:
    """Require two configs to differ only in the frozen ``use_shape`` switch."""
    if not isinstance(shape_aware, ShapeMixConfig) or not isinstance(
        count_only, ShapeMixConfig
    ):
        raise TypeError("Both ablation configurations must be ShapeMixConfig instances.")
    if shape_aware.use_shape is not True or count_only.use_shape is not False:
        raise ValueError("Ablation configs must be ordered as shape-aware then count-only.")
    mismatched = [
        field.name
        for field in fields(shape_aware)
        if field.name != "use_shape"
        and getattr(shape_aware, field.name) != getattr(count_only, field.name)
    ]
    if mismatched:
        raise ValueError(
            "Shape-aware and count-only configs differ outside use_shape: "
            f"{', '.join(mismatched)}."
        )
