"""Differentiable AutoDock-Vina affinity-map interpolation.

The upstream ``grid`` object stores ``n_voxels + 1`` samples.  Its physical
origin is ``center - 0.5 * n_voxels * spacing`` and the flat map-file order is
x, then y, then z.  This module keeps that convention while exposing a small,
pure-Python map container and ChainRules-compatible interpolation/scoring
rules.  Map generation remains an imperative operation in the upstream Vina
binding; the differentiable boundary starts with the values supplied here.
"""

from __future__ import annotations

import glob
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .core import _rows
from .protocol import NonDifferentiablePoint, UnsupportedWrt, ZERO, rules


class GridBoundaryError(ValueError):
    """A query is outside the recorded affinity-map box."""


_XS_TYPE_NAMES = (
    "C_H", "C_P", "N_P", "N_D", "N_A", "N_DA", "O_P", "O_D", "O_A",
    "O_DA", "S_P", "P_P", "F_H", "Cl_H", "Br_H", "I_H", "Si", "At",
    "Met_D", "C_H_CG0", "C_P_CG0", "G0", "C_H_CG1", "C_P_CG1", "G1",
    "C_H_CG2", "C_P_CG2", "G2", "C_H_CG3", "C_P_CG3", "G3", "W",
)
_XS_NAME_TO_TYPE = {name.lower(): i for i, name in enumerate(_XS_TYPE_NAMES)}
# Vina's cache uses the ordinary C_H/C_P map for all macrocycle closure
# variants.  Keep those aliases available when loading or querying maps.
_XS_BASE_TYPE = {
    19: 0, 20: 1, 22: 0, 23: 1, 25: 0, 26: 1, 28: 0, 29: 1,
    21: 21, 24: 24, 27: 27, 30: 30,
}


def _finite_vector(value: Any, name: str) -> tuple[float, float, float]:
    if isinstance(value, (str, bytes)):
        raise TypeError(f"{name} must be a length-3 real sequence")
    try:
        values = list(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be a length-3 real sequence") from exc
    if len(values) != 3:
        raise ValueError(f"{name} must have length 3")
    out = []
    for item in values:
        if isinstance(item, bool):
            raise TypeError(f"{name} must contain real numbers")
        try:
            number = float(item)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{name} must contain real numbers") from exc
        if not math.isfinite(number):
            raise ValueError(f"{name} must be finite")
        out.append(number)
    return tuple(out)  # type: ignore[return-value]


def _spacing(value: Any) -> tuple[float, float, float]:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        result = (float(value),) * 3
    else:
        result = _finite_vector(value, "spacing")
    if any(not math.isfinite(item) or item <= 0 for item in result):
        raise ValueError("spacing must be positive and finite")
    return result


def _same_vector(left: Sequence[float], right: Sequence[float], *, tolerance: float = 1e-10) -> bool:
    return all(abs(float(left[i]) - float(right[i])) <= tolerance * max(1.0, abs(float(left[i])), abs(float(right[i]))) for i in range(3))


def _shape3(values: Any) -> tuple[int, int, int]:
    shape = getattr(values, "shape", None)
    if shape is not None:
        try:
            result = tuple(int(item) for item in shape)
        except (TypeError, ValueError) as exc:
            raise TypeError("grid values must be a three-dimensional array") from exc
        if len(result) != 3:
            raise ValueError("grid values must have shape (nx, ny, nz)")
    else:
        if isinstance(values, (str, bytes)):
            raise TypeError("grid values must be a three-dimensional array")
        try:
            x = list(values)
        except TypeError as exc:
            raise TypeError("grid values must be a three-dimensional array") from exc
        if not x:
            raise ValueError("grid values cannot be empty")
        try:
            y = [list(row) for row in x]
        except TypeError as exc:
            raise TypeError("grid values must be a three-dimensional array") from exc
        if not y or not y[0]:
            raise ValueError("grid values cannot have an empty dimension")
        try:
            z = [list(column) for column in y[0]]
        except TypeError as exc:
            raise TypeError("grid values must have shape (nx, ny, nz)") from exc
        # ``z`` contains one entry for each Y column in the first X row;
        # each entry is itself the Z column.  The third dimension therefore
        # comes from the length of one column, not from the number of columns.
        result = (len(y), len(y[0]), len(z[0]))
        if any(len(row) != result[1] for row in y) or any(
            len(column) != result[2] for row in y for column in row
        ):
            raise ValueError("grid values must be rectangular")
    if any(item < 2 for item in result):
        raise ValueError("each grid dimension must contain at least two samples")
    return result  # type: ignore[return-value]


def _as_type(value: Any) -> int:
    if isinstance(value, bool):
        raise TypeError("atom_types must contain integer XS_TYPE values or names")
    if isinstance(value, int):
        if value < 0 or value >= len(_XS_TYPE_NAMES):
            raise ValueError(f"atom type must be in [0, {len(_XS_TYPE_NAMES)})")
        return value
    if isinstance(value, str):
        key = value.strip().lower()
        if key in _XS_NAME_TO_TYPE:
            return _XS_NAME_TO_TYPE[key]
    raise TypeError("atom_types must contain integer XS_TYPE values or names")


def _type_key(value: Any) -> int:
    return _as_type(value)


def _grid_value(values: Any, i: int, j: int, k: int) -> float:
    try:
        result = values[i][j][k]
    except (IndexError, KeyError, TypeError):
        result = values[i, j, k]
    try:
        number = float(result)
    except (TypeError, ValueError) as exc:
        raise TypeError("grid values must contain real numbers") from exc
    if not math.isfinite(number):
        raise ValueError("grid values must be finite")
    return number


def _nested_zeros(shape: tuple[int, int, int]) -> list[list[list[float]]]:
    return [
        [[0.0 for _ in range(shape[2])] for _ in range(shape[1])]
        for _ in range(shape[0])
    ]


def _zeros_like(values: Any, shape: tuple[int, int, int]) -> Any:
    try:
        import numpy as np  # optional dependency

        if hasattr(values, "shape"):
            return np.zeros_like(values, dtype=float)
    except ImportError:
        pass
    return _nested_zeros(shape)


def _add_grid_value(values: Any, i: int, j: int, k: int, amount: float) -> None:
    try:
        values[i, j, k] += amount
        return
    except (TypeError, IndexError):
        pass
    values[i][j][k] += amount


def _copy_grid_tangent(values: Any, shape: tuple[int, int, int]) -> Any:
    if hasattr(values, "copy"):
        try:
            copied = values.copy()
            # Force validation and a floating interpretation at query time.
            _shape3(copied)
            return copied
        except (AttributeError, TypeError, ValueError):
            pass
    return values


@dataclass(frozen=True)
class AffinityGrid:
    """One atom-type map and its physical grid provenance.

    ``center`` follows Vina's map header.  ``box_size`` is optional metadata;
    when omitted it is derived as ``(shape - 1) * spacing``.
    """

    values: Any
    center: Sequence[float]
    spacing: float | Sequence[float]
    box_size: Sequence[float] | None = None
    source: str | None = None

    def __post_init__(self) -> None:
        shape = _shape3(self.values)
        center = _finite_vector(self.center, "center")
        spacing = _spacing(self.spacing)
        derived = tuple((shape[i] - 1) * spacing[i] for i in range(3))
        if self.box_size is None:
            box = derived
        else:
            box = _finite_vector(self.box_size, "box_size")
            if any(item <= 0 for item in box):
                raise ValueError("box_size must be positive")
            if any(abs(box[i] - derived[i]) > 1e-8 * max(1.0, abs(derived[i])) for i in range(3)):
                raise ValueError("box_size must equal (grid_shape - 1) * spacing")
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "spacing", spacing)
        object.__setattr__(self, "box_size", box)

    @property
    def shape(self) -> tuple[int, int, int]:
        return _shape3(self.values)

    @property
    def origin(self) -> tuple[float, float, float]:
        return tuple(self.center[i] - 0.5 * self.box_size[i] for i in range(3))  # type: ignore[index]

    @property
    def provenance(self) -> dict[str, Any]:
        return {
            "center": self.center,
            "box_size": self.box_size,
            "spacing": self.spacing,
            "shape": self.shape,
            "source": self.source,
        }


@dataclass(frozen=True)
class AffinityMaps:
    """A family of X-Score atom-type maps sharing one grid geometry."""

    values: Mapping[Any, Any]
    center: Sequence[float]
    spacing: float | Sequence[float]
    box_size: Sequence[float] | None = None
    source: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.values, Mapping) or not self.values:
            raise ValueError("values must be a non-empty mapping of atom types to grids")
        center = _finite_vector(self.center, "center")
        spacing = _spacing(self.spacing)
        normalized: dict[int, Any] = {}
        shape: tuple[int, int, int] | None = None
        for key, values in self.values.items():
            atom_type = _type_key(key)
            if isinstance(values, AffinityGrid):
                if not _same_vector(values.center, center) or not _same_vector(values.spacing, spacing):
                    raise ValueError("AffinityGrid geometry does not match the map family")
                if self.box_size is not None and not _same_vector(values.box_size, _finite_vector(self.box_size, "box_size")):
                    raise ValueError("AffinityGrid box_size does not match the map family")
                values = values.values
            current_shape = _shape3(values)
            if shape is None:
                shape = current_shape
            elif current_shape != shape:
                raise ValueError("all affinity maps must have identical shapes")
            if atom_type in normalized:
                raise ValueError(f"duplicate affinity map atom type {atom_type}")
            normalized[atom_type] = values
        assert shape is not None
        derived = tuple((shape[i] - 1) * spacing[i] for i in range(3))
        if self.box_size is None:
            box = derived
        else:
            box = _finite_vector(self.box_size, "box_size")
            if any(item <= 0 for item in box):
                raise ValueError("box_size must be positive")
            if any(abs(box[i] - derived[i]) > 1e-8 * max(1.0, abs(derived[i])) for i in range(3)):
                raise ValueError("box_size must equal (grid_shape - 1) * spacing")
        object.__setattr__(self, "values", normalized)
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "spacing", spacing)
        object.__setattr__(self, "box_size", box)

    @property
    def shape(self) -> tuple[int, int, int]:
        return _shape3(next(iter(self.values.values())))

    @property
    def origin(self) -> tuple[float, float, float]:
        return tuple(self.center[i] - 0.5 * self.box_size[i] for i in range(3))  # type: ignore[index]

    @property
    def provenance(self) -> dict[str, Any]:
        return {
            "center": self.center,
            "box_size": self.box_size,
            "spacing": self.spacing,
            "shape": self.shape,
            "atom_types": tuple(sorted(self.values)),
            "source": self.source,
        }

    def grid(self, atom_type: Any) -> Any:
        key = _as_type(atom_type)
        base = _XS_BASE_TYPE.get(key, key)
        if base in self.values:
            return self.values[base]
        if key in self.values:
            return self.values[key]
        name = _XS_TYPE_NAMES[key]
        raise KeyError(f"no affinity map is present for atom type {name}")

    def with_values(self, values: Mapping[Any, Any]) -> "AffinityMaps":
        return AffinityMaps(values, self.center, self.spacing, self.box_size, self.source)


def _map_family(maps: AffinityMaps | Mapping[Any, Any], *, center: Any = None, spacing: Any = None) -> AffinityMaps:
    if isinstance(maps, AffinityMaps):
        if center is not None or spacing is not None:
            raise TypeError("center/spacing must not be supplied with AffinityMaps")
        return maps
    if center is None or spacing is None:
        raise TypeError("a raw map mapping requires center and spacing")
    return AffinityMaps(maps, center=center, spacing=spacing)


def _query_point(point: Sequence[float], grid: AffinityGrid | AffinityMaps) -> tuple[float, float, float]:
    return _finite_vector(point, "coordinates")


def _trilinear_one(
    values: Any,
    point: Sequence[float],
    geometry: AffinityGrid | AffinityMaps,
    *,
    differentiate: bool,
    return_weights: bool = False,
) -> tuple[float, tuple[float, float, float], tuple[tuple[int, int, int, float], ...]]:
    shape = _shape3(values)
    center = geometry.center
    spacing = geometry.spacing
    box_size = geometry.box_size
    assert box_size is not None
    origin = tuple(center[i] - 0.5 * box_size[i] for i in range(3))
    coordinates = _query_point(point, geometry)
    logical = tuple((coordinates[i] - origin[i]) / spacing[i] for i in range(3))
    eps = 1e-12
    if any(value < -eps or value > shape[i] - 1 + eps for i, value in enumerate(logical)):
        raise GridBoundaryError(
            f"coordinate {coordinates!r} is outside affinity-map box "
            f"[{origin!r}, {tuple(origin[i] + box_size[i] for i in range(3))!r}]"
        )
    logical = tuple(min(max(logical[i], 0.0), shape[i] - 1.0) for i in range(3))
    if differentiate:
        for value in logical:
            if abs(value - round(value)) <= eps:
                raise NonDifferentiablePoint("grid-cell boundary has no unique coordinate derivative")
    base = tuple(min(int(math.floor(logical[i])), shape[i] - 2) for i in range(3))
    frac = tuple(logical[i] - base[i] for i in range(3))
    x, y, z = frac
    mx, my, mz = 1.0 - x, 1.0 - y, 1.0 - z
    weights = (
        (base[0], base[1], base[2], mx * my * mz),
        (base[0] + 1, base[1], base[2], x * my * mz),
        (base[0], base[1] + 1, base[2], mx * y * mz),
        (base[0] + 1, base[1] + 1, base[2], x * y * mz),
        (base[0], base[1], base[2] + 1, mx * my * z),
        (base[0] + 1, base[1], base[2] + 1, x * my * z),
        (base[0], base[1] + 1, base[2] + 1, mx * y * z),
        (base[0] + 1, base[1] + 1, base[2] + 1, x * y * z),
    )
    value = sum(_grid_value(values, i, j, k) * weight for i, j, k, weight in weights)
    f000, f100, f010, f110, f001, f101, f011, f111 = (
        _grid_value(values, *indices[:3]) for indices in weights
    )
    dlogical = (
        (f100 - f000) * my * mz + (f110 - f010) * y * mz + (f101 - f001) * my * z + (f111 - f011) * y * z,
        (f010 - f000) * mx * mz + (f110 - f100) * x * mz + (f011 - f001) * mx * z + (f111 - f101) * x * z,
        (f001 - f000) * mx * my + (f101 - f100) * x * my + (f011 - f010) * mx * y + (f111 - f110) * x * y,
    )
    gradient = tuple(dlogical[i] / spacing[i] for i in range(3))
    return value, gradient, weights if return_weights else ()


def _curl_value(value: float, clip_value: float | None) -> tuple[float, float]:
    """Return Vina's soft upper-bound curl and its derivative.

    ``cache::eval`` calls ``curl(e, 1000)`` for map scoring.  Negative map
    values are left untouched; positive values are mapped to
    ``1000*e/(1000+e)``.  ``clip_value=None`` exposes the raw interpolation.
    """
    if clip_value is None or value <= 0.0:
        return value, 1.0
    if not math.isfinite(clip_value) or clip_value <= 0.0:
        raise ValueError("clip_value must be positive and finite, or None")
    factor = (clip_value / (clip_value + value)) ** 2
    return clip_value * value / (clip_value + value), factor


def _coordinate_rows(value: Any) -> list[list[float]]:
    rows, _ = _rows(value, tangent=True)
    return rows


def _restore_rows(original: Any, values: list[list[float]]) -> Any:
    if isinstance(original, tuple):
        return tuple(tuple(row) for row in values)
    try:
        import numpy as np

        if hasattr(original, "shape"):
            return np.asarray(values, dtype=float).reshape(original.shape)
    except ImportError:
        pass
    return values


def _restore_vector(original: Any, values: tuple[float, float, float]) -> Any:
    if isinstance(original, tuple):
        return values
    try:
        import numpy as np

        if hasattr(original, "shape"):
            return np.asarray(values, dtype=float).reshape(original.shape)
    except ImportError:
        pass
    return list(values)


def _restore_map_mapping(original: Any, family: AffinityMaps, values: Mapping[int, Any]) -> Any:
    """Restore a raw mapping's original (possibly named) keys in a pullback."""
    if isinstance(original, AffinityMaps):
        return family.with_values(values)
    restored: dict[Any, Any] = {}
    for original_key in original:
        key = _as_type(original_key)
        canonical = _XS_BASE_TYPE.get(key, key)
        if canonical in values:
            restored[original_key] = values[canonical]
        elif key in values:
            restored[original_key] = values[key]
    return restored


def interpolate_grid(
    grid_values: Any,
    coordinates: Any,
    center: Sequence[float] | None = None,
    spacing: float | Sequence[float] | None = None,
    box_size: Sequence[float] | None = None,
) -> tuple[float, ...]:
    """Trilinearly interpolate one grid at every ``(N, 3)`` coordinate."""
    if isinstance(grid_values, AffinityGrid):
        if center is not None or spacing is not None or box_size is not None:
            raise TypeError("grid provenance is already recorded in AffinityGrid")
        geometry = grid_values
        grid_values = geometry.values
    else:
        if center is None or spacing is None:
            raise TypeError("center and spacing are required for raw grid values")
        geometry = AffinityGrid(grid_values, center=center, spacing=spacing, box_size=box_size)
    rows = _coordinate_rows(coordinates)
    return tuple(_trilinear_one(grid_values, row, geometry, differentiate=False)[0] for row in rows)


trilinear_interpolate = interpolate_grid
interpolate = interpolate_grid
interpolate_map = interpolate_grid


def _map_values_for_tangent(tangent: Any, primal: AffinityMaps) -> AffinityMaps:
    if isinstance(tangent, AffinityMaps):
        if not (_same_vector(tangent.center, primal.center) and _same_vector(tangent.spacing, primal.spacing) and tangent.shape == primal.shape):
            raise ValueError("map tangent geometry must match primal AffinityMaps")
        raw = tangent.values
    elif isinstance(tangent, Mapping):
        raw = tangent
    else:
        raise TypeError("maps tangent must be an AffinityMaps or atom-type mapping")
    normalized = {_as_type(key): value for key, value in raw.items()}
    filled: dict[int, Any] = {}
    for key, values in primal.values.items():
        tangent_values = normalized.get(key)
        if tangent_values is None:
            aliases = [candidate for candidate in normalized if _XS_BASE_TYPE.get(candidate, candidate) == key]
            tangent_values = normalized[aliases[0]] if aliases else _zeros_like(values, _shape3(values))
        if _shape3(tangent_values) != _shape3(values):
            raise ValueError("map tangent geometry must match primal AffinityMaps")
        filled[key] = tangent_values
    return AffinityMaps(filled, primal.center, primal.spacing, primal.box_size)


def interpolate_maps(maps: AffinityMaps | Mapping[Any, Any], coordinates: Any, atom_types: Any, center: Any = None, spacing: Any = None) -> tuple[float, ...]:
    """Interpolate the atom-type map selected by each atom type."""
    family = _map_family(maps, center=center, spacing=spacing)
    rows = _coordinate_rows(coordinates)
    types = tuple(_as_type(value) for value in atom_types)
    if len(types) != len(rows):
        raise ValueError("atom_types must have one type per coordinate")
    return tuple(_trilinear_one(family.grid(atom_type), row, family, differentiate=False)[0] for row, atom_type in zip(rows, types))


def _interpolate_maps_linearisation(
    family: AffinityMaps,
    rows: list[list[float]],
    types: tuple[int, ...],
    *,
    differentiate: bool,
    with_map_gradient: bool,
) -> tuple[tuple[float, ...], list[list[float]], AffinityMaps | None]:
    outputs: list[float] = []
    coordinate_gradient = [[0.0, 0.0, 0.0] for _ in rows]
    map_gradient: dict[int, Any] | None = None
    if with_map_gradient:
        map_gradient = {key: _zeros_like(values, _shape3(values)) for key, values in family.values.items()}
    for row_index, (row, atom_type) in enumerate(zip(rows, types)):
        grid = family.grid(atom_type)
        sample, gradient, weights = _trilinear_one(grid, row, family, differentiate=differentiate, return_weights=True)
        outputs.append(sample)
        if differentiate:
            coordinate_gradient[row_index] = list(gradient)
        if map_gradient is not None:
            key = _XS_BASE_TYPE.get(atom_type, atom_type)
            for i, j, k, weight in weights:
                _add_grid_value(map_gradient[key], i, j, k, weight)
    return tuple(outputs), coordinate_gradient, family.with_values(map_gradient) if map_gradient is not None else None


@rules.jvp_for(interpolate_maps)
def _interpolate_maps_jvp(tangents: dict[str, Any], maps: AffinityMaps | Mapping[Any, Any], coordinates: Any, atom_types: Any, center: Any = None, spacing: Any = None) -> tuple[tuple[float, ...], Any]:
    supported = {"maps", "coordinates"}
    unsupported = set(tangents) - supported
    if unsupported:
        raise UnsupportedWrt(interpolate_maps, unsupported, supported=supported)
    family = _map_family(maps, center=center, spacing=spacing)
    rows = _coordinate_rows(coordinates)
    types = tuple(_as_type(value) for value in atom_types)
    if len(types) != len(rows):
        raise ValueError("atom_types must have one type per coordinate")
    coordinate_tangent = tangents.get("coordinates", ZERO)
    map_tangent = tangents.get("maps", ZERO)
    values, coordinate_gradient, _ = _interpolate_maps_linearisation(family, rows, types, differentiate=coordinate_tangent is not ZERO, with_map_gradient=False)
    directional = 0.0
    if coordinate_tangent is not ZERO:
        tangent_rows = _coordinate_rows(coordinate_tangent)
        if len(tangent_rows) != len(rows):
            raise ValueError("coordinates tangent must have shape (N, 3)")
        directional += sum(coordinate_gradient[i][k] * tangent_rows[i][k] for i in range(len(rows)) for k in range(3))
    if map_tangent is not ZERO:
        tangent_family = _map_values_for_tangent(map_tangent, family)
        directional += sum(_trilinear_one(tangent_family.grid(atom_type), row, tangent_family, differentiate=False)[0] for row, atom_type in zip(rows, types))
    return values, directional


@rules.vjp_for(interpolate_maps)
def _interpolate_maps_vjp(wrt: tuple[str, ...], maps: AffinityMaps | Mapping[Any, Any], coordinates: Any, atom_types: Any, center: Any = None, spacing: Any = None) -> tuple[tuple[float, ...], Any]:
    supported = {"maps", "coordinates"}
    unsupported = set(wrt) - supported
    if unsupported:
        raise UnsupportedWrt(interpolate_maps, unsupported, supported=supported)
    family = _map_family(maps, center=center, spacing=spacing)
    original_coordinates = coordinates
    rows = _coordinate_rows(coordinates)
    types = tuple(_as_type(value) for value in atom_types)
    if len(types) != len(rows):
        raise ValueError("atom_types must have one type per coordinate")
    values, coordinate_gradient, map_gradient = _interpolate_maps_linearisation(family, rows, types, differentiate="coordinates" in wrt, with_map_gradient="maps" in wrt)

    def pullback(cotangent: Any) -> dict[str, Any]:
        if cotangent is ZERO:
            return {name: ZERO for name in wrt}
        if len(rows) == 1 and isinstance(cotangent, (int, float)) and not isinstance(cotangent, bool):
            cotangents = [cotangent]
        else:
            try:
                cotangents = list(cotangent)
            except TypeError as exc:
                raise TypeError("interpolate_maps cotangent must have one value per coordinate") from exc
        if len(cotangents) != len(rows):
            raise ValueError("interpolate_maps cotangent must have one value per coordinate")
        result: dict[str, Any] = {}
        if map_gradient is not None:
            scaled = {key: _zeros_like(grid, _shape3(grid)) for key, grid in map_gradient.values.items()}
            for row, (point, atom_type) in enumerate(zip(rows, types)):
                _, _, weights = _trilinear_one(family.grid(atom_type), point, family, differentiate=False, return_weights=True)
                key = _XS_BASE_TYPE.get(atom_type, atom_type)
                for i, j, k, weight in weights:
                    _add_grid_value(scaled[key], i, j, k, float(cotangents[row]) * weight)
            result["maps"] = _restore_map_mapping(maps, family, scaled)
        if "coordinates" in wrt:
            result["coordinates"] = _restore_rows(original_coordinates, [[float(cotangents[i]) * component for component in coordinate_gradient[i]] for i in range(len(rows))])
        return result

    return values, pullback


def _rotation_matrix(rotation: Sequence[float]) -> tuple[tuple[float, ...], ...]:
    rx, ry, rz = rotation
    theta = math.sqrt(rx * rx + ry * ry + rz * rz)
    a = ((0.0, -rz, ry), (rz, 0.0, -rx), (-ry, rx, 0.0))
    a2 = tuple(tuple(sum(a[i][k] * a[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    if theta < 1e-8:
        first = 1.0 - theta * theta / 6.0 + theta**4 / 120.0
        second = 0.5 - theta * theta / 24.0 + theta**4 / 720.0
    else:
        first = math.sin(theta) / theta
        second = (1.0 - math.cos(theta)) / (theta * theta)
    return tuple(tuple((1.0 if i == j else 0.0) + first * a[i][j] + second * a2[i][j] for j in range(3)) for i in range(3))


def _right_jacobian(rotation: Sequence[float]) -> tuple[tuple[float, ...], ...]:
    rx, ry, rz = rotation
    theta = math.sqrt(rx * rx + ry * ry + rz * rz)
    a = ((0.0, -rz, ry), (rz, 0.0, -rx), (-ry, rx, 0.0))
    a2 = tuple(tuple(sum(a[i][k] * a[k][j] for k in range(3)) for j in range(3)) for i in range(3))
    if theta < 1e-6:
        c1 = 0.5 - theta * theta / 24.0 + theta**4 / 720.0
        c2 = 1.0 / 6.0 - theta * theta / 120.0 + theta**4 / 5040.0
    else:
        c1 = (1.0 - math.cos(theta)) / (theta * theta)
        c2 = (theta - math.sin(theta)) / (theta**3)
    return tuple(tuple((1.0 if i == j else 0.0) - c1 * a[i][j] + c2 * a2[i][j] for j in range(3)) for i in range(3))


def _transform_point(point: Sequence[float], translation: Sequence[float], rotation: Sequence[float]) -> tuple[float, float, float]:
    matrix = _rotation_matrix(rotation)
    return tuple(sum(matrix[i][j] * point[j] for j in range(3)) + translation[i] for i in range(3))  # type: ignore[return-value]


def _transform_jacobians(point: Sequence[float], rotation: Sequence[float]) -> tuple[tuple[float, float, float], ...]:
    matrix = _rotation_matrix(rotation)
    jac = _right_jacobian(rotation)
    result = []
    for parameter in range(3):
        omega = tuple(jac[i][parameter] for i in range(3))
        cross = (omega[1] * point[2] - omega[2] * point[1], omega[2] * point[0] - omega[0] * point[2], omega[0] * point[1] - omega[1] * point[0])
        result.append(tuple(sum(matrix[i][j] * cross[j] for j in range(3)) for i in range(3)))
    return tuple(result)


def transform_pose(coordinates: Any, translation: Any = (0.0, 0.0, 0.0), rotation: Any = (0.0, 0.0, 0.0)) -> tuple[tuple[float, float, float], ...]:
    """Apply a Rodrigues rotation vector followed by a Cartesian translation."""
    rows = _coordinate_rows(coordinates)
    shift = _finite_vector(translation, "translation")
    angle = _finite_vector(rotation, "rotation")
    return tuple(_transform_point(row, shift, angle) for row in rows)


@rules.jvp_for(transform_pose)
def _transform_pose_jvp(tangents: dict[str, Any], coordinates: Any, translation: Any = (0.0, 0.0, 0.0), rotation: Any = (0.0, 0.0, 0.0)) -> tuple[tuple[float, float, float], ...]:
    supported = {"coordinates", "translation", "rotation"}
    unsupported = set(tangents) - supported
    if unsupported:
        raise UnsupportedWrt(transform_pose, unsupported, supported=supported)
    rows = _coordinate_rows(coordinates)
    shift = _finite_vector(translation, "translation")
    angle = _finite_vector(rotation, "rotation")
    matrix = _rotation_matrix(angle)
    coordinate_tangent = _coordinate_rows(tangents["coordinates"]) if tangents.get("coordinates", ZERO) is not ZERO else None
    if coordinate_tangent is not None and len(coordinate_tangent) != len(rows):
        raise ValueError("coordinates tangent must have shape (N, 3)")
    translation_tangent = _finite_vector(tangents["translation"], "translation tangent") if tangents.get("translation", ZERO) is not ZERO else (0.0, 0.0, 0.0)
    rotation_tangent = _finite_vector(tangents["rotation"], "rotation tangent") if tangents.get("rotation", ZERO) is not ZERO else (0.0, 0.0, 0.0)
    outputs = []
    for index, row in enumerate(rows):
        tangent = [sum(matrix[i][j] * coordinate_tangent[index][j] for j in range(3)) if coordinate_tangent is not None else 0.0 for i in range(3)]
        tangent = [tangent[i] + translation_tangent[i] for i in range(3)]
        jacobians = _transform_jacobians(row, angle)
        for parameter in range(3):
            for component in range(3):
                tangent[component] += jacobians[parameter][component] * rotation_tangent[parameter]
        outputs.append(tuple(tangent))
    return transform_pose(coordinates, translation, rotation), tuple(outputs)


@rules.vjp_for(transform_pose)
def _transform_pose_vjp(wrt: tuple[str, ...], coordinates: Any, translation: Any = (0.0, 0.0, 0.0), rotation: Any = (0.0, 0.0, 0.0)) -> tuple[tuple[tuple[float, float, float], ...], Any]:
    supported = {"coordinates", "translation", "rotation"}
    unsupported = set(wrt) - supported
    if unsupported:
        raise UnsupportedWrt(transform_pose, unsupported, supported=supported)
    rows = _coordinate_rows(coordinates)
    shift = _finite_vector(translation, "translation")
    angle = _finite_vector(rotation, "rotation")
    matrix = _rotation_matrix(angle)
    outputs = transform_pose(coordinates, translation, rotation)

    def pullback(cotangent: Any) -> dict[str, Any]:
        if cotangent is ZERO:
            return {name: ZERO for name in wrt}
        try:
            cotangent_rows = [list(row) for row in cotangent]
        except TypeError as exc:
            raise TypeError("transform_pose cotangent must have shape (N, 3)") from exc
        if len(cotangent_rows) != len(rows) or any(len(row) != 3 for row in cotangent_rows):
            raise ValueError("transform_pose cotangent must have shape (N, 3)")
        result: dict[str, Any] = {}
        if "coordinates" in wrt:
            result["coordinates"] = _restore_rows(coordinates, [[sum(matrix[j][i] * cotangent_rows[row][j] for j in range(3)) for i in range(3)] for row in range(len(rows))])
        if "translation" in wrt:
            result["translation"] = _restore_vector(translation, tuple(sum(row[i] for row in cotangent_rows) for i in range(3)))
        if "rotation" in wrt:
            gradients = [0.0, 0.0, 0.0]
            for row, cotangent_row in zip(rows, cotangent_rows):
                for parameter, jacobian in enumerate(_transform_jacobians(row, angle)):
                    gradients[parameter] += sum(cotangent_row[i] * jacobian[i] for i in range(3))
            result["rotation"] = _restore_vector(rotation, tuple(gradients))
        return result

    return outputs, pullback


def _score_map_linearisation(
    family: AffinityMaps,
    rows: list[list[float]],
    types: tuple[int, ...],
    translation: tuple[float, float, float],
    rotation: tuple[float, float, float],
    *,
    clip_value: float | None,
    active: frozenset[str],
) -> tuple[float, list[list[float]], tuple[float, float, float], tuple[float, float, float], AffinityMaps | None]:
    if clip_value is not None and (not math.isfinite(clip_value) or clip_value <= 0.0):
        raise ValueError("clip_value must be positive and finite, or None")
    derivative_position = bool(active & {"coordinates", "translation", "rotation"})
    matrix = _rotation_matrix(rotation)
    transformed = [_transform_point(row, translation, rotation) for row in rows]
    local_gradient = [[0.0, 0.0, 0.0] for _ in rows]
    translation_gradient = [0.0, 0.0, 0.0]
    rotation_gradient = [0.0, 0.0, 0.0]
    map_gradients: dict[int, Any] | None = None
    if "maps" in active:
        map_gradients = {key: _zeros_like(values, _shape3(values)) for key, values in family.values.items()}
    value = 0.0
    for index, (row, atom_type, point) in enumerate(zip(rows, types, transformed)):
        grid = family.grid(atom_type)
        raw_sample, raw_gradient, weights = _trilinear_one(grid, point, family, differentiate=derivative_position, return_weights=True)
        sample, curl_factor = _curl_value(raw_sample, clip_value)
        gradient = tuple(component * curl_factor for component in raw_gradient)
        value += sample
        if map_gradients is not None:
            key = _XS_BASE_TYPE.get(atom_type, atom_type)
            if key not in map_gradients:
                # A closure type may resolve through its base map; this branch
                # is only reachable for an invalid/mismatched tangent family.
                key = atom_type
            for i, j, k, weight in weights:
                _add_grid_value(map_gradients[key], i, j, k, curl_factor * weight)
        if derivative_position:
            for component in range(3):
                translation_gradient[component] += gradient[component]
            jacobians = _transform_jacobians(row, rotation)
            for parameter in range(3):
                rotation_gradient[parameter] += sum(gradient[i] * jacobians[parameter][i] for i in range(3))
            for local_component in range(3):
                local_gradient[index][local_component] += sum(matrix[world][local_component] * gradient[world] for world in range(3))
    map_result = family.with_values(map_gradients) if map_gradients is not None else None
    return value, local_gradient, tuple(translation_gradient), tuple(rotation_gradient), map_result


def score_affinity_maps(
    maps: AffinityMaps | Mapping[Any, Any],
    coordinates: Any,
    atom_types: Any,
    *,
    translation: Any = (0.0, 0.0, 0.0),
    rotation: Any = (0.0, 0.0, 0.0),
    clip_value: float | None = 1000.0,
    center: Any = None,
    spacing: Any = None,
) -> float:
    """Sum interpolated atom-type affinity values after a rigid pose transform."""
    family = _map_family(maps, center=center, spacing=spacing)
    rows = _coordinate_rows(coordinates)
    types = tuple(_as_type(value) for value in atom_types)
    if len(types) != len(rows):
        raise ValueError("atom_types must have one type per coordinate")
    shift = _finite_vector(translation, "translation")
    angle = _finite_vector(rotation, "rotation")
    return _score_map_linearisation(family, rows, types, shift, angle, clip_value=clip_value, active=frozenset())[0]


score_maps = score_affinity_maps
grid_score = score_affinity_maps
score_grid = score_affinity_maps
map_score = score_affinity_maps


def _pose_vector(value: Any) -> tuple[float, float, float, float, float, float]:
    if isinstance(value, (str, bytes)):
        raise TypeError("pose must be a length-6 real sequence (translation xyz, rotation xyz)")
    try:
        values = list(value)
    except TypeError as exc:
        raise TypeError("pose must be a length-6 real sequence (translation xyz, rotation xyz)") from exc
    if len(values) != 6:
        raise ValueError("pose must have length 6 (translation xyz, rotation xyz)")
    converted = tuple(float(item) for item in values)
    if any(not math.isfinite(item) for item in converted):
        raise ValueError("pose must be finite")
    return converted  # type: ignore[return-value]


def score_pose(
    maps: AffinityMaps | Mapping[Any, Any],
    coordinates: Any,
    atom_types: Any,
    *,
    pose: Any = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    clip_value: float | None = 1000.0,
    center: Any = None,
    spacing: Any = None,
) -> float:
    """Map score using one six-vector ``(translation, rotation)`` pose."""
    pose_values = _pose_vector(pose)
    return score_affinity_maps(
        maps,
        coordinates,
        atom_types,
        translation=pose_values[:3],
        rotation=pose_values[3:],
        clip_value=clip_value,
        center=center,
        spacing=spacing,
    )


@rules.jvp_for(score_pose)
def _score_pose_jvp(tangents: dict[str, Any], maps: AffinityMaps | Mapping[Any, Any], coordinates: Any, atom_types: Any, *, pose: Any = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), clip_value: float | None = 1000.0, center: Any = None, spacing: Any = None) -> tuple[float, Any]:
    supported = {"maps", "coordinates", "pose"}
    unsupported = set(tangents) - supported
    if unsupported:
        raise UnsupportedWrt(score_pose, unsupported, supported=supported)
    pose_values = _pose_vector(pose)
    inner_tangents: dict[str, Any] = {name: tangents[name] for name in ("maps", "coordinates") if name in tangents}
    if tangents.get("pose", ZERO) is not ZERO:
        tangent_pose = _pose_vector(tangents["pose"])
        inner_tangents["translation"] = tangent_pose[:3]
        inner_tangents["rotation"] = tangent_pose[3:]
    return _score_affinity_maps_jvp(
        inner_tangents,
        maps,
        coordinates,
        atom_types,
        translation=pose_values[:3],
        rotation=pose_values[3:],
        clip_value=clip_value,
        center=center,
        spacing=spacing,
    )


@rules.vjp_for(score_pose)
def _score_pose_vjp(wrt: tuple[str, ...], maps: AffinityMaps | Mapping[Any, Any], coordinates: Any, atom_types: Any, *, pose: Any = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), clip_value: float | None = 1000.0, center: Any = None, spacing: Any = None) -> tuple[float, Any]:
    supported = {"maps", "coordinates", "pose"}
    unsupported = set(wrt) - supported
    if unsupported:
        raise UnsupportedWrt(score_pose, unsupported, supported=supported)
    pose_values = _pose_vector(pose)
    inner_wrt = tuple(name for name in ("maps", "coordinates") if name in wrt)
    if "pose" in wrt:
        inner_wrt += ("translation", "rotation")
    value, inner_pullback = _score_affinity_maps_vjp(
        inner_wrt,
        maps,
        coordinates,
        atom_types,
        translation=pose_values[:3],
        rotation=pose_values[3:],
        clip_value=clip_value,
        center=center,
        spacing=spacing,
    )

    def pullback(cotangent: Any) -> dict[str, Any]:
        if cotangent is ZERO:
            return {name: ZERO for name in wrt}
        inner = inner_pullback(cotangent)
        result: dict[str, Any] = {}
        if "maps" in wrt:
            result["maps"] = inner["maps"]
        if "coordinates" in wrt:
            result["coordinates"] = inner["coordinates"]
        if "pose" in wrt:
            result["pose"] = tuple(inner["translation"]) + tuple(inner["rotation"])
        return result

    return value, pullback


pose_score = score_pose
GridMap = AffinityGrid
MapFamily = AffinityMaps


@rules.jvp_for(interpolate_grid)
def _interpolate_grid_jvp(tangents: dict[str, Any], grid_values: Any, coordinates: Any, center: Sequence[float] | None = None, spacing: float | Sequence[float] | None = None, box_size: Sequence[float] | None = None) -> tuple[tuple[float, ...], Any]:
    supported = {"grid_values", "coordinates"}
    unsupported = set(tangents) - supported
    if unsupported:
        raise UnsupportedWrt(interpolate_grid, unsupported, supported=supported)
    if isinstance(grid_values, AffinityGrid):
        if center is not None or spacing is not None or box_size is not None:
            raise TypeError("grid provenance is already recorded in AffinityGrid")
        geometry = grid_values
        grid_values = geometry.values
    else:
        if center is None or spacing is None:
            raise TypeError("center and spacing are required for raw grid values")
        geometry = AffinityGrid(grid_values, center=center, spacing=spacing, box_size=box_size)
    rows = _coordinate_rows(coordinates)
    coordinate_tangent = tangents.get("coordinates", ZERO)
    grid_tangent = tangents.get("grid_values", ZERO)
    if coordinate_tangent is ZERO and grid_tangent is ZERO:
        return interpolate_grid(grid_values, coordinates, center=geometry.center, spacing=geometry.spacing, box_size=geometry.box_size), ZERO
    tangent_rows = _coordinate_rows(coordinate_tangent) if coordinate_tangent is not ZERO else None
    if tangent_rows is not None and len(tangent_rows) != len(rows):
        raise ValueError("coordinates tangent must have shape (N, 3)")
    if grid_tangent is not ZERO:
        if isinstance(grid_tangent, AffinityGrid):
            if not (_same_vector(grid_tangent.center, geometry.center) and _same_vector(grid_tangent.spacing, geometry.spacing) and grid_tangent.shape == geometry.shape):
                raise ValueError("grid tangent geometry must match the primal AffinityGrid")
            grid_tangent = grid_tangent.values
        _shape3(grid_tangent)
    outputs = []
    for index, row in enumerate(rows):
        sample, gradient, weights = _trilinear_one(grid_values, row, geometry, differentiate=coordinate_tangent is not ZERO, return_weights=True)
        directional = 0.0
        if coordinate_tangent is not ZERO:
            directional += sum(gradient[k] * tangent_rows[index][k] for k in range(3))  # type: ignore[index]
        if grid_tangent is not ZERO:
            directional += sum(_grid_value(grid_tangent, i, j, k) * weight for i, j, k, weight in weights)
        outputs.append(directional)
    return tuple(_trilinear_one(grid_values, row, geometry, differentiate=False)[0] for row in rows), tuple(outputs)


@rules.vjp_for(interpolate_grid)
def _interpolate_grid_vjp(wrt: tuple[str, ...], grid_values: Any, coordinates: Any, center: Sequence[float] | None = None, spacing: float | Sequence[float] | None = None, box_size: Sequence[float] | None = None) -> tuple[tuple[float, ...], Any]:
    supported = {"grid_values", "coordinates"}
    unsupported = set(wrt) - supported
    if unsupported:
        raise UnsupportedWrt(interpolate_grid, unsupported, supported=supported)
    if isinstance(grid_values, AffinityGrid):
        if center is not None or spacing is not None or box_size is not None:
            raise TypeError("grid provenance is already recorded in AffinityGrid")
        geometry = grid_values
        grid_values = geometry.values
    else:
        if center is None or spacing is None:
            raise TypeError("center and spacing are required for raw grid values")
        geometry = AffinityGrid(grid_values, center=center, spacing=spacing, box_size=box_size)
    original_coordinates = coordinates
    rows = _coordinate_rows(coordinates)
    active = frozenset(wrt)
    values = tuple(_trilinear_one(grid_values, row, geometry, differentiate="coordinates" in active)[0] for row in rows)
    gradient_values = _zeros_like(grid_values, geometry.shape) if "grid_values" in active else None
    coordinate_gradient = [[0.0, 0.0, 0.0] for _ in rows] if "coordinates" in active else None
    for row_index, row in enumerate(rows):
        sample, gradient, weights = _trilinear_one(grid_values, row, geometry, differentiate="coordinates" in active, return_weights=True)
        if gradient_values is not None:
            for i, j, k, weight in weights:
                _add_grid_value(gradient_values, i, j, k, weight)
        if coordinate_gradient is not None:
            coordinate_gradient[row_index] = list(gradient)

    def pullback(cotangent: Any) -> dict[str, Any]:
        if cotangent is ZERO:
            return {name: ZERO for name in wrt}
        if len(rows) == 1 and isinstance(cotangent, (int, float)) and not isinstance(cotangent, bool):
            cotangents = [cotangent]
        else:
            try:
                cotangents = list(cotangent)
            except TypeError as exc:
                raise TypeError("interpolate_grid cotangent must have one value per coordinate") from exc
        if len(cotangents) != len(rows):
            raise ValueError("interpolate_grid cotangent must have one value per coordinate")
        result: dict[str, Any] = {}
        if gradient_values is not None:
            scaled = _zeros_like(grid_values, geometry.shape)
            for row, coordinate in enumerate(rows):
                _, _, weights = _trilinear_one(grid_values, coordinate, geometry, differentiate=False, return_weights=True)
                factor = float(cotangents[row])
                for i, j, k, weight in weights:
                    _add_grid_value(scaled, i, j, k, factor * weight)
            result["grid_values"] = scaled
        if coordinate_gradient is not None:
            result["coordinates"] = _restore_rows(original_coordinates, [[float(cotangents[i]) * component for component in coordinate_gradient[i]] for i in range(len(rows))])
        return result

    return values, pullback


@rules.jvp_for(score_affinity_maps)
def _score_affinity_maps_jvp(tangents: dict[str, Any], maps: AffinityMaps | Mapping[Any, Any], coordinates: Any, atom_types: Any, *, translation: Any = (0.0, 0.0, 0.0), rotation: Any = (0.0, 0.0, 0.0), clip_value: float | None = 1000.0, center: Any = None, spacing: Any = None) -> tuple[float, Any]:
    supported = {"maps", "coordinates", "translation", "rotation"}
    unsupported = set(tangents) - supported
    if unsupported:
        raise UnsupportedWrt(score_affinity_maps, unsupported, supported=supported)
    family = _map_family(maps, center=center, spacing=spacing)
    rows = _coordinate_rows(coordinates)
    types = tuple(_as_type(value) for value in atom_types)
    if len(types) != len(rows):
        raise ValueError("atom_types must have one type per coordinate")
    shift = _finite_vector(translation, "translation")
    angle = _finite_vector(rotation, "rotation")
    active = frozenset(name for name, tangent in tangents.items() if tangent is not ZERO)
    value, local_gradient, translation_gradient, rotation_gradient, _ = _score_map_linearisation(family, rows, types, shift, angle, clip_value=clip_value, active=active)
    directional = 0.0
    if tangents.get("maps", ZERO) is not ZERO:
        tangent_family = _map_values_for_tangent(tangents["maps"], family)
        transformed = transform_pose(rows, shift, angle)
        for atom_type, point in zip(types, transformed):
            raw, _, _ = _trilinear_one(family.grid(atom_type), point, family, differentiate=False, return_weights=True)
            tangent_raw, _, _ = _trilinear_one(tangent_family.grid(atom_type), point, tangent_family, differentiate=False, return_weights=True)
            directional += _curl_value(raw, clip_value)[1] * tangent_raw
    if tangents.get("coordinates", ZERO) is not ZERO:
        tangent_rows = _coordinate_rows(tangents["coordinates"])
        if len(tangent_rows) != len(rows):
            raise ValueError("coordinates tangent must have shape (N, 3)")
        directional += sum(local_gradient[i][k] * tangent_rows[i][k] for i in range(len(rows)) for k in range(3))
    if tangents.get("translation", ZERO) is not ZERO:
        tangent_shift = _finite_vector(tangents["translation"], "translation tangent")
        directional += sum(translation_gradient[k] * tangent_shift[k] for k in range(3))
    if tangents.get("rotation", ZERO) is not ZERO:
        tangent_angle = _finite_vector(tangents["rotation"], "rotation tangent")
        directional += sum(rotation_gradient[k] * tangent_angle[k] for k in range(3))
    return value, directional


@rules.vjp_for(score_affinity_maps)
def _score_affinity_maps_vjp(wrt: tuple[str, ...], maps: AffinityMaps | Mapping[Any, Any], coordinates: Any, atom_types: Any, *, translation: Any = (0.0, 0.0, 0.0), rotation: Any = (0.0, 0.0, 0.0), clip_value: float | None = 1000.0, center: Any = None, spacing: Any = None) -> tuple[float, Any]:
    supported = {"maps", "coordinates", "translation", "rotation"}
    unsupported = set(wrt) - supported
    if unsupported:
        raise UnsupportedWrt(score_affinity_maps, unsupported, supported=supported)
    family = _map_family(maps, center=center, spacing=spacing)
    rows = _coordinate_rows(coordinates)
    types = tuple(_as_type(value) for value in atom_types)
    if len(types) != len(rows):
        raise ValueError("atom_types must have one type per coordinate")
    shift = _finite_vector(translation, "translation")
    angle = _finite_vector(rotation, "rotation")
    original_coordinates = coordinates
    value, local_gradient, translation_gradient, rotation_gradient, map_gradient = _score_map_linearisation(family, rows, types, shift, angle, clip_value=clip_value, active=frozenset(wrt))

    def pullback(cotangent: Any) -> dict[str, Any]:
        if cotangent is ZERO:
            return {name: ZERO for name in wrt}
        factor = float(cotangent)
        result: dict[str, Any] = {}
        if "maps" in wrt:
            assert map_gradient is not None
            scaled = {key: _zeros_like(grid, _shape3(grid)) for key, grid in map_gradient.values.items()}
            for key, grid in map_gradient.values.items():
                source = scaled[key]
                for i in range(map_gradient.shape[0]):
                    for j in range(map_gradient.shape[1]):
                        for k in range(map_gradient.shape[2]):
                            value_at = _grid_value(grid, i, j, k)
                            _add_grid_value(source, i, j, k, factor * value_at)
            result["maps"] = _restore_map_mapping(maps, family, scaled)
        if "coordinates" in wrt:
            result["coordinates"] = _restore_rows(original_coordinates, [[factor * component for component in row] for row in local_gradient])
        if "translation" in wrt:
            result["translation"] = _restore_vector(translation, tuple(factor * component for component in translation_gradient))
        if "rotation" in wrt:
            result["rotation"] = _restore_vector(rotation, tuple(factor * component for component in rotation_gradient))
        return result

    return value, pullback


def load_maps(map_prefix: str | Path) -> AffinityMaps:
    """Load a Vina ``<prefix>.<atom-type>.map`` family and its provenance."""
    prefix = str(map_prefix)
    paths = sorted(Path(path) for path in glob.glob(prefix + ".*.map"))
    if not paths:
        raise FileNotFoundError(f"no affinity maps found with prefix {prefix!r}")
    values: dict[int, Any] = {}
    center: tuple[float, float, float] | None = None
    spacing: tuple[float, float, float] | None = None
    box_size: tuple[float, float, float] | None = None
    for path in paths:
        suffix = path.name[: -len(".map")].rsplit(".", 1)[-1]
        try:
            atom_type = _as_type(suffix)
        except (TypeError, ValueError):
            # AD4 maps use names such as C/OA and are not an XS/Vina map
            # family.  Ignore them rather than silently exposing wrong types.
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        if len(lines) < 7:
            raise ValueError(f"map file {path} has an incomplete header")
        headers = {}
        for line in lines[:6]:
            fields = line.split(maxsplit=1)
            if len(fields) == 2:
                headers[fields[0].upper()] = fields[1]
        try:
            step = float(headers["SPACING"])
            dimensions = tuple(int(item) for item in headers["NELEMENTS"].split())
            map_center = tuple(float(item) for item in headers["CENTER"].split())
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"map file {path} has invalid Vina provenance header") from exc
        if len(dimensions) != 3 or len(map_center) != 3 or step <= 0:
            raise ValueError(f"map file {path} has invalid grid dimensions/provenance")
        if any(item < 1 for item in dimensions):
            raise ValueError(f"map file {path} must have at least two samples per axis")
        shape = tuple(item + 1 for item in dimensions)
        data: list[float] = []
        for line in lines[6:]:
            if line.strip():
                try:
                    value = float(line.split()[0])
                except ValueError as exc:
                    raise ValueError(f"map file {path} contains a non-numeric value") from exc
                if not math.isfinite(value):
                    raise ValueError(f"map file {path} contains a non-finite value")
                data.append(value)
        expected = shape[0] * shape[1] * shape[2]
        if len(data) != expected:
            raise ValueError(f"map file {path} contains {len(data)} values; expected {expected}")
        array = [[[0.0 for _ in range(shape[2])] for _ in range(shape[1])] for _ in range(shape[0])]
        cursor = 0
        for z in range(shape[2]):
            for y in range(shape[1]):
                for x in range(shape[0]):
                    array[x][y][z] = data[cursor]
                    cursor += 1
        current_center = tuple(map_center)
        current_spacing = (step,) * 3
        current_box = tuple(dimensions[i] * step for i in range(3))
        if center is None:
            center, spacing, box_size = current_center, current_spacing, current_box
        elif not (_same_vector(current_center, center) and _same_vector(current_spacing, spacing) and _same_vector(current_box, box_size)):
            raise ValueError("all affinity maps must share center, spacing and dimensions")
        base_type = _XS_BASE_TYPE.get(atom_type, atom_type)
        values[base_type] = array
    if not values:
        raise ValueError(f"no X-Score/Vina affinity maps found with prefix {prefix!r}")
    return AffinityMaps(values, center=center, spacing=spacing, box_size=box_size, source=prefix)


__all__ = [
    "AffinityGrid", "AffinityMaps", "GridMap", "MapFamily", "GridBoundaryError", "load_maps",
    "interpolate_grid", "trilinear_interpolate", "interpolate", "interpolate_map", "interpolate_maps",
    "transform_pose", "score_affinity_maps", "score_maps", "grid_score", "score_grid", "map_score", "pose_score",
]
