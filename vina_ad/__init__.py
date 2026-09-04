"""Differentiable AutoDock-Vina sidecar with explicit ChainRules rules."""

from .core import DEFAULT_VINA_WEIGHTS, energy, score, score_coordinates
from .maps import (
    AffinityGrid,
    AffinityMaps,
    GridBoundaryError,
    grid_score,
    interpolate,
    interpolate_grid,
    interpolate_map,
    interpolate_maps,
    load_maps,
    pose_score,
    score_pose,
    score_affinity_maps,
    score_grid,
    map_score,
    score_maps,
    transform_pose,
    trilinear_interpolate,
)
from .protocol import NonDifferentiablePoint, RuleNotFound, UnsupportedWrt, ZERO
from .protocol import grad, jvp, value_and_grad, vjp

__version__ = "0.1.0"
__all__ = [
    "DEFAULT_VINA_WEIGHTS",
    "score_coordinates",
    "score",
    "energy",
    "jvp",
    "vjp",
    "grad",
    "value_and_grad",
    "ZERO",
    "RuleNotFound",
    "UnsupportedWrt",
    "NonDifferentiablePoint",
    "AffinityGrid",
    "AffinityMaps",
    "GridBoundaryError",
    "load_maps",
    "interpolate_grid",
    "trilinear_interpolate",
    "interpolate",
    "interpolate_map",
    "interpolate_maps",
    "transform_pose",
    "score_affinity_maps",
    "score_maps",
    "grid_score",
    "score_grid",
    "map_score",
    "pose_score",
    "score_pose",
]
