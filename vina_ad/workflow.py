"""Workflow adapters with an explicit boundary around the compiled Vina API.

``run_demo`` is retained as a deterministic toy diagnostic.  It is deliberately
not described as a real docking workflow.  ``run_official_workflow`` loads the
official example PDBQT files, invokes the real ``vina.Vina`` binding when it is
installed, and runs the restricted sidecar through public AD interfaces on the
sourced receptor/ligand pair set.  Supplied-map interpolation is covered by
`vina_ad.score_affinity_maps`; full map generation and search-loop derivatives
remain out of scope.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from . import DEFAULT_VINA_WEIGHTS, grad, jvp, score_coordinates, value_and_grad


_PDBQT_TO_XS = {
    "C": 0,
    "A": 1,
    "N": 2,
    "NA": 4,
    "OA": 8,
    "SA": 8,
    "HD": 3,
    "O": 6,
    "S": 10,
    "F": 12,
    "CL": 13,
    "BR": 14,
    "I": 15,
}


def _read_pdbqt(path: Path) -> tuple[list[tuple[float, float, float]], list[int]]:
    coordinates: list[tuple[float, float, float]] = []
    atom_types: list[int] = []
    for line in path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        fields = line.split()
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            pdbqt_type = fields[-1].upper()
        except (ValueError, IndexError) as exc:
            raise ValueError(f"invalid PDBQT atom line in {path}") from exc
        if pdbqt_type not in _PDBQT_TO_XS:
            raise ValueError(f"unsupported PDBQT atom type {pdbqt_type!r} in {path}")
        coordinates.append(xyz)
        atom_types.append(_PDBQT_TO_XS[pdbqt_type])
    if not coordinates:
        raise ValueError(f"no ATOM/HETATM records found in {path}")
    return coordinates, atom_types


def _read_torsion_count(path: Path) -> float:
    """Read the fixed ``TORSDOF`` state from an official ligand file."""
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields and fields[0].upper() == "TORSDOF" and len(fields) == 2:
            try:
                value = float(fields[1])
            except ValueError as exc:
                raise ValueError(f"invalid TORSDOF record in {path}") from exc
            if value < 0:
                raise ValueError(f"invalid TORSDOF record in {path}")
            return value
    return 0.0


def run_demo() -> dict[str, float | str]:
    """Run a labelled, install-independent toy pair diagnostic."""
    coordinates = ((0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (0.0, 4.0, 0.0))
    atom_types = (0, 0, 0)
    value = score_coordinates(coordinates, atom_types)
    gradients = grad(score_coordinates, coordinates, atom_types, wrt="coordinates")["coordinates"]
    norm = sum(component * component for row in gradients for component in row) ** 0.5
    return {"status": "toy-diagnostic", "score": value, "gradient_l2": norm, "n_atoms": float(len(coordinates))}


def run_official_workflow(
    receptor: str | Path | None = None,
    ligand: str | Path | None = None,
) -> dict[str, Any]:
    """Run the sourced official Python workflow and a public AD replay.

    The default files are resolved from the source checkout.  Installed wheels
    intentionally do not bundle the immutable upstream snapshot, so callers
    should pass explicit PDBQT paths in an installed environment.  If the real
    binding or paths are unavailable, the returned status is ``deferred`` with
    a reason rather than a synthetic completion.
    """
    if receptor is None:
        receptor = Path(__file__).resolve().parents[1] / "upstream/example/python_scripting/1iep_receptor.pdbqt"
    if ligand is None:
        ligand = Path(__file__).resolve().parents[1] / "upstream/example/python_scripting/1iep_ligand.pdbqt"
    receptor_path, ligand_path = Path(receptor), Path(ligand)
    if not receptor_path.exists() or not ligand_path.exists():
        return {"status": "deferred", "reason": "official PDBQT inputs are unavailable after installation"}
    try:
        from vina import Vina
    except (ImportError, ModuleNotFoundError) as exc:
        return {"status": "deferred", "reason": f"real vina binding unavailable: {exc}"}
    try:
        vina = Vina(sf_name="vina", cpu=1, seed=1, verbosity=0)
        vina.set_receptor(str(receptor_path))
        vina.set_ligand_from_file(str(ligand_path))
        vina.compute_vina_maps(center=[15.190, 53.903, 16.917], box_size=[20, 20, 20])
        real_score = float(vina.score()[0])
        # Execute one bounded local-minimisation iteration from the sourced
        # official pose.  This supplies a reproducible workflow-iteration
        # metric without invoking the expensive stochastic global search.
        optimized_score = float(vina.optimize(max_steps=1)[0])
    except Exception as exc:  # binding errors are an explicit deferred result
        return {"status": "deferred", "reason": f"real Vina workflow failed: {exc}"}
    receptor_coordinates, receptor_types = _read_pdbqt(receptor_path)
    ligand_coordinates, ligand_types = _read_pdbqt(ligand_path)
    torsion_count = _read_torsion_count(ligand_path)

    # The installed sidecar cannot claim to reproduce Vina's grid scorer.  It
    # can run the faithful SF_VINA pair kernel on the complete sourced
    # receptor/ligand cross topology.  The two files' atoms remain separate so
    # the active direction below is a rigid translation of the ligand.
    receptor_count = len(receptor_coordinates)
    pair_coordinates = tuple(receptor_coordinates + ligand_coordinates)
    pair_types = tuple(receptor_types + ligand_types)
    pairs = tuple(
        (i, receptor_count + j)
        for i in range(receptor_count)
        for j in range(len(ligand_coordinates))
    )
    ligand_direction = tuple((0.0, 0.0, 0.0) for _ in receptor_coordinates) + tuple(
        (1.0, 0.0, 0.0) for _ in ligand_coordinates
    )
    ad_value, gradients = value_and_grad(
        score_coordinates,
        pair_coordinates,
        pair_types,
        pairs=pairs,
        torsion_count=torsion_count,
        wrt=("coordinates", "weights"),
    )
    coordinate_gradient = gradients["coordinates"]
    weight_gradient = gradients["weights"]
    directional = sum(
        row[0] * direction[0] + row[1] * direction[1] + row[2] * direction[2]
        for row, direction in zip(coordinate_gradient, ligand_direction)
    )
    jvp_value, jvp_direction = jvp(
        score_coordinates,
        pair_coordinates,
        pair_types,
        pairs=pairs,
        torsion_count=torsion_count,
        tangents={"coordinates": ligand_direction},
    )
    step = 1e-5
    plus = tuple(
        tuple(value + step * delta for value, delta in zip(row, direction))
        for row, direction in zip(pair_coordinates, ligand_direction)
    )
    minus = tuple(
        tuple(value - step * delta for value, delta in zip(row, direction))
        for row, direction in zip(pair_coordinates, ligand_direction)
    )
    finite_difference = (
        score_coordinates(plus, pair_types, pairs=pairs, torsion_count=torsion_count)
        - score_coordinates(minus, pair_types, pairs=pairs, torsion_count=torsion_count)
    ) / (2.0 * step)
    coordinate_gradient_l2 = sum(
        component * component for row in coordinate_gradient for component in row
    ) ** 0.5
    weight_gradient_l2 = sum(component * component for component in weight_gradient) ** 0.5
    derivative_abs_error = abs(float(jvp_direction) - finite_difference)
    duality_abs_error = abs(float(jvp_direction) - directional)
    comparable = len(receptor_coordinates) == 1 and len(ligand_coordinates) == 1
    result = {
        # Completion means both the sourced real-Vina call and the public AD
        # replay ran.  ``deviation_bound`` is only meaningful for the existing
        # one-atom oracle fixture; multi-atom grid-vs-pair scores are reported
        # quantitatively without pretending they are the same model.
        "status": "completed",
        "real_vina_score": real_score,
        "real_vina_optimized_score": optimized_score,
        "restricted_pair_score": ad_value,
        "absolute_deviation": abs(real_score - ad_value),
        "deviation_bound": 0.05 if comparable else None,
        "receptor_atoms": float(receptor_count),
        "ligand_atoms": float(len(ligand_coordinates)),
        "interaction_pairs": float(len(pairs)),
        "workflow_scale": {
            "receptor_atoms": float(receptor_count),
            "ligand_atoms": float(len(ligand_coordinates)),
            "interaction_pairs": float(len(pairs)),
        },
        "iterations": {
            "vina_score": 1.0,
            "vina_optimize": 1.0,
            "ad_value_and_grad": 1.0,
            "ad_jvp": 1.0,
            "finite_difference": 2.0,
        },
        "torsion_count": torsion_count,
        "ad_primal": ad_value,
        "jvp_primal": jvp_value,
        "ad_coordinate_gradient_l2": coordinate_gradient_l2,
        "ad_weight_gradient_l2": weight_gradient_l2,
        "ad_directional_derivative": float(jvp_direction),
        "gradient_directional_derivative": directional,
        "finite_difference_directional_derivative": finite_difference,
        "derivative_abs_error": derivative_abs_error,
        "duality_abs_error": duality_abs_error,
        "finite_difference_step": step,
        "vina_score_evaluations": 1.0,
        "vina_optimize_max_steps": 1.0,
        "vina_optimize_evaluations": 1.0,
        "ad_value_and_grad_evaluations": 1.0,
        "ad_jvp_evaluations": 1.0,
        "finite_difference_evaluations": 2.0,
        "remaining_coverage": (
            "supplied-map grid interpolation and pose-transform derivatives are "
            "implemented; full Vina search-loop, map-generation, and stochastic "
            "docking derivatives remain deferred; AD also covers fixed SF_VINA "
            "receptor/ligand cross pairs"
        ),
    }
    if not comparable:
        result["reason"] = (
            "real Vina and the sourced restricted AD replay completed; their "
            "primal scores are different models"
        )
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receptor", type=Path, default=None)
    parser.add_argument("--ligand", type=Path, default=None)
    args = parser.parse_args()
    result = run_official_workflow(args.receptor, args.ligand)
    if result["status"] == "deferred":
        print("status=deferred reason={reason}".format(**result))
    else:
        print("status=completed " + json.dumps(result, sort_keys=True))
