import math
from pathlib import Path

import pytest

import vina_ad


def _linear_grid(shape=(4, 4, 4), center=(0.0, 0.0, 0.0), spacing=1.0):
    origin = tuple(center[i] - 0.5 * (shape[i] - 1) * spacing for i in range(3))
    return [
        [
            [2.0 * (origin[0] + i * spacing) + 3.0 * (origin[1] + j * spacing) + 5.0 * (origin[2] + k * spacing) for k in range(shape[2])]
            for j in range(shape[1])
        ]
        for i in range(shape[0])
    ]


def test_trilinear_interpolation_and_coordinate_map_derivative():
    grid = _linear_grid()
    point = ((-0.23, 0.17, 0.41),)
    assert vina_ad.interpolate_grid(grid, point, center=(0.0, 0.0, 0.0), spacing=1.0) == pytest.approx((2 * point[0][0] + 3 * point[0][1] + 5 * point[0][2],))
    value, pullback = vina_ad.vjp(vina_ad.interpolate_grid, grid, point, center=(0.0, 0.0, 0.0), spacing=1.0, wrt="coordinates")
    assert value[0] == pytest.approx(2 * point[0][0] + 3 * point[0][1] + 5 * point[0][2])
    assert pullback(1.0)["coordinates"][0] == pytest.approx([2.0, 3.0, 5.0])


def test_map_value_jvp_and_vjp_are_trilinear_weights():
    grid = _linear_grid()
    point = ((-0.23, 0.17, 0.41),)
    tangent = [[[1.0 for _ in range(4)] for _ in range(4)] for _ in range(4)]
    _, direction = vina_ad.jvp(vina_ad.interpolate_grid, grid, point, center=(0, 0, 0), spacing=1, tangents={"grid_values": tangent})
    assert direction == pytest.approx((1.0,))
    _, pullback = vina_ad.vjp(vina_ad.interpolate_grid, grid, point, center=(0, 0, 0), spacing=1, wrt="grid_values")
    gradient = pullback(1.0)["grid_values"]
    assert sum(sum(sum(row) for row in plane) for plane in gradient) == pytest.approx(1.0)


def test_affinity_maps_pose_chain_rule_matches_finite_difference():
    grid = _linear_grid()
    maps = vina_ad.AffinityMaps({"C_H": grid}, center=(0, 0, 0), spacing=1)
    coordinates = ((-0.23, 0.17, 0.41),)
    translation = (0.13, -0.11, 0.07)
    rotation = (0.06, -0.04, 0.03)
    value, pullback = vina_ad.vjp(vina_ad.score_affinity_maps, maps, coordinates, (0,), translation=translation, rotation=rotation, wrt=("coordinates", "translation", "rotation"))
    gradients = pullback(1.0)
    h = 1e-6
    for name, base, expected in (("translation", translation, gradients["translation"]), ("rotation", rotation, gradients["rotation"])):
        for component in range(3):
            plus = list(base)
            minus = list(base)
            plus[component] += h
            minus[component] -= h
            plus_value = vina_ad.score_affinity_maps(maps, coordinates, (0,), **{name: plus, "rotation" if name == "translation" else "translation": rotation if name == "translation" else translation})
            minus_value = vina_ad.score_affinity_maps(maps, coordinates, (0,), **{name: minus, "rotation" if name == "translation" else "translation": rotation if name == "translation" else translation})
            assert expected[component] == pytest.approx((plus_value - minus_value) / (2 * h), abs=2e-5)
    assert math.isfinite(value)


def test_six_vector_pose_wrapper_and_aliases():
    grid = _linear_grid()
    maps = vina_ad.AffinityMaps({0: grid}, center=(0, 0, 0), spacing=1)
    coordinates = ((-0.23, 0.17, 0.41),)
    pose = (0.13, -0.11, 0.07, 0.06, -0.04, 0.03)
    expected = vina_ad.score_affinity_maps(maps, coordinates, (0,), translation=pose[:3], rotation=pose[3:])
    assert vina_ad.score_pose(maps, coordinates, (0,), pose=pose) == pytest.approx(expected)
    assert vina_ad.pose_score is vina_ad.score_pose
    _, pullback = vina_ad.vjp(vina_ad.score_pose, maps, coordinates, (0,), pose=pose, wrt="pose")
    gradient = pullback(1.0)["pose"]
    assert len(gradient) == 6 and all(math.isfinite(item) for item in gradient)


def test_grid_boundaries_are_explicit_and_map_only_derivative_is_defined():
    grid = _linear_grid()
    with pytest.raises(vina_ad.GridBoundaryError):
        vina_ad.interpolate_grid(grid, ((2.0, 0.0, 0.0),), center=(0, 0, 0), spacing=1)
    with pytest.raises(vina_ad.NonDifferentiablePoint):
        vina_ad.vjp(vina_ad.interpolate_grid, grid, ((-0.5, 0.17, 0.41),), center=(0, 0, 0), spacing=1, wrt="coordinates")
    maps = vina_ad.AffinityMaps({0: grid}, center=(0, 0, 0), spacing=1)
    value, pullback = vina_ad.vjp(vina_ad.score_affinity_maps, maps, ((-0.5, 0.17, 0.41),), (0,), wrt="maps")
    assert math.isfinite(value)
    assert isinstance(pullback(1.0)["maps"], vina_ad.AffinityMaps)


def test_load_maps_records_vina_header_and_matches_map_score(tmp_path: Path):
    vina = pytest.importorskip("vina")
    receptor = tmp_path / "receptor.pdbqt"
    ligand = tmp_path / "ligand.pdbqt"
    def atom(serial, x):
        return f"ATOM  {serial:5d}  C   LIG A   1    {x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00     0.000 C\n"
    receptor.write_text(atom(1, 0.0))
    ligand.write_text("ROOT\n" + atom(1, 0.73) + "ENDROOT\nTORSDOF 0\n")
    prefix = tmp_path / "maps"
    upstream = vina.Vina(sf_name="vina", cpu=1, seed=1, verbosity=0, no_refine=True)
    upstream.set_receptor(str(receptor))
    upstream.set_ligand_from_file(str(ligand))
    upstream.compute_vina_maps(center=[0.0, 0.0, 0.0], box_size=[4.0, 4.0, 4.0], spacing=1.0)
    upstream.write_maps(str(prefix), overwrite=True)
    maps = vina_ad.load_maps(prefix)
    assert maps.center == (0.0, 0.0, 0.0)
    assert maps.spacing == (1.0, 1.0, 1.0)
    expected = float(upstream.score()[0])
    actual = vina_ad.score_affinity_maps(maps, ((0.73, 0.0, 0.0),), (0,))
    # Vina writes maps with four significant decimal digits, so one
    # interpolation can differ from the in-memory map by a few 1e-4.
    assert actual == pytest.approx(expected, abs=2e-3)
