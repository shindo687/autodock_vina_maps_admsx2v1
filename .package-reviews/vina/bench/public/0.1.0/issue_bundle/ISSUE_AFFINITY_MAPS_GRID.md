<!-- generated from ISSUE_DRAFT_AFFINITY_MAPS_GRID.md; edit the source draft and regenerate -->
<!-- classification: upstream_parity_gap; submit_enabled: true -->
# Add differentiable affinity maps and grid interpolation for pose optimization

## Summary

The pinned upstream binding can precompute affinity maps
(`compute_vina_maps`) and load map families (`load_maps`); grid quantities are
the substrate of both Vina's own scoring and the grid-based pose optimization
workflows in this corpus. The reviewed sidecar has no map or interpolation API
at all, so any task whose derivative crosses an atom-type grid score — pose
optimization on atomic grids, receptor-ensemble and cross-docking
conformation sensitivities — has nothing to differentiate.

**Classification:** `upstream_parity_gap`

## Versions and scope

- Upstream software: AutoDock Vina, pinned commit `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645` (Python binding 1.2.x series)
- AD package: `https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj` @ `955ec33c019cdcdbc4c0aeaa0c5394accb9de95c` (`vina_ad` 0.1.0)
- Capability review: `.package-reviews/vina/bench/public/0.1.0/TASK_CAPABILITY_REVIEW.md`
- Bench: public (`ad-software-public-benchmark` @ `96c49e5a033b311416f18d8eba5d14570714e1ec`)

## Related tasks and papers

| Task | What it needs | Paper |
|---|---|---|
| [`vina.atomic_grid_score_pose_optimization`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L66) | interpolation of atom-type grid values + pose transform derivatives | [arXiv:1710.07400](https://doi.org/10.48550/arXiv.1710.07400) |
| [`vina.receptor_ensemble_sensitivity`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L253) | score sensitivity to receptor conformation via the map/model path | [Proteins 2020](https://doi.org/10.1002/prot.25899) |
| [`vina.cross_docking_receptor_sensitivity`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L270) | receptor/ligand conformation derivatives for cross-docking | [Mol Sim 2014](https://doi.org/10.1080/08927022.2014.917300) |

## Evidence of the gap

### Upstream operation

- Map precomputation and loading are public binding methods in the pinned
  snapshot:
  [`build/python/vina/vina.py#L232`](https://github.com/ccsb-scripps/AutoDock-Vina/blob/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645/build/python/vina/vina.py#L232)
  (`compute_vina_maps(center, box_size, spacing=0.375, ...)`),
  [`#L268`](https://github.com/ccsb-scripps/AutoDock-Vina/blob/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645/build/python/vina/vina.py#L268)
  (`load_maps`), plus `write_maps` at `#L280`. The grid machinery itself is in
  the same snapshot (`src/lib/grid.cpp`, `src/lib/array3d.h`, `src/lib/cache.cpp`).

### AD boundary

- The sidecar exposes no map/grid surface at all; its score entry point takes
  explicit pair lists and fixed atom types:
  [`vina_ad/core.py#L352`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj/-/blob/main/vina_ad/core.py#L352).
- Map-related binding methods are classed out of scope:
  [`vina_ad/support_status.md#L10`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj/-/blob/main/vina_ad/support_status.md#L10)
  (`map, pose and energy I/O methods ... not suitable/deferred`).

### Minimal reproduction

```python
# upstream (pinned binding) can build and reload maps:
#   Vina().compute_vina_maps(center, box_size); Vina().load_maps(prefix)
# sidecar @ 955ec33c — no grid layer:
import vina_ad
print([n for n in dir(vina_ad) if "map" in n.lower() or "grid" in n.lower()])
# -> []  : nothing to interpolate or differentiate through
```

Expected: a differentiable grid-interpolation score (maps + trilinear atom-type
interpolation + pose-transform composition). Observed: absent; only explicit
pair replay exists.

## Expected capability

- Load (or receive) atom-type affinity maps with recorded provenance
  (center/box/spacing), expose trilinear interpolation of the grid value per
  atom, and register JVP/VJP over grid values, atom coordinates, and the
  transforming pose parameters composed above it.
- Map-generation itself (`compute_vina_maps`) can stay an imperative upstream
  call; the differentiable contract begins at interpolation/composition.

## Acceptance criteria

- Interpolated score matches upstream map-based scoring at primal level on a
  published example.
- Gradients w.r.t. coordinates match independent finite differences of the
  interpolation plus map values.
- Cell-boundary/outside-box behavior raises `NonDifferentiablePoint` or an
  explicit boundary error, never a silent zero.
- Existing tests remain passing.

## Non-goals

- Differentiable compute_vina_maps generation (upstream precompute stays a
  discrete step).
- AD4/Vinardo term families and per-term decomposition (tracked by the full
  scoring-surface parity issue).
- Search-loop derivatives and DOF front-ends (draft `POSE_DOF_SEARCH_FRONTEND`).
