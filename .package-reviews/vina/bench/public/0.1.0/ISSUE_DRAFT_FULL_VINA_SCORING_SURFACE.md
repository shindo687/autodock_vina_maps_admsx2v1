# Differentiable full Vina scoring surface: SF_VINA recombination, per-term decomposition, and AD4/Vinardo families

## Summary

The reviewed sidecar exposes a differentiable path only for a caller-assembled
replay of the six pair potentials (`score_coordinates` with explicit `pairs`,
`atom_types`, fixed `torsion_count`). The pinned upstream scorer exposes the
complete scoring surface — full SF_VINA recombination over precomputed
interactions, the AD4 and Vinardo scoring families via `sf_name`, weight
reconfiguration via `set_weights`, and per-term weighted accumulation — none of
which has a registered ChainRules path here. This blocks every task that needs
derivatives of the real Vina score (weight calibration, pose-vs-affinity
tuning, flexible/AD4 docking, free-energy surrogates, cross-docking
sensitivity).

**Classification:** `upstream_parity_gap`

## Versions and scope

- Upstream software: AutoDock Vina, pinned commit `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645` (Python binding 1.2.x series)
- AD package: `https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj` @ `955ec33c019cdcdbc4c0aeaa0c5394accb9de95c` (`vina_ad` 0.1.0)
- Capability review: `.package-reviews/vina/bench/public/0.1.0/TASK_CAPABILITY_REVIEW.md`
- Bench: public (`ad-software-public-benchmark` @ `96c49e5a033b311416f18d8eba5d14570714e1ec`)

## Related tasks and papers

| Task | What it needs | Paper |
|---|---|---|
| [`vina.score_weight_affinity_calibration`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L168) | per-term energies as regression features + weight VJP | [BPC 2018](https://doi.org/10.1016/j.bpc.2018.05.010) |
| [`vina.pose_vs_affinity_weight_optimization`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L185) | weight VJP over full score for pose/affinity losses | [Comput Biol Chem 2016](https://doi.org/10.1016/j.compbiolchem.2016.04.005) |
| [`vina.target_family_robust_score_tuning`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L202) | weight VJP at family scale with full-term coverage | [JCAMD 2009](https://doi.org/10.1007/s10822-009-9276-1) |
| [`vina.multiobjective_score_training`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L219) | term-vector features into a multi-task loss | [Chem Sci 2023](https://doi.org/10.1039/d3sc02044d) |
| [`vina.vina_term_free_energy_surrogate`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L236) | interpretable term vectors as surrogate inputs | [J Cheminform 2021](https://doi.org/10.1186/s13321-021-00536-w) |
| [`vina.end_to_end_pose_parameter_optimization`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L49) | full-score objective under 6+k pose variables | [Brief Bioinform 2022](https://doi.org/10.1093/bib/bbac520) |
| [`vina.native_bfgs_pose_gradient_local_search`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L83) | native-score gradients inside Vina's refinement | [J Comput Chem 2010](https://doi.org/10.1002/jcc.21334) |
| [`vina.flexible_receptor_sidechain_gradient_docking`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L117) | score derivatives under selective flexibility | [JCAMD 2008](https://doi.org/10.1007/s10822-007-9148-5) |
| [`vina.ad4_energy_gradient_alternative_scoring`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L151) | AD4 electrostatic/desolvation/H-bond term derivatives | [JCIM 2021](https://doi.org/10.1021/acs.jcim.1c00203) |
| [`vina.cross_docking_receptor_sensitivity`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L270) | full-score derivatives vs receptor/ligand conformations | [Mol Sim 2014](https://doi.org/10.1080/08927022.2014.917300) |
| [`vina.peptide_flexible_pose_gradient`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L389) | full-score gradients for peptide internal coordinates | [Protein Sci 2002](https://doi.org/10.1110/ps.0202302) |

## Evidence of the gap

### Upstream operation

- Scoring families are selectable in the pinned binding:
  [`build/python/vina/vina.py#L19`](https://github.com/ccsb-scripps/AutoDock-Vina/blob/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645/build/python/vina/vina.py#L19)
  (`sf_name='vina'|'vinardo'|'ad4'`),
  [`#L407`](https://github.com/ccsb-scripps/AutoDock-Vina/blob/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645/build/python/vina/vina.py#L407)
  (`score()`), with weight reconfiguration at
  [`#L209`](https://github.com/ccsb-scripps/AutoDock-Vina/blob/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645/build/python/vina/vina.py#L209)
  (`set_weights`).
- The scorer accumulates per-term weighted potentials:
  [`src/lib/scoring_function.h#L123`](https://github.com/ccsb-scripps/AutoDock-Vina/blob/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645/src/lib/scoring_function.h#L123)
  (`acc += m_weights[i] * m_potentials[i]->eval(...)`), weights exposed at
  [`#L159`](https://github.com/ccsb-scripps/AutoDock-Vina/blob/3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645/src/lib/scoring_function.h#L159)
  — the per-term decomposition is a thin slicing of this accumulation.
- AD4 terms live in the same snapshot (`src/lib/ad4cache.{h,cpp}`).

### AD boundary

- The sidecar's only differentiable entry point is the caller-assembled pair
  replay:
  [`vina_ad/core.py#L352`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj/-/blob/main/vina_ad/core.py#L352)
  (`score_coordinates`; "coordinates and weights are the active differentiable
  inputs", `atom_types`/`pairs`/`torsion_count` fixed state),
  with exactly six SF_VINA pair terms in
  [`core.py#L219`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj/-/blob/main/vina_ad/core.py#L219).
- The full scorer is explicitly deferred:
  [`support_status.md#L11`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj/-/blob/main/vina_ad/support_status.md#L11)
  (`vina.Vina.score ... deferred; restricted pair replay ... explicitly
  separate`). The sidecar kernel deliberately does not ship the compiled
  `vina_wrapper` (`vina_ad/requirements.md#L15-L16`).
- Term sums are computed internally for the weights VJP
  (`core.py#L308-L346`) but never returned: no per-term output API exists.

### Minimal reproduction

```python
# vina_ad 0.1.0 @ 955ec33c
import vina_ad
print([n for n in dir(vina_ad) if any(k in n.lower() for k in ("term", "ad4", "vinardo"))])
# -> []  : no scoring-family selection, no per-term decomposition

from vina_ad import score_coordinates, value_and_grad, jvp
coords = ((0.0, 0.0, 0.0), (2.9, 0.0, 0.0))          # two carbons, 2.9 A
s = score_coordinates(coords, (0, 0), pairs=((0, 1),))  # only reruns this pair assembly
g = value_and_grad(score_coordinates, wrt="weights")(coords, (0, 0), pairs=((0, 1),))
# -> build succeeds, but there is no API for
#    Vina("ad4").score(), per-term vectors, or recombined native scores
#    (core.py L422: supported wrt keys are exactly {"coordinates", "weights"})
```

Expected: a registered path whose forward reproduces the upstream scorer for
the selected family, with per-term outputs and weight/coordinate JVP/VJP.
Observed: pair replay only; everything else deferred or unregistered.

## Expected capability

- A `score`-family API that selects the scoring family (`vina`/`vinardo`/`ad4`)
  as an explicit contract input, with primal parity against the pinned native
  scorer under identical inputs.
- Per-term decomposition (vector per pose, same length as the family's weight
  vector) as a public differentiably-composed output, so regression/calibration
  tasks consume features without reimplementing the scorer.
- JVP/VJP preserved for `coordinates` and `weights` (existing pair rules should
  remain a sub-path); `RuleNotFound`/`NonDifferentiablePoint` boundaries kept
  at the documented non-smooth points.

## Acceptance criteria

- Weighted-sum reconstruction of the term vectors equals the composed score
  (term consistency check).
- Primal parity vs the binding `vina.Vina.score()` on a published example
  within the tolerances used by the existing oracle tests.
- Weight and coordinate gradients match independent finite differences on the
  new surface (direction, duality, and zero checks).
- Existing tests remain passing.

## Non-goals

- Differentiable maps/grid interpolation (tracked by the separate maps parity
  issue).
- Driving or differentiating through Vina's internal search loop (draft
  `POSE_DOF_SEARCH_FRONTEND`).
- Hydration/metal/RNA scoring extensions and mutation ΔΔG front-ends (draft
  issues, upstream lacks the forward surface).