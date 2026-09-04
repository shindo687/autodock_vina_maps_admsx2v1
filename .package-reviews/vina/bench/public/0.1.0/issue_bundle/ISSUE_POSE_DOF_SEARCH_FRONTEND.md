<!-- generated from ISSUE_DRAFT_POSE_DOF_SEARCH_FRONTEND.md; edit the source draft and regenerate -->
<!-- classification: scientific_extension; submit_enabled: false -->
# Pose DOF front-end and search composition (6+k DOF, flexible sidechains, macrocycle/peptide closure, line-search instrumentation)

Draft only — not submitted to GitLab.

## Summary

Tasks that optimize or sample over pose variables beyond atomic coordinates —
translation/rotation/torsion pose vectors, selective-flexibility receptor
sidechains, macrocycle ring closure, peptide backbone closure, and
instrumentation of Vina's internal line search — currently have no
differentiable front-end in `vina_ad`. The pair replay provides coordinate and
weight derivatives, and a caller *can* build transform layers on top, but the
sidecar offers no DOF→coordinates layer, no closure constraints, and no path
into the upstream search loop. Neither the upstream binding nor the sidecar
promises such a front-end as an API, so this is recorded as an extension
candidate rather than an upstream-parity gap.

**Classification:** `scientific_extension`

## Versions and scope

- Upstream software: AutoDock Vina, pinned commit `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645`
- AD package: `https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj` @ `955ec33c019cdcdbc4c0aeaa0c5394accb9de95c` (`vina_ad` 0.1.0)
- Capability review: `.package-reviews/vina/bench/public/0.1.0/TASK_CAPABILITY_REVIEW.md`

## Related tasks and papers

| Task | What it needs | Paper |
|---|---|---|
| [`vina.end_to_end_pose_parameter_optimization`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L49) | 6+k pose-vector front-end | [Brief Bioinform 2022](https://doi.org/10.1093/bib/bbac520) |
| [`vina.native_bfgs_pose_gradient_local_search`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L83) | gradient supply into Vina's refinement | [J Comput Chem 2010](https://doi.org/10.1002/jcc.21334) |
| [`vina.gradient_consistent_line_search_sampling`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L100) | line-search instrumentation | [ChemRxiv 2021](https://doi.org/10.26434/chemrxiv.15004371/v1) |
| [`vina.flexible_receptor_sidechain_gradient_docking`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L117) | flex-residue torsion DOF | [JCAMD 2008](https://doi.org/10.1007/s10822-007-9148-5) |
| [`vina.macrocycle_torsion_gradient_sampling`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L134) | ring-closure constraints | [JCAMD 2024](https://doi.org/10.1007/s10822-024-00574-0) |
| [`vina.peptide_flexible_pose_gradient`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L389) | backbone closure + torsion DOF | [Protein Sci 2002](https://doi.org/10.1110/ps.0202302) |

## Gap sketch

- `vina_ad` active inputs are `coordinates` and `weights` only
  (core.py L422); there is no pose-parameter or internal-coordinate layer.
- `vina.Vina.optimize`/`dock` are classed not suitable/deferred
  (support_status.md L12) and the upstream search loop is a compiled private
  implementation, so no registered path can observe per-step gradients.
- Closure constraints (macrocycle/peptide) and selective-flexibility DOF have
  no representation.

## Candidate design (if ever enabled)

- A `pose` front-end: rigid transform + torsion vector → coordinates, with
  JVP/VJP composed over `score_coordinates`.
- Closure: constrained torsion coordinates with an explicit violation term or
  Jacobian-projected updates.
- Search composition: expose per-step score/gradient hooks by re-driving a
  local BFGS/line-search in Python against the sidecar surface (upstream's
  loop itself is not instrumentable without patching).

## Non-goals

- Changing upstream Vina internals; claiming derivatives through the native
  compiled search loop; full paper-scale docking pipelines.
