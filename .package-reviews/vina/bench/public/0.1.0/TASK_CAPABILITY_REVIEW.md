# TASK_CAPABILITY_REVIEW — AutoDock Vina / bench: public / vina_ad 0.1.0

## Review identity

| Item | Value |
|---|---|
| Benchmark corpus (explicit) | **public** (`bench=public`); the private corpus was intentionally **not reviewed** |
| Benchmark repository | `https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark.git`, branch `main`, commit `96c49e5a033b311416f18d8eba5d14570714e1ec` |
| Task ledger | `vina/vina-task-ledger.json` (same repo/commit), 22 tasks, all `existing_ad_acceleration` |
| Package under review | `https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj`, branch `main`, commit `955ec33c019cdcdbc4c0aeaa0c5394accb9de95c` (batch-accepted round-3 commit) |
| Package version | `0.1.0` (declared in `pyproject.toml`) |
| Upstream software | AutoDock Vina, fixed upstream commit `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645` (snapshot vendored under `upstream/`, unmodified; Python binding 1.2.x series) |
| Runtime used to freeze evidence | macOS/CPython 3.14 local read-only probe (see §Evidence probes); package self-tests and fix receipts from the remote build machines (`vina_ad/receipts/`) |

## Method note

Each of the 22 ledger tasks was mapped to a compact **task contract** (scientific
output, forward path, active inputs, derivative products). The ledger's
`task`/`implementation` was the first hypothesis. For the 7 papers whose PDFs
are part of this bench checkout, the abstract/methods of the relevant paper
were consulted; for the other 15 papers only the ledger record and the
paper manifest metadata were available. No full paper reproduction was
attempted for any task. Package capability was checked by reading the package
surface at the pinned commit, and a minimal bounded probe was run against the
actual `vina_ad` checkout (see below). The report therefore distinguishes
capability verdicts from reproduction claims via a separate
`paper_workflow_status` column.

### Evidence probes (run at review time, package commit 955ec33c)

- Import surface: `vina_ad` exports only
  `score_coordinates`/`score`/`energy`, `jvp`, `vjp`, `grad`, `value_and_grad`,
  plus constants/protocol types (`DEFAULT_VINA_WEIGHTS`, `ZERO`,
  `NonDifferentiablePoint`, `RuleNotFound`, `UnsupportedWrt`, `core`,
  `protocol`). No map, grid, per-term, AD4, optimization, or pose-parameter
  API exists.
- Minimal score: two carbons at 2.9 Å, `score_coordinates(...)` → `1.871479`
  (pair replay works with explicit pairs/types).
- Weight VJP: `vjp(score_coordinates, coords, types, pairs=..., wrt="weights")`
  → pullback over the 7 SF_VINA weights returns term-sum gradients
  `[0.03916, 0.02231, 0.81, 1.0, 0.0, 0.0, -0.0]` (weight path works; term
  sums exist internally in `core.py::_linearisation` but are not public API).

## Totals

| Status | Count |
|---|---|
| can_implement | 2 |
| partially_supported | 12 |
| cannot_implement | 8 |
| uncertain | 0 |

## Task matrix

| task_id | bench | scientific output | required derivative/API path | package evidence | capability status | confidence | paper workflow status | blocker |
|---|---|---|---|---|---|---|---|---|
| `vina.ligand_pose_coordinate_gradient_refinement` | public | Refined pose via gradient of Vina pair interaction score w.r.t. ligand Cartesian coordinates | coords JVP/VJP of the pair score, composed into an external line search/quasi-Newton loop | `core.py#L352` score_coordinates; `core.py#L410` jvp; `support_status.md#L5-L7`; minimal probe above | can_implement | high | workflow_partially_checked | — |
| `vina.end_to_end_pose_parameter_optimization` | public | Joint 6+k pose/torsion optimization against a Vina-like docking objective | pose→coords transform (external) composed with coords VJP; full-score parity first | coords/weights path works; but only the restricted 6-pair replay, no native full-score recombination | partially_supported | high | not_assessed | B1 full SF_VINA recombination absent |
| `vina.atomic_grid_score_pose_optimization` | public | Grid-interpolated atom-type score derivatives over rigid-body/torsion vars | maps/interpolation layer; not pair-based | no maps API in package exports; `support_status.md#L10` map methods not suitable/deferred | cannot_implement | high | not_assessed | B2 affine maps/grid layer missing |
| `vina.native_bfgs_pose_gradient_local_search` | public | Standard Vina local refinement with score gradients, optionally flexible receptor residues | gradient supply into Vina's internal BFGS; flexible-receptor DOF | coords grads exist for a pair model; `support_status.md#L12` optimize/deferred; flexible receptor absent | partially_supported | medium | not_assessed | B3 native search not driven by sidecar; B1 |
| `vina.gradient_consistent_line_search_sampling` | public | Improve sampling by changing the line search inside Vina's local refinement | instrumentation of upstream search loop | package has no search loop at all; upstream does not expose one as API | cannot_implement | medium | not_assessed | B3 search-machinery instrumentation absent |
| `vina.flexible_receptor_sidechain_gradient_docking` | public | Ligand + receptor side-chain torsion sensitivities through selective flexibility | flex-residue DOF extension + score grads over complex coords | pair replay could carry complex pairs, but no flex-residue front-end, and full score missing | partially_supported | medium | not_assessed | B3; B1 |
| `vina.macrocycle_torsion_gradient_sampling` | public | Gradient through ring torsions under closure constraints | ring-closure constrained DOF layer | not present; `upstream` does macrocycle sampling natively but sidecar has no interface | cannot_implement | medium | not_assessed | B4 ring-closure constraints absent |
| `vina.ad4_energy_gradient_alternative_scoring` | public | AD4 force-field score (electrostatics, desolvation, H-bond) derivatives | AD4 term rules | AD4 absent: only SF_VINA 6-term replay (`core.py#L20` weights default) | cannot_implement | high | not_assessed | B1-alt AD4/Vinardo scoring families absent |
| `vina.score_weight_affinity_calibration` | public | Fit linear weights over per-term energies to predict ΔG | per-term decomposition + weight VJP into regression | weight VJP implemented and tested (`support_status.md#L7`; probe); per-term output not public | partially_supported | high | workflow_partially_checked | B5 per-term decomposition not exposed |
| `vina.pose_vs_affinity_weight_optimization` | public | Weight tuning for pose prediction vs affinity ranking | weight VJP + full score over pose sets | weight path ok; full score + pose machinery external | partially_supported | medium | not_assessed | B1; B5 |
| `vina.target_family_robust_score_tuning` | public | Target-family-wide robust weight optimization | weight VJP at scale over family | weight path ok; full-term coverage for family fairness missing | partially_supported | medium | not_assessed | B1; B5 |
| `vina.multiobjective_score_training` | public | Balanced scoring (affinity, ranking, screening) training | term-vector features + weight/composite grads | features not exposed (probe: no term API) | partially_supported | medium | not_assessed | B5; B1 |
| `vina.vina_term_free_energy_surrogate` | public | Differentiable surrogate built from interpretable Vina terms | per-term vectors as regression features | absent (internal `term_sums` in `core.py#L308-L346` only) | partially_supported | high | not_assessed | B5 |
| `vina.receptor_ensemble_sensitivity` | public | Score sensitivity to receptor conformation ensemble | receptor-conformation dependent score path | pair replay has no receptor-field/maps model; grid term missing | partially_supported | medium | not_assessed | B2; B1 |
| `vina.cross_docking_receptor_sensitivity` | public | Derivatives of cross-docking score vs receptor conformation | receptor coords sensitivity via maps or grid scoring | ligand-side coords grads exist; receptor-side absent | partially_supported | medium | not_assessed | B1; B2 |
| `vina.ligand_substituent_score_sensitivity` | public | Local score sensitivity to continuous substituent geometry | ligand coords grads (geometry changes) | covered by the pair replay coords path (probe path) | can_implement | medium | workflow_partially_checked | — (atom-type changes are discrete; out of scope) |
| `vina.protein_mutation_affinity_sensitivity` | public | ΔΔG by differentiating vs residue identity | derivative w.r.t. discrete residue identity | no discrete-identity derivative; upstream has none either | cannot_implement | high | not_assessed | B6 discrete identity derivative |
| `vina.water_placement_gradient_docking` | public | Pose + ordered-water placement joint optimization | explicit-water scoring terms | no hydrated scoring in SF_VINA; upstream lacks it | cannot_implement | high | not_assessed | B7 water terms outside upstream surface |
| `vina.coordinated_water_interaction_sensitivity` | public | Water-mediated contact sensitivity | coordinated-water terms | absent upstream and in sidecar | cannot_implement | high | not_assessed | B7 |
| `vina.metal_coordination_gradient_docking` | public | Metal-coordination-geometry-sensitive scoring | metal atom types / coordination terms | Vina atom typing has no metals; paper uses its own scorer | cannot_implement | high | not_assessed | B7 |
| `vina.rna_ligand_pose_gradient` | public | RNA pocket pose gradients | RNA-specific scoring terms | paper modifies AD4 itself (ledger: other_implementation) | cannot_implement | high | not_assessed | B7 |
| `vina.peptide_flexible_pose_gradient` | public | Flexible peptide docking, many backbone/sidechain DOF | backbone DOF front-end + score grads | coords grads exist; backbone closure/torsion layer and full score absent | partially_supported | medium | not_assessed | B4; B1 |

## Blockers grouped by type

| Blocker | Type | Tasks | Detail |
|---|---|---|---|
| B1 full SF_VINA recombination | missing API (scoring surface) | 2, 4, 6, 8 (AD4 variant), 10, 11, 12, 13, 14, 15, 22 | `score_coordinates` replays only the 6 explicit pair terms assembled by the caller; the compiled scorer's recombination (maps-based interactions, AD4/Vinardo families) has no public differentiable path |
| B2 affine maps/grid layer | missing API | 3, 14, 15 | no compute/maps/interpolate API; grid-scoring tasks have nothing to differentiate |
| B3 upstream search machinery | unsupported control flow | 4, 5, 6, 22 | Vina's internal BFGS/line search/selective-flexibility loop is compiled; the sidecar does not drive or instrument it |
| B4 ring/backbone closure constraints | unsupported control flow | 7, 22 | no macrocycle/peptide closure-constrained DOF layer |
| B5 per-term decomposition | missing API | 9, 12, 13 | term sums exist internally (weights VJP) but are not exposed as callable outputs |
| B6 discrete identity derivative | unsupported control flow | 17 | derivative w.r.t. residue identity is discrete; not an AD contract |
| B7 beyond-upstream scoring surfaces | external solver/scientific scope | 18, 19, 20, 21 | water/metal/RNA terms do not exist in SF_VINA; papers rely on own/modified scorers |

## Issue candidates — not submitted here

(this review does not submit issues; grouping and submission follow the
ad-task-issue-publisher skill)

1. **Parity — differentiable full Vina scoring surface** (SF_VINA recombination, per-term decomposition, AD4/Vinardo families) — tasks 2, 4, 6, 8, 9, 10, 11, 12, 13, 15, 22. Upstream `Vina(sf_name=...)`/`set_weights`/`score` forward exists; the sidecar exposes only the 6-term pair replay. Minimal reproduction: `import vina_ad` has no AD4/maps/term API; `RuleNotFound` on `vina.Vina.score`.
2. **Parity — differentiable affinity maps / grid interpolation** — tasks 3, 14, 15. Upstream `vina.Vina.compute_vina_maps`/`load_maps` forward; sidecar has no map layer to differentiate.
3. **Extension draft — pose DOF front-end and search composition** (6+k DOF, flex sidechains, macrocycle/peptide closure, line-search instrumentation) — tasks 4 (partial), 5, 6 (partial), 7, 22 (partial).
4. **Extension draft — hydration / metal / RNA scoring extensions** — tasks 18, 19, 20, 21.
5. **Extension draft — mutation ΔΔG and analogue descriptor front-end** — task 17 (and the discrete part of 16).

## Assumptions, gaps, unreviewed material

- The **private** benchmark corpus (`flyingwagner/ad-software-private-benchmark`)
  was not accessed and not reviewed; no fallback copy was used.
- 15 of 22 papers have no PDF in this public checkout (manifest marks their
  downloads unavailable); their contracts rely on the ledger record only.
  This is recorded per row via `paper_workflow_status=not_assessed`.
- No full paper reproduction was attempted. Probe evidence covers only the
  package API mechanics listed above.
- The compiled `vina_wrapper` was not exercised from this Mac; the sidecar's
  own requirement that the restricted kernel drop it (see
  `vina_ad/requirements.md`) was taken as documented behavior.
- Package self-tests and machine receipts from the batch build
  (`955ec33c` round-3) were used for test-level evidence; they were not
  re-run inside this review.