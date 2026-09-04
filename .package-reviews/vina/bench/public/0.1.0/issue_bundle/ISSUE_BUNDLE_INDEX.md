# AutoDock Vina consolidated issue bundle

> Generated locally. Remote submission is a separate guarded action.

Package project: `matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj`
Package commit: `955ec33c019cdcdbc4c0aeaa0c5394accb9de95c`
Problem project/commit: `` / ``
Capability review: `.package-reviews/vina/bench/public/0.1.0/TASK_CAPABILITY_REVIEW.md`

## Submission queue

Only `upstream_parity_gap` groups with `submit_enabled: true` are eligible.

| Group | Classification | Enabled | Tasks | Rendered body | Submission |
|---|---|---:|---:|---|---|
| `full_vina_scoring_surface` | `upstream_parity_gap` | yes | 11 | [ISSUE_FULL_VINA_SCORING_SURFACE.md](ISSUE_FULL_VINA_SCORING_SURFACE.md) | not submitted |
| `affinity_maps_grid` | `upstream_parity_gap` | yes | 3 | [ISSUE_AFFINITY_MAPS_GRID.md](ISSUE_AFFINITY_MAPS_GRID.md) | not submitted |
| `pose_dof_search_frontend` | `scientific_extension` | no | 6 | [ISSUE_POSE_DOF_SEARCH_FRONTEND.md](ISSUE_POSE_DOF_SEARCH_FRONTEND.md) | not submitted |
| `hydration_metal_rna_extensions` | `scientific_extension` | no | 4 | [ISSUE_HYDRATION_METAL_RNA_EXTENSIONS.md](ISSUE_HYDRATION_METAL_RNA_EXTENSIONS.md) | not submitted |
| `mutation_delta_delta_g` | `scientific_extension` | no | 2 | [ISSUE_MUTATION_DELTA_DELTA_G.md](ISSUE_MUTATION_DELTA_DELTA_G.md) | not submitted |

Enabled parity issues: **2**

## Task mapping

### `full_vina_scoring_surface` (upstream_parity_gap)

- `vina.score_weight_affinity_calibration`
- `vina.pose_vs_affinity_weight_optimization`
- `vina.target_family_robust_score_tuning`
- `vina.multiobjective_score_training`
- `vina.vina_term_free_energy_surrogate`
- `vina.end_to_end_pose_parameter_optimization`
- `vina.native_bfgs_pose_gradient_local_search`
- `vina.flexible_receptor_sidechain_gradient_docking`
- `vina.ad4_energy_gradient_alternative_scoring`
- `vina.cross_docking_receptor_sensitivity`
- `vina.peptide_flexible_pose_gradient`

### `affinity_maps_grid` (upstream_parity_gap)

- `vina.atomic_grid_score_pose_optimization`
- `vina.receptor_ensemble_sensitivity`
- `vina.cross_docking_receptor_sensitivity`

### `pose_dof_search_frontend` (scientific_extension)

- `vina.end_to_end_pose_parameter_optimization`
- `vina.native_bfgs_pose_gradient_local_search`
- `vina.gradient_consistent_line_search_sampling`
- `vina.flexible_receptor_sidechain_gradient_docking`
- `vina.macrocycle_torsion_gradient_sampling`
- `vina.peptide_flexible_pose_gradient`

### `hydration_metal_rna_extensions` (scientific_extension)

- `vina.water_placement_gradient_docking`
- `vina.coordinated_water_interaction_sensitivity`
- `vina.metal_coordination_gradient_docking`
- `vina.rna_ligand_pose_gradient`

### `mutation_delta_delta_g` (scientific_extension)

- `vina.protein_mutation_affinity_sensitivity`
- `vina.ligand_substituent_score_sensitivity`
