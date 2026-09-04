<!-- generated from ISSUE_DRAFT_HYDRATION_METAL_RNA_EXTENSIONS.md; edit the source draft and regenerate -->
<!-- classification: scientific_extension; submit_enabled: false -->
# Hydration / metal / RNA scoring extensions

Draft only — not submitted to GitLab.

## Summary

Four ledger tasks require scoring surfaces that the pinned upstream SF_VINA
force field does not provide: explicit ordered-water placement terms,
metal-coordination geometry terms, and RNA-specific pocket scoring. The papers
build or modify their own scorers for these. Because the upstream software has
no corresponding forward API, these cannot be upstream-parity gaps and are
kept as non-submitted extension drafts.

**Classification:** `scientific_extension`

## Versions and scope

- Upstream software: AutoDock Vina, pinned commit `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645`
- AD package: `https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj` @ `955ec33c019cdcdbc4c0aeaa0c5394accb9de95c` (`vina_ad` 0.1.0)
- Capability review: `.package-reviews/vina/bench/public/0.1.0/TASK_CAPABILITY_REVIEW.md`

## Related tasks and papers

| Task | What it needs | Paper |
|---|---|---|
| [`vina.water_placement_gradient_docking`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L321) | ordered-water DOF + hydrated score | [J Med Chem 2008](https://doi.org/10.1021/jm8006239) |
| [`vina.coordinated_water_interaction_sensitivity`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L338) | water-mediated contact contributions | [PLoS Comput Biol 2020](https://doi.org/10.1371/journal.pcbi.1008103) |
| [`vina.metal_coordination_gradient_docking`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L355) | metal-donor geometry terms | [JCIM 2015](https://doi.org/10.1021/ci500647f) |
| [`vina.rna_ligand_pose_gradient`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L372) | RNA pocket scoring/topology | [Methods 2021](https://doi.org/10.1016/j.ymeth.2021.01.009) |

## Notes

- Guard against overreach: these surfaces would complement the existing
  SF_VINA coverage, not replace upstream scoring.
- Any future work should pin the paper-specific scoring formulas and
  validation protocols in its own spec, mirroring the review's philosophy.
