<!-- generated from ISSUE_DRAFT_MUTATION_DELTA_DELTA_G.md; edit the source draft and regenerate -->
<!-- classification: scientific_extension; submit_enabled: false -->
# Mutation ΔΔG and discrete-descriptor derivative front-end

Draft only — not submitted to GitLab.

## Summary

`vina.protein_mutation_affinity_sensitivity` asks for derivatives of a
Vina-compatible complex score with respect to residue identity / local
side-chain parameters. Residue identity is a discrete input; the sidecar's
active inputs are `coordinates` and `weights` only, and no upstream Vina API
provides a mutation-derivative contract. The adjacent
`vina.ligand_substituent_score_sensitivity` task is supported for continuous
geometry via coordinate gradients, but its discrete part (atom-type changes
along an analogue series) has no derivative either. This is therefore recorded
as a non-submitted extension draft for a possible ΔΔG / analogue front-end.

**Classification:** `scientific_extension`

## Versions and scope

- Upstream software: AutoDock Vina, pinned commit `3c65c0b3e6c2c1d183f6a175ecb65e3c5ba91645`
- AD package: `https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/autodock-vina_v3_adms_xj` @ `955ec33c019cdcdbc4c0aeaa0c5394accb9de95c` (`vina_ad` 0.1.0)
- Capability review: `.package-reviews/vina/bench/public/0.1.0/TASK_CAPABILITY_REVIEW.md`

## Related tasks and papers

| Task | What it needs | Paper |
|---|---|---|
| [`vina.protein_mutation_affinity_sensitivity`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L304) | derivative w.r.t. residue identity/local side-chain params | [bioRxiv 2025](https://doi.org/10.1101/2025.07.19.665665) |
| [`vina.ligand_substituent_score_sensitivity`](https://git.gewu-lab.ai/matrixlab/sci-ad-shadows/ad-software-public-benchmark/-/blob/main/vina/vina-task-ledger.json#L287) | geometry derivatives (supported); atom-type changes (discrete) | [J Biomol Struct Dyn 2020](https://doi.org/10.1080/07391102.2020.1792987) |

## Candidate design (if ever enabled)

- ΔΔG as paired-score differences over explicitly enumerated variants, with
  coordinate/weight gradients from the existing path where inputs stay
  continuous; discrete side-chain choices enumerated, not differentiated.
- Continuous atom-level descriptors (radii, well depths, electronegativity
  proxies) as optional active inputs beside `coordinates`/`weights`.

## Non-goals

- Differentiating through residue identity as if it were a continuous
  variable; mutation scanning tooling; force-field reparameterization.
