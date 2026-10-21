# Curated public-release safety audit — 2026-10-21

**Scope:** isolated v1.5 generic engine and local web UI, no prior git
history. This snapshot was created by copying selected code and documentation
from a larger private research workspace; it was not a wholesale upload.

## Included

- Structural engine: core `v15_engine/` (parsing, geometric gates, scoring), excluding private-fixture validation runner.
- Local Flask UI: `v15_web/`, bounded in-memory upload/JSON API, 3D viewer.
- Public-domain RCSB 1UBQ structure and provenance SHA256 receipt.
- 3Dmol.js vendor distribution and BSD-3-Clause license/attribution notice.
- Public-only tests, setup metadata, documentation and GitHub CI.

## Deliberately excluded

- All raw and processed lab structural files, coordinates, cryo-EM/tomo data,
  live annotation databases, private model predictions and derived visualizations.
- All original `data/`, `output/`, `releases/`, `reviews/`, `runs/`
  folders and original Git history.
- Research-specific SLC34A2/NK2R/DRG/Cage test fixtures or outputs, model designs,
  internal reviewers' transcripts and development logs.
- HPC/Slurm deployment recipes, machine aliases, personal paths, local authentication,
  cloud credentials and hidden assistant configuration directories.
- Untracked v16–v18 experimental work and unpublished benchmarks.
- No secrets or access tokens were copied knowingly.

## Checks

- Source sizes and file inventory logged locally before publishing.
- Python syntax compilation, Node JS syntax, and public test suite.
- Scanned staged text for secrets, email/identity strings, private absolute
  paths, internal targets and embedded external URLs. Audit must be rerun
  after changes.
- Public test cases use only 1UBQ and synthetic sequences; no test-time download.

## Blocking publication issues

1. **Intellectual property / copyright ownership**: obtain confirmation from
   the software rights holder(s), collaborators and institutional requirements.
2. **License**: add an approved OSI license only after item 1. Until then this
   repository is a private release candidate, NOT yet open source.
3. **Usability**: verify clean install/CI on an independent machine, documentation,
   and the public repo's version/release workflow.
4. **Adoption**: new project has no verified external stars, downloads or community
   maintenance. Applications must not suggest otherwise.

Do not change GitHub visibility or apply as a qualifying public repository until
items 1 and 2 are resolved.
