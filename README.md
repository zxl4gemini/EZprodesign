# EZprodesign

**Interpretable protein surface accessibility and exploratory epitope review**

[![MIT License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Tests](https://github.com/zxl4gemini/EZprodesign/actions/workflows/tests.yml/badge.svg)](https://github.com/zxl4gemini/EZprodesign/actions/workflows/tests.yml)

EZprodesign is a **local-first open-source structural biology workbench**.
It combines CPU-only protein structure parsing, solvent-exposure and steric
approach heuristics, candidate surface patch ranking, interactive 3Dmol.js
visualization, scientific uncertainty warnings, and exportable review data.
Researchers can explore soluble or membrane-associated proteins on their
own computer without a cloud account or GPU.

**Scientific nonclaims:** This is *not* an antibody sequence generator,
experimentally validated epitope predictor, docking simulator, or affinity
predictor. Its numeric values are heuristic priorities, **not binding
probabilities**. Membrane orientation, glycan occupancy, and detergent-micelle
boundaries are not independently verified. Read [scientific scope](docs/SCIENTIFIC_SCOPE.md).

## Browser quick start — Python 3.10+

```bash
git clone https://github.com/zxl4gemini/EZprodesign.git
cd EZprodesign
python -m venv .venv
source .venv/bin/activate           # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -e ".[web]"
python -m v15_web                  # optional: --port 18786
```

Open **http://127.0.0.1:18785** in Chrome or Edge and click the bundled
**1UBQ ubiquitin** example, or upload your own PDB, mmCIF or FASTA file.
The UI is currently Chinese-language; English documentation is available.
Uploaded structure data are analyzed in local memory and are not saved
server-side by default. **Never expose the local Flask server to the Internet.**

## Command-line quick start

```bash
python -m v15_engine examples/1UBQ.pdb --environment soluble --summary
python -m v15_engine examples/1UBQ.pdb --environment soluble > result.json
python -m v15_engine examples/1UBQ.pdb --environment bilayer --binder Fab --summary
```

Accepted extensions: `.pdb`, `.cif`, `.mmcif`, `.fa`, `.fasta`, `.faa`
(maximum 15 MiB). Environments: `soluble`, `detergent` (default), `bilayer`.
Binder-envelope proxies: `VHH` (default) or `Fab`. A FASTA-only analysis
returns `needs_structure` and **no artificial 3D hotspots**.

Python API:

```python
from pathlib import Path
from v15_engine import analyze_bytes

source = Path("examples/1UBQ.pdb")
result = analyze_bytes(source.name, source.read_bytes(), environment="soluble")
print(result["status"], len(result["patches"]))
```

## Reproducibility and tests

```bash
python -m pip install -e ".[web,dev]"
python -m unittest discover -s tests -v
python -m pytest -q tests
node --check v15_web/static/app.js     # Node.js optional
```

The supplied `examples/1UBQ.pdb` is public RCSB/wwPDB structure data
(CC0-1.0), with checksum and attribution in `examples/1UBQ.receipt.json`.
It verifies execution, **not experimental antibody binding**.
Continuous integration runs public-only regression tests across Python versions.

## Repository layout

| Path | Purpose |
|---|---|
| `v15_engine/` | Deterministic structural engine and CLI |
| `v15_web/` | Loopback Flask UI and 3Dmol viewer |
| `examples/` | Public 1UBQ demo and provenance |
| `tests/` | Public-only regression tests |
| `docs/SCIENTIFIC_SCOPE.md` | Scientific assumptions and limitations |
| `OSS_RELEASE_AUDIT.md` | Published-code scope and exclusions |

## Community and licenses

EZprodesign is a **new research software project**; it has no claimed
downstream adoption, experimentally validated binding accuracy or hit rate.
We welcome reproducible issues, public-data benchmarks, scientific controls,
translations and tested PRs. See [CONTRIBUTING.md](CONTRIBUTING.md) and
[SECURITY.md](SECURITY.md).

Generic project software is [MIT-licensed](LICENSE), © 2026 EZprodesign
contributors. Bundled 3Dmol.js retains its **BSD-3-Clause** license.
1UBQ public data are **wwPDB CC0-1.0**. See [LICENSING.md](LICENSING.md).
Unpublished laboratory experiments and private designs are excluded.
