# EZprodesign — structural epitope accessibility explorer

**Research software preview | Engine v1.5 | GitHub release not yet approved**

EZprodesign is an interpretable, local-first workbench for exploring **candidate
accessible surface patches** on protein structures. It supports uploaded
PDB/mmCIF files and FASTA input, offers geometric exposure and approach checks,
glycosylation-motif cautions, a 3Dmol.js visualization, and a JSON/CLI engine.
This repo is a deliberately minimal **generic software extraction** from a
larger private lab research workspace.

**Important:** This is *not* an antibody generator, docking predictor, binding
affinity predictor, experimental epitope validation or a calibrated score.
No predicted pLDDT is invented from crystallographic B-factors.
A membrane orientation not independently supplied and verified remains unknown.
Results are hypotheses for human expert review. See [scientific limitations](docs/SCIENTIFIC_SCOPE.md).

## Quick start (Python 3.10+)

From this repository's root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[web,dev]'
python -m v15_engine examples/1UBQ.pdb --environment soluble > /tmp/ezprodesign_demo.json
python -m v15_web
# open http://127.0.0.1:18785
```

Python core library:

```python
from v15_engine import analyze_bytes
raw = open("examples/1UBQ.pdb", "rb").read()
report = analyze_bytes("1UBQ.pdb", raw, environment="soluble", binder_type="VHH")
print(report["status"], len(report["patches"]))
```

Browser: Open localhost:18785, click the included **1UBQ** structure or
upload your own PDB/mmCIF/FASTA. The service binds to 127.0.0.1 only.
Upload bytes are processed in memory and are not saved by default; results
are sent back to your browser. Do not deploy this development server to a
public interface.

## Demo and verifiability

- `examples/1UBQ.pdb`: the public-domain RCSB PDB 1UBQ structure with a
  SHA-256 provenance receipt and attribution in `examples/1UBQ.receipt.json`.
- `tests/test_public_release.py`: synthetic and public-structure tests;
  never needs internal structures or laboratory data.
- `python -m unittest discover -s tests -v` or `pytest -q tests`
- `node --check v15_web/static/app.js` (Node optional)
- Scientific assumptions and limitations: `docs/SCIENTIFIC_SCOPE.md`
- Security considerations: `SECURITY.md`

## Code layout

| Path | Role |
|---|---|
| `v15_engine/` | Data-independent Python structure parsing/scoring engine |
| `v15_web/` | Local Flask browser workbench |
| `v15_web/vendor/` | Unmodified 3Dmol.js redistribution and upstream license notices |
| `examples/` | Public RCSB test structure and provenance receipt |
| `tests/` | Public safe contracts and API smoke tests |

## Status and contribution

This is an initial extracted package, not evidence of ecosystem adoption.
No verified binding specificity, precision, recall, or utility beyond
documented tests is claimed. We invite real reproducibility feedback and
independent benchmarking once the repository is publicly released.

**Licensing gate:** This private staging version deliberately does not
assert that institution-developed code may be publicly relicensed.
The rights-holder must clear IP/third-party terms and add an approved
open-source license **before changing the repository to Public**.
See `LICENSING.md`.
