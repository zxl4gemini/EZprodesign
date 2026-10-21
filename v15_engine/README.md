# Structure analysis engine (v1.5)

The Python module `v15_engine` implements CPU-only, deterministic and
local structural surface exploration for uploaded PDB/mmCIF or FASTA:

```python
from v15_engine.core import analyze_bytes, UploadValidationError
report = analyze_bytes("1UBQ.pdb", open("examples/1UBQ.pdb","rb").read(),
                       target_chain="auto", environment="soluble",
                       binder_type="VHH", settings=None)
```

The CLI is `python -m v15_engine examples/1UBQ.pdb --environment soluble`.
No network, disk writes or deep learning inference occur in the engine.
The web layer enforces a 16 MiB request bound and runs at localhost only.

Input formats: PDB, mmCIF, FASTA; compressed input is rejected. Only
the first structural model is scored; non-target heavy atoms may act as
occluders. Structures use observed residues only. FASTA returns
`needs_structure` and no geometry, and pLDDT is **never** inferred from
arbitrary B factors. Valid environments are `soluble`, `detergent`,
`bilayer`; binders are `VHH` and `Fab`. Without a independently verified
membrane orientation, bilayer access remains uncertified. Detergent mode
uses a hydropathy proxy, **not** an atomic micelle model.

These are heuristic targetability and geometric plausibility checks,
not verified epitopes, affinity predictions, or binder designs.
See ../docs/SCIENTIFIC_SCOPE.md for limitations and ../README.md for
reproduction and test commands.
