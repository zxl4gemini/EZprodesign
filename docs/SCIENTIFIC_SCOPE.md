# Scientific scope and nonclaims

The core currently computes heuristic solvent/surface geometry, residue
patch grouping, steric approach directions and uncertainty labels on uploaded
coordinate models. It reports glycosylation **sequons**, not observed glycans;
a possible motif is not proven occupied. The local hydropathy context is **not**
a DDM/GDN micelle simulation and the tool does **not** verify bilayer or
membrane orientation. PDB B factors do not establish AlphaFold pLDDT.
FASTA inputs do not receive invented 3D coordinates; they return
`needs_structure` and no geometric patches.

`DESIGN_EXPLORATION` marks an unconstrained research hypothesis, not a
binding probability, docking hit or validated epitope. `HOLD_REVIEW` means a
geometry or context caveat requires human inspection. Protein flexibility,
oligomeric partners, detergents, glycans, membrane curvature and antibodies
may change accessibility. The method has **not** been demonstrated to predict
binding affinity, binder success, epitope specificity, FDR, precision or
recall. Outputs should not be used as medical or diagnostic evidence.

Reproducible baseline: publicly available 1UBQ PDB with SHA-256 receipt.
Independent validation on structurally diverse membrane proteins and negative
controls, verified conformations, ligand/oligomer contexts, real binders and
matched blinded experimental outcomes remains future work.
