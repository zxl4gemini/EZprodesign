# Contributing to EZprodesign

We welcome reproducible scientific bug reports and focused code changes.
Share only **public** test structures (or synthetic FASTA), Python version,
reproduction command, expected and observed output. Do not submit private
research structures, unpublished binders, patient data, passwords or HPC logs.

Before sending a pull request, run:

```bash
python -m pip install -e ".[web,dev]"
python -m unittest discover -s tests -v
python -m pytest -q tests
```

A scientific change must state its limitations and add public-data tests;
software execution is not experimental validation. Good contributions
include cross-platform installation, accessibility, translation, performance,
negative controls and deterministic scientific regression checks.
