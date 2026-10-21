"""CLI writes only to stdout unless the caller explicitly redirects it."""
import argparse
import json
from pathlib import Path
import sys

from .core import MAX_BYTES, UploadValidationError, analyze_bytes


def main():
    parser = argparse.ArgumentParser(description="CPU protein surface exploration; JSON to stdout")
    parser.add_argument("input", type=Path)
    parser.add_argument("--chain", default="auto")
    parser.add_argument("--environment", choices=("detergent", "soluble", "bilayer"), default="detergent")
    parser.add_argument("--binder", choices=("VHH", "Fab"), default="VHH")
    parser.add_argument("--sphere-points", type=int, default=64)
    parser.add_argument("--summary", action="store_true", help="print a concise human-readable research summary instead of JSON")
    args = parser.parse_args()
    try:
        with args.input.open("rb") as stream:
            raw = stream.read(MAX_BYTES+1)
        result = analyze_bytes(args.input.name, raw, args.chain, args.environment,
                               args.binder, {"sphere_points":args.sphere_points})
    except (OSError, UploadValidationError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    if args.summary:
        print("EZprodesign | exploratory geometry; NOT validated binding or affinity")
        print(f"Input: {result['source']['filename']}")
        print(f"Status: {result['status']} | Environment: {result['environment']['mode']} | Binder proxy: {args.binder}")
        print(f"Residues: {result['target']['residue_count']} | Candidate patches: {len(result['patches'])}")
        for patch in result['patches']:
            hot = ",".join(str(n) for n in patch['hotspots'])
            print(f"{patch['id']}: {patch['status']} | priority {patch['score']:.3f} | hotspot positions {hot}")
        if result['status'] == 'needs_structure':
            print("No 3D geometry: supply a structure to assess candidate patches.")
        if result.get('warnings'):
            print("Warnings: " + ", ".join(result['warnings']))
    else:
        print(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
