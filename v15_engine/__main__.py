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
    args = parser.parse_args()
    try:
        with args.input.open("rb") as stream:
            raw = stream.read(MAX_BYTES+1)
        result = analyze_bytes(args.input.name, raw, args.chain, args.environment,
                               args.binder, {"sphere_points":args.sphere_points})
    except (OSError, UploadValidationError) as exc:
        print(str(exc), file=sys.stderr)
        return 2
    print(json.dumps(result, ensure_ascii=False, allow_nan=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
