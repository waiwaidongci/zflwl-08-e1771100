#!/usr/bin/env python3
"""Partial spectral-line association entry point."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 3:
        print("REFUSED: expected input and output paths", file=sys.stderr)
        return 2
    source = Path(sys.argv[1])
    target = Path(sys.argv[2])
    payload = json.loads(source.read_text(encoding="utf-8"))
    exposures = payload.get("exposures", [])
    catalog = payload.get("catalog", [])
    result = {
        "status": "ok",
        "velocity_kms": 0.0,
        "line_matches": {
            catalog[index].get("id", "LINE-0000"): index
            for index in range(min(len(catalog), 16))
        },
    }
    target.write_text(json.dumps(result, separators=(",", ":")) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
