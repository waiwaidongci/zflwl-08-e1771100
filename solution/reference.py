#!/usr/bin/env python3
"""Reference solver for Doppler spectral-line association."""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path
from typing import Any


SPEED_OF_LIGHT_KM_S = 299792.458
MATCH_TOLERANCE_NM = 0.001
TOP_KEYS = {"schema", "catalog", "exposures"}
LINE_ID = re.compile(r"LINE-[0-9]{4}\Z")
FRAME_ID = re.compile(r"FRAME-[0-9]{3}\Z")


class Refused(ValueError):
    pass


def _number(value: Any, label: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise Refused(f"{label} is not numeric")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0.0):
        raise Refused(f"{label} is outside its domain")
    return result


def parse_payload(path: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise Refused("input JSON is unreadable") from exc
    if not isinstance(payload, dict) or set(payload) != TOP_KEYS:
        raise Refused("top-level keys are invalid")
    if payload.get("schema") != "astro-lines-v1":
        raise Refused("schema is unsupported")

    raw_catalog = payload.get("catalog")
    if not isinstance(raw_catalog, list) or not 32 <= len(raw_catalog) <= 80:
        raise Refused("catalog size is outside 32 to 80")
    catalog: list[dict[str, Any]] = []
    line_ids: set[str] = set()
    wavelengths: set[float] = set()
    for item in raw_catalog:
        if not isinstance(item, dict) or set(item) != {"id", "wavelength_nm"}:
            raise Refused("catalog entry keys are invalid")
        line_id = item.get("id")
        if not isinstance(line_id, str) or LINE_ID.fullmatch(line_id) is None:
            raise Refused("catalog line id is invalid")
        if line_id in line_ids:
            raise Refused("catalog line ids must be unique")
        line_ids.add(line_id)
        wavelength = _number(
            item.get("wavelength_nm"), "catalog wavelength", positive=True
        )
        if wavelength in wavelengths:
            raise Refused("catalog wavelengths must be unique")
        wavelengths.add(wavelength)
        catalog.append({"id": line_id, "wavelength_nm": wavelength})

    raw_exposures = payload.get("exposures")
    if not isinstance(raw_exposures, list) or not 6 <= len(raw_exposures) <= 12:
        raise Refused("exposure count is outside 6 to 12")
    exposures: list[dict[str, Any]] = []
    frame_ids: set[str] = set()
    for item in raw_exposures:
        if not isinstance(item, dict) or set(item) != {
            "frame_id",
            "pixel_positions",
            "calibration",
        }:
            raise Refused("exposure object keys are invalid")
        frame_id = item.get("frame_id")
        if not isinstance(frame_id, str) or FRAME_ID.fullmatch(frame_id) is None:
            raise Refused("frame id is invalid")
        if frame_id in frame_ids:
            raise Refused("frame ids must be unique")
        frame_ids.add(frame_id)
        calibration = item.get("calibration")
        if not isinstance(calibration, dict) or set(calibration) != {
            "pixel_zero",
            "nm_per_pixel",
        }:
            raise Refused("calibration keys are invalid")
        pixel_zero = _number(calibration.get("pixel_zero"), "pixel_zero")
        dispersion = _number(
            calibration.get("nm_per_pixel"), "nm_per_pixel", positive=True
        )
        raw_positions = item.get("pixel_positions")
        if not isinstance(raw_positions, list) or not 16 <= len(raw_positions) <= 28:
            raise Refused("pixel_positions must contain 16 to 28 peaks")
        positions = [
            _number(value, "pixel position", positive=True) for value in raw_positions
        ]
        if any(positions[index] >= positions[index + 1] for index in range(len(positions) - 1)):
            raise Refused("pixel positions must be strictly increasing")
        exposures.append(
            {
                "frame_id": frame_id,
                "pixel_positions": positions,
                "calibration": {
                    "pixel_zero": pixel_zero,
                    "nm_per_pixel": dispersion,
                },
            }
        )
    return catalog, exposures


def match_at_velocity(
    observed: list[float],
    catalog: list[dict[str, Any]],
    velocity_kms: float,
) -> tuple[int, tuple[int, ...], dict[str, int]]:
    factor = 1.0 + velocity_kms / SPEED_OF_LIGHT_KM_S
    expected = sorted(
        (line["wavelength_nm"] * factor, index)
        for index, line in enumerate(catalog)
    )
    used: set[int] = set()
    observed_matches: list[int] = []
    by_line: dict[str, int] = {}
    for observed_index, wavelength in enumerate(observed):
        # 只保留容差区间内的未用谱线，候选集合按最近波长收敛。
        candidates = [
            (abs(wavelength - shifted), line_index)
            for shifted, line_index in expected
            if abs(wavelength - shifted) <= MATCH_TOLERANCE_NM
            and line_index not in used
        ]
        if not candidates:
            observed_matches.append(-1)
            continue
        _, line_index = min(candidates)
        used.add(line_index)
        observed_matches.append(line_index)
        by_line[catalog[line_index]["id"]] = observed_index
    support = sum(index >= 0 for index in observed_matches)
    return support, tuple(observed_matches), dict(sorted(by_line.items()))


def recover_exposure(
    exposure: dict[str, Any], catalog: list[dict[str, Any]]
) -> dict[str, Any]:
    calibration = exposure["calibration"]
    observed = [
        calibration["pixel_zero"] + pixel * calibration["nm_per_pixel"]
        for pixel in exposure["pixel_positions"]
    ]
    hypotheses: dict[float, float] = {}
    for wavelength in observed:
        for line in catalog:
            velocity = (
                wavelength / line["wavelength_nm"] - 1.0
            ) * SPEED_OF_LIGHT_KM_S
            if -300.0 - 1e-9 <= velocity <= 300.0 + 1e-9:
                velocity = min(300.0, max(-300.0, velocity))
                hypotheses[round(velocity, 9)] = velocity

    if not hypotheses:
        raise Refused("no velocity candidate lies inside the public range")
    scored: list[tuple[int, float, tuple[int, ...], dict[str, int]]] = []
    for velocity in sorted(hypotheses.values()):
        support, mapping, by_line = match_at_velocity(observed, catalog, velocity)
        scored.append((support, velocity, mapping, by_line))
    best_support = max(item[0] for item in scored)
    best = [item for item in scored if item[0] == best_support]
    distinct_mappings = {item[2] for item in best}
    if best_support != len(observed):
        raise Refused("no velocity covers every observed peak")
    if len(distinct_mappings) > 1:
        return {"candidates": len(distinct_mappings), "status": "ambiguous"}
    _, velocity, mapping, by_line = best[0]
    matched_pairs = [
        (catalog[line_index]["wavelength_nm"], observed[observed_index])
        for observed_index, line_index in enumerate(mapping)
        if line_index >= 0
    ]
    factor = sum(rest * shifted for rest, shifted in matched_pairs) / sum(
        rest * rest for rest, _ in matched_pairs
    )
    velocity = (factor - 1.0) * SPEED_OF_LIGHT_KM_S
    return {
        "status": "ok",
        "velocity_kms": round(velocity, 6),
        "line_matches": by_line,
    }


def main() -> int:
    if len(sys.argv) != 3:
        print("REFUSED: expected input and output paths", file=sys.stderr)
        return 2
    try:
        catalog, exposures = parse_payload(Path(sys.argv[1]))
        results = [
            recover_exposure(exposure, catalog) for exposure in exposures
        ]
    except Refused as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2
    ambiguous = {
        exposure["frame_id"]: result["candidates"]
        for exposure, result in zip(exposures, results)
        if result["status"] == "ambiguous"
    }
    if ambiguous:
        output: dict[str, Any] = {
            "status": "ambiguous",
            "candidates": dict(sorted(ambiguous.items())),
        }
    else:
        output = {
            "status": "ok",
            "frames": {
                exposure["frame_id"]: {
                    "velocity_kms": result["velocity_kms"],
                    "line_matches": result["line_matches"],
                }
                for exposure, result in zip(exposures, results)
            },
        }
    # 结果按规范编码写出紧凑 JSON，保证归档字节稳定。
    Path(sys.argv[2]).write_text(
        json.dumps(output, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
