from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path


ENTRY = os.environ.get("LINE_LINKER", "/app/link_lines.py")

SAFE_ENV = {
    "PATH": "/usr/bin:/bin",
    "HOME": "/tmp",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PYTHONDONTWRITEBYTECODE": "1",
    "PYTHONHASHSEED": "0",
}
C_KM_S = 299792.458


def _run(payload: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-S", "-P", "-B", ENTRY, str(payload), str(output)],
        text=True,
        capture_output=True,
        env=SAFE_ENV,
        timeout=30,
        check=False,
    )


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload), encoding="utf-8")


def _catalog() -> list[dict[str, object]]:
    return [
        {
            "id": f"LINE-{index:04d}",
            "wavelength_nm": 420.0 + 3.71 * index + 0.19 * (index % 7),
        }
        for index in range(64)
    ]


def _payload(
    exposure_count: int = 8,
    *,
    first_velocity: float = -182.4,
) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    catalog = _catalog()
    velocities = [first_velocity, 51.25, 244.0, 117.75, 12.5, 88.125, -64.5, 199.875]
    exposures: list[dict[str, object]] = []
    expected: dict[str, dict[str, object]] = {}
    pixel_zero = 380.0
    dispersion = 0.12
    for exposure_index in range(exposure_count):
        frame_id = f"FRAME-{exposure_index:03d}"
        velocity = velocities[exposure_index]
        selected = [(exposure_index * 11 + step * 17) % len(catalog) for step in range(24)]
        factor = 1.0 + velocity / C_KM_S
        observed_values = [
            catalog[index]["wavelength_nm"] * factor for index in selected
        ]
        ordered = sorted(
            (value, index) for value, index in zip(observed_values, selected)
        )
        positions = [
            round((value - pixel_zero) / dispersion, 9) for value, _ in ordered
        ]
        line_matches = {
            catalog[index]["id"]: observed_index
            for observed_index, (_, index) in enumerate(ordered)
        }
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
        expected[frame_id] = {
            "velocity_kms": round(velocity, 6),
            "line_matches": line_matches,
        }
    return {
        "schema": "astro-lines-v1",
        "catalog": catalog,
        "exposures": exposures,
    }, expected


def test_six_frame_spectral_subsets_are_identified(tmp_path: Path) -> None:
    """题面：
    「`exposures` 包含 6 到 12 帧。每帧只能有 `frame_id`、`pixel_positions`、`calibration`；」
    「每一帧只取一个速度值，观察到的每个峰位都要分到不同的目录谱线。选择能覆盖全部峰位」
    「若没有任何速度能在容差内覆盖全部峰位，不写结果文件，标准错误输出一行」
    「目录谱线号，值是该谱线对应的峰位下标，并且每个峰位恰好出现一次。」
    「输入结构、类型或值域错误时也按拒收处理，不产生输出文件。生成文件删去多余空格，」
    「键按字母升序，尾部保留单个换行。不得读取评测资料，也不能内置具体帧答案。」
    """

    payload, expected = _payload(6)
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert completed.returncode == 0
    assert completed.stdout == ""
    assert completed.stderr == ""
    assert result["status"] == "ok"
    assert result["frames"] == expected, "CASE-STDLIB-01 CASE-FRAME-01 CASE-FRAME-MATCH-01 CASE-MAP-01"


def test_eight_frame_catalog_association_is_complete(tmp_path: Path) -> None:
    """题面：
    「`frame_id` 采用 `FRAME-[0-9]{3}`，六帧以内不能出现重复编号。`calibration` 只能有 `pixel_zero` 与」
    """

    payload, expected = _payload(8)
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert completed.returncode == 0
    assert set(result["frames"]) == set(expected), "CASE-INPUT-02 CASE-CATALOG-02 CASE-OUTPUT-01"
    assert all(
        len(result["frames"][frame]["line_matches"]) == 24
        for frame in expected
    )


def test_doppler_velocity_is_recovered_from_many_lines(tmp_path: Path) -> None:
    """题面：
    「每帧只包含 `velocity_kms` 和 `line_matches`。速度保留六位小数；`line_matches` 的键是」
    """

    payload, expected = _payload(7)
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert completed.returncode == 0
    assert abs(result["frames"]["FRAME-000"]["velocity_kms"] + 182.4) < 1e-5
    assert abs(result["frames"]["FRAME-006"]["velocity_kms"] + 64.5) < 1e-5
    assert result["frames"]["FRAME-002"]["line_matches"] == expected["FRAME-002"]["line_matches"], "CASE-CATALOG-01 CASE-FRAME-02 CASE-CAL-02 CASE-OUTPUT-02"


def test_negative_velocity_limit_is_accepted(tmp_path: Path) -> None:
    """题面：
    「`pixel_positions` 是 16 到 28 个有限数，必须严格递增且不小于 0。实际观测波长等于」
    """

    payload, expected = _payload(6, first_velocity=-300.0)
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert completed.returncode == 0
    assert result["frames"]["FRAME-000"]["velocity_kms"] == -300.0
    assert result["frames"]["FRAME-000"]["line_matches"] == expected["FRAME-000"]["line_matches"], "CASE-SCHEMA-01 CASE-CAL-01 CASE-VELOCITY-01 CASE-FRAME-03"


def test_zero_pixel_position_is_accepted(tmp_path: Path) -> None:
    """题面：
    「`pixel_positions` 是 16 到 28 个有限数，必须严格递增且不小于 0。」
    """

    # CASE-PIXEL-02
    payload, expected = _payload(6)
    for exposure in payload["exposures"]:
        first = exposure["pixel_positions"][0]
        exposure["calibration"]["pixel_zero"] += (
            first * exposure["calibration"]["nm_per_pixel"]
        )
        exposure["pixel_positions"] = [
            position - first for position in exposure["pixel_positions"]
        ]
        assert abs(exposure["pixel_positions"][0]) <= 1e-12
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert completed.returncode == 0
    assert completed.stdout == ""
    assert completed.stderr == ""
    assert result["status"] == "ok"
    assert result["frames"] == expected


def test_repeated_line_patterns_report_frame_ambiguity(tmp_path: Path) -> None:
    """题面：
    「`schema` 必须为 `"astro-lines-v1"`。」
    """

    base = [430.0 + 8.3 * index for index in range(16)]
    shifted = [value * 1.0005 for value in base]
    catalog = [
        {"id": f"LINE-{index:04d}", "wavelength_nm": value}
        for index, value in enumerate(base + shifted)
    ]
    exposures = [
        {
            "frame_id": f"FRAME-{index:03d}",
            "pixel_positions": [round((value - 300.0) / 0.1, 6) for value in base],
            "calibration": {"pixel_zero": 300.0, "nm_per_pixel": 0.1},
        }
        for index in range(6)
    ]
    payload = {
        "schema": "astro-lines-v1",
        "catalog": catalog,
        "exposures": exposures,
    }
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    result = json.loads(output.read_text(encoding="utf-8"))
    assert completed.returncode == 0
    assert set(result) == {"status", "candidates"}
    assert result["status"] == "ambiguous"
    assert set(result["candidates"]) == {f"FRAME-{index:03d}" for index in range(6)}
    assert all(value >= 2 for value in result["candidates"].values()), "CASE-AMBIG-01 CASE-CANDIDATES-01"


def test_extra_catalog_key_is_refused(tmp_path: Path) -> None:
    """题面：
    「JSON 顶层只允许 `schema`、`catalog`、`exposures` 三个键，其他键一律视为坏件。」
    「速度候选只取 -300 到 300 km/s，两个端点也允许。」
    """

    payload, _ = _payload(6)
    payload["observer"] = "roof"
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    assert completed.returncode == 2
    assert completed.stdout == ""
    assert completed.stderr.startswith("REFUSED: ")
    assert not output.exists()


def test_success_output_is_compact_sorted_and_newline_terminated(tmp_path: Path) -> None:
    """题面：
    「`python3 /app/link_lines.py <input.json> <output.json>`，成功时不向标准输出或标准错误」
    """

    payload, _ = _payload(6)
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    raw = output.read_bytes()
    parsed = json.loads(raw)
    canonical = json.dumps(parsed, sort_keys=True, separators=(",", ":")).encode() + b"\n"
    assert completed.returncode == 0
    assert raw == canonical


def test_positive_velocity_limit_exceeded_is_refused(tmp_path: Path) -> None:
    """题面：
    「`REFUSED: <原因>`，退出码为 2。」
    """

    payload, _ = _payload(6, first_velocity=301.0)
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    assert completed.returncode == 2
    assert completed.stderr.startswith("REFUSED: ")
    assert not output.exists(), "CASE-PIXEL-01 CASE-PIXEL-02 CASE-SAFE-01"


def test_equal_pixel_positions_are_refused(tmp_path: Path) -> None:
    """题面：
    「写内容。运行只使用镜像内 Python 标准库。」
    """

    # CASE-PIXEL-03
    catalog = [
        {
            "id": f"LINE-{index:04d}",
            "wavelength_nm": 425.0
            if index == 0
            else 425.0006
            if index == 1
            else 428.0 + index
            if index < 16
            else 500.0 + index,
        }
        for index in range(32)
    ]
    positions = [425.0, 425.0] + [430.0 + index for index in range(14)]
    payload = {
        "schema": "astro-lines-v1",
        "catalog": catalog,
        "exposures": [
            {
                "frame_id": f"FRAME-{index:03d}",
                "pixel_positions": positions,
                "calibration": {"pixel_zero": 0.0, "nm_per_pixel": 1.0},
            }
            for index in range(6)
        ],
    }
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    report = tmp_path / "result-report.json"
    summary = tmp_path / "result-summary.json"
    _write_json(source, payload)
    completed = _run(source, output)
    assert completed.returncode == 2
    assert completed.stdout == ""
    assert completed.stderr.startswith("REFUSED: ")
    assert not output.exists()
    assert not report.exists()
    assert not summary.exists()


def test_zero_dispersion_is_refused(tmp_path: Path) -> None:
    """题面：
    「`nm_per_pixel`，二者均为有限数，后者必须为正。像素读数按」
    「的速度；若有多个速度给出同样最多的有效对应，则结果只能包含 `status` 与 `candidates`，」
    「歧义结果里的 `candidates` 是一个对象，键只列出现并列最优的帧号，值为该帧的候选数量。」
    """

    # CASE-CAL-03
    catalog = [
        {
            "id": f"LINE-{index:04d}",
            "wavelength_nm": 500.0 + 0.00006 * index
            if index < 16
            else 600.0 + index,
        }
        for index in range(32)
    ]
    exposures = [
        {
            "frame_id": f"FRAME-{index:03d}",
            "pixel_positions": [float(pixel) for pixel in range(1, 17)],
            "calibration": {"pixel_zero": 500.0, "nm_per_pixel": 0.0},
        }
        for index in range(6)
    ]
    payload = {
        "schema": "astro-lines-v1",
        "catalog": catalog,
        "exposures": exposures,
    }
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    report = tmp_path / "result-report.json"
    journal = tmp_path / "result-journal.jsonl"
    _write_json(source, payload)
    completed = _run(source, output)
    assert completed.returncode == 2
    assert completed.stdout == ""
    assert completed.stderr.startswith("REFUSED: ")
    assert not output.exists()
    assert not report.exists()
    assert not journal.exists()


def test_gate_maps_every_peak_to_a_distinct_line(tmp_path: Path) -> None:
    """题面：
    「`catalog` 包含 32 到 80 条谱线。每条只能有 `id`、`wavelength_nm`，编号匹配」
    """

    payload, _ = _payload(6)
    source = tmp_path / "input.json"
    output = tmp_path / "result.json"
    _write_json(source, payload)
    completed = _run(source, output)
    result = json.loads(output.read_text(encoding="utf-8"))
    first = result["frames"]["FRAME-000"]["line_matches"]
    assert completed.returncode == 0
    assert completed.stdout == ""
    assert completed.stderr == ""
    assert len(first) == 24
    assert sorted(first.values()) == list(range(24)), "CASE-MAP-02"
