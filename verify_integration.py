"""Run numerical controls, formal checks and timed model calculations.

Uses adjacent source modules, the pinned numerical dependencies, and an
explicitly selected mathlib installation for formal checks.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PINS = {"python-flint": "0.8.0", "numpy": "2.4.4", "scipy": "1.17.1", "mpmath": "1.3.0"}
MODULES = ["MomentEnclosure", "MomentConvolution", "MomentCircuit", "MomentBridge"]
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}


def save(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")


def numerical(output):
    import mpmath as mp
    from flint import ctx
    from integration_moment_controls import (
        run_leaf_controls, run_derivative_controls, run_nested_controls)
    import integration_query_controls as query_controls

    ctx.prec = 128
    mp.mp.dps = 80
    records = {"leaf": run_leaf_controls(),
               "derivatives": run_derivative_controls(),
               "nested": run_nested_controls()}
    save(output / "moment_controls.json", records)
    result = output / "query_controls.json"
    query_controls.main(result)
    data = json.loads(result.read_text())
    if data["status"] != "ALL_CONTROLS_PASSED":
        raise RuntimeError("Query controls did not pass")
    save(output / "query_controls.json", data)
    two = two_moment(output)
    return {"passed": True, "moment_controls": "moment_controls.json",
            "query_controls": "query_controls.json",
            "two_moment_controls": two,
            "limits": "Finite numerical comparisons, not proofs of the reference quadrature or the implementation"}


def two_moment(output, *, final_replay_only=False):
    from integration_two_moment_controls import run_controls
    name = "two_moment_final_replay.json" if final_replay_only else "two_moment_controls.json"
    record = run_controls(output / name, final_replay_only=final_replay_only)
    if record["status"] != "ALL_REQUESTED_CONTROLS_PASSED":
        raise RuntimeError("Two-component controls/replay failed")
    return {"passed": True, "result": name,
            "scope": "Four final recorded meshes, no timing requirement" if final_replay_only else
                     "Finite controls and four first-mesh endpoint replays"}


def proofs(output, mathlib):
    if mathlib is None:
        raise ValueError("--mathlib is required for proof verification")
    mathlib = mathlib.resolve(strict=True)

    def probe(command):
        return subprocess.run(command, cwd=mathlib, capture_output=True,
                              text=True, timeout=60, check=True).stdout.strip()

    version = probe(["lake", "env", "lean", "--version"])
    mathlib_toolchain = (mathlib/"lean-toolchain").read_text().strip()
    if "version 4.34.0" not in version or mathlib_toolchain != "leanprover/lean4:v4.34.0":
        raise RuntimeError("Lean/mathlib versions differ from the required pins")
    lean = Path(probe(["lake", "env", "lean", "--print-prefix"])) / "bin/lean"
    import_path = probe(["lake", "env", "printenv", "LEAN_PATH"])
    stage = output / "lean"
    stage.mkdir(exist_ok=False)
    env = os.environ.copy()
    env["LEAN_PATH"] = str(stage) + os.pathsep + import_path
    record = {"Lean": version, "mathlib_toolchain": mathlib_toolchain, "modules": [],
              "limits": "Required Lean toolchain and supplied mathlib installation; no independent kernel or Python/Arb verification"}
    for module in MODULES:
        source = HERE / (module + ".lean")
        text = source.read_text()
        targets = re.findall(r"^#print axioms ([\w.]+)\s*$", text, re.M)
        if not targets or re.search(r"^\s*(?:axiom|opaque|unsafe)\s|\b(?:sorry|admit|native_decide)\b", text, re.M):
            raise RuntimeError("Forbidden proof shortcut or absent targets")
        staged = stage / source.name
        shutil.copyfile(source, staged)
        command = [str(lean), "--trust=0", "--root=" + str(stage),
                   "-o", str(stage / (module + ".olean")), str(staged)]
        cp = subprocess.run(command, cwd=mathlib, env=env, capture_output=True,
                            text=True, timeout=180)
        log = stage / (module + ".log")
        log.write_text(cp.stdout + cp.stderr)
        observed = {name: [v.strip() for v in values.split(",") if v.strip()]
                    for name, values in re.findall(
                        r"'([^']+)' depends on axioms:\s*\[([^]]*)\]", log.read_text())}
        for name in re.findall(r"'([^']+)' does not depend on any axioms", log.read_text()):
            observed[name] = []
        passed = (cp.returncode == 0 and set(observed) == set(targets)
                  and all(set(a) <= ALLOWED_AXIOMS for a in observed.values())
                  and (stage / (module + ".olean")).is_file())
        record["modules"].append({"module": module, "exit_code": cp.returncode, "axioms": observed, "passed": passed, "log": str(log.relative_to(output))})
        save(output / "proofs.json", record)
        if not passed:
            raise RuntimeError("Strict proof/axiom check failed: " + module)
    record["passed"] = True
    save(output / "proofs.json", record)
    return record


def experiments(output):
    from integration_queries import QUERY_CASES
    from integration_scaling import CASES
    jobs = [("moment_scaling", "integration_moment_scaling.py", [c[0] for c in CASES], 60,
             "TOLERANCE_REACHED_PENDING_PROCESS_BUDGET_CHECK"),
            ("two_moment_scaling", "integration_two_moment_scaling.py", [c[0] for c in CASES], 60,
             "TARGET_REACHED_PENDING_CUMULATIVE_PROCESS_BUDGET"),
            ("queries", "integration_queries.py", QUERY_CASES, 180,
             "QUERY_TOLERANCE_REACHED_PENDING_BUDGET_CHECK"),
            ("ordinary_fft", "integration_point_queries.py", QUERY_CASES, 180,
             "ALL_ELEVEN_MESHES_FINISHED")]
    record = {"runs": [], "passed": True,
              "scope": "Timed model calculations; recorded wall times are not expected to be identical. No GuBPI execution is included."}
    for group, filename, cases, budget, expected_status in jobs:
        directory = output / group
        directory.mkdir()
        for case in cases:
            dest = directory / (case + ".json")
            command = [sys.executable, str(HERE / filename), "--worker", case,
                       "--output", str(dest)]
            started = time.perf_counter()
            with (directory / (case + ".log")).open("w") as log:
                try:
                    cp = subprocess.run(command, cwd=output, stdout=log,
                                        stderr=subprocess.STDOUT, timeout=budget)
                    code, timed_out = cp.returncode, False
                except subprocess.TimeoutExpired:
                    code, timed_out = None, True
            elapsed = time.perf_counter() - started
            data = json.loads(dest.read_text()) if dest.exists() else {}
            passed = code == 0 and not timed_out and data.get("status") == expected_status and elapsed <= budget
            row = {"group": group, "case": case, "seconds": elapsed,
                   "budget_seconds": budget, "exit_code": code, "timed_out": timed_out,
                   "worker_status": data.get("status"), "passed": passed,
                   "result": str(dest.relative_to(output))}
            record["runs"].append(row)
            record["passed"] = record["passed"] and passed
            save(output / "experiments.json", record)
            print(json.dumps(row), flush=True)
    if not record["passed"]:
        raise RuntimeError("One or more experiment targets or executions failed; inspect retained results")
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["numerical", "proofs", "all", "experiments",
                                           "two-moment", "replay-final"], default="all")
    parser.add_argument("--mathlib", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sys.flags.optimize:
        raise RuntimeError("Do not disable assertion-based numerical controls with -O")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    record = {"started_at": datetime.now(timezone.utc).isoformat(), "status": "RUNNING", "mode": args.mode, "Python": sys.version, "dependencies": {}, "limitations": "Finite numerical controls and local formal obligations; not whole-program verification"}
    try:
        if args.mode != "proofs":
            for package, pin in PINS.items():
                version = importlib.metadata.version(package)
                record["dependencies"][package] = version
                if version != pin:
                    raise RuntimeError(f"Version mismatch: {package} {version}, expected {pin}")
            if sys.version_info[:2] != (3, 12):
                raise RuntimeError("Use the recorded Python 3.12 series for this verification")
        if args.mode in ["all", "numerical"]:
            record["numerical"] = numerical(output)
        if args.mode in ["all", "proofs"]:
            record["proofs"] = proofs(output, args.mathlib)
        if args.mode == "experiments":
            record["experiments"] = experiments(output)
        if args.mode in ["two-moment", "replay-final"]:
            record["two_moment"] = two_moment(output, final_replay_only=args.mode == "replay-final")
        record["status"] = "PASSED_REQUESTED_SCOPE"
        returncode = 0
    except Exception as exc:
        record["status"] = "FAILED"
        record["error"] = type(exc).__name__ + ": " + str(exc)
        returncode = 1
    record["completed_at"] = datetime.now(timezone.utc).isoformat()
    save(output / "verification.json", record)
    print(json.dumps({"status": record["status"], "mode": args.mode,
                      "error": record.get("error")}), flush=True)
    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
