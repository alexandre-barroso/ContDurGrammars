"""Run two-moment normalizer calculations on four fixed synthetic models.

A time budget is a computational target, not a hardware-independent correctness
criterion. The recorded model trees specify the calculation inputs.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

sys.dont_write_bytecode = True
from integration_two_moments import partition_two
import flint
from flint import ctx
from integration_bounds import ball, fraction_string
from integration_scaling import CASES, build_case, describe

HERE = Path(__file__).resolve().parent
DEADLINE = 60


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def peak_bytes():
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw if platform.system() == "Darwin" else raw*1024 if platform.system() == "Linux" else None


def sources():
    names = ["integration_two_moment_scaling.py", "integration_two_moments.py",
             "TWO_MOMENT.md", "integration_bounds.py", "integration_moments.py",
             "integration_scaling.py", "published_timing_parameters.json",
             "results/moment_comparison.json"]
    return {name: sha(HERE/name) for name in names}


def save(path, data):
    path = path.resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2)+"\n")
    temporary.replace(path)


def worker(case, path):
    if path.exists():
        raise FileExistsError("Use a fresh worker output file")
    ctx.prec = 128
    node = build_case(case)
    tree = describe(node)
    tree_hash = hashlib.sha256(json.dumps(tree, sort_keys=True).encode()).hexdigest()
    historical = json.loads((HERE/"results/moment_comparison.json").read_text())["cases"][case]["three"]
    if historical["tree_sha256"] != tree_hash:
        raise AssertionError("frozen tree mismatch")
    initial_sources = sources()
    started = time.perf_counter()
    result = {"case": case, "status": "RUNNING", "dimension": len(node.leaves),
              "tree": tree, "tree_sha256": tree_hash, "attempts": [],
              "started_at": datetime.now(timezone.utc).isoformat(),
              "method": "two_component_mass_and_signed_first_enclosures",
              "precision_bits": ctx.prec, "python_flint": flint.__version__,
              "sources": initial_sources,
              "memory_note": "ru_maxrss bytes on macOS, KiB converted to bytes on Linux, unavailable on other platforms; observed after completed meshes; timeout records are lower bounds"}
    save(path, result)
    for j in range(15):
        h = F(20, 2**j)
        result["in_progress_mesh_ms"] = str(h)
        save(path, result)
        before = time.perf_counter()
        lower, upper, tail, width, count = partition_two(node, h, F(800))
        reached = bool(width.is_finite() and width <= ball(F(1,1000)))
        row = {"mesh_ms": str(h), "seconds": time.perf_counter()-before,
               "cumulative_worker_seconds": time.perf_counter()-started,
               "Z_lower_exact_rational": fraction_string(lower),
               "Z_upper_exact_rational": fraction_string(upper),
               "tail_upper_exact_rational": fraction_string(tail),
               "log_width_upper_exact_rational": fraction_string(width) if width.is_finite() else None,
               "log_width_approx": float(width) if width.is_finite() else None,
               "target_met_by_implemented_bounds": reached,
               "coefficients_at_root": count,
               "peak_resident_bytes_after_completed_mesh": peak_bytes()}
        result["attempts"].append(row)
        del result["in_progress_mesh_ms"]
        if reached:
            result["status"] = "TARGET_REACHED_PENDING_CUMULATIVE_PROCESS_BUDGET"
        save(path, result)
        print(json.dumps({"case": case, "mesh_ms": str(h), "log_width": row["log_width_approx"],
                          "cumulative_worker_seconds": row["cumulative_worker_seconds"]}), flush=True)
        if reached:
            break
    else:
        result["status"] = "MAXIMUM_MESH_WITHOUT_TARGET"
    if sources() != initial_sources:
        raise AssertionError("source changed during worker")
    result["source_unchanged_during_worker"] = True
    save(path, result)
    return bool(result["attempts"] and result["attempts"][-1]["target_met_by_implemented_bounds"])


def main(output):
    if output.exists():
        raise FileExistsError("Output exists; choose a new result directory")
    output.mkdir(parents=True)
    initial_sources = sources()
    result = {"started_at": datetime.now(timezone.utc).isoformat(), "status": "RUNNING", "cases": [],
              "process_budget_seconds_per_case": DEADLINE, "sources": initial_sources,
              "platform": platform.platform(), "machine": platform.machine(), "python": sys.version,
              "interpretation": "Calculation on fixed models. Missing a time budget is not a mathematical counterexample."}
    save(output/"results.json", result)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    for case, *_ in CASES:
        path = (output/(case+".json")).resolve()
        before = time.perf_counter()
        timed_out = False
        try:
            process = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", case,
                                      "--output", str(path)], timeout=DEADLINE, check=False, env=env)
            exit_code = process.returncode
        except subprocess.TimeoutExpired:
            # subprocess.run kills and waits for this direct child; no subprocesses in worker.
            timed_out, exit_code = True, None
        elapsed = time.perf_counter()-before
        record = json.loads(path.read_text()) if path.exists() else {"case": case, "attempts": []}
        record.update(process_elapsed_seconds=elapsed, process_exit_code=exit_code,
                      hard_timeout=timed_out, child_reaped=True)
        numerical_target = bool(record.get("attempts") and record["attempts"][-1]["target_met_by_implemented_bounds"])
        record["numerical_target_reached"] = numerical_target
        if timed_out or elapsed > DEADLINE:
            record["status"] = "TIME_BUDGET_EXCEEDED"
        elif exit_code != 0:
            record["status"] = "PROCESS_ERROR"
        elif numerical_target:
            record["status"] = "FEASIBILITY_TARGET_MET"
        else:
            record["status"] = "MAXIMUM_MESH_WITHOUT_TARGET"
        save(path, record)
        result["cases"].append({"case": case, "status": record["status"], "seconds": elapsed,
                                "numerical_target_reached": numerical_target, "result": path.name})
        save(output/"results.json", result)
    if sources() != initial_sources:
        raise AssertionError("source changed during experiment")
    result.update(status="ALL_FOUR_FROZEN_CASES_HANDLED", source_unchanged_during_run=True,
                  completed_at=datetime.now(timezone.utc).isoformat())
    save(output/"results.json", result)
    return 0 if all(row["status"] == "FEASIBILITY_TARGET_MET" for row in result["cases"]) else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", choices=[c[0] for c in CASES])
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if sys.flags.optimize:
        raise RuntimeError("Do not disable numerical assertions with -O")
    raise SystemExit((0 if worker(args.worker, args.output) else 1) if args.worker else main(args.output))
