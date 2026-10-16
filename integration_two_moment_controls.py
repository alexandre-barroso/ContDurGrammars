"""Finite controls and exact endpoint replay for the two-component normalizer.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import time

from integration_two_moments import (TwoMessage, leaf_two, merge_two,
                                     factor_two, enclose_two, partition_two)
from integration_bounds import Potential, Node, ball, fraction_string
from flint import arb, ctx
import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent


def mm(value):
    q = F(value)
    return mp.mpf(q.numerator)/q.denominator


def ends(value):
    return tuple(F(fraction_string(x)) for x in (value.lower(), value.upper()))


def reference_factor(terms, x):
    cost = sum((mm(a if x >= mm(t) else b)*(x-mm(t))**2 for t,a,b in terms), mp.mpf(0))
    return mp.exp(-cost)


def assert_contains_number(enclosure, value, slack=mp.mpf(0)):
    lo, hi = (mm(x) for x in ends(enclosure))
    assert lo-slack <= value <= hi+slack, (str(lo), str(value), str(hi))


def check_leaves():
    fixtures = [
        ((), F(-2), F(3), F(0)),
        (((F(1,7), F(2,3), F(5,4)),), F(-1), F(2), F(1,2)),
        (((F(-2,3), F(0), F(2)), (F(4,3), F(3), F(0))), F(-2), F(2), F(1,4)),
        (((F(1000000), F(2,7), F(3,7)),), F(999999), F(1000002), F(2000001,2)),
        (((F(1,2), F(7), F(11)),), F(499,1000), F(501,1000), F(1,2)),
        (((0.1, 0.3, 0.7),), F(0), F(1), F(1,2)),
    ]
    records = []
    for terms, lo, hi, center in fixtures:
        potential = Potential(terms)
        actual = leaf_two(potential, lo, hi, center)
        cuts = [mm(lo), *(mm(t) for t in potential.knots if lo < t < hi), mm(hi)]
        refs = [mp.quad(lambda x: (x-mm(center))**j*reference_factor(potential.terms,x), cuts)
                for j in (0,1)]
        for enclosure, value in zip(actual, refs):
            assert_contains_number(enclosure, value, mp.mpf('1e-65')*max(1,abs(value)))
        records.append({"terms": [[str(x) for x in t] for t in potential.terms],
                        "bounds": [str(lo),str(hi)], "center": str(center),
                        "reference_moments": [str(x) for x in refs],
                        "enclosures": [[str(x) for x in ends(v)] for v in actual]})
    empty = leaf_two(Potential(), F(1), F(1), F(1))
    assert all(ends(x) == (F(0),F(0)) for x in empty)
    return {"cases": records, "zero_length_leaf": "exact zero", "quadrature_digits": mp.mp.dps,
            "reference_tolerance": "1e-65*max(1,abs(reference)); not certified quadrature"}


def check_merge():
    aa = [[F(1),F(0),F(2)],[F(-1,4),F(0),F(1,2)]]
    bb = [[F(0),F(3),F(0)],[F(0),F(-1,2),F(0)]]
    a = TwoMessage([[ball(q) for q in row] for row in aa], 1, F(2))
    b = TwoMessage([[ball(q) for q in row] for row in bb], 1, F(-3))
    c = merge_two(a,b)
    expected = [[F(0)]*5 for _ in range(2)]
    for i in range(3):
        for j in range(3):
            expected[0][i+j] += aa[0][i]*bb[0][j]
            expected[1][i+j] += aa[1][i]*bb[0][j]+aa[0][i]*bb[1][j]
    for row, exact in zip(c.moments, expected):
        assert len(row) == 5
        for value,q in zip(row, exact):
            lo, hi = ends(value)
            assert lo <= q <= hi
    assert c.leaves == 2 and c.origin == F(-1)
    return {"signed_merge_exact": True, "preserved_length": 5,
            "trailing_zero_cells_preserved": True, "expected": [[str(q) for q in r] for r in expected]}


def check_factors():
    pots = [Potential(), Potential(((F(1,3),F(2,7),F(3,5)),)),
            Potential(((F(-2,5),F(0),F(4)),(F(7,5),F(3),F(0)))),
            Potential(((F(1),F(1,100),F(9)),))]
    atoms = [(F(-1),F(1,4)), (F(-1,4),F(1,2)), (F(3,4),F(1,4))]
    m0 = sum((p for _,p in atoms), F(0))
    m1 = sum((u*p for u,p in atoms), F(0))
    records = []
    for center in (F(-1,3), F(1,2), F(2)):
        for potential in pots:
            actual = factor_two((ball(m0),ball(m1)), potential, center, F(1))
            refs = [sum((mm(p)*mm(u)**j*reference_factor(potential.terms,mm(center+u))
                         for u,p in atoms), mp.mpf(0)) for j in (0,1)]
            for enclosure, value in zip(actual,refs):
                assert_contains_number(enclosure,value,mp.mpf('1e-65'))
            zero = factor_two((arb(0),arb(0)),potential,center,F(1))
            assert all(ends(x) == (F(0),F(0)) for x in zero)
            records.append({"center":str(center), "terms":[[str(x) for x in t] for t in potential.terms],
                            "reference":[str(x) for x in refs],
                            "enclosures":[[str(x) for x in ends(v)] for v in actual]})
    return {"cases":records, "all_zero_measure_updates_exact":True,
            "input_signed_first":str(m1), "scope":"finite atomic integral checks, not quadrature of a density"}


def nested_model():
    leaves = [Node(Potential(((F(t),F(a),F(b)),))) for t,a,b in
              [(F(2,5),F(2,7),F(1,3)),(F(7,5),F(1,4),F(3,8)),(F(6,5),F(2,9),F(1,5))]]
    pair = Node(Potential(((F(7,5),F(1,11),F(2,13)),)), tuple(leaves[:2]))
    return Node(Potential(((F(16,5),F(1,17),F(2,19)),)), (pair,leaves[2]))


def independent_tensor(order):
    p,w = np.polynomial.legendre.leggauss(order)
    p,w = (p+1)/2,w/2
    indices = np.array(list(itertools.product(range(order),repeat=3)))
    offsets = p[indices]
    weights = np.prod(w[indices],axis=1)
    result = np.zeros((2,4))
    def penalty(x,t,a,b):
        return np.where(x>=t,a,b)*(x-t)**2
    for cells in itertools.product(range(2),repeat=3):
        x = offsets+np.array(cells)
        cost = (penalty(x[:,0],2/5,2/7,1/3)+penalty(x[:,1],7/5,1/4,3/8)
                +penalty(x[:,2],6/5,2/9,1/5)+penalty(x[:,0]+x[:,1],7/5,1/11,2/13)
                +penalty(x.sum(axis=1),16/5,1/17,2/19))
        k = sum(cells)
        weighted = weights*np.exp(-cost)
        result[0,k] += weighted.sum()
        result[1,k] += np.dot(weighted,x.sum(axis=1)-(k+1.5))
    return result


def check_nested():
    message = enclose_two(nested_model(),F(1),F(0),F(2),{})
    a,b = independent_tensor(20),independent_tensor(32)
    delta = float(np.max(np.abs(a-b)))
    assert delta < 2e-4, delta
    for j in range(2):
        for k,value in enumerate(message.moments[j]):
            lo,hi = map(float,ends(value))
            assert lo <= a[j,k] <= hi and lo <= b[j,k] <= hi, (j,k,lo,a[j,k],b[j,k],hi)
    return {"cartesian_cells":8, "moment_values_checked":8, "orders":[20,32],
            "maximum_reference_change":delta, "scope":"independent tensor estimates; not certified oracle values"}


def load_records():
    """Validate recorded calculations against the local model definitions."""
    from integration_scaling import CASES, build_case, describe
    records = json.loads((HERE / "results/moment_comparison.json").read_text())
    expected = [row[0] for row in CASES]
    assert list(records["cases"]) == expected
    assert records["design"] == {
        "cube_upper": "800", "initial_mesh": "20", "maximum_mesh_index": 14,
        "precision_bits": 128, "log_width_target": "1/1000",
        "process_budget_seconds": 60}
    common = []
    for case, methods in records["cases"].items():
        tree = describe(build_case(case))
        dimension = len(build_case(case).leaves)
        for method in ("two", "three"):
            record = methods[method]
            assert record["tree"] == tree
            assert record["precision_bits"] == 128
            assert record["status"] == "FEASIBILITY_TARGET_MET"
            assert 0 < record["process_seconds"] <= 60
            attempts = record["attempts"]
            assert 1 <= len(attempts) <= 15
            for index, attempt in enumerate(attempts):
                mesh = F(20, 2**index)
                assert F(attempt["mesh_ms"]) == mesh
                assert 0 < F(attempt["Z_lower_exact_rational"]) <= F(attempt["Z_upper_exact_rational"])
                assert F(attempt["tail_upper_exact_rational"]) >= 0
                width = F(attempt["log_width_upper_exact_rational"])
                assert width >= 0
                assert attempt["target_met_by_implemented_bounds"] == (width <= F(1, 1000))
                assert attempt["coefficients_at_root"] == dimension * (int(F(800)/mesh)-1)+1
                assert attempt["seconds"] >= 0
                assert attempt["target_met_by_implemented_bounds"] == (index == len(attempts)-1)
        two, three = methods["two"]["attempts"], methods["three"]["attempts"]
        assert F(two[-1]["mesh_ms"])*4 == F(three[-1]["mesh_ms"])
        by_mesh = {a["mesh_ms"]: a for a in two}
        for a3 in three:
            a2 = by_mesh[a3["mesh_ms"]]
            assert a2["tail_upper_exact_rational"] == a3["tail_upper_exact_rational"]
            assert F(a2["Z_lower_exact_rational"]) <= F(a3["Z_lower_exact_rational"])
            assert F(a3["Z_upper_exact_rational"]) <= F(a2["Z_upper_exact_rational"])
            common.append({"case": case, "mesh": a3["mesh_ms"],
                "three_inside_two": True,
                "two_over_three_width": float(F(a2["log_width_upper_exact_rational"])
                                              / F(a3["log_width_upper_exact_rational"]))})
    assert len(common) == 19
    return records, {"cases": expected, "common_meshes": common,
        "scope": "Integrity/consistency of recorded calculations, not recomputation of every mesh or authentication against deliberate metadata rewriting"}


def replay_records(which="first", progress=None):
    """Recompute one exact-endpoint mesh for each case, without a time criterion.

    `first` is the small portable smoke test. `last` recomputes the four final
    meshes. Neither mode reruns adaptive timing or requires completion in 60s.
    """
    from integration_scaling import build_case
    if which not in {"first", "last"}:
        raise ValueError("Choose first or last recorded mesh")
    records, audit = load_records()
    results = []
    for case, methods in records["cases"].items():
        expected = methods["two"]["attempts"][0 if which == "first" else -1]
        before = time.perf_counter()
        lo, hi, tail, width, count = partition_two(build_case(case), F(expected["mesh_ms"]), F(800))
        actual = {"Z_lower_exact_rational": fraction_string(lo),
                  "Z_upper_exact_rational": fraction_string(hi),
                  "tail_upper_exact_rational": fraction_string(tail),
                  "log_width_upper_exact_rational": fraction_string(width),
                  "coefficients_at_root": count}
        matches = {name: value == expected[name] for name, value in actual.items()}
        row = {"case": case, "mesh": expected["mesh_ms"], "matches": matches,
               "actual": actual, "validation_seconds": time.perf_counter()-before,
               "timing_role": "Verification time only; not an experiment timing measurement"}
        results.append(row)
        if progress:
            progress(results)
        print(json.dumps({"case": case, "mesh": row["mesh"],
                          "exact_endpoint_replay_passed": all(matches.values()),
                          "validation_seconds": row["validation_seconds"]}), flush=True)
        assert all(matches.values()), (case, matches)
    return {"mesh_selection": which, "cases": results, "record_audit": audit,
            "time_limit_is_a_correctness_requirement": False,
            "scope": "Exact rational endpoint reproduction in the pinned numerical environment; timings are not compared"}


def write_tables(directory):
    """Generate parameter and normalizer-summary tables from local records."""
    records, _ = load_records()
    labels = ["Balanced 8", "Balanced 16", "Left branching 16", "Stiff balanced 16"]
    lines = [r"\begin{tabular}{lrrrr}\toprule",
             r"Case & $h_2$ (ms) & $h_3$ (ms) & Two (s) & Three (s)\\\midrule"]
    for label, methods in zip(labels, records["cases"].values()):
        two, three = methods["two"], methods["three"]
        h2, h3 = (float(F(v["attempts"][-1]["mesh_ms"])) for v in (two, three))
        lines.append(f"{label} & {h2:g} & {h3:g} & {two['process_seconds']:.2f} & {three['process_seconds']:.2f} " + r"\\")
    lines.append(r"\bottomrule\end{tabular}")
    tables = {"two_moment_table.tex": "\n".join(lines)+"\n"}
    parameters = json.loads((HERE / "published_timing_parameters.json").read_text())
    weights = parameters["weights_in_column_order_AU_to_BB"]
    lines = [r"\begin{tabular}{crrr}\toprule",
             r"$j$ & $t_j$ & Raw $w_{j,+}$ & Raw $w_{j,-}$\\\midrule"]
    for i, key in enumerate(["mora", "syllable", "grid", "hemimetron"]):
        vals = [parameters["targets"][key], *weights[2*i:2*i+2]]
        lines.append(str(i+1) + " & " + " & ".join(format(float(F(x)), ".17g") for x in vals) + " " + r"\\")
    lines.append(r"\bottomrule\end{tabular}")
    tables["parameter_table.tex"] = "\n".join(lines)+"\n"
    records = []
    for name, contents in tables.items():
        dest = Path(directory) / name
        with dest.open("x") as stream:
            stream.write(contents)
        records.append({"path": name})
    return records


def run_controls(output, *, final_replay_only=False):
    """Save finite controls plus portable records replay to an unoccupied JSON path."""
    import sys
    if sys.flags.optimize:
        raise RuntimeError("Do not disable assertion-based controls with -O")
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError("Output exists; choose a new file")
    output.parent.mkdir(parents=True, exist_ok=True)
    ctx.prec, mp.mp.dps = 128, 80
    names = ["integration_two_moments.py", "integration_two_moment_controls.py",
             "integration_bounds.py", "integration_moments.py", "integration_scaling.py",
             "published_timing_parameters.json", "results/moment_comparison.json", "TWO_MOMENT.md"]
    source_contents = lambda: {name: (HERE/name).read_bytes() for name in names}
    initial_sources = source_contents()
    before = time.perf_counter()
    record = {"started_at": datetime.now(timezone.utc).isoformat(), "status": "RUNNING",
              "precision_bits": 128, "sources": names,
              "scope": "Finite mathematical fixtures and exact endpoint replay; not a formal primitive proof or a timed performance experiment"}
    def save():
        temporary = output.with_suffix(output.suffix + ".tmp")
        temporary.write_text(json.dumps(record, indent=2)+"\n")
        temporary.replace(output)
    def progress(rows):
        record["partial_replay"] = rows
        save()
    save()
    try:
        if not final_replay_only:
            record.update(leaves=check_leaves(), merge=check_merge(),
                          factors=check_factors(), nested=check_nested())
            save()
        record["record_replay"] = replay_records("last" if final_replay_only else "first", progress)
        record.pop("partial_replay", None)
        record["generated_tables"] = write_tables(output.parent)
        assert source_contents() == initial_sources, "sources changed during verification"
        record["status"] = "ALL_REQUESTED_CONTROLS_PASSED"
    except Exception as exc:
        record.update(status="CONTROL_FAILED", error=type(exc).__name__ + ": " + str(exc))
        raise
    finally:
        record["seconds"] = time.perf_counter()-before
        save()
    return record


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--final-replay-only", action="store_true")
    args = parser.parse_args()
    run_controls(args.output, final_replay_only=args.final_replay_only)
