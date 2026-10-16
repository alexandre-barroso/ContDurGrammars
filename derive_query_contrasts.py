"""Derive model contrasts and check stored records using exact rationals.

Uses the Python standard library. It checks stored intervals and model identities;
it does not run the enclosing integrator or build the formal modules.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys
from typing import Any

F = Fraction
CASES = ("HH0H", "HHLH", "HLLH", "balanced_8", "balanced_16",
         "left_branching_16", "stiff_balanced_16")
NAMES = ("P_A", "P_B", "P_AB", "P_A_given_B", "Delta")
SYNTHETIC = CASES[3:]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


@dataclass(frozen=True)
class Interval:
    lower: F
    upper: F

    def __post_init__(self) -> None:
        require(self.lower <= self.upper, "Reversed interval")

    @classmethod
    def read(cls, data: dict[str, Any]) -> Interval:
        return cls(F(data["lower"]), F(data["upper"]))

    def contains(self, other: Interval) -> bool:
        return self.lower <= other.lower <= other.upper <= self.upper

    def subtract(self, other: Interval) -> Interval:
        return Interval(self.lower - other.upper, self.upper - other.lower)

    def encode(self) -> dict[str, str]:
        return {"lower": str(self.lower), "upper": str(self.upper)}


def decimal_bound(value: F, places: int = 7, *, upper: bool = False) -> str:
    """Round an exact rational down/up without a floating-point conversion."""
    require(isinstance(places, int) and places >= 0, "Nonnegative precision required")
    scale = 10 ** places
    numerator = value.numerator * scale
    integer = numerator // value.denominator
    if upper and numerator % value.denominator:
        integer += 1
    sign = "-" if integer < 0 else ""
    integer = abs(integer)
    if places == 0:
        return sign + str(integer)
    return f"{sign}{integer // scale}.{integer % scale:0{places}d}"


def ordered_leaves(tree: dict[str, Any]) -> list[dict[str, Any]]:
    if not tree["children"]:
        return [tree]
    return [leaf for child in tree["children"] for leaf in ordered_leaves(child)]


def check_stiffness(base: dict[str, Any], stiff: dict[str, Any]) -> int:
    require(base["name"] == stiff["name"], "Changed node name")
    require(len(base["children"]) == len(stiff["children"]), "Changed branching")
    if not base["children"]:
        require(base == stiff, "Changed unary potential in stiffness contrast")
        return 0
    require(len(base["terms"]) == len(stiff["terms"]), "Changed number of cost terms")
    for t1, t2 in zip(base["terms"], stiff["terms"]):
        target, above, below = map(F, t1)
        require(tuple(map(F, t2)) == (target, 4 * above, 4 * below),
                "Stiffness contrast is not exactly a fourfold internal reweighting")
    return 1 + sum(check_stiffness(a, b)
                   for a, b in zip(base["children"], stiff["children"]))


def probability_hull(numerator: Interval, denominator: Interval) -> Interval:
    require(numerator.lower >= 0 and denominator.lower > 0,
            "Invalid nonnegative-mass/positive-denominator contract")
    return Interval(max(F(0), numerator.lower / denominator.upper),
                    min(F(1), numerator.upper / denominator.lower))


def derive(root: Path) -> tuple[dict[str, Any], str]:
    inputs: list[str] = []

    def read(relative: str) -> dict[str, Any]:
        path = root / relative
        inputs.append(relative)
        return json.loads(path.read_text())

    records = {case: read(f"results/queries/{case}.json")
               for case in CASES}
    fft = {case: read(f"results/point_queries/{case}.json")
           for case in CASES}
    final: dict[str, dict[str, Interval]] = {}
    attempt_count = ratio_count = width_count = fft_count = 0
    for case, record in records.items():
        require(record["status"] == "QUERY_TOLERANCE_MET", f"Unresolved query: {case}")
        require(fft[case]["tree"] == record["tree"], f"FFT model differs: {case}")
        require(len(fft[case]["attempts"]) == 11, f"Unexpected FFT mesh count: {case}")
        require(fft[case]["status"] == "ALL_ELEVEN_MESHES_FINISHED", f"Unfinished FFT: {case}")
        for attempt in record["attempts"]:
            attempt_count += 1
            masses = {k: Interval.read(v) for k, v in attempt["integrals"].items()}
            probabilities = {k: Interval.read(attempt["probabilities"][k]) for k in NAMES}
            rational_ratios: dict[str, Interval] = {}
            for key, numerator, denominator in (
                ("P_A", "A", "Z"), ("P_B", "B", "Z"),
                ("P_AB", "AB", "Z"), ("P_A_given_B", "AB", "B"),
            ):
                expected = probability_hull(masses[numerator], masses[denominator])
                require(probabilities[key].contains(expected), f"Ratio hull not enclosed: {case}/{key}")
                rational_ratios[key] = expected
                ratio_count += 1
            # Check the rational query obtained from the stored mass bounds.
            # Serializing each probability endpoint may add its own outward
            # rounding, so subtracting those exported probability hulls is a
            # different operation from subtracting the Arb balls before endpoint export.
            expected_delta = rational_ratios["P_A_given_B"].subtract(rational_ratios["P_A"])
            require(probabilities["Delta"].contains(expected_delta), f"Delta hull not enclosed: {case}")
            upper_widths = []
            for key in NAMES:
                value = probabilities[key]
                stored_width = F(attempt["probabilities"][key]["width_upper"])
                require(stored_width >= value.upper - value.lower,
                        f"Reported width is not an upper bound: {case}/{key}")
                upper_widths.append(stored_width)
            require(attempt["all_widths_at_most_0_002"] == all(w <= F(1, 500) for w in upper_widths),
                    f"Stopping flag inconsistent: {case}")
        last = record["attempts"][-1]
        require(last["all_widths_at_most_0_002"], f"Final query misses target: {case}")
        final[case] = {k: Interval.read(last["probabilities"][k]) for k in NAMES}
        for key in NAMES:
            width_count += 1
            point = F(fft[case]["attempts"][-1]["probabilities"][key])
            require(final[case][key].contains(Interval(point, point)), f"Final FFT point outside interval: {case}/{key}")
            fft_count += 1

    b, l, s = (records[n] for n in ("balanced_16", "left_branching_16", "stiff_balanced_16"))
    require(ordered_leaves(b["tree"]) == ordered_leaves(l["tree"]) == ordered_leaves(s["tree"]),
            "Model contrasts do not share the same ordered unary potentials")
    require(b["events"] == l["events"] == s["events"], "Model contrasts do not share events")
    require(b["tree"] != l["tree"], "No tree change in hierarchy contrast")
    require(len(ordered_leaves(b["tree"])) == 16, "Contrasts are not sixteen-leaf models")
    internal_stiffness_checks = check_stiffness(b["tree"], s["tree"])

    records = read("results/moment_comparison.json")
    normalizers: dict[str, Any] = {}
    common_meshes = 0
    for case in SYNTHETIC:
        pair = records["cases"][case]
        require(pair["two"]["tree"] == pair["three"]["tree"], f"Normalizer models differ: {case}")
        by_mesh = {method: {F(a["mesh_ms"]): a for a in pair[method]["attempts"]}
                   for method in ("two", "three")}
        common = sorted(set(by_mesh["two"]) & set(by_mesh["three"]))
        for mesh in common:
            a, c = (by_mesh[m][mesh] for m in ("two", "three"))
            enclosing = Interval(F(a["Z_lower_exact_rational"]), F(a["Z_upper_exact_rational"]))
            enclosed = Interval(F(c["Z_lower_exact_rational"]), F(c["Z_upper_exact_rational"]))
            require(enclosing.contains(enclosed), f"Observed interval containment fails: {case}/{mesh}")
            common_meshes += 1
        last = {m: pair[m]["attempts"][-1] for m in ("two", "three")}
        for m in last:
            require(F(last[m]["log_width_upper_exact_rational"]) <= F(1, 1000),
                    f"Final normalizer misses target: {case}/{m}")
        mesh_ratio = F(last["three"]["mesh_ms"]) / F(last["two"]["mesh_ms"])
        require(mesh_ratio == 4, f"Mesh ratio not four: {case}")
        normalizers[case] = {"common_meshes": len(common), "final_mesh_ratio_three_over_two": str(mesh_ratio)}

    contrasts = []
    for case, tex_name, meaning in (
        ("left_branching_16", "L16", "group-constraint hierarchy"),
        ("stiff_balanced_16", "S16", "fourfold internal group-cost stiffness"),
    ):
        difference = final[case]["Delta"].subtract(final["balanced_16"]["Delta"])
        require(difference.lower > 0, f"Model difference unresolved: {case}")
        fft_difference = (F(fft[case]["attempts"][-1]["probabilities"]["Delta"])
                          - F(fft["balanced_16"]["attempts"][-1]["probabilities"]["Delta"]))
        require(difference.contains(Interval(fft_difference, fft_difference)),
                f"FFT model difference outside derived interval: {case}")
        contrasts.append({"case": case, "baseline": "balanced_16", "quantity": "Delta_case_minus_Delta_baseline",
                          "controlled_change": meaning, "exact_interval": difference.encode(),
                          "decimal_outward_7": {"lower": decimal_bound(difference.lower),
                                                "upper": decimal_bound(difference.upper, upper=True)},
                          "strictly_positive": True, "final_fft_difference_inside": True,
                          "tex_name": tex_name})
    rows = [r"\begin{tabular}{llc}\toprule",
            r"Change & Prediction contrast & Enclosure\\\midrule"]
    for contrast in contrasts:
        name = contrast["tex_name"]
        label = "Hierarchy" if name == "L16" else "Stiffness"
        dec = contrast["decimal_outward_7"]
        rows.append(label + r" & $\Delta_{\mathrm{" + name + r"}}-\Delta_{\mathrm{B16}}$ & $["
                    + dec["lower"] + ", " + dec["upper"] + r"]$ \\")
    rows.append(r"\bottomrule\end{tabular}")
    report = {"status": "PASSED_STORED_RECORD_AND_DERIVATION_CHECKS",
              "scope": "Exact rational record checks and derived model contrasts; does not execute integration, proof builds, timing measurements, or speech evaluation",
              "query_attempts_checked": attempt_count, "ratio_hulls_checked": ratio_count,
              "delta_hulls_checked": attempt_count, "final_width_targets_checked": width_count,
              "final_fft_points_inside_stored_intervals": fft_count,
              "shared_ordered_unaries_and_events": True,
              "internal_nodes_checked_for_exact_fourfold_stiffness": internal_stiffness_checks,
              "common_normalizer_mesh_containments": common_meshes,
              "normalizers": normalizers, "contrasts": contrasts,
              "inputs": inputs}
    return report, "\n".join(rows) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path, required=True, help="A new output directory")
    args = parser.parse_args()
    output = args.output.resolve()
    try:
        output.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        parser.error(f"Cannot create fresh output directory: {exc}")
    try:
        report, table = derive(args.root.resolve(strict=True))
        (output / "query_contrasts.tex").write_text(table)
        (output / "derivation_checks.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({k: report[k] for k in ("status", "final_width_targets_checked",
                                               "common_normalizer_mesh_containments")}))
        return 0
    except Exception as exc:
        (output / "derivation_checks.json").write_text(json.dumps({
            "status": "FAILED", "error": f"{type(exc).__name__}: {exc}"}, indent=2) + "\n")
        print(f"FAILED: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
