"""Standard-library tests for exact-rational result analysis.
"""
from __future__ import annotations
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch
import derive_query_contrasts as d

ROOT = Path(__file__).resolve().parent

class DerivedAnalysisTests(unittest.TestCase):
    def test_outward_decimal_positive(self):
        self.assertEqual(d.decimal_bound(F(1, 3)), "0.3333333")
        self.assertEqual(d.decimal_bound(F(1, 3), upper=True), "0.3333334")

    def test_outward_decimal_negative(self):
        self.assertEqual(d.decimal_bound(F(-1, 3)), "-0.3333334")
        self.assertEqual(d.decimal_bound(F(-1, 3), upper=True), "-0.3333333")

    def test_exact_decimal_and_zero(self):
        self.assertEqual(d.decimal_bound(F(-1, 8)), "-0.1250000")
        self.assertEqual(d.decimal_bound(F(0)), "0.0000000")
        self.assertEqual(d.decimal_bound(F(1, 2), 0, upper=True), "1")

    def test_interval_subtraction(self):
        out = d.Interval(F(2), F(4)).subtract(d.Interval(F(1), F(3)))
        self.assertEqual(out, d.Interval(F(-1), F(3)))

    def test_reversed_interval_rejected(self):
        with self.assertRaises(ValueError):
            d.Interval(F(1), F(0))

    def test_probability_contract(self):
        self.assertEqual(d.probability_hull(d.Interval(F(1), F(2)), d.Interval(F(4), F(5))),
                         d.Interval(F(1, 5), F(1, 2)))
        with self.assertRaises(ValueError):
            d.probability_hull(d.Interval(F(0), F(1)), d.Interval(F(0), F(5)))

    def test_stiffness_check_rejects_changed_target(self):
        leaf = {"name": "x1", "terms": [["2", "1", "3"]], "children": []}
        base = {"name": "M", "terms": [["4", "1", "2"]], "children": [leaf]}
        stiff = deepcopy(base)
        stiff["terms"] = [["4", "4", "8"]]
        self.assertEqual(d.check_stiffness(base, stiff), 1)
        stiff["terms"][0][0] = "5"
        with self.assertRaises(ValueError):
            d.check_stiffness(base, stiff)

    def test_recorded_intervals_and_table(self):
        report, table = d.derive(ROOT)
        self.assertEqual(report["final_width_targets_checked"], 35)
        self.assertEqual(report["common_normalizer_mesh_containments"], 19)
        self.assertEqual(report["internal_nodes_checked_for_exact_fourfold_stiffness"], 15)
        self.assertIn("[0.0015193, 0.0037408]", table)
        self.assertIn("[0.0012885, 0.0032420]", table)
        self.assertEqual(table, (ROOT / "results/derived/query_contrasts.tex").read_text())

    def test_altered_comparison_event_rejected(self):
        # Patch only the file read, not a value in the supplied input files.
        original = Path.read_text
        def read(path, *args, **kwargs):
            text = original(path, *args, **kwargs)
            if path == ROOT / "results/queries/stiff_balanced_16.json":
                value = json.loads(text)
                event = value["events"]["A"]
                name = next(iter(event))
                event[name][0] = "0"
                return json.dumps(value)
            return text
        with patch.object(Path, "read_text", read):
            with self.assertRaisesRegex(ValueError, "share events"):
                d.derive(ROOT)

if __name__ == "__main__":
    unittest.main(verbosity=2)
