"""Standard-library tests for deterministic table and plot-data generation."""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import unittest
from make_integration_artifacts import build_tables

ROOT=Path(__file__).resolve().parent

class ArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tables, cls.refinement, cls.checks = build_tables(ROOT)

    def test_five_tables(self):
        self.assertEqual(set(self.tables), {'parameter_table.tex','normalizer_table.tex',
            'query_table.tex','worked_query.tex','query_contrasts.tex'})

    def test_tables_match_local_outputs(self):
        for name,text in self.tables.items():
            with self.subTest(name=name):
                self.assertEqual(text,(ROOT/'results/tables'/name).read_text())

    def test_refinement_count_and_target(self):
        data=self.refinement
        self.assertEqual(F(data['target']),F(1,1000))
        self.assertEqual(sum(len(m['points']) for c in data['cases'].values() for m in c.values()),46)

    def test_stopping_widths_and_meshes(self):
        for case,methods in self.refinement['cases'].items():
            for method,record in methods.items():
                with self.subTest(case=case,method=method):
                    points=record['points']
                    for i,p in enumerate(points):
                        self.assertEqual(F(p['mesh_ms']),F(20,2**i))
                    self.assertLessEqual(F(points[-1]['log_width_upper_exact_rational']),F(1,1000))
                    self.assertTrue(all(F(p['log_width_upper_exact_rational'])>F(1,1000) for p in points[:-1]))

    def test_query_records_checked_before_export(self):
        self.assertEqual(self.checks['status'],'PASSED_STORED_RECORD_AND_DERIVATION_CHECKS')
        self.assertEqual(self.checks['ratio_hulls_checked'],116)
        self.assertEqual(self.checks['common_normalizer_mesh_containments'],19)

if __name__=='__main__':
    unittest.main(verbosity=2)
