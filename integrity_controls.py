"""Regression controls for payload integrity and process-budget enforcement.

All data are generated fixtures. No external executable is started.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import reproduce_integration as runner

FIXTURE_PARENT = None


class ReleaseControls(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=FIXTURE_PARENT)
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.root = self.base / 'release'
        self.root.mkdir()
        self.raw = b'Generated regression fixture, not research data.\n'
        (self.root / 'payload.txt').write_bytes(self.raw)
        self.row = {'path': 'payload.txt', 'bytes': len(self.raw)}
        self.write_metadata([self.row])

    def write_metadata(self, rows):
        raw = (json.dumps({'files': rows}) + '\n').encode()
        (self.root / 'manifest.json').write_bytes(raw)

    def audit(self, clean=False):
        with patch.object(runner, 'HERE', self.root):
            return runner.verify_inventory(clean=clean)

    def test_valid_payload(self):
        self.assertEqual(self.audit(clean=True), 1)

    def test_changed_payload(self):
        (self.root / 'payload.txt').write_text('Changed fixture')
        with self.assertRaises(ValueError):
            self.audit()

    def test_missing_payload(self):
        (self.root / 'payload.txt').unlink()
        with self.assertRaises(ValueError):
            self.audit()

    def test_truncated_inventory(self):
        (self.root / 'manifest.json').write_text('{"files": []}')
        with self.assertRaises(ValueError):
            self.audit()

    def test_unlisted_file(self):
        (self.root / 'unexpected.txt').write_text('Generated extra fixture')
        with self.assertRaises(ValueError):
            self.audit()

    def test_invalid_inventory(self):
        (self.root / 'manifest.json').write_text('invalid')
        with self.assertRaises(ValueError):
            self.audit()

    def test_wrong_byte_count(self):
        self.write_metadata([dict(self.row, bytes=1)])
        with self.assertRaises(ValueError):
            self.audit()

    def test_duplicate_entries(self):
        self.write_metadata([self.row, self.row])
        with self.assertRaises(ValueError):
            self.audit()

    def test_path_traversal(self):
        (self.base / 'outside.txt').write_bytes(self.raw)
        self.write_metadata([dict(self.row, path='../outside.txt')])
        with self.assertRaises(ValueError):
            self.audit()

    def test_symlinked_payload(self):
        outside = self.base / 'outside.txt'
        outside.write_bytes(self.raw)
        (self.root / 'payload.txt').unlink()
        (self.root / 'payload.txt').symlink_to(outside)
        with self.assertRaises(ValueError):
            self.audit()

    def test_symlinked_parent(self):
        (self.root / 'nested').symlink_to(self.base, target_is_directory=True)
        self.write_metadata([self.row, dict(self.row, path='nested/release/payload.txt')])
        with self.assertRaises(ValueError):
            self.audit()

    def test_runtime_cache_requires_working_audit(self):
        cache = self.root / 'runs'
        cache.mkdir()
        (cache / 'local.txt').write_text('Synthetic runtime output')
        self.assertEqual(self.audit(), 1)
        with self.assertRaises(ValueError):
            self.audit(clean=True)

    def test_empty_runtime_directory_rejected_when_clean(self):
        (self.root / 'runs').mkdir()
        with self.assertRaises(ValueError):
            self.audit(clean=True)

    def test_unreadable_subtree_rejected(self):
        def unreadable_walk(*args, **kwargs):
            error = PermissionError('Synthetic unreadable directory')
            if kwargs.get('onerror'):
                kwargs['onerror'](error)
            return iter(())
        with patch.object(runner.os, 'walk', unreadable_walk):
            with self.assertRaises(PermissionError):
                self.audit(clean=True)

    def fake_run(self, seconds, encoding='pdf', code=0, terminated=True,
                 timeout=False, first_wide=False):
        dummy = self.base / 'not_executable.txt'
        dummy.write_text('Never executed.')
        output = self.base / 'fake_process'
        output.mkdir()
        invocations = []

        class FakeProcess:
            pid = 999999999

            def __init__(self, command, **kwargs):
                invocations.append(command)
                self.waits = 0
                wide = first_wide and len(invocations) == 1
                kwargs['stdout'].write('Normalizing constant: 1 <= Z <= ' + ('2' if wide else '1') + '\n')

            def wait(self, timeout=None):
                self.waits += 1
                if timeout_flag and self.waits == 1:
                    raise subprocess.TimeoutExpired('synthetic process', timeout)
                return code

            def poll(self):
                return code if terminated else None

        timeout_flag = timeout
        clock = iter(seconds)
        args = SimpleNamespace(dotnet=dummy, gubpi_dll=dummy, case='HH0H', encoding=encoding)
        with (patch.object(runner.subprocess, 'Popen', FakeProcess),
              patch.object(runner.time, 'perf_counter', lambda: next(clock)),
              patch.object(runner.os, 'killpg')):
            return runner.external_gubpi(args, output)

    def test_on_time_narrow_output(self):
        self.assertTrue(self.fake_run([0, 59.9])['target_met'])

    def test_exact_deadline_is_inclusive(self):
        self.assertTrue(self.fake_run([0, 60])['target_met'])

    def test_late_narrow_output_rejected(self):
        result = self.fake_run([0, 61])
        self.assertFalse(result['target_met'])
        self.assertTrue(result['attempts'][0]['reported_target_met'])
        self.assertFalse(result['attempts'][0]['within_time_budget'])
        self.assertEqual(result['attempts'][0]['status'], 'TIME_BUDGET_EXCEEDED')

    def test_cumulative_budget_across_settings(self):
        result = self.fake_run([0, 35, 35, 61], first_wide=True)
        self.assertEqual(len(result['attempts']), 2)
        self.assertFalse(result['target_met'])
        self.assertEqual(result['attempts'][-1]['cumulative_process_seconds'], 61)

    def test_nonzero_exit_rejected(self):
        self.assertFalse(self.fake_run([0, 1], code=1)['target_met'])

    def test_unconfirmed_termination_rejected(self):
        self.assertFalse(self.fake_run([0, 1], terminated=False)['target_met'])

    def test_timeout_rejected(self):
        self.assertFalse(self.fake_run([0, 61], timeout=True)['target_met'])

    def test_gaussian_proposal_path(self):
        self.assertTrue(self.fake_run([0, 1], encoding='gaussian')['target_met'])


def main():
    global FIXTURE_PARENT
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    FIXTURE_PARENT = args.output.resolve()
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(ReleaseControls))
    record = {'checked_at': datetime.now(timezone.utc).isoformat(),
              'tests_run': result.testsRun, 'failures': len(result.failures),
              'errors': len(result.errors), 'passed': result.wasSuccessful(),
              'scope': 'Generated integrity fixtures and simulated process/clock; no external inference or scientific experiment'}
    (args.output / 'integrity_controls.json').write_text(json.dumps(record, indent=2) + '\n')
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    raise SystemExit(main())
