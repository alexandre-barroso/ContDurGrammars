"""Run the fifteen file-inventory fixtures using the Python standard library.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import unittest
import integrity_controls

TESTS = (
    "test_valid_payload", "test_changed_payload", "test_missing_payload",
    "test_truncated_manifest_without_updated_sums",
    "test_truncated_manifest_with_updated_sums", "test_unlisted_file",
    "test_corrupted_checksum_file", "test_wrong_byte_count", "test_duplicate_entries",
    "test_path_traversal", "test_symlinked_payload", "test_symlinked_parent",
    "test_runtime_cache_requires_working_audit",
    "test_empty_runtime_directory_rejected_when_clean", "test_unreadable_subtree_rejected",
)

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    integrity_controls.FIXTURE_PARENT = args.output.resolve()
    suite = unittest.TestSuite(integrity_controls.ReleaseControls(name) for name in TESTS)
    with (args.output / "integrity_tests.log").open("w") as log:
        result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    record = {"checked_at": datetime.now(timezone.utc).isoformat(),
              "status": "PASSED_INVENTORY_FIXTURES" if result.wasSuccessful() else "FAILED",
              "tests_run": result.testsRun, "failures": len(result.failures),
              "errors": len(result.errors), "test_names": list(TESTS),
              "not_run": "Eight process-budget fixtures require numerical dependencies",
              "scope": "Generated payload-integrity fixtures only; no numerical or formal replay"}
    (args.output / "integrity_tests.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps({k: record[k] for k in ("status", "tests_run", "failures", "errors")}))
    return 0 if result.wasSuccessful() else 1

if __name__ == "__main__":
    raise SystemExit(main())
