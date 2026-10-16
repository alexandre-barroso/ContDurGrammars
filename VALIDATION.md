# Validation

Validation has separate scopes. Passing one level does not imply that the next
has run or passed.

## File integrity

```sh
python3 reproduce_integration.py --mode audit --output runs/audit_check
python3 integrity_checks.py --output runs/inventory_tests
```

manifest.json lists the project files by relative path and byte count. The audit checks completeness, file sizes and safe paths. Runtime directories such as .git, .venv, .lake, __pycache__ and runs are excluded in a working checkout; --clean-tree requires a clean distribution. The inventory check does not prove numerical correctness or detect every possible content change.

The integrity command runs fourteen filesystem fixtures. The full integrity_controls.py suite also exercises process-budget handling. The speech application's validation compares split assignments, explicit model inputs and exact probability endpoints; its deadline evaluation checks every feasibility and minimality inequality.

## Stored-record analysis and export

```sh
python3 -m unittest test_derived_analysis test_artifacts test_repository -v
python3 reproduce_integration.py --mode analyze --output runs/record_checks
python3 make_integration_artifacts.py --output runs/table_checks
```

Tests cover rational endpoint rounding, subtraction, denominator contracts,
controlled model identities, altered-input rejection, table regeneration,
local links, module syntax, and local file-path consistency. The full record
analysis checks ratio hulls, difference hulls, final width flags, final FFT
membership, and common-mesh normalizer containment. It assumes the stored
integral enclosures are valid; it does not recompute those integrals.

## Numerical execution

`reproduce_integration.py --mode check` runs finite primitive, event, and
nested-model controls plus first-mesh two-moment replay in the pinned numerical
environment. `verify_integration.py --mode replay-final` recomputes four final
two-moment meshes. Assertion controls must not be disabled with `-O`.

Finite numerical reference checks are not proofs of their reference quadrature.
Endpoint replay compares exact rational values. It does not require reproducing
recorded timing and memory observations.

## Formal execution

`reproduce_integration.py --mode proofs --mathlib .lake/packages/mathlib` builds
the four local modules and checks eighteen declaration targets. The
`--output` destination must not already exist. [FORMAL_MAP.md](FORMAL_MAP.md) gives the assumptions and
axiom policy. The mathematical contracts do not formalize Python or Arb.

## Reports

Each command records its actual scope and status in its selected output
directory and returns nonzero on failure. `validation_checks.json` describes
the checks run on the supplied project files, including numerical and formal
checks that were not executed. Numerical, formal, and timed experiments must
not be inferred to have run from a successful integrity or stored-record check.
