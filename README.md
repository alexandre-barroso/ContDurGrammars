# Repository for verification

Certified probability queries for continuous duration grammars.

This project computes enclosing intervals for normalizers, event probabilities,
and conditional probabilities in continuous duration models. Costs act on
individual nonnegative durations and on nested or disjoint duration sums.
Three grid-centred moment arrays carry integral information through the hierarchy.

The project contains the numerical implementation, exact model parameters,
recorded calculations, result-analysis tools, finite numerical controls, and
four Lean modules for local analytic and interval-arithmetic obligations.

## Inspect the recorded results

These commands use only Python's standard library. Run them from this directory:

```sh
python3 reproduce_integration.py --mode audit --output runs/audit
python3 reproduce_integration.py --mode analyze --output runs/analysis
python3 reproduce_integration.py --mode artifacts --output runs/tables
python3 integrity_checks.py --output runs/integrity
python3 -m unittest test_derived_analysis test_artifacts test_repository -v
```

`audit` checks the file inventory, sizes, SHA-256 values, and checksum index.
`analyze` checks exact-rational probability bounds and derives model contrasts.
`artifacts` creates five LaTeX table bodies and exact refinement data.
These checks do not perform numerical integration or compile proofs.
Each command that takes `--output` requires a new destination. Generated
outputs belong under `runs/`; the recorded inputs under `results/` are read-only.

## Run numerical controls

The numerical environment uses Python 3.12, python-flint 0.8.0, NumPy 2.4.4,
SciPy 1.17.1, and mpmath 1.3.0. The recorded calculations used Python 3.12.3.

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-integration.txt
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  reproduce_integration.py --mode check --output runs/numerical
```

The controls cover Gaussian leaf moments, signed convolution, derivative bounds,
nested factors, clipped and empty events, exterior events, and the two-moment
normalizer. They also replay the first stored two-moment mesh for each of four
synthetic models. Do not run assertion-based controls with `python -O`.

For the four final two-moment meshes, without a timing deadline:

```sh
.venv/bin/python verify_integration.py --mode replay-final --output runs/final_meshes
```

## Run calculations

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  reproduce_integration.py --mode experiments --output runs/experiments
```

This runs four three-moment normalizers, four two-moment normalizers, seven
conditional queries, and seven ordinary FFT comparisons. Jobs run sequentially.
All attempted meshes and unresolved outcomes are retained. A missed time budget
is reported as a failed accuracy request, not silently counted as a success.

A single conditional-query calculation is also available:

```sh
.venv/bin/python integration_queries.py --worker HH0H --output runs/HH0H.json
```

[EXPERIMENTS.md](EXPERIMENTS.md) defines the cases, tolerances, and budgets.
[RESULTS.md](RESULTS.md) explains every family of stored records.

## Generate refinement plots

```sh
.venv/bin/python -m pip install -r requirements-figures.txt
.venv/bin/python reproduce_integration.py --mode figures --output runs/figures
```

This adds one PNG plot per synthetic model to the table and data outputs.
No special font or document template is required.

## Check the formal modules

Use Lean/Lake 4.34.0 and mathlib at commit
`5ed2965256430c3649e86755f9576b54eca72435`. With that dependency installed at
`.lake/packages/mathlib`, run:

```sh
python3 reproduce_integration.py --mode proofs \
  --mathlib .lake/packages/mathlib --output runs/proofs
```

The command checks the toolchain identity, builds the four local modules with
`--trust=0`, and checks the printed axiom dependencies. Proof-only checks do not
require the numerical Python packages. The precise scope is in
[FORMAL_MAP.md](FORMAL_MAP.md).

## Documentation

| File | Contents |
| --- | --- |
| [METHOD.md](METHOD.md) | Model, message geometry, nonlinear updates, and probability contracts |
| [PARAMETERS.md](PARAMETERS.md) | Units, exact parameter encoding, source attribution, and model construction |
| [EXPERIMENTS.md](EXPERIMENTS.md) | Cases, event definitions, stopping rules, and execution commands |
| [RESULTS.md](RESULTS.md) | Record schema, outcomes, and exact-rational model contrasts |
| [TWO_MOMENT.md](TWO_MOMENT.md) | Support-only comparator and endpoint replay |
| [FORMAL_MAP.md](FORMAL_MAP.md) | Eighteen local declarations and their assumptions |
| [GUBPI.md](GUBPI.md) | Optional probabilistic-program baseline and retained attempts |
| [VALIDATION.md](VALIDATION.md) | Integrity, record analysis, numerical tests, and formal-check boundaries |
