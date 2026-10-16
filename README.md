# Certified syllable-duration inference for computational prosody

This repository implements certified inference for computational prosody. Given a supplied duration grammar and timing budget, it computes syllable-duration probabilities with bounded numerical error and verifies timing decisions against explicit probability requirements. It contains the speech annotations, fitted models, inference routines, experiment inputs, numerical results and evaluation scripts needed to run the experiments independently.

The method encloses normalizers, event probabilities and conditional probabilities. Three grid-centred moment arrays carry integral information through a hierarchy of nonnegative durations. The mathematical contracts, limits and formal obligations are documented alongside the code.

## Start with the speech experiments

The [speech application](speech/README.md) includes all annotation files, licenses, fixed split assignments, fitting scripts, model inputs, probability certificates and evaluation results. The [deadline planner](speech/deadlines/README.md) compares interval-certified timing decisions with ordinary numerical approximation. The annotation data needed for these experiments are included in `speech/data/`.

On 836 held-out words from two CMU ARCTIC voices, the shared timing grammar achieves a Brier score of 0.0847, compared with 0.1050 for independent durations and 0.1039 for proportional allocation. Every probability enclosure has width at most 0.002. For the 95% model-probability deadline task, certification resolves 787 feasible and 49 infeasible inputs. Ordinary 10 ms integration makes 18 decision errors, eight of which survive a successive-mesh stability check. Fine ordinary integration agrees with every certified decision; certificates establish feasibility and minimality under the supplied model. These experiments measure timing prediction and numerical decision reliability, not perceptual synthesis quality.

## Repository guide

The components below support the corresponding parts of *Certified syllable-duration inference for computational prosody*. All commands operate on files within this repository; a manuscript checkout is not needed.

| Manuscript component | Repository contents |
| --- | --- |
| Duration grammar and inference contract (Sections 2 and 4) | [METHOD.md](METHOD.md), `integration_bounds.py`, `integration_moments.py` and `integration_queries.py` |
| Recorded-speech prediction (Section 3) | `speech/data/`, `speech/split_assignments.json`, extraction and fitting scripts, exact model inputs, certificates and scoring outputs in `speech/` |
| Certified speech-planning decisions (Section 3.5) | Ordinary controls, interval certificates, complete decision records and independent checks in `speech/deadlines/` |
| Numerical evaluation (Section 5.1) | Normalizer, conditional-query, ordinary FFT and probabilistic-program results in `results/`; experiment drivers and controls in the repository root |
| Reproducibility (Section 5.2) | Dependency files, `reproduce_integration.py`, validation scripts, inventory and unit tests |
| Enclosure proofs and moment comparison (Appendices A and B) | Four `Moment*.lean` modules, their assumptions in [FORMAL_MAP.md](FORMAL_MAP.md), and the two-moment implementation and controls |
| Hausa timing example (Appendix C) | `published_timing_parameters.json`, [PARAMETERS.md](PARAMETERS.md), and model construction in `integration_bounds.py` |

Results are stored as ordinary files in `results/`, `speech/certificates/` and `speech/deadlines/`. Read or run them directly after cloning the repository.

## Python environment

Set up Python 3.12:

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-integration.txt
```

Inspect the supplied speech results:

```sh
.venv/bin/python speech/validate.py
.venv/bin/python speech/evaluate.py
.venv/bin/python speech/deadlines/evaluate.py
.venv/bin/python speech/deadlines/independent_checks.py
```

Regenerate the speech calculations:

```sh
.venv/bin/python speech/prepare_data.py
.venv/bin/python speech/fit_models.py
.venv/bin/python speech/certify_predictions.py --workers 4
.venv/bin/python speech/deadlines/run.py --ordinary
.venv/bin/python speech/deadlines/run.py --certify --workers 4
```

The certificate drivers reuse matching completed records after comparing their explicit model and event inputs. Use empty certificate directories for fresh timing measurements. Run the evaluation commands again after regenerating calculations. The data protocol and reproduction details are in the speech directory.

## Inspect the general numerical results

These commands use Python's standard library and do not perform integration or compile proofs:

```sh
python3 reproduce_integration.py --mode audit --output runs/audit
python3 reproduce_integration.py --mode analyze --output runs/analysis
python3 reproduce_integration.py --mode artifacts --output runs/tables
python3 integrity_checks.py --output runs/integrity
python3 -m unittest test_derived_analysis test_artifacts test_repository -v
```

The inventory check verifies relative file paths and byte counts. Record analysis checks explicit model identities, exact-rational bounds and model contrasts. Artifact generation creates the LaTeX tables and refinement data. Commands taking `--output` require a destination that does not yet exist; generated outputs belong under `runs/`.

## Numerical controls and calculations

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  reproduce_integration.py --mode check --output runs/numerical
.venv/bin/python verify_integration.py --mode replay-final --output runs/final_meshes
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python \
  reproduce_integration.py --mode experiments --output runs/experiments
```

Controls cover Gaussian leaf moments, signed convolution, derivative bounds, nested factors, clipped and empty events, exterior events and the two-moment normalizer. Experiments include four three-moment normalizers, four two-moment normalizers, seven conditional queries and seven ordinary FFT comparisons. All attempts and unresolved outcomes are retained. Missing a time budget is a failed accuracy request, not a success. Do not disable assertion checks with python -O.

The numerical environment uses Python 3.12, python-flint 0.8.0, NumPy 2.4.4, SciPy 1.17.1 and mpmath 1.3.0. The normalizer and conditional-query timings were recorded with Python 3.12.3; speech computations used Python 3.12.14.

## Figures and formal modules

```sh
.venv/bin/python -m pip install -r requirements-figures.txt
.venv/bin/python reproduce_integration.py --mode figures --output runs/figures
```

To check the four Lean modules, use Lean/Lake 4.34.0 with a compatible mathlib installation at .lake/packages/mathlib:

```sh
python3 reproduce_integration.py --mode proofs \
  --mathlib .lake/packages/mathlib --output runs/proofs
```

The command checks the Lean toolchain, builds the local modules with --trust=0 and checks their printed axiom dependencies. These obligations do not formalize the Python implementation or Arb.

## Documentation

| File | Contents |
| --- | --- |
| [speech/README.md](speech/README.md) | Recorded-speech inputs, models, predictions and reproduction |
| [speech/deadlines/README.md](speech/deadlines/README.md) | Certified planning decisions and numerical comparisons |
| [METHOD.md](METHOD.md) | Model, message geometry, nonlinear updates and probability contracts |
| [PARAMETERS.md](PARAMETERS.md) | Units, exact parameters, source attribution and model construction |
| [EXPERIMENTS.md](EXPERIMENTS.md) | Cases, events, tolerances and budgets |
| [RESULTS.md](RESULTS.md) | Numerical records and model contrasts |
| [TWO_MOMENT.md](TWO_MOMENT.md) | Support-only comparator and endpoint replay |
| [FORMAL_MAP.md](FORMAL_MAP.md) | Eighteen local declarations and their assumptions |
| [GUBPI.md](GUBPI.md) | Probabilistic-program baseline and retained attempts |
| [VALIDATION.md](VALIDATION.md) | File, record, numerical and formal checks |
