# Calculations and stopping rules

## Cases

| Record key | Short label | Structure |
| --- | --- | --- |
| `HH0H`, `HHLH`, `HLLH` | Same | Two- or three-coordinate duration models |
| `balanced_8` | B8 | Balanced, eight leaves |
| `balanced_16` | B16 | Balanced, sixteen leaves |
| `left_branching_16` | L16 | Left branching, sixteen leaves |
| `stiff_balanced_16` | S16 | B16 hierarchy with fourfold internal weights |

[PARAMETERS.md](PARAMETERS.md) specifies every input. Each enclosing calculation
uses `[0,800]^n`, exact rational inputs, 128-bit Arb precision, and meshes
`h=20/2^j`. Normalizer drivers allow `j=0,...,14`.

## Normalizers

Two- and three-moment normalizers request log-width at most `0.001` with a
60-second end-to-end budget per model. The budget includes interpreter startup
and all completed refinements. The two-moment comparator replaces stored
second moments with their support bound; [TWO_MOMENT.md](TWO_MOMENT.md)
defines the difference.

```sh
.venv/bin/python integration_moment_scaling.py --worker balanced_8 --output runs/three_B8.json
.venv/bin/python integration_two_moment_scaling.py --worker balanced_8 --output runs/two_B8.json
```

Worker commands perform refinements and write attempt records; they do not
independently enforce the enclosing process budget. The full runner launches
workers with time limits and checks both status and elapsed time.

Range-only normalizers use either scalar factor extrema or a quadratic-grid
recurrence. They are separate algorithms, named by their operation rather than
by a development sequence:

```sh
.venv/bin/python reproduce_integration.py --mode range --output runs/range
```

Unresolved cases cause a nonzero exit while retaining completed attempts.

## Conditional queries

For each model, `A` states that its first leaf exceeds its mora target and `B`
states the same for its second leaf. Four integrals yield the five outputs
listed in [METHOD.md](METHOD.md). Every output requests width at most `0.002`.
A query has a 180-second budget for the complete four-integral case.

Ordinary FFT queries use midpoint quadrature, fractional occupancy of clipped
cells, and eleven meshes through `h=5/256`. Tiny negative convolution
coefficients are clipped and the removed mass is recorded. Successive-mesh
agreement is a diagnostic, not a full integration-error bound.

## Recorded environment and interpretation

The stored measurements use Python 3.12.3 on macOS 15.5, arm64, with
single-thread backend settings. They are individual elapsed-time observations,
not replicated hardware-independent performance estimates. `process_seconds`
and `process_elapsed_seconds` denote elapsed wall-clock time despite the field
names. A completed-mesh resident-memory sample is not an upper bound on a
subsequent unfinished refinement's peak.

Correctness replay does not require matching these timings. The ordinary FFT
and enclosing methods have different stopping rules and guarantees; dividing
their timing columns is not a matched-accuracy speed comparison. The
controlled two-/three-moment cost comparison concerns normalizers, not the
conditional-query suite.
