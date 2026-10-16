# Probabilistic-program baseline

`gubpi_encoding.py` constructs SPCF input programs and parses the numerical
outputs for a broader probabilistic-inference baseline. All recorded programs,
encoding checks, stdout logs, and attempt summaries are in `results/gubpi/`.
The enclosing duration integrator does not depend on this baseline.

## Encodings and records

| Directory | Contents |
| --- | --- |
| `results/gubpi/uniform/` | Box and polytope arithmetic encodings with uniform proposals |
| `results/gubpi/gaussian_boxes/` | Arithmetic box encodings with Gaussian proposals |
| `results/gubpi/gaussian_pdf/` | Gaussian-PDF branch encodings with linear guards |

The uniform and Gaussian-box directories group execution batches by their
recorded start time on 26 September 2026. Every batch retains its attempts,
including parser errors, nonfinite outputs, and timeouts. Batch labels identify
recorded executions, not interchangeable algorithms.

The Gaussian-PDF encoding produces finite outputs for the three small models.
Their best completed log-widths are approximately `1.4081`, `3.7800`, and
`4.0140`, above `0.001`. The eight-leaf case has no finite bound before timeout;
termination was not confirmed in its bounded post-kill wait. That case is not a
clean inference-time measurement.

Gaussian constants and parsed endpoints are binary64. Sampled high-precision
factor checks are not a proof of uniform equivalence to the exact-rational
continuous model. The stored .NET runtime is 7.0.20, x86-64 on macOS arm64
translation; the moment calculations use native Arb. Neither a general speed
claim nor a model-equivalence theorem follows from these comparisons.

## Optional execution

GuBPI and .NET are optional installed dependencies, not required to inspect the records or run duration
queries. Put the runtime and compiled application in `.deps/dotnet` and
`.deps/gubpi`, respectively, to use the example below.

```sh
.venv/bin/python reproduce_integration.py --mode gubpi-inputs --output runs/gubpi_inputs
.venv/bin/python reproduce_integration.py --mode gubpi --case HH0H --encoding pdf \
  --dotnet .deps/dotnet/dotnet --gubpi-dll .deps/gubpi/GuBPI.dll --output runs/gubpi_HH0H
```

Input generation writes explicit programs and finite encoding checks. Inference
adds process execution with a cumulative 60-second budget and bounded teardown.
Failed targets return nonzero and retain output. Executable names and sizes alone
do not certify a complete transitive build identity.
