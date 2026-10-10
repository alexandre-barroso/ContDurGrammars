# Recorded results

All inputs for result analysis are local. Exact rational strings are authoritative
for interval endpoints; decimal table endpoints are rounded outward.

| Location | Contents |
| --- | --- |
| `results/moment_comparison.json` | Two- and three-moment normalizer trees, meshes, endpoints, timing, and memory |
| `results/normalizers/three_moments/` | Complete three-moment case records and suite summary |
| `results/normalizers/range_scalar/` | Scalar factor-range normalizer attempts, including timeouts |
| `results/normalizers/range_grid/` | Grid-recurrence range normalizer attempts, including timeouts |
| `results/queries/` | Seven enclosing query calculations with all attempted meshes |
| `results/point_queries/` | Eleven ordinary FFT meshes for each query |
| `results/gubpi/` | Probabilistic-program inputs, encoding checks, outputs, and attempt summaries |
| `results/derived/` | Exact-rational analysis of the stored query and normalizer intervals |
| `results/tables/` | Five generated LaTeX table bodies and refinement coordinates |

`results.json` inside an experiment directory summarizes the corresponding
case files. A case's `attempts` array retains completed meshes and their outcomes;
unfinished work may also be identified by an in-progress field. Status flags
must be read alongside elapsed time, process exit, and timeout information.

## Normalizer comparison

At log-width target `0.001`, both moment methods meet the requested accuracy on
all four synthetic models. The two-moment stopping meshes are four times finer.
All nineteen common-mesh three-moment intervals lie inside their two-moment
counterparts in these records. This is an observed comparison, not a general
interval-nesting theorem. Both constructions admit quadratic exact-real
mass-width analysis.

## Query records

Each query case contains its exact tree, events, cube, precision, and attempt
array. An attempt contains four `integrals`, five `probabilities`, upper bounds
on widths, and its stopping flag. The final five widths are at most `0.002`
for all seven models. All final ordinary FFT estimates lie in their respective
stored intervals.

For HH0H, the recorded enclosure of `Delta` is
`[-0.0128211, -0.0115736]` after outward decimal rounding. Removing the group
factor makes the density factorize and gives `Delta=0` exactly. The enclosing
result quantifies a prediction of this fixed model, not a measured speech effect.

## Controlled model contrasts

`derive_query_contrasts.py` checks that B16, L16, and S16 have identical ordered
unary potentials and events. B16 and S16 share the same hierarchy and targets;
all fifteen internal nodes have exactly fourfold directional weights in S16.
The script subtracts the stored exact endpoints to obtain:

| Difference | Outward-rounded enclosure |
| --- | --- |
| `Delta(L16) - Delta(B16)` | `[0.0015193, 0.0037408]` |
| `Delta(S16) - Delta(B16)` | `[0.0012885, 0.0032420]` |

Both intervals exclude zero. Conditional on the input enclosures, numerical
error within their bounds cannot reverse either ordering. The point estimates
give the same orderings. This does not imply a general monotonic relationship
between stiffness and dependence or an empirical difference between speakers.

The analysis checks 29 query attempts, 116 rational ratio hulls, 29 difference
hulls, 35 final width targets, 35 final FFT inclusions, and 19 common-mesh
normalizer containments. It does not recompute the enclosing integrals.

## Recompute derived outputs

```sh
python3 derive_query_contrasts.py --output runs/contrasts
python3 make_integration_artifacts.py --output runs/export
```

The exporters use only standard-library arithmetic unless `--plots` is selected.
Output reports record their input and output hashes. Timing and memory columns
remain descriptive recorded measurements; table generation does not retime jobs.
