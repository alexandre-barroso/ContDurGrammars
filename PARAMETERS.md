# Parameters and model construction

`published_timing_parameters.json` supplies all numerical parameters used by
`integration_bounds.hayes_models()`. Computation reads this local file only.
Durations and thresholds are in milliseconds; divided directional weights
are in inverse squared milliseconds. No observed duration tokens are used.

## Attribution and interpretation

The parameter values are attributed to Bruce Hayes and Russell G. Schuh,
*Metrical structure and sung rhythm of the Hausa rajaz* (2019), duration
analysis in Section 8. The local metadata retains the parameter cell locators:
weights `AU3:BB3`, mora target `BC3`, syllable target `AR3`, grid target `AS3`,
and group target `AT3` on the `Main` sheet.

The numerical model here is continuous with product Lebesgue measure. The
source analysis used a finite duration grid; changing the reference space
changes the model, so its fitted probabilities are not reproduced by assumption.
The directional names in the printed definitions and parameter formulas differ;
the local fields label the reconstructed directions explicitly as `above` and
`below`.

## Exact encoding

Rational strings encode the binary64 parameter values exactly. Multiply each
stored raw weight by `1/10000` before the count scaling below. Displayed decimal
values in generated tables are descriptive rounded summaries, not alternative
inputs. All computations use the rational strings.

For leaf `i`, with mora count `n_i` and grid-column count `c_i`, the cost is

```text
psi(n_i*t_mora, a_mora/n_i, b_mora/n_i; x_i)
+ psi(t_syllable, a_syllable, b_syllable; x_i)
+ psi(c_i*t_grid, a_grid/c_i, b_grid/c_i; x_i).
```

Here the `a` and `b` weights have already been divided by 10000, and `psi` is
defined in [METHOD.md](METHOD.md). The group cost uses the hemimetron target
and its divided directional weights.

## Three small models

| Label | Active coordinates | Mora counts | Grid-column counts |
| --- | --- | --- | --- |
| HH0H | 2 | `(2,2)` | `(2,2)` |
| HHLH | 3 | `(2,2,1)` | `(1,2,1)` |
| HLLH | 3 | `(2,1,1)` | `(2,1,1)` |

The constructor also provides `LHLH` and `LLHH` as coordinate permutations of
HLLH. The query suite uses the three representative labels above.

Mora counts and grid-column counts are separate inputs. The labels retain the
duration-analysis convention that excludes the long final syllable from the
measured group. The empty metrical position in HH0H is not a zero-valued random
variable. The hemimetron label is a measurement convention, not a claim of a
phonological constituent.

## Synthetic hierarchies

`integration_scaling.build_case()` constructs balanced eight- and sixteen-leaf
models, a left-branching sixteen-leaf model, and a stiffer balanced sixteen-leaf
model. Leaves alternate the heavy and light unary costs of HLLH. For a node
containing `M` moras, the target is `M*t_hemimetron/4` and each directional weight
is `4/M` times the divided root weight. The stiff model multiplies all internal
weights by four. These are numerical test structures, not attested prosodic trees.
