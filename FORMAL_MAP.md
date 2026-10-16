# Local formal obligations

The four Lean modules contain eighteen checked-declaration targets. Their
analytic assumptions and primitive interval contracts are explicit; the modules
do not verify the entire numerical executable.

## MomentEnclosure

- `integral_error_bound`
- `mass_propagation`
- `first_moment_propagation`
- `second_moment_propagation`
- `first_support_bound`
- `second_support_bounds`
- `factor_range_mass`
- `factor_range_second`

These declarations express integral residual bounds and their consequences for
factor updates, support intersections, and positive-moment range bounds. There
is no independent positive range identity for a signed first moment.

## MomentConvolution

- `first_moment_product`
- `second_moment_product`
- `tail_dominated_by_envelope`
- `full_integral_enclosure`

These cover product-measure moment identities and extension of cube integrals
by an exterior envelope, under the stated integrability and domination hypotheses.

## MomentCircuit

- `mul_contains`
- `factor_moments_contains`
- `trace_contains`
- `ratio_enclosure`

These cover finite signed interval calculations and positive-denominator ratios
under their primitive and analytic contracts.

## MomentBridge

- `factor_contract_from_integrals`
- `enclosed_factor_from_integrals`

These connect integral hypotheses with the finite enclosing factor calculation.
The concrete aggregate-set geometry is described in [METHOD.md](METHOD.md).

## Build and scope

`verify_integration.py --mode proofs` requires Lean 4.34.0 and a mathlib
installation using the same toolchain, builds each local module with
`--trust=0`, and checks the `#print axioms` output. Only `propext`,
`Classical.choice`, and `Quot.sound` are accepted. Missing targets, unsupported
proof shortcuts, failed builds, or unexpected axioms cause a nonzero exit.
Logs and declaration checks are stored in the selected output directory.

A successful build establishes the displayed local statements under their
hypotheses. It is not a single theorem about the whole-tree executable, the
width recurrence, Python, Arb, or an independently implemented proof kernel.
