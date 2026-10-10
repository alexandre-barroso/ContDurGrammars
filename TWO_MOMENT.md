# Two-moment normalizer

`integration_two_moments.py` retains mass and signed first grid-centred moment,
using the same Gaussian primitives, model constructors, mesh, and precision as
the three-moment implementation. It does not calculate or store a second-moment
array or evaluate the third Gaussian primitive.

At support radius `r` and mass upper bound `U0`, its candidate factor updates are

```text
I0' = g0*I0 + g1*I1 + [-K2*r^2*U0/2, K2*r^2*U0/2]
I1' = g0*I1         + [-K1*r^2*U0,   K1*r^2*U0].
```

Valid mass-range and signed-support intersections follow. The first remainder
uses `M2 <= r^2*U0`; the second uses the bound on `g'`. These updates admit
quadratic exact-real mass-width refinement. The comparison therefore measures
enclosure sharpness and computational cost, not a uniquely higher order for
three moments.

## Stored comparison

`results/moment_comparison.json` records the exact trees and all completed
normalizer meshes for B8, B16, L16, and S16. Both methods meet log-width `0.001`
within the recorded 60-second budget. The two-moment method needs two additional
halvings in each case. There are nineteen common-mesh interval containments.
These are observations about the specified algorithms and cases, not a proof
that every two-component algorithm must be inferior.

## Controls and replay

`integration_two_moment_controls.py` checks six leaf cases, a signed merge,
twelve finite-atomic factor updates, and a nested three-coordinate calculation.
It validates the recorded trees, mesh schedule, support counts, stopping flags,
and exact endpoints. Finite reference quadrature checks are not independent
certificates of the reference quadrature itself.

```sh
.venv/bin/python verify_integration.py --mode two-moment --output runs/two_controls
.venv/bin/python verify_integration.py --mode replay-final --output runs/two_final
```

The first command includes the first mesh of each stored two-moment case.
The second replays all four final meshes. Replay compares exact endpoint strings
without imposing a timing deadline. The conditional-query implementation uses
three moments; no matched query-level two-moment cost advantage is asserted.
