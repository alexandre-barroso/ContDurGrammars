# Method and inference contract

## Model

Durations are distinct coordinates `x_i >= 0`, measured in milliseconds. A finite
rooted tree specifies nested or disjoint groups. For an internal node `v`, `S_v`
is the sum of its descendant durations. The unnormalized density is

```text
f(x) = product_i exp(-V_i(x_i)) * product_v exp(-V_v(S_v(x))).
```

Every cost is a finite sum of squared hinges:

```text
psi(t, a, b; s) = a * max(s-t, 0)^2 + b * max(t-s, 0)^2,   a,b >= 0.
```

The reference measure is product Lebesgue measure. All sum loadings are one.
General overlapping groups cannot be inserted without changing the factorization.
Each unary factor `q_i = exp(-V_i)` must be integrable on `[0, infinity)`; a
strictly positive positive-tail quadratic coefficient suffices. Since each
internal factor lies in `(0,1]`, the product of the unary factors is an integrable
upper envelope. Joint propriety alone is not enough for this envelope construction.

Targets, weights, event boundaries, and mesh coordinates denote exact rationals.
A supplied binary floating-point number denotes its exact binary value, not an
inferred decimal intention. [PARAMETERS.md](PARAMETERS.md) defines the examples.

## Events and requested outputs

A supported event is a product of leaf intervals. Intervals may be unrestricted,
empty, or partially aligned with the mesh. For product events `A` and `B`, four
integrals over the whole domain, `A`, `B`, and their intersection give

```text
P(A), P(B), P(A intersection B), P(A | B), Delta = P(A | B) - P(A).
```

A conditional probability requires positive conditioning mass. A computation
returns enclosing intervals or an unresolved accuracy request. Failed intersections,
nonfinite values, and nonpositive denominator lower bounds are not certificates.
Finite unions can be split into separate requests; there is no compressed
representation for arbitrary measurable events.

## Grid geometry

Choose `L > 0`, `h > 0`, and integer `N = L/h`. Leaf cells retain their grid
indices and centres, including when an event clips part of a cell. Empty cells
keep zero entries at their indices. For a node with `m` leaves, array index `k`
represents the disjoint union of event-clipped leaf boxes whose indices sum to `k`.
Its centre and support radius are

```text
c(k) = (k + m/2) * h,      r = m*h/2.
```

There are `m*(N-1)+1` array positions rather than `N^m` Cartesian boxes. The
aggregate set is not the entire inverse image of `[c-r,c+r]` under the duration
sum: using that larger strip as the integration region would overlap and
miscount mass. The larger interval is used only to bound factors and derivatives.
Recentring clipped cells would break the additive geometry.

## Three-component messages

At an aggregate set, write `d = S-c`. The message stores intervals enclosing

```text
M0 = integral 1 dnu,   M1 = integral d dnu,   M2 = integral d^2 dnu,
```

where `nu` is the positive measure including all subtree factors already applied.
These are unnormalized integrals, not a fitted distribution for the sum. `M1`
is signed. Each leaf splits retained cell pieces at every cost knot and computes
the three moments with Gaussian/error-function primitives; constant branches
use polynomial integration.

For two disjoint child messages `A` and `B`, centres and offsets add. With `*`
denoting convolution of array indices, the pre-factor moments are

```text
M0 = A0*B0
M1 = A1*B0 + A0*B1
M2 = A2*B0 + 2*A1*B1 + A0*B2.
```

Signed interval arithmetic is required, especially for `2*A1*B1`. Higher arity
uses a fixed binary grouping. Different groupings can change overestimation
without changing these exact identities.

## Nonlinear factor transfer

For factor `g(S)`, let `g0=g(c)`, `g1=g'(c)`, `|g'| <= K1`, and let `g'` be
`K2`-Lipschitz on `[c-r,c+r]`. Let `U2` be a nonnegative upper bound on the true
second moment. The post-factor moments obey

```text
M0' in g0*M0 + g1*M1 + [-K2*U2/2, K2*U2/2]
M1' in g0*M1 + g1*M2 + [-K2*r*U2/2, K2*r*U2/2]
M2' in g0*M2         + [-K1*r*U2,   K1*r*U2].
```

The first two remainders follow from the first-order Taylor bound; the last
uses the first-derivative bound. On this support, `|d|^3 <= r*d^2`, so a third
moment need not be stored. Remainders use an upper bound on `M2`, never its
midpoint.

Valid support intersections impose `M0 >= 0`, `|M1| <= r*M0`, and
`0 <= M2 <= r^2*M0`. A factor range `ell <= g <= u` also bounds the positive
moments by `ell*M0 <= M0' <= u*M0` and `ell*M2 <= M2' <= u*M2`.
There is no independent positive range rule for signed `M1`. Intersections
must retain the true moments; an empty intersection is a failure.

On a quadratic branch, `V(s)=alpha*s^2-2*beta*s+gamma`, so

```text
g'(s)  = -(2*alpha*s-2*beta)*g(s)
g''(s) = ((2*alpha*s-2*beta)^2-2*alpha)*g(s).
```

Branchwise bounds extend across knots because the squared-hinge cost and its
first derivative are continuous. A second derivative at the knot itself is
not required. These operations are in `integration_moments.py`.

## Full-domain and conditional enclosures

Summing root mass intervals encloses the event integral inside the cube. An
upper bound on the omitted mass is obtained from the product envelope:

```text
T >= product_i integral_0^infinity q_i - product_i integral_0^L q_i.
```

If the cube event mass is in `[l,u]`, the full event mass is in
`[max(0,l), u+T]`. For nested events with numerator mass `[lN,uN]` and
denominator mass `[lD,uD]`, where `lN>=0` and `lD>0`,

```text
P(N | D) in [lN/uD, uN/lD] intersection [0,1].
```

Endpoint operations must themselves be enclosing. No independence of numerical
errors is assumed. Interval subtraction bounds `Delta` even though its terms
share integrals.

## Accuracy and implementation assumptions

A fixed-tree exact-real width analysis and a finite-precision stopping decision
are different matters. At fixed cube, with bounded factor derivatives and exact
leaf moments, the mass-width analysis permits quadratic refinement. Both the
two- and three-moment constructions admit that order; retaining `M2` can sharpen
remainders but is not uniquely responsible for the order.

Actual computations test their returned widths. Mesh refinement alone cannot
remove a fixed excessive exterior allowance. Rounding may require increased
precision. Conditioning on very small mass can demand impractical relative
accuracy. Array length alone is not a complete runtime or memory bound.

Containment assumes correct Gaussian primitives, signed ball and polynomial
operations, endpoint extraction, cache identities, and runtime behavior. The
local formal obligations are listed in [FORMAL_MAP.md](FORMAL_MAP.md); the
finite numerical controls are described in [VALIDATION.md](VALIDATION.md).
