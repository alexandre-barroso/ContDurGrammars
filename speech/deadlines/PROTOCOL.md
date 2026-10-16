# Certified syllable-deadline audit

This numerical evaluation uses the fitted models and test inputs from the conditional-allocation experiment. It tests model-level planning decisions, not independent predictive generalization.

## Speech planning task

A duration-planning component receives the fitted linguistic timing model and supplied word-duration budget T. It must return the earliest deadline d on a 10 ms grid, with 0 < d <= T, for which the model probability that the first syllable ends by d is at least 0.95. If no admissible deadline meets this condition, it must return INFEASIBLE. The 10 ms clock and 0.95 requirement are declared engineering choices for this demonstration, not claimed speech-industry standards or perceptual thresholds. The practical output is a scheduling decision: reserve that many frames, or request a revised word budget. No waveform or listening-test result is implied.

All 836 held-out word inputs are included. Parameters, feature fits and word budgets are specified in the parent test_requests.json and fitted_models.json files. Internal observed syllable boundaries are irrelevant to this numerical audit and never determine a deadline. All cases and the stated settings are reported.

## Ordinary controls

1. The repository's ordinary midpoint/FFT integrator, at fixed integration meshes 20, 10, 5, 2.5 and 1.25 ms. All five results are reported. The 10 ms mesh is the primary ordinary proposal. Event-cell fractions and clipping diagnostics are retained exactly as in that implementation. The integration cube is [0,1600]^2, as in the midpoint prediction task.
2. A mesh-stability heuristic: stop at the first successive pair from the same sequence whose selected deadline is unchanged and whose CDF estimates at that deadline and its predecessor change by at most 0.001. For INFEASIBLE, check the last admissible deadline. If none qualifies, report NO_STABLE_PAIR; do not silently substitute a result. This is a heuristic, not an error bound.
3. A high-accuracy ordinary control exploiting the quadratic special case: complete the square and integrate the truncated Gaussian with SciPy quad (epsabs 1e-12, epsrel 1e-10, limit 150), followed by binary search over the same deadline grid. Report its estimated quadrature errors and all results, including if it agrees throughout. It is a point calculation with error estimates, not an outward-rounded certificate.

Binary search uses the same >= 0.95 predicate and upper deadline floor(T/10)*10 for all controls. A deadline error is a different frame from the certified earliest deadline; premature and excessive allocations are counted separately. An incorrect INFEASIBLE decision or a feasible decision on a certified infeasible input is also counted.

## Certified planner

Use the repository's enclosing integration modules. Start from the ordinary 10 ms proposal. On meshes 20/2^j ms, j=0,...,10, compute enclosing CDF ratios with the unary exterior bound, 128-bit Arb and exact serialized [0,1] intersection. Accept a frame only if its CDF lower bound is >= 0.95 and its predecessor's upper bound is < 0.95. Reject all admissible frames only if the CDF upper bound at the last frame is < 0.95. Monotonicity of a CDF makes these certificates sufficient; this applies the probability enclosures described in the paper.

Move a demonstrably too-early or too-late proposal by one frame and check again. Refine ambiguous comparisons. Stop after the mesh limit or a 180-second allowance checked after completed computations; preserve unresolved cases, failures and every attempt. Report exact boundary intervals so correctness does not depend on treating the ordinary quadratic control as ground truth.

## Reporting and limits

Report every input, control, failure, abstention, selected frame, interval and observed elapsed time. Summarize disagreement counts and certified coverage; identify any illustrated disagreement by the first occurrence in the input ordering. If no control makes a wrong decision, report that outcome. Do not claim certification is faster or that ordinary approximation necessarily fails. The intended added value is a verifiable numerical feasibility and minimality decision under a supplied speech model. A model-level 95% requirement is not a guarantee of 95% empirical speech coverage, perceptual naturalness or correct corpus annotation.
