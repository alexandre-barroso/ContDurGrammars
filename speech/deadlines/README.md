# Verified syllable deadlines for a speech-duration planner

The planner takes a fitted grammar, linguistic inputs and a word budget T. It returns the earliest deadline d on a 10 ms clock satisfying P(first syllable ends by d | inputs, T) >= 0.95, with d <= T. If no frame qualifies, the result requests a revised word budget. This is an executable numerical planning component; the clock and probability requirement are design choices, not perceptual standards.

All 836 test inputs from the recorded-speech task are included. The same fitted models and supplied word budgets are used. This is a numerical decision evaluation, not an independent predictive test.

## Results

Exact interval comparisons resolve all 836 decisions: 787 have a feasible earliest deadline, and 49 have no feasible frame within the supplied budget. No request is unresolved or omitted.

| Ordinary computation | Wrong frame | Wrong feasibility | Total wrong decisions |
|---|---:|---:|---:|
| Midpoint/FFT, 20 ms | 86 | 3 | 89 |
| Midpoint/FFT, 10 ms | 17 | 1 | 18 |
| Midpoint/FFT, 5 ms | 7 | 1 | 8 |
| Midpoint/FFT, 2.5 ms | 2 | 0 | 2 |
| Midpoint/FFT, 1.25 ms | 0 | 0 | 0 |
| Successive-mesh stability | 7 | 1 | 8 |
| High-accuracy quadratic quadrature | 0 | 0 | 0 |

Wrong feasibility includes either a falsely accepted budget or a false infeasibility declaration. At 10 ms, the errors are 17 premature deadlines and one falsely accepted budget. All methods use the same 10 ms output clock; the integration mesh is a separate quantity. The stability heuristic requires an unchanged decision and CDF changes <= 0.001 at both decisive endpoints. It returns a decision on all inputs, but eight decisions are wrong.

The first primary disagreement is Chinese, SLT/arctic_a0120/8. Its phone sequence is ch ay n | iy z, stress 0/1 and supplied total 512.5 ms. Ordinary 10 and 5 ms integration assign probabilities 0.950366 and 0.950011 to a 340 ms deadline, so both accept it. The enclosing computation proves F(340) is in [0.9498569, 0.9499287] and F(350) is in [0.9714059, 0.9714794]. It therefore reserves 350 ms, one additional frame. The sole false-feasibility case, little in BDL/arctic_b0340/2, has last admissible frame 240 ms; its probability upper bound is 0.9499527, so the budget requires revision.

The four-worker certified batch took 139.92 seconds; median accumulated worker time was 0.445 seconds and the maximum was 48.504 seconds. All ordinary controls together took 2.75 seconds in one serial run. These are elapsed measurements on a shared arm64 macOS machine, not a controlled speed benchmark. Certification is more expensive; fine ordinary integration succeeds on every input here.

## Why the certificate adds value

The stopping rule verifies both feasibility and minimality: the selected frame's CDF lower bound is >= 0.95, and the preceding frame's upper bound is < 0.95. For an infeasible budget, the upper bound at its last frame is < 0.95. CDF monotonicity then establishes the decision. The decision rule applies the model-probability enclosures described in the paper. A stable ordinary estimate cannot itself establish these inequalities; the observed wrong decisions illustrate the difference.

This guarantees the numerical decision under the supplied grammar and the paper's primitive-operation contract. It does not establish empirical 95% coverage, better audible synthesis, or correctness of automatic speech annotations. The violations in the illustrated cases are small in probability but change a declared discrete planning decision. Whether an application needs this strict contract is a separate engineering choice. No claim is made that certification is necessary to obtain a correct point answer, or that it is faster.

## Reproduction

With the repository requirements-integration.txt installed, run from the parent application directory:

~~~sh
python deadlines/run.py --ordinary
python deadlines/run.py --certify --workers 4
python deadlines/evaluate.py
python deadlines/independent_checks.py
~~~

The certificate driver resumes matching completed models. For fresh timing, move deadlines/certificates/ and its certification_run.json aside first. Every test request is included. Models are compared directly to their stored trees before a certificate is reused. Both point and enclosing computations use the modules in the repository root.

ordinary.json retains every point decision, queried CDF, timing and clipping diagnostic. certificates/ retains every interval, integral, mesh attempt and decision. certification_run.json records the execution environment and outcomes. evaluation.json and decisions.json include every comparison and exact endpoint verification. All 5,407 independently computed CDF values checked lie inside their corresponding enclosures. independent_checks.json additionally records 50- and 80-digit calculations for the two illustrated cases, plus direct two-dimensional midpoint sums agreeing with the FFT controls within 1e-13. High-precision agreement is a validation check; only the outward-rounded interval computation supplies the numerical guarantee.
