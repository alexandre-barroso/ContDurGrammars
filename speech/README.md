# Recorded speech timing experiments

This application evaluates conditional syllable-duration allocation and certified deadline planning on CMU ARCTIC speech annotations. Inputs are phones, lexical stress and a supplied word duration used as a soft timing target. Internal syllable boundaries are withheld from prediction and planning. The numerical routines are the Python modules in the repository root.

The experiments use automatic phone labels and lexicon-derived syllable structure. They evaluate a component of duration-controlled synthesis; they do not generate speech or measure perceptual quality.

## Results

All 836 eligible test words are included: 418 per voice and 203 shared test prompts. Both readings of a prompt belong to the same split.

| Model | SLT Brier | BDL Brier | Equal-voice Brier | Equal-voice log loss |
|---|---:|---:|---:|---:|
| Grammar with shared word target | 0.0917 | 0.0777 | 0.0847 | 0.2752 |
| Independent unary durations | 0.1110 | 0.0991 | 0.1050 | 0.3429 |
| Proportional allocation | 0.1156 | 0.0922 | 0.1039 | 0.3299 |

Every probability enclosure has width at most 0.002. Exact endpoint scoring encloses the pooled Brier score in [0.0845104, 0.0848983]. The paired prompt-cluster bootstrap gives 95% intervals of [-0.02939, -0.01184] for the grammar-minus-independent difference and [-0.02842, -0.00987] for grammar-minus-proportional, using 2,000 replicates and seed 20261016. Both voices and all words within a sampled prompt stay together. Numerical enclosures and sampling uncertainty have different meanings.

The [deadline experiment](deadlines/README.md) chooses the earliest 10 ms syllable deadline meeting a 95% model-probability requirement. Certification resolves 787 feasible and 49 infeasible inputs. Ordinary 10 ms integration makes 18 decision errors; eight survive a successive-mesh stability check. The 1.25 ms mesh and high-accuracy quadratic control agree with every certified decision.

## Inspect and reproduce

Install the repository's requirements-integration.txt in Python 3.12. From this directory:

```sh
python validate.py
python evaluate.py
python deadlines/evaluate.py
python deadlines/independent_checks.py
```

To regenerate the data extraction, fits and probability certificates:

```sh
python prepare_data.py
python fit_models.py
python certify_predictions.py --workers 4
python evaluate.py
python validate.py
```

The included data/ directory contains all annotation files needed for these commands, with source notices and licenses. [DATA_SOURCES.json](DATA_SOURCES.json) records the archive URLs. fetch_data.py can download those archives and extract the annotation, transcript and license files. Full audio is not required for the reported experiments.

Matching certificates are reused after direct comparison of their model and event inputs. For a fresh timing measurement, use an empty certificates/ directory. The recorded four-worker probability batch took 152.95 seconds, with median worker time 0.70 seconds, on a shared arm64 macOS machine. Different numerical-library builds may slightly change the fitted binary64 values; the supplied exact rational inputs specify the reported models.

## Data and outputs

- [PROTOCOL.md](PROTOCOL.md) defines the task, model, fitting, split and scoring.
- [split_assignments.json](split_assignments.json) assigns all 1,132 prompt identifiers to training, validation or test. The assignment is shared by the voices.
- prepared_words.json contains every eligible word, exact millisecond durations, phones, stress, source indices and split. data_audit.json records counts and exclusions.
- fitted_models.json contains training fits, validation scores and the selected 40 ms shared-target tolerance.
- test_requests.json contains every test input, exact rational targets and weights, baseline probabilities and the independent quadratic control. Internal durations are included for scoring and provenance; validate.py checks that changing them while holding the budget fixed cannot change model inputs.
- certificates/ contains every mesh attempt, model tree, event, integral and probability endpoint, exterior allowance and elapsed time. certification_run.json records the batch environment and outcomes.
- test_predictions.json, evaluation.json and validation.json contain scored predictions, aggregate comparisons, bootstrap settings and validation results.
- deadlines/ contains the planning controls, all decisions and their interval certificates.

The strict midpoint event scores 14 exact reference ties as zero. Log loss uses a numerical probability floor of 10^-15; no reported value required clipping. Grammar point scores use enclosure midpoints, while Brier bounds use exact rational endpoints. Ratio endpoints are intersected with [0,1] after serialization to avoid outward expansion when reboxing a bounded probability in midpoint-radius arithmetic.

An independent calculation completes the square and integrates a truncated bivariate Gaussian by one-dimensional quadrature. Its probability lies inside every enclosure. It checks the application calculation; it is not a certificate or evidence that a general integrator is necessary for the quadratic special case. Predictive improvement belongs to the fitted timing model, while certification bounds numerical error in evaluating that model.

## Data attribution

John Kominek and Alan W. Black (2004), “The CMU Arctic Speech Databases,” Proceedings of the 5th ISCA Speech Synthesis Workshop, pp. 223–224. Primary paper: https://www.cs.cmu.edu/~awb/papers/ssw5/arctic.pdf. Release website: http://festvox.org/cmu_arctic/.

CMU ARCTIC 0.95 SLT and BDL annotations are distributed under each voice's COPYING notice. Source annotation files and notices are unmodified. Derived data, fits and evaluation outputs are provided with the application.
