# Conditional syllable-duration evaluation

The task is duration allocation in two-syllable English words under a supplied word-duration budget. Inputs are syllable phone sequences, lexical stress, position and total word duration. The output is the probability that the first syllable extends beyond the word midpoint. The internal boundary is hidden at prediction time. The task concerns conditional allocation, not text-only prediction of total word duration.

## Data and split

CMU ARCTIC 0.95 SLT and BDL provide Festival syllable structure and recorded-speech phone end times. The parser checks the phone sequences and times against the release's lab files and ignores source_end, dur_factor, pitch and generated or acoustic features. Labels are automatic reference annotations, not manual syllable boundaries. Source notices and licenses accompany the data.

Eligible words have exactly two nonempty syllables, positive phone durations and contiguous non-pause phones. All exclusions are recorded. split_assignments.json gives the fixed assignment of utterance identifiers to training, validation and test; both voices' readings of a prompt share its split. Fits are voice-specific. Test outcomes do not inform feature selection, regularization, residual scales or shared-factor tolerance.

## Models

The unary model is ridge regression for duration in milliseconds, with penalty 10 and an unpenalized intercept. Features are phoneme counts, lexical stress, syllable position and their interaction. Predicted targets have a 5 ms floor. Training residual root-mean-square errors within stress-position classes supply scales, with a 10 ms floor. The nonnegative duration density has symmetric quadratic costs.

The shared factor penalizes the difference between the sum of syllable durations and the supplied word duration. Its tolerance is selected from 5, 10, 20 and 40 ms by validation Brier score, weighting voices equally. Targets and weights are exported as exact rational encodings of binary64 fitted inputs.

The independent baseline removes the shared factor from the same unary model. Proportional allocation uses the fitted target ratio plus a training-fitted Gaussian residual mean and scale for the first-syllable fraction, truncated to [0,1]. Both return probabilities under the same scoring rule.

## Evaluation and inference

The primary metric is Brier score and the secondary metric is log loss. Report each voice and their equally weighted mean. Paired prompt-cluster bootstrap intervals use 2,000 replicates and seed 20261016, retaining both voices and all words of a sampled prompt together. These intervals concern sampling variation in this cohort; they do not estimate new-speaker generalization.

The enclosing computation uses 128-bit Arb, the product-event routine and its unary exterior allowance. Meshes 20/2^j ms on [0,1600]^2 are refined until probability width is at most 0.002, with all attempts and failures retained. Only enclosing computations supply a certificate. Point quadrature is used for fitting and independent numerical checks.

The planning experiment in deadlines/ uses the same fitted models and test inputs. It evaluates numerical feasibility and minimality decisions, not a second independent predictive test. Neither experiment establishes perceptual synthesis gains or generalization beyond the tested corpus and voices.
