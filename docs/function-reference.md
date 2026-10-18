# Function reference

This page inventories the named functions defined in the repository's active R source files. Names beginning with `.` are implementation helpers; they are documented here for maintainers and are not intended as stable public interfaces. Functions nested inside other functions are implementation details of their enclosing function.

## Data generation — `create_multivariate_survival_data.R`

- `incidence_to_hazards(eval_time, per_year_incidence, RR)`: Converts baseline incidence and a risk ratio to exponential event hazards. `eval_time` is retained for interface compatibility; the current calculation uses the annual incidence inputs.
- `simulate_correlated_survival_data(...)`: Simulates a cohort with treatment, primary outcome, negative-control outcome, censoring, and optional baseline covariate `X`. See the function signature and inline comments for the full parameter definitions.

## Counterfactual survival estimation — `counterfactual_survival_estimators.R`

- `CFsurvival.nuisance.options(...)`: Builds the nuisance-model configuration list passed to `CFsurvival::CFsurvival` in the simulation.
- `.estimate.uniform.quantile(...)`: Estimates the critical value for a uniform confidence band.
- `.risk.ratio(...)`: Constructs risk-ratio contrasts and their influence-function-based uncertainty; used by the oracle and EM-based estimators.

## Survival Super Learner — `survival_super_learner.R`

- `survSuperLearner2(...)`: Fits event and censoring survival Super Learners, allowing separate covariate sets, prediction rows, weights, and cluster-aware folds.
- `.checkInputs(...)`: Checks and normalizes inputs to `survSuperLearner2`.
- `survSuperLearner.control(...)`: Creates fitting-control options, including initial weights and saved learner fits.
- `survSuperLearner.CV.control(...)`: Creates cross-validation controls.
- `.createLibrary(...)`: Normalizes a learner library specification.
- `.survCVFolds2(...)`: Creates cross-validation folds, respecting IDs when supplied.
- `.surviterativeSL2(...)`: Performs iterative ensemble fitting for event and censoring survival predictions.
- `.survcomputeCoef_stream(...)`: Computes ensemble coefficients using streamed risk calculations.

## EM-based NCI adjustment — `nci_adjustment_em.R`

- `AIPWCC_EIF_em(...)`: Estimates counterfactual survival contrasts with the EM-augmented AIPWCC influence-function estimator. Supports cross-fitting, nuisance libraries, and EM updates. This is the adjusted estimator called by `nci_sim.R`.
- `survSL.coxph.ps(...)` / `predict.survSL.coxph.ps(...)`: Cox proportional-hazards learner and prediction method using a piecewise-spline baseline.
- `.CFsurvival.nuis.opt(...)`: Creates nuisance-model options for the NCI-adjusted estimator.
- `.pSSL_event_only(...)`: Extracts/predicts event survival from a fitted event-only Super Learner.
- `.compute_aug_likelihood(...)`: Evaluates the likelihood contributions for augmented training records.
- `.update_aug_weights(...)`: Updates latent-record weights from the current survival fit.
- `.fit_SY_em(...)`: Fits the conditional primary-event survival model through EM updates.
- `.augment_fold(...)`: Builds augmented training data within a cross-fitting fold.
- `.expand_test(...)`: Expands held-out records over the negative-control evaluation grid.
- `.compute_posterior_fN_2(...)`: Computes posterior negative-control event probabilities using the current survival fits (vectorized implementation).
- `.get.survival.AIPWCC(...)`: Computes treatment-specific survival estimates and influence-function terms for the NCI-adjusted estimator.

## Oracle comparator — `oracle_nci_adjustment.R`

- `AIPWCC_newEIF_oracle(...)`: Computes the oracle NCI-adjusted estimator using known data-generating quantities; used as a benchmark in the simulation.
- `.compute_posterior_fN(...)`: Computes posterior negative-control event probabilities for the oracle implementation.
- `.get.survival.AIPWCC.oracle(...)`: Computes oracle survival estimates and influence-function contributions.

## Where to find argument details

The estimator functions include substantial inline documentation and comments. Their signatures are the authoritative source for defaults; this reference describes each function's role and the simulation's call graph.
