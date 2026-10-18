# NCI-Adjust

Simulation code for evaluating negative-control adjustment of treatment effects in right-censored survival data. The repository contains the data-generating mechanism, counterfactual survival estimators, an EM-based NCI-adjusted estimator, an oracle comparator, and the `SimEngine` study driver.

## Repository map

| File | Purpose |
| --- | --- |
| `nci_sim.R` | Main `SimEngine` simulation specification and analysis pipeline. |
| `create_multivariate_survival_data.R` | Generates treatment, primary-event, negative-control-event, and censoring data. |
| `counterfactual_survival_estimators.R` | Nuisance-option defaults and risk-ratio inference utilities used by the simulation estimators. |
| `survival_super_learner.R` | Super Learner implementation for conditional event and censoring survival functions. |
| `nci_adjustment_em.R` | EM-based NCI-adjusted AIPWCC estimator and supporting learners. |
| `oracle_nci_adjustment.R` | Oracle comparator that uses known data-generating distributions. |
| `docs/function-reference.md` | Function-by-function guide to the source files. |

## Simulation workflow

The driver runs four analyses on each generated data set: unadjusted Kaplan–Meier, covariate-adjusted `CFsurvival`, oracle NCI adjustment, and cross-fitted EM-based NCI adjustment. It returns replicate-level estimates, confidence limits, variances, bias, relative efficiency, and coverage.

Run from the repository root so the relative source paths resolve:

```r
source("nci_sim.R")
```

`nci_sim.R` is configured for the Fred Hutch Rhino SLURM environment. It loads `SimEngine` and several survival/learning packages and currently uses `/home/eashby/R_lib` as the cluster library path. Install the required packages in your environment and adjust that path before running elsewhere. The default grid, replication count, seed, core count, and SLURM configuration are set in the driver.

## Reproducibility notes

- The simulation sources the estimator and data-generation files explicitly from the repository root. Keep these filenames synchronized if the layout changes.
- The oracle estimator calls the same data generator to construct its large reference sample; its `dgm_params$n` is intentionally much larger than the analysis sample.
- `CFsurvival_funcs.R` and `survSuperLearner_funcs.R` were renamed to descriptive filenames; function names used by the source code remain available.
- See the function reference for the implemented entry points, arguments, and internal helpers.

## Citation and manuscript use

This repository accompanies a manuscript submission. Add the manuscript citation, software version/commit, and any package citation requirements here once the manuscript metadata are final. The estimator code notes its adaptation from the CFsurvival package in the source file.
