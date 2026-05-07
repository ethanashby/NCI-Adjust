# =============================================================================
# sandwich_variance.r
#
# Sandwich variance correction for the AIPWCC estimator, following the
# stacked estimating equations approach of Lee et al. (SPARCC, 2025).
#
# The problem:
#   var(phi_obs) / n underestimates the true variance of the AIPWCC estimator
#   because it treats hat_f_{N|A,X} as fixed. In reality, the augmented
#   dataset weights AND the posterior f_{N|Y,DeltaY,A,X} both depend on
#   hat_alpha (the parameters of f_{N|A,X}), so estimation error in hat_alpha
#   propagates into phi_obs.
#
# The fix (stacked M-estimation / sandwich):
#   Stack the score equations for alpha (f_{N|A,X} model) with the estimating
#   equation for psi (the counterfactual survival). The corrected influence
#   function for participant i is:
#
#     phi_tilde_i = (phi_obs_i - psi_hat) + A21 %*% A11_inv %*% S_alpha_i
#
#   where:
#     S_alpha_i  = score of log f_{N|A,X} for participant i
#     A11        = -E[d S_alpha_i / d alpha^T]  (Fisher info of N model)
#     A21        = -E[d phi_obs_i / d alpha^T]  (sensitivity of AIPWCC IF
#                   to perturbations in alpha; computed by finite differences)
#
#   The corrected variance is then Var(phi_tilde_i) / n.
#
# Structure:
#   1. fit_fN_parametric()    - fits f_{N|A,X} parametrically (Cox or exp),
#                               returns alpha_hat and score function
#   2. compute_score_fN()     - computes S_alpha_i for each participant
#   3. compute_A11()          - estimates bread A11 = Fisher info of N model
#   4. compute_A21_numderiv() - estimates A21 by finite differences
#   5. sandwich_IF_correction()- applies the correction to IF.vals
#   6. sandwich_logRR_var()   - computes corrected log-RR SE
#
# Integration with AIPWCC_CFsurvival:
#   After the main estimation loop, call sandwich_IF_correction() using the
#   fitted FN_fit object (from the marginal f_{N|A,X} fit), the original
#   data, and the IF.vals matrices from .get.survival.AIPWCC.
# =============================================================================


# =============================================================================
# 1. fit_fN_parametric()
#
# Fits f_{N|A,X} on the FULL dataset (marginal, not cross-fitted) using
# either a Cox or exponential parametric model.
# Returns the fitted object and a function to evaluate the score S_alpha_i.
#
# Arguments:
#   time_N    - observed N times (length n)
#   DeltaN_eff- recoded event indicator (admin cens = 1)
#   A         - treatment (length n)
#   W         - covariate data frame (n x p)
#   model     - "cox" or "exponential"
# =============================================================================

fit_fN_parametric <- function(time_N, DeltaN_eff, A, W, model = "cox") {
  
  dat <- data.frame(time_N = time_N, DeltaN_eff = DeltaN_eff, A = A, W)
  
  if (model == "cox") {
    # Fit Cox model: alpha = regression coefficients
    # Score = partial likelihood score
    
    fit <- survSL.coxph(time = dat$time_N,
                        event = dat$DeltaN_eff,
                        X = bind_cols("A"=dat$A, W),
                        newX = bind_cols("A"=dat$A, W),
                        new.times = sort(unique(dat$time_N)),
                        obsWeights = rep(1, nrow(dat)))
    
    # Predicted baseline hazard at each observed time
    bh   <- survival::basehaz(fit$fit$object, centered = FALSE)
    
    survest_N <- function(t, x, a) {
      x <- as.data.frame(x)
      colnames(x) <- colnames(W)
      if (nrow(x) == 1L) { x <- rbind(x, x); drop1 <- TRUE } else drop1 <- FALSE
      newdat <- bind_cols(A = a, x)
      out <- predict.survSL.coxph(fit$fit, newX=newdat, new.times=t)
      if (drop1) out <- out[1,, drop=FALSE]
      as.vector(out)
    }
    
    return(list(
      fit       = fit,
      model     = "cox",
      alpha_hat = coef(fit$fit$object),
      bh        = bh,
      dat       = dat,
      survest_N = survest_N
    ))
    
  } else if (model == "exponential") {
    # Fit exponential AFT: log(T) ~ A + W, i.e. hazard = exp(-X beta)
    
    fit <- survSL.expreg(time = dat$time_N,
                        event = dat$DeltaN_eff,
                        X = bind_cols("A"=dat$A, W),
                        newX = bind_cols("A"=dat$A, W),
                        new.times = sort(unique(dat$time_N)),
                        obsWeights = rep(1, nrow(dat)))
    
    survest_N <- function(t, x, a) {
      x <- as.data.frame(x)
      colnames(x) <- colnames(W)
      if (nrow(x) == 1L) { x <- rbind(x, x); drop1 <- TRUE } else drop1 <- FALSE
      newdat <- bind_cols(A = a, x)
      out <- predict.survSL.expreg(fit$fit, newX=newdat, new.times=t)
      if (drop1) out <- out[1,, drop=FALSE]
      as.vector(out)
    }
    
    return(list(
      fit       = fit,
      model     = "exponential",
      alpha_hat = coef(fit$fit$reg.object),
      scale     = fit$fit$reg.object$scale,
      dat       = dat,
      survest_N = survest_N
    ))
  } else if(model=="weibull"){
  
    fit <- survSL.weibreg(time = dat$time_N,
                         event = dat$DeltaN_eff,
                         X = bind_cols("A"=dat$A, W),
                         newX = bind_cols("A"=dat$A, W),
                         new.times = sort(unique(dat$time_N)),
                         obsWeights = rep(1, nrow(dat)))
    
    survest_N <- function(t, x, a) {
      x <- as.data.frame(x)
      colnames(x) <- colnames(W)
      if (nrow(x) == 1L) { x <- rbind(x, x); drop1 <- TRUE } else drop1 <- FALSE
      newdat <- bind_cols(A = a, x)
      out <- predict.survSL.weibreg(fit$fit, newX=newdat, new.times=t)
      if (drop1) out <- out[1,, drop=FALSE]
      as.vector(out)
    }
    
    return(list(
      fit       = fit,
      model     = "weibull",
      alpha_hat = coef(fit$fit$reg.object),
      scale     = fit$fit$reg.object$scale,
      dat       = dat,
      survest_N = survest_N
    ))
  } else {
    stop("model must be 'cox', 'weibull', or 'exponential'")
  }
}


# =============================================================================
# 2. compute_score_fN()
#
# Computes the score S_alpha_i = d log f_{N,DeltaN | A,X}(N_i, DeltaN_i) / d alpha
# for each participant i. Returns an (n x p_alpha) matrix.
#
# For the Cox model: uses the Schoenfeld residuals (partial likelihood score
# contributions). survival::residuals(fit, type="score") gives these directly.
#
# For the exponential model: closed-form score of the log-likelihood.
# =============================================================================

compute_score_fN <- function(fN_fit) {
  
  if (fN_fit$model == "cox") {
    # residuals(..., type="score") returns the n x p score matrix
    # where p = number of covariates
    S <- residuals(fN_fit$fit$fit$object, type = "score")
    if (is.vector(S)) S <- matrix(S, ncol = 1)
    return(S)   # n x p_alpha matrix
    
  } else if (fN_fit$model == "exponential") {
    dat     <- fN_fit$dat
    fit     <- fN_fit$fit$fit$reg.object
    n       <- nrow(dat)
    X       <- model.matrix(fit)          # n x p design matrix
    beta    <- coef(fit)
    sigma   <- fit$scale                  # = 1 for exponential
    mu_i    <- as.vector(X %*% beta)      # log(mean) for each obs
    t_i     <- dat$time_N
    d_i     <- dat$DeltaN_eff
    # Log-likelihood for exponential: d_i * (log(t_i) - mu_i) - t_i * exp(-mu_i)
    # Score w.r.t. beta: (d_i - t_i * exp(-mu_i)) * (-X_i)  [AFT parameterisation]
    resid_i <- d_i - t_i * exp(-mu_i)    # n x 1
    S       <- sweep(-X, 1, resid_i, FUN = "*")   # n x p
    return(S)
  }
}

# =============================================================================
# 3. compute_A11()
#
# Estimates A11 = -E[d S_alpha_i / d alpha^T] = Fisher information matrix
# of the f_{N|A,X} model. For Cox: use the observed information (negative
# Hessian of the partial log-likelihood). For exponential: use the observed
# information of the log-likelihood.
# =============================================================================

compute_A11 <- function(fN_fit) {
  
  if (fN_fit$model == "cox") {
    # vcov(coxph) = (observed information)^{-1}, so A11 = solve(vcov(fit))
    A11 <- Matrix::solve(vcov(fN_fit$fit$fit$object))
    return(A11)   # p_alpha x p_alpha
    
  } else if (fN_fit$model == "exponential") {
    A11 <- Matrix::solve(vcov(fN_fit$fit$fit$reg.object))
    return(A11)
  }
}


# =============================================================================
# 4. compute_A21_numderiv()
#
# Estimates A21 = -E[d phi_obs_i(alpha) / d alpha^T] by numerical
# differentiation. For each direction j in alpha, we perturb alpha_hat
# by +/- epsilon, refit the augmented dataset and recompute phi_obs, then
# take the central difference.
#
# This is the most expensive step: requires 2 * p_alpha full passes through
# the augmented dataset pipeline (rebuild aug_train, refit S_Y, recompute
# phi_obs). To make this tractable:
#   - Use a FIXED eval.times grid (no refitting SuperLearner for S_Y)
#   - Only recompute the posterior f_N and the AIPWCC IF values
#   - Hold S_Y hat fixed (the sensitivity of phi_obs to alpha through S_Y
#     is second-order when cross-fitted; the dominant first-order term comes
#     through the posterior and the augmented weights)
#
# Arguments:
#   fN_fit        - output of fit_fN_parametric()
#   phi_obs_0     - (n x K) matrix of baseline IF values at alpha_hat
#   time, event, time_N, DeltaN_eff, A, W, admin_t - original data
#   S.hats_0, S.hats_1 - (n_aug x K) matrices of S_Y predictions (held fixed)
#   GY0, GY1, GN0, GN1 - (n x K) censoring predictions (held fixed)
#   g.hats        - propensity scores (held fixed)
#   eval.times, fit.times, ids - grid and id vectors
#   treat         - original treatment vector
#   epsilon       - perturbation size (default 1e-4)
#   fit.treat     - c(0,1) or subset
#   verbose       - logical
#
# Returns A21: (1 x p_alpha) row vector (for a scalar psi at fixed t0)
# Extended to (K x p_alpha) for a vector of time points.
# =============================================================================

compute_A21_numderiv <- function(
    fN_fit,
    phi_obs_0,          # n x K baseline IF values
    time, event, time_N, DeltaN_eff, A_orig, W, admin_t,
    S.hats_0, S.hats_1, # aug-aligned, held fixed
    GY0, GY1, GN0, GN1,
    g.hats,
    eval.times, fit.times, ids,
    treat,
    W_names,
    epsilon = 1e-4,
    verbose = TRUE
) {
  
  n       <- length(time)
  K       <- length(fit.times)
  p_alpha <- length(fN_fit$alpha_hat)
  
  A21 <- matrix(0, nrow = K, ncol = p_alpha)  # K x p_alpha
  
  for (j in seq_len(p_alpha)) {
    
    if (verbose) message(sprintf("  A21: perturbing alpha[%d] / %d", j, p_alpha))
    
    # ---- Perturb alpha up and down ----
    alpha_up   <- fN_fit$alpha_hat
    alpha_down <- fN_fit$alpha_hat
    alpha_up[j]   <- alpha_up[j]   + epsilon
    alpha_down[j] <- alpha_down[j] - epsilon
    
    phi_up   <- .recompute_phi_perturbed(
      fN_fit, alpha_up, time, event, time_N, DeltaN_eff, A_orig,
      W, admin_t, S.hats_0, S.hats_1, GY0, GY1, GN0, GN1,
      g.hats, eval.times, fit.times, ids, treat, W_names
    )
    
    phi_down <- .recompute_phi_perturbed(
      fN_fit, alpha_down, time, event, time_N, DeltaN_eff, A_orig,
      W, admin_t, S.hats_0, S.hats_1, GY0, GY1, GN0, GN1,
      g.hats, eval.times, fit.times, ids, treat, W_names
    )
    
    # Central difference: d phi_obs / d alpha_j ≈ (phi_up - phi_down) / (2 eps)
    dphi_dalpha_j <- (phi_up - phi_down) / (2 * epsilon)   # n x K
    
    # A21[k, j] = -mean_i(d phi_obs_i_k / d alpha_j)
    A21[, j] <- -colMeans(dphi_dalpha_j, na.rm = TRUE)
  }
  
  return(A21)  # K x p_alpha
}


# =============================================================================
# Internal helper: recompute phi_obs with a perturbed alpha
#
# Holds S_Y fixed; only rebuilds:
#   (a) The survest_N function using the perturbed alpha
#   (b) The posterior f_N via Bayes rule
#   (c) The AIPWCC IF values
#
# This is the key approximation: the dominant sensitivity of phi_obs to alpha
# flows through the posterior (which reweights phi^F over N grid points)
# and through the augmented training weights (which affect S_Y indirectly).
# We hold S_Y fixed to avoid refitting SuperLearner on each perturbation.
# =============================================================================

.recompute_phi_perturbed <- function(
    fN_fit, alpha_perturbed,
    time, event, time_N, DeltaN_eff, A_orig,
    W, admin_t, S.hats_0, S.hats_1,
    GY0, GY1, GN0, GN1, g.hats,
    eval.times, fit.times, ids, treat, W_names
) {
  
  n <- length(time)
  l <- length(eval.times)
  
  # ------------------------------------------------------------------
  # Build perturbed survest_N function
  # ------------------------------------------------------------------
  survest_N_pert <- .make_survest_N_perturbed(fN_fit, alpha_perturbed, W_names)
  
  # ------------------------------------------------------------------
  # Recompute posterior f_N on full dataset using perturbed survest_N
  # We use all participants as both "train" and "test" since we are
  # doing a marginal (non-cross-fitted) perturbation analysis
  # ------------------------------------------------------------------
  dat_all <- data.frame(
    id = seq_len(n), Y = time, DeltaY = event,
    N = time_N, DeltaN_eff = DeltaN_eff,
    A = treat, admin_t = admin_t
  )
  dat_all <- cbind(dat_all, W)
  
  # aug_test: expand all participants over eval.times grid
  aug_all_list <- vector("list", n)
  for (i in seq_len(n)) {
    ri      <- dat_all[i, ]
    times_i <- eval.times[eval.times <= ri$admin_t]
    l_i     <- length(times_i)
    base    <- ri[rep(1L, l_i), c("id","Y","DeltaY","A","admin_t", W_names), drop=FALSE]
    base$N      <- times_i
    base$DeltaN <- c(rep(1L, l_i - 1L), 0L)
    aug_all_list[[i]] <- base
  }
  aug_all <- do.call(rbind, aug_all_list)
  
  # Recompute posteriors using perturbed survest_N but fixed S_Y
  post_pert <- .compute_posterior_fN(
    dat_te    = dat_all,
    aug_te    = aug_all,
    SY0       = S.hats_0,
    SY1       = S.hats_1,
    survest_N = survest_N_pert,
    fit.times = fit.times,
    eval.times = eval.times,
    W_names   = W_names
  )
  
  # ------------------------------------------------------------------
  # Recompute IF values for A=0 curve
  # ------------------------------------------------------------------
  phi_0_pert <- .get.survival.AIPWCC(
    Y = time, DeltaY = event, N = time_N, DeltaN = DeltaN_eff,
    A = 1L - treat, admin_t = admin_t,
    fit.times = fit.times, eval.times = eval.times, ids = ids,
    S.hats         = S.hats_0,
    GY.hats_a0     = GY0, GY.hats_alt  = GY1,
    GN.hats_a0     = GN0, GN.hats_alt  = GN1,
    g.hats         = 1 - g.hats,
    fN_a0          = post_pert$fN0,
    fN_alt         = post_pert$fN1,
    isotonize      = FALSE,
    verbose        = FALSE
  )$IF.vals
  
  # ------------------------------------------------------------------
  # Recompute IF values for A=1 curve
  # ------------------------------------------------------------------
  phi_1_pert <- .get.survival.AIPWCC(
    Y = time, DeltaY = event, N = time_N, DeltaN = DeltaN_eff,
    A = treat, admin_t = admin_t,
    fit.times = fit.times, eval.times = eval.times, ids = ids,
    S.hats         = S.hats_1,
    GY.hats_a0     = GY1, GY.hats_alt  = GY0,
    GN.hats_a0     = GN1, GN.hats_alt  = GN0,
    g.hats         = g.hats,
    fN_a0          = post_pert$fN1,
    fN_alt         = post_pert$fN0,
    isotonize      = FALSE,
    verbose        = FALSE
  )$IF.vals
  
  # Return a list for both arms
  list(phi_0 = phi_0_pert, phi_1 = phi_1_pert)
}


# Helper: build a survest_N function from perturbed alpha
.make_survest_N_perturbed <- function(fN_fit, alpha_perturbed, W_names) {
  
  if (fN_fit$model == "cox") {
    
    fit_orig <- fN_fit$fit
    dat_orig <- fN_fit$dat
    bh       <- fN_fit$bh
    
    # Cumulative baseline hazard at each time
    H0_fun <- stepfun(bh$time, c(0, bh$hazard), right = FALSE)
    
    survest_N_pert <- function(t, x, a) {
      x  <- as.data.frame(x)
      colnames(x) <- W_names
      xmat <- as.matrix(cbind(A = a, x))
      lp   <- as.vector(xmat %*% alpha_perturbed)
      H0_t <- H0_fun(t)
      S    <- exp(-outer(H0_t, exp(lp)))  # length(t) x nrow(x)
      as.vector(t(S))  # return as vector: rows = x-observations, cols = times
      # Actually survest_N returns for ONE x, multiple t → length(t) vector
      # Redo: for single x, return S for each time
      as.vector(exp(-H0_t * exp(lp[1])))
    }
    
    return(survest_N_pert)
    
  } else if (fN_fit$model == "exponential") {
    
    fit_orig <- fN_fit$fit
    
    survest_N_pert <- function(t, x, a) {
      x    <- as.data.frame(x)
      colnames(x) <- W_names
      xmat <- as.matrix(cbind(1, A = a, x))  # add intercept
      mu   <- as.vector(xmat %*% alpha_perturbed)  # log(mean)
      rate <- exp(-mu)
      as.vector(exp(-rate * t))
    }
    
    return(survest_N_pert)
  }
}


# =============================================================================
# 5. sandwich_IF_correction()
#
# Applies the sandwich correction to the IF matrices for both arms.
# Returns corrected IF.vals.center for arm 0 and arm 1.
#
# Corrected IF for arm a:
#   phi_tilde_i^(a) = phi_obs_i^(a) + correction_i^(a)
# where:
#   correction_i^(a) = A21^(a) %*% A11_inv %*% S_alpha_i
#
# A21^(a) is computed separately for each arm because phi_obs depends on
# the treatment arm through the propensity weighting and the posterior.
#
# Arguments:
#   IF.vals.0, IF.vals.1  - (n x K) uncentered IF matrices from .get.survival.AIPWCC
#   A21_0, A21_1          - (K x p_alpha) sensitivity matrices for each arm
#   A11_inv               - (p_alpha x p_alpha) inverse Fisher info
#   S_alpha               - (n x p_alpha) score matrix from compute_score_fN()
#
# Returns list with corrected IF.vals.center.0 and IF.vals.center.1
# =============================================================================

sandwich_IF_correction <- function(
    IF.vals.0, IF.vals.1,
    A21_0, A21_1,
    A11_inv,
    S_alpha
) {
  n <- nrow(IF.vals.0)
  K <- ncol(IF.vals.0)
  
  # Correction matrix: n x K
  # correction[i, k] = A21[k, ] %*% A11_inv %*% S_alpha[i, ]
  #                  = (A21 %*% A11_inv %*% t(S_alpha))[k, i]
  # = t(S_alpha %*% t(A21 %*% A11_inv))
  
  # Precompute A21 %*% A11_inv: (K x p_alpha) %*% (p_alpha x p_alpha) = K x p_alpha
  C0 <- A21_0 %*% A11_inv    # K x p_alpha
  C1 <- A21_1 %*% A11_inv    # K x p_alpha
  
  # Correction for each arm: n x K
  # correction[i, k] = C[k, ] %*% S_alpha[i, ] = (S_alpha %*% t(C))[i, k]
  correction_0 <- S_alpha %*% t(C0)   # n x K
  correction_1 <- S_alpha %*% t(C1)   # n x K
  
  # Corrected uncentered IF
  IF_corrected_0 <- IF.vals.0 + correction_0
  IF_corrected_1 <- IF.vals.1 + correction_1
  
  # Center
  surv_0 <- colMeans(IF.vals.0, na.rm = TRUE)
  surv_1 <- colMeans(IF.vals.1, na.rm = TRUE)
  
  IF_center_corrected_0 <- sweep(IF_corrected_0, 2, surv_0, "-")
  IF_center_corrected_1 <- sweep(IF_corrected_1, 2, surv_1, "-")
  
  list(
    IF.vals.center.0.corrected = IF_center_corrected_0,
    IF.vals.center.1.corrected = IF_center_corrected_1,
    correction_0 = correction_0,
    correction_1 = correction_1
  )
}


# =============================================================================
# 6. sandwich_logRR_var()
#
# Computes the sandwich-corrected SE for log(RR) at each time point,
# and constructs CIs.
#
# Arguments:
#   surv.0, surv.1                   - point estimates (length K)
#   IF.vals.center.0, IF.vals.center.1 - corrected centered IFs (n x K)
#   fit.times                         - time grid (length K)
#   conf.level                        - default 0.95
#
# Returns a data frame with columns: time, logRR, se.logRR, lower, upper, RR
# =============================================================================

sandwich_logRR_var <- function(
    surv.0, surv.1,
    IF.vals.center.0.corrected,
    IF.vals.center.1.corrected,
    fit.times,
    conf.level = 0.95
) {
  n   <- nrow(IF.vals.center.0.corrected)
  z   <- qnorm((1 + conf.level) / 2)
  
  # Delta method IF for log(RR) = log(1 - S_1) - log(1 - S_0)
  # d/d S_a [log(1 - S_a)] = -1/(1 - S_a)
  risk_0 <- pmax(1 - surv.0, 1e-8)
  risk_1 <- pmax(1 - surv.1, 1e-8)
  
  # log-RR IF for participant i at time k:
  # psi_i_k = -IF_1[i,k] / risk_1[k]  +  IF_0[i,k] / risk_0[k]
  # (negative because d(log risk)/d(surv) = -1/risk)
  psi <- sweep(-IF.vals.center.1.corrected, 2, risk_1, "/") +
    sweep( IF.vals.center.0.corrected, 2, risk_0, "/")
  
  se_logRR <- apply(psi, 2, sd, na.rm = TRUE) / sqrt(n)
  logRR    <- log(risk_1) - log(risk_0)
  
  data.frame(
    time     = fit.times,
    logRR    = logRR,
    RR       = exp(logRR),
    se.logRR = se_logRR,
    lower    = exp(logRR - z * se_logRR),
    upper    = exp(logRR + z * se_logRR)
  )
}


# =============================================================================
# 7. Master wrapper: compute_sandwich_variance()
#
# Wraps the full sandwich correction pipeline.
# Call this after AIPWCC_CFsurvival() has returned its result.
#
# Arguments:
#   result       - output of AIPWCC_CFsurvival()
#   fN_fit       - output of fit_fN_parametric() on full data
#   time, event, time_N, DeltaN_eff, A_orig (= treat), W, admin_t
#   nuis         - nuis list from result$nuisance
#   fit.times    - time grid
#   eval.times   - N grid
#   treat        - original treatment vector
#   epsilon      - finite difference step (default 1e-4)
#   verbose      - logical
#
# Returns the sandwich-corrected log-RR data frame and corrected IF matrices.
# =============================================================================

compute_sandwich_variance <- function(
    result,
    fN_fit,
    time, event, time_N, DeltaN_eff, treat, W, admin_t,
    nuis,
    fit.times,
    eval.times,
    epsilon = 1e-4,
    verbose = TRUE
) {
  
  W_names <- colnames(W)
  n       <- length(time)
  
  if (verbose) message("Sandwich correction: computing score S_alpha...")
  S_alpha  <- compute_score_fN(fN_fit)          # n x p_alpha
  
  if (verbose) message("Sandwich correction: computing A11 (Fisher info)...")
  A11      <- compute_A11(fN_fit)               # p_alpha x p_alpha
  A11_inv  <- solve(A11)
  
  # Retrieve IF.vals (uncentered) from result
  # result$IF.vals.0 and result$IF.vals.1 are centered; we need uncentered.
  # Reconstruct: uncentered = centered + surv (point estimate)
  surv.0 <- result$surv.df$surv[result$surv.df$trt == 0 & result$surv.df$time > 0]
  surv.1 <- result$surv.df$surv[result$surv.df$trt == 1 & result$surv.df$time > 0]
  
  IF.vals.0 <- sweep(result$IF.vals.0, 2, surv.0, "+")  # uncentered
  IF.vals.1 <- sweep(result$IF.vals.1, 2, surv.1, "+")  # uncentered
  
  if (verbose) message("Sandwich correction: computing A21 by finite differences...")
  message("  (This requires 2 x p_alpha repasses through the posterior/IF pipeline)")
  
  # For A21, we need to perturb alpha and see how phi_obs changes.
  # We do this separately for each arm by calling compute_A21_numderiv.
  # The function returns perturbed phi for BOTH arms simultaneously to avoid
  # redundant computation.
  
  p_alpha <- length(fN_fit$alpha_hat)
  K       <- length(fit.times)
  A21_0   <- matrix(0, K, p_alpha)
  A21_1   <- matrix(0, K, p_alpha)
  
  for (j in seq_len(p_alpha)) {
    if (verbose) message(sprintf("  A21: direction %d / %d", j, p_alpha))
    
    alpha_up   <- fN_fit$alpha_hat; alpha_up[j]   <- alpha_up[j]   + epsilon
    alpha_down <- fN_fit$alpha_hat; alpha_down[j] <- alpha_down[j] - epsilon
    
    res_up   <- .recompute_phi_perturbed(
      fN_fit, alpha_up, time, event, time_N, DeltaN_eff, treat,
      W, admin_t, nuis$event.pred.0, nuis$event.pred.1,
      nuis$censY.pred.0, nuis$censY.pred.1,
      nuis$censN.pred.0, nuis$censN.pred.1,
      nuis$prop.pred, eval.times, fit.times, nuis$id, treat, W_names
    )
    
    res_down <- .recompute_phi_perturbed(
      fN_fit, alpha_down, time, event, time_N, DeltaN_eff, treat,
      W, admin_t, nuis$event.pred.0, nuis$event.pred.1,
      nuis$censY.pred.0, nuis$censY.pred.1,
      nuis$censN.pred.0, nuis$censN.pred.1,
      nuis$prop.pred, eval.times, fit.times, nuis$id, treat, W_names
    )
    
    dphi_0 <- (res_up$phi_0 - res_down$phi_0) / (2 * epsilon)   # n x K
    dphi_1 <- (res_up$phi_1 - res_down$phi_1) / (2 * epsilon)
    
    A21_0[, j] <- -colMeans(dphi_0, na.rm = TRUE)
    A21_1[, j] <- -colMeans(dphi_1, na.rm = TRUE)
  }
  
  if (verbose) message("Sandwich correction: applying correction to IF matrices...")
  corrected <- sandwich_IF_correction(
    IF.vals.0, IF.vals.1, A21_0, A21_1, A11_inv, S_alpha
  )
  
  if (verbose) message("Sandwich correction: computing log-RR with corrected variance...")
  rr_df <- sandwich_logRR_var(
    surv.0      = surv.0,
    surv.1      = surv.1,
    IF.vals.center.0.corrected = corrected$IF.vals.center.0.corrected,
    IF.vals.center.1.corrected = corrected$IF.vals.center.1.corrected,
    fit.times   = fit.times,
    conf.level  = 0.95
  )
  
  list(
    risk.ratio.sandwich    = rr_df,
    IF.vals.center.0.sw    = corrected$IF.vals.center.0.corrected,
    IF.vals.center.1.sw    = corrected$IF.vals.center.1.corrected,
    A21_0                  = A21_0,
    A21_1                  = A21_1,
    A11_inv                = A11_inv,
    S_alpha                = S_alpha
  )
}


# =============================================================================
# Example usage (after running AIPWCC_CFsurvival):
#
#  # Step 1: fit marginal f_N model on full data
#  fN_fit <- fit_fN_parametric(
#    time_N     = data$N1,
#    DeltaN_eff = DeltaN_eff,   # pre-computed in AIPWCC_CFsurvival
#    A          = data$A,
#    W          = data.frame(X = data$X),
#    model      = "cox"
#  )
#
#  # Step 2: run main estimator (using fN_fit$survest_N in the fold loop)
#  est <- AIPWCC_CFsurvival(...)
#
#  # Step 3: compute sandwich variance correction
#  sw <- compute_sandwich_variance(
#    result     = est,
#    fN_fit     = fN_fit,
#    time       = data$Y,
#    event      = data$DeltaY,
#    time_N     = data$N1,
#    DeltaN_eff = DeltaN_eff,
#    treat      = data$A,
#    W          = data.frame(X = data$X),
#    admin_t    = admin_t_vec,
#    nuis       = est$nuisance,
#    fit.times  = est$fit.times,
#    eval.times = nuis_opts$eval.times,
#    epsilon    = 1e-4,
#    verbose    = TRUE
#  )
#
#  # Sandwich-corrected log-RR confidence interval:
#  sw$risk.ratio.sandwich
# =============================================================================