# =============================================================================
# AIPWCC_reverse_factorization.r
#
# Implements the AIPWCC estimator for counterfactual survival using the
# REVERSE factorization of nuisances:
#
#   f_{Y,N | A,X} = f_{Y | A,X} x f_{N | A,X, tilde_Y, Delta_Y}
#
# where the second factor decomposes as:
#
#   f_{N | A,X,tilde_Y,Delta_Y} = Delta_Y  * f^1(N | A,X,Y  = tilde_Y)
#                                + (1-Delta_Y) * f^0(N | A,X,Y  > tilde_Y)
#
# S_{0Y}(y | a, x, n) is recovered by Bayes-rule inversion ON THE HAZARD
# SCALE, guaranteeing monotonicity in y by construction. Smoothness in n
# is controlled by regularizing the likelihood ratio r(y,n,a,x).
#
# Variance: plain empirical variance of the centered influence function
# scores — no sandwich correction for the N model needed, because both
# nuisance components are estimated from i.i.d. observed data with no
# augmentation step.
#
# =============================================================================
#
# NOTATION GUIDE
#   tilde_Y  / time    : observed Y time  (min of true Y and censoring time)
#   Delta_Y  / event   : 1 = Y observed, 0 = censored
#   tilde_N  / time_N  : observed N time
#   Delta_N  / event_N : 1 = N observed, 0 = censored
#   a0                 : target treatment arm (scalar, 0 or 1)
#   fit.times          : grid of K time points at which to estimate S(t|a0)
#   eval.times         : grid of l N-values at which to evaluate r(y,n,a,x)
#
# WORKFLOW
#   Step 1  Fit S_{Y|A,X}        -> marginal survival, all obs
#   Step 2  Fit f^1_{N|A,X,Y}   -> N model in Delta_Y=1 stratum
#   Step 3  Fit f^0_{N|A,X,Y>y} -> N model in Delta_Y=0 stratum
#   Step 4  Compute r(y,n,a,x)  -> likelihood ratio, regularised + clipped
#   Step 5  Build S_{0Y}(y|a,x,n) -> product-limit form, monotone by construction
#   Step 6  Compute phi^F and Q_0 -> full-data IF and its conditional expectation
#   Step 7  Assemble AIPWCC IF   -> three-term EIF for observed data
#   Step 8  Point estimate + CI  -> colMeans of IF + empirical variance
#
# =============================================================================


# ===========================================================================
# MASTER FUNCTION
# ===========================================================================

AIPWCC_reverse <- function(
    time,           # n-vector: observed Y times
    event,          # n-vector: Delta_Y (1=event, 0=censored)
    time_N,         # n-vector: observed N times
    event_N,        # n-vector: Delta_N (1=event, 0=censored)
    admin_t,        # scalar or n-vector: administrative censoring time
    treat,          # n-vector: treatment indicator (0/1)
    W,              # n x p data.frame of baseline covariates
    fit.times,      # K-vector: times at which to estimate S(t|a0)
    eval.times,     # l-vector: N grid for r(y,n,a,x) evaluation
    fit.treat       = c(0, 1),
    nuisance.options = list(),   # passed to .reverse.nuisance.options()
    r.clip          = c(0.05, 20),  # [lower, upper] clip for likelihood ratio r
    r.smooth.bw     = NULL,         # bandwidth for kernel smoothing r in n;
    # NULL = use model smoothness only
    conf.level      = 0.95,
    verbose         = TRUE
) {
  
  # r.clip=c(0.05, 20)
  # r.smooth.bw = 5
  # time = data$Y;
  # event = data$DeltaY;
  # time_N = data$N1;
  # event_N = data$DeltaN1;
  # admin_t = round(365*L$admin_date);
  # treat = data$A;
  # W = data.frame(X=data$X);
  # id=1:nrow(data);
  # fit.times=sort(unique(c(data$Y[data$Y>0])));
  # fit.treat=c(0,1);
  # conf.band=FALSE;
  # conf.level=.95;
  # contrasts = c("risk.ratio");
  # verbose=TRUE;
  # nuisance.options = list(
  #   cross.fit = TRUE,
  #   V = 5,
  #   folds = NULL,
  #   Y.df = 3,
  #   event.SL.library = c("survSL.km",
  #                        "survSL.coxph",
  #                        "survSL.expreg"),
  #   cens.SL.library = c("survSL.km",
  #                       "survSL.coxph",
  #                       "survSL.expreg"),
  #   prop.SL.library = NULL,
  #   prop.trunc = 0,
  #   prop.pred = prop.fit,
  #   survSL.control = list(initWeightAlg = "survSL.km", verbose = FALSE),
  #   survSL.cvControl = list(V=5)
  # )
  # eval.times = round(seq(min(data$N1[data$N1>0]), max(data$N1[data$N1>0]), length.out=100))
  
  nuis <- do.call(".reverse.nuisance.options", nuisance.options)
  
  n       <- length(time)
  W       <- as.data.frame(W)
  W_names <- colnames(W)
  fit.times  <- sort(fit.times[fit.times > 0])
  eval.times <- sort(eval.times[eval.times > 0])
  K <- length(fit.times)
  l <- length(eval.times)
  
  if (length(admin_t) == 1L) admin_t <- rep(admin_t, n)
  
  # DeltaN_eff: administrative censoring of N is treated as observed
  # (N is known to be >= admin_t so there is no missing-data uncertainty)
  DeltaN <- as.integer(event_N == 1L)
  DeltaN_eff <- as.integer(event_N == 1L | time_N >= admin_t)
  
  # ------------------------------------------------------------------
  # Cross-fitting folds (stratified on Delta_Y)
  # ------------------------------------------------------------------
  if (!nuis$cross.fit) {
    nuis$V     <- 1L
    nuis$folds <- rep(1L, n)
  }
  if (is.null(nuis$folds)) {
    nuis$folds <- integer(n)
    for (ev in c(0L, 1L)) {
      idx <- which(event == ev)
      nuis$folds[idx] <- sample(rep(seq_len(nuis$V), length.out = length(idx)))
    }
  }
  
  # Storage for cross-fitted nuisances (all in original-data row order)
  # SY[i, k]      = S_{Y|A,X}(fit.times[k] | A_i, X_i)  under a0
  # SY_alt[i, k]  = same under 1-a0
  # GY[i, k]      = G_C^Y(fit.times[k] | A_i, X_i)      under a0
  # GN[i, k]      = G_C^N(fit.times[k] | A_i, X_i)
  # pi_hat[i]     = P(A=a0 | X_i)
  # R0[i, l, k]   = r(fit.times[k], eval.times[l] | A=a0, X_i)
  #                 NOTE: this is an n x l x K array — see Step 4
  
  cf <- list(
    SY_a0  = matrix(NA_real_, n, K),   # S_Y(t | a0, X)
    SY_alt = matrix(NA_real_, n, K),   # S_Y(t | 1-a0, X)
    GY_a0  = matrix(NA_real_, n, K),   # G_C^Y(t | a0, X)
    GY_alt = matrix(NA_real_, n, K),
    GN     = matrix(NA_real_, n, K),   # G_C^N(t | A, X) — arm-specific
    pi_hat = rep(NA_real_, n),
    # r is stored as n x l x K; initialized after we know n, l, K
    R0      = array(NA_real_, dim = c(n, l, K)),
    f10     = matrix(NA_real_, n, l),
    f00     = matrix(NA_real_, n, l),
    R1      = array(NA_real_, dim = c(n, l, K)),
    f11     = matrix(NA_real_, n, l),
    f01     = matrix(NA_real_, n, l)
  )
  
  # ===========================================================================
  # CROSS-FITTING LOOP
  # ===========================================================================
  
  for (v in seq_len(nuis$V)) {
    
    if (verbose) message(sprintf("Fold %d / %d", v, nuis$V))
    
    tr <- which(nuis$folds != v)
    te <- which(nuis$folds == v)
    
    dat_tr <- .make_df(time, event, time_N, DeltaN, treat, admin_t, W, tr)
    dat_te <- .make_df(time, event, time_N, DeltaN, treat, admin_t, W, te)
    n_te   <- nrow(dat_te)
    W_tr   <- dat_tr[, W_names, drop = FALSE]
    W_te   <- dat_te[, W_names, drop = FALSE]
    
    # ---------------------------------------------------------------
    # STEP 1: S_{Y | A, X}  — marginal survival of Y given treatment
    #         and covariates, NO conditioning on N.
    #         Fit on training fold, predict for test fold under both arms.
    #
    # MONOTONE by default
    # ---------------------------------------------------------------
    if (verbose) message("  [1] S_{Y|A,X}")
    
    SY_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr$Y,
      event            = dat_tr$DeltaY,
      X                = cbind(A = dat_tr$A, W_tr),
      newX             = as.data.frame(rbind(
        cbind(A = 0L, W_te),
        cbind(A = 1L, W_te)
      )),
      new.times        = fit.times,
      event.SL.library = nuis$event.SL.library,
      cens.SL.library  = nuis$cens.SL.library,
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl
    )
    
    SY_pred <- SY_fit$event.SL.predict   # (2 n_te) x K
    SY_0    <- SY_pred[seq_len(n_te),          , drop = FALSE]   # A=0 block
    SY_1    <- SY_pred[seq_len(n_te) + n_te,   , drop = FALSE]   # A=1 block
    SY_0    <- pmax(SY_0, 1e-8)
    SY_1    <- pmax(SY_1, 1e-8)
    
    # store in original-data order: rows indexed by te
    # a0 = 1 means SY_a0 gets the A=1 block; handled at assembly below
    cf$SY_a0[te,  ] <- SY_1    # will be swapped per-arm at assembly
    cf$SY_alt[te, ] <- SY_0
    
    # ---------------------------------------------------------------
    # STEP 2: G_C^Y  — censoring survival for Y process
    # ---------------------------------------------------------------
    if (verbose) message("  [2] G_C^Y")
    
    rcY <- as.integer(dat_tr$DeltaY == 0L & dat_tr$Y < dat_tr$admin_t)
    
    GY_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr$Y,
      event            = rcY,
      X                = cbind(A = dat_tr$A, W_tr),
      newX             = as.data.frame(rbind(
        cbind(A = 0L, W_te),
        cbind(A = 1L, W_te)
      )),
      new.times        = fit.times,
      event.SL.library = nuis$cens.SL.library,
      cens.SL.library  = nuis$cens.SL.library,
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl
    )
    
    GY_pred       <- GY_fit$event.SL.predict
    cf$GY_a0[te,] <- pmax(GY_pred[seq_len(n_te) + n_te, , drop=FALSE], 1e-4)  # A=1
    cf$GY_alt[te,]<- pmax(GY_pred[seq_len(n_te),         , drop=FALSE], 1e-4)  # A=0
    
    # ---------------------------------------------------------------
    # STEP 3: G_C^N  — censoring survival for N process
    # ---------------------------------------------------------------
    if (verbose) message("  [3] G_C^N")
    
    rcN <- as.integer(dat_tr$DeltaN == 0L & dat_tr$N < dat_tr$admin_t)
    
    GN_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr$N,
      event            = rcN,
      X                = cbind(A = dat_tr$A, W_tr),
      newX             = as.data.frame(cbind(A = dat_te$A, W_te)),
      new.times        = fit.times,
      event.SL.library = nuis$cens.SL.library,
      cens.SL.library  = nuis$cens.SL.library,
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl
    )
    
    cf$GN[te, ] <- pmax(GN_fit$event.SL.predict, 1e-4)
    
    # ---------------------------------------------------------------
    # STEP 4: Propensity score
    # ---------------------------------------------------------------
    if (verbose) message("  [4] Propensity score")
    
    if(!is.null(nuis$prop.pred)){ 
    cf$pi_hat[te]<-nuis$prop.pred[te]
    }else{
    cf$pi_hat[te] <- .estimate.propensity(
      A        = dat_tr$A,
      W        = W_tr,
      newW     = W_te,
      SL.library = nuis$prop.SL.library,
      fit.treat  = fit.treat,
      prop.trunc = nuis$prop.trunc,
      save.fit   = FALSE,
      verbose    = FALSE
    )$prop.pred
    }
    # ---------------------------------------------------------------
    # STEP 5: Likelihood ratio r(y, n, a, x)
    #
    # r(y, n, a, x) = f^1_N(n | a, x, Y=y) / f^0_N(n | a, x, Y>y)
    #
    # HOW: We fit two separate survival models for N:
    #   Model f^1: fit on rows with DeltaY=1 (Y observed), using tilde_Y
    #              as an additional covariate. This estimates the density
    #              of N conditional on (A, X, Y=tilde_Y).
    #   Model f^0: fit on rows with DeltaY=0 (Y censored), using tilde_Y
    #              as an additional covariate. This estimates the density
    #              of N conditional on (A, X, Y > tilde_Y).
    #
    # WHY TWO MODELS: The decomposition
    #   f_{N|A,X,tY,dY} = dY * f^1 + (1-dY) * f^0
    # means f^1 and f^0 are densities in DIFFERENT subpopulations.
    # Fitting a single model with dY as a covariate would conflate them.
    #
    # SMOOTHNESS IN n: Both models include N as a prediction target
    # (it is the time variable in the survival model). The predicted
    # hazard increments dLambda^1 and dLambda^0 are smooth functions
    # of n by construction when using Cox with spline(tilde_Y) or a
    # parametric AFT model. We further apply:
    #   (a) Kernel smoothing of r in n if r.smooth.bw is supplied
    #   (b) Hard clipping of r to r.clip = [lower, upper]
    # to prevent a 1-day change in N from producing wild variation.
    #
    # MONOTONICITY IN y: r(y, n, a, x) is evaluated at each fit.times[k]
    # and multiplied INTO the marginal hazard lambda_{Y|A,X}(y|a,x).
    # The product-limit form in Step 6 guarantees monotone S_{0Y} in y
    # regardless of how r varies with y.
    # ---------------------------------------------------------------
    if (verbose) message("  [5] Likelihood ratio r(y, n, a, x)")
    
    R_0 <- .compute_r_2(
      dat_tr    = dat_tr,
      dat_te    = dat_te,
      a0 = 0, #FIX 3
      W_names   = W_names,
      fit.times = fit.times,
      eval.times = eval.times,
      r.clip    = r.clip,
      r.smooth.bw = r.smooth.bw,
      nuis      = nuis,
      verbose   = verbose
    ) #n_te x k x l
    
    R_1 <- .compute_r_2(
      dat_tr    = dat_tr,
      dat_te    = dat_te,
      a0 = 1, #FIX 3
      W_names   = W_names,
      fit.times = fit.times,
      eval.times = eval.times,
      r.clip    = r.clip,
      r.smooth.bw = r.smooth.bw,
      nuis      = nuis,
      verbose   = verbose
    ) #n_te x k x l
    
    cf$R0[te, , ] <- R_0$R
    cf$f00[te, ] <- R_0$f0
    cf$f10[te, ] <- R_0$f1
    
    cf$R1[te, , ] <- R_1$R
    cf$f01[te, ] <- R_1$f0
    cf$f11[te, ] <- R_1$f1
    
  } # end cross-fitting loop
  
  # ===========================================================================
  # STEP 6: Build S_{0Y}(y | a, x, n)  via product-limit Bayes inversion
  #
  # For each participant i, arm a0, N-grid point n_j, and time t_k:
  #
  #   S_{0Y}(t_k | a0, X_i, n_j)
  #     = prod_{m=1}^{k} max(0, 1 - dLambda_{Y|A,X}(t_m | a0, X_i) * r(t_m, n_j, a0, X_i))
  #
  # where dLambda_{Y|A,X}(t_m | a0, X_i) = 1 - S_{Y|A,X}(t_m) / S_{Y|A,X}(t_{m-1})
  #                                        (discrete hazard increment)
  #
  
  # MONOTONICITY GUARANTEE:
  #   Each factor max(0, 1 - h_m) is in [0,1] because:
  #     - dLambda >= 0 (survival model always gives non-increasing S)
  #     - r >= r.clip[1] > 0 by clipping in Step 5
  #     - the max(0, ...) clips the rare case where r is large enough
  #       to push the conditional hazard above 1
  #   A cumulative product of [0,1] factors is non-increasing in k.
  #   => S_{0Y}(t_k | a0, X_i, n_j) is monotone non-increasing in k BY CONSTRUCTION.
  #   No isotonic regression or post-hoc fix is needed.
  #
  # SMOOTHNESS IN n:
  #   r(t_m, n_j, a0, X_i) is smooth in n_j by the survival model choice
  #   (and optional kernel smoothing). The product-limit form inherits this
  #   smoothness: a small change in n_j produces a small change in each factor
  #   and hence a small change in the product.
  # ===========================================================================
  
  if (verbose) message("Building S_{0Y}(y|a,x,n) via product-limit inversion...")
  
  # dLambda_{Y|A,X}(t_k | a0, X_i): n x K matrix of discrete hazard increments
  # dLambda[i,k] = 1 - SY_a0[i,k] / SY_a0[i,k-1], with SY_a0[i,0] := 1
  S1Y_list <- .build_S0Y(
    SY    = cf$SY_a0,     # n x K: S_{Y|A=a0,X}
    R        = cf$R1,         # n x l x K: r(t_k, n_j, a0, X_i)
    fit.times  = fit.times,
    eval.times = eval.times
  )
  
  S0Y_list <- .build_S0Y(
    SY    = cf$SY_alt,     # n x K: S_{Y|A=a0,X}
    R        = cf$R0,         # n x l x K: r(t_k, n_j, a0, X_i)
    fit.times  = fit.times,
    eval.times = eval.times
  )
  # S0Y_list$S0Y: n x l x K array
  # S0Y_list$dLambda: n x K matrix of marginal hazard increments
  
  # ===========================================================================
  # STEP 7: Compute phi^F and Q_0, then assemble AIPWCC influence function
  #
  # phi^F(tilde_y, delta_y, a, x, n) for target arm a0:
  #
  #   phi^F_i(n) = S_{0Y}(t0 | a0, X_i, n)
  #              - [A_i = a0]/pi(a0|X_i) * S_{0Y}(t0|a0,X_i,n)
  #                * { I(tY_i<=t0, dY_i=1) / [S_{0Y}(tY_i|a0,X_i,n) * G_C^Y(tY_i|a0,X_i)]
  #                    - int_0^{min(t0,tY_i)} dLambda_{0Y}(u|a0,X_i,n) / [S_{0Y}(u|a0,X_i,n)*G_C^Y(u|a0,X_i)] }
  #
  # Q_0(tN_i | tY_i, dY_i, A_i, X_i) = E[phi^F(N) | N >= tN_i, tY_i, dY_i, A_i, X_i]
  #   = sum_{j: n_j >= tN_i} phi^F_i(n_j) * f_{N|...}(n_j | ...) / P(N >= tN_i | ...)
  #   where f_{N|...} is the posterior mixing density (Delta_Y=1: f^1, else f^0)
  #
  # AIPWCC IF (three terms):
  #   term1 = [DeltaN_i / G_C^N(tN_i | A_i,X_i)] * phi^F_i(tN_i)
  #   term2 = [(1-DeltaN_i) / G_C^N(tN_i | A_i,X_i)] * Q_0(tN_i | ...)
  #   term3 = -int_0^{tN_i} Q_0(u | ...) / G_C^N(u | A_i,X_i) * dLambda_C^N(u | A_i,X_i) du
  # ===========================================================================
  
  if (verbose) message("Computing AIPWCC influence functions...")
  
  result <- list(fit.times = fit.times, surv.df = data.frame())
  
  for (a0 in fit.treat) {
    
    if (verbose) message(sprintf("  Arm a0 = %d", a0))
    
    # For arm a0: SY_a0 is S_{Y|A=a0,X}, pi_a0 = P(A=a0|X)
    # The cross-fitting stored SY_a0 as the A=1 block; swap if a0=0
    SY_use <- if (a0 == 1L) cf$SY_a0 else cf$SY_alt
    pi_use <- if (a0 == 1L) cf$pi_hat else 1 - cf$pi_hat
    GY_use <- if (a0 == 1L) cf$GY_a0 else cf$GY_alt
    SY_N <- if(a0==1L) S1Y_list$S0Y else S0Y_list$S0Y
    dLambaY_N <- if(a0==1L) S1Y_list$dLambda_cond else S0Y_list$dLambda_cond
    f1 <- if(a0==1L) cf$f11 else cf$f10
    f0 <- if(a0==1L) cf$f01 else cf$f00
    
    IF_mat <- .compute_AIPWCC_IF(
      time       = time,
      event      = event,
      time_N     = time_N,
      DeltaN     = DeltaN_eff,
      A          = treat,
      a0         = a0,
      admin_t    = admin_t,
      fit.times  = fit.times,
      eval.times = eval.times,
      SY         = SY_use,          # n x K: S_{Y|A=a0,X}
      S0Y        = SY_N,    # n x l x K: S_{0Y}(t|a0,X,n)
      dLambda0Y  = dLambaY_N, # n x l x K: dLambda_{0Y}(t|a0,X,n)
      GY         = GY_use,          # n x K: G_C^Y(t|a0,X)
      GN         = cf$GN,           # n x K: G_C^N(t|A,X)
      pi_hat     = pi_use,          # n: P(A=a0|X)
      #fN_obs     = cf$R,            # n x l x K: r used to build posterior weights
      # NOTE: the posterior f_{N|tY,dY,A,X} is reconstructed inside
      # .compute_AIPWCC_IF from the two component models' predictions
      fN1_pred   = f1,     # n x l: f^1(n_j | A_i, X_i, Y=tY_i)
      fN0_pred   = f0,     # n x l: f^0(n_j | A_i, X_i, Y>tY_i)
      verbose    = TRUE
    )
    
    surv_hat <- colMeans(IF_mat, na.rm = TRUE)
    IF_center <- sweep(IF_mat, 2, surv_hat, "-")
    
    # ---------------------------------------------------------------
    # STEP 8: Variance — plain empirical variance of centered IF scores
    #
    # WHY NO SANDWICH CORRECTION NEEDED:
    #   Both nuisance components (S_{Y|A,X} and f_{N|A,X,tY,dY}) are
    #   estimated by solving i.i.d. estimating equations on the OBSERVED
    #   data, with no augmentation step and no data-dependent weights.
    #   Under cross-fitting, the estimation error in both nuisances
    #   contributes only at rate o_p(n^{-1/2}) to the EIF, so the
    #   empirical variance of phi^obs_i is a valid variance estimator.
    #   This is the key advantage of the reverse factorization.
    # ---------------------------------------------------------------
    var_hat <- apply(IF_center, 2, var, na.rm = TRUE) / n
    se_hat  <- sqrt(var_hat)
    
    z <- qnorm((1 + conf.level) / 2)
    
    # Logit-scale CI (better coverage near 0 and 1)
    logit_s  <- log(surv_hat / (1 - surv_hat))
    se_logit <- se_hat / (surv_hat * (1 - surv_hat))
    lo_logit <- logit_s - z * se_logit
    hi_logit <- logit_s + z * se_logit
    
    df <- data.frame(
      time         = c(0, fit.times),
      trt          = a0,
      surv         = c(1, pmin(1, pmax(0, surv_hat))),
      se           = c(0, se_hat),
      ptwise.lower = c(1, plogis(lo_logit)),
      ptwise.upper = c(1, plogis(hi_logit))
    )
    
    result$surv.df        <- rbind(result$surv.df, df)
    result[[paste0("IF.vals.", a0)]] <- IF_center
    
  }
  
  result$nuisance <- cf
  result
}


# ===========================================================================
# HELPER: .compute_r()
#
# Fits f^1 and f^0 models for N and returns the likelihood ratio
# r(y, n, a, x) = f^1(n | a, x, Y=y) / f^0(n | a, x, Y>y)
# evaluated on the grid (fit.times x eval.times) for each test participant.
#
# INPUTS
#   dat_tr, dat_te : training and test fold data frames
#   fit.times      : K-vector, y-grid
#   eval.times     : l-vector, n-grid
#   r.clip         : c(lower, upper) hard clip for r values
#   r.smooth.bw    : optional kernel bandwidth for smoothing r in n
#
# OUTPUT
#   Array of dimension (n_te x l x K):
#     R[i, j, k] = r(fit.times[k], eval.times[j], A_i=a0, X_i)
#
# SMOOTHNESS IN n (design notes):
#   f^1 and f^0 are survival models for N. The predicted SURVIVAL functions
#   S^1_N(n | a, x, y) and S^0_N(n | a, x, y) are smooth in n when using
#   Cox with spline(tilde_Y) or a parametric model. The density increments
#   dLambda^1 and dLambda^0 are differences of the survival function on the
#   eval.times grid, so they inherit the smoothness.
#
#   Additional smoothing via kernel (r.smooth.bw): if supplied, for each
#   participant i and time t_k, we smooth r(t_k, ., a0, X_i) over n using
#   a Gaussian kernel with bandwidth r.smooth.bw. This prevents a 1-day
#   change in n from producing discontinuous jumps in r.
#
#   Hard clip (r.clip): ensures r stays in [lower, upper] regardless of
#   model predictions. Prevents extreme S_{0Y} values from a single outlier
#   N observation.
# ===========================================================================

.compute_r_2 <- function(
    dat_tr, dat_te, W_names,
    a0,
    fit.times, eval.times, cuts=NULL,
    r.clip, r.smooth.bw, nuis,
    verbose = FALSE
) {
  
  #a0=0
  n_te <- nrow(dat_te)
  K    <- length(fit.times)
  l    <- length(eval.times)
  W_te <- dat_te[, W_names, drop = FALSE]
  
  # ---------------------------------------------------------------
  # Fit f^1: survival model for N in the Delta_Y = 1 stratum
  # Covariates: (A, X, tilde_Y) — tilde_Y enters as a covariate
  # because we condition on Y = tilde_Y.
  # ---------------------------------------------------------------
  tr1 <- dat_tr[dat_tr$DeltaY == 1L, ]
  
  if (nrow(tr1) < 10L) {
    warning("Fewer than 10 DeltaY=1 obs in training fold; r set to 1.")
    f1_obs <- matrix(1/l, nrow = n_te, ncol = l)
    f0_obs <- matrix(1/l, nrow = n_te, ncol = l)
    R_arr  <- array(1,    dim  = c(n_te, l, K))
    return(list(R = R_arr, f1 = f1_obs, f0 = f0_obs))
  }
  
  # newX for f^1: for each test participant i and each fit.time t_k,
  # we want f^1_N(n | A_i, X_i, Y = t_k). So we expand the test set
  # over fit.times: n_te * K rows with Y = fit.times[k] per block.
  newX_f1_grid <- do.call(rbind, lapply(fit.times, function(tk) {
    cbind(A = a0, Y = tk, W_te) #A=a0 not dat_te$A
  }))   # (n_te * K) x (2 + p), column named "Y" to match coxph formula
  
  # newX for observed Y_i — used for f^1(n | a0, X_i, Y = tY_i) weights
  newX_f1_obs <- cbind(A = a0, Y = dat_te$Y, W_te)   # n_te rows
  
  fit1 <- survival::coxph(
    survival::Surv(N, DeltaN) ~ A + splines::ns(Y, df = nuis$Y.df) + .,
    data = cbind(
      tr1[, c("N", "DeltaN", "A", "Y"), drop = FALSE],
      tr1[, W_names, drop = FALSE]
    ),
    x = TRUE
  )
  # Predict survival of N on eval.times grid for each row of newX_f1
  # Result: (n_te * K) x l matrix — row indexed by (i, k), col by j
  SN1_mat <- .predict_cox_surv_v2(fit1, newX_f1_grid, eval.times)
  
  # Predict survival of N on eval.times grid for data (observed Y)
  # Result: (n_te x l) matrix
  surv_d1_yi <- .predict_cox_surv_v2(fit1, newX_f1_obs, eval.times)
  
  # ---------------------------------------------------------------
  # Fit f^0: survival model for N in the Delta_Y = 0 stratum
  # Covariates: (A, X, tilde_Y) — tilde_Y enters as a covariate
  # because we condition on Y > tilde_Y (i.e. Y was not yet observed).
  # ---------------------------------------------------------------
  make_pooled_coarse <- function(df, cuts) {
    # cuts define interval endpoints: e.g. c(0, 2, 4, 6, ...)
    
    n <- nrow(df)
    J <- length(cuts) - 1
    
    # Expand
    out <- df[rep(1:n, each = J), ]
    
    out$tstart <- rep(cuts[-length(cuts)], times = n)
    out$tstop  <- rep(cuts[-1], times = n)
    
    # At-risk: alive at start of interval
    out$at_risk <- as.integer(out$Y > out$tstart)
    
    # Event occurs in interval
    out$event <- as.integer(
      out$DeltaY == 1 &
        out$Y > out$tstart &
        out$Y <= out$tstop
    )
    
    # Drop intervals after failure/censoring
    out <- out[out$at_risk == 1, ]
    
    return(out)
  }
  
  if(is.null(cuts)){cuts <- fit.times}
  
  pool_tr<-make_pooled_coarse(dat_tr, cuts=unique(c(0, cuts)))
  pool_te<-make_pooled_coarse(dat_te, cuts=unique(c(0, cuts)))
  
  # newX for observed Y_i — used for f^1(n | a0, X_i, Y = tY_i) weights
  newX_f1_obs <- cbind(A = a0, W_te)   # n_te rows
  
  gridpoints_y <- unique(cuts)
  
  S0N1 <- vector(mode="list", length=length(gridpoints_y))
  surv_d0_yi <- matrix(NA, nrow=n_te, ncol=l)

  for(i in 1:c(length(gridpoints_y))){
    #print(i)
    k=gridpoints_y[i]
    
    tmp_haz <- coxph(survival::Surv(N, DeltaN) ~ A + .,
          pool_tr %>% 
            filter(tstop == k) %>% 
            dplyr::select(c("N", "DeltaN", "A", W_names)),
          x = TRUE)
    
    S0N1[[i]]<-.predict_cox_surv_v2(tmp_haz, 
                            newX = newX_f1_obs, 
                            eval.times = eval.times) #each is a n_te x l matrix
   
    #FLAG
    lower = ifelse(i==1, 0, gridpoints_y[(i-1)])
    upper = k
    indx<-which(dat_te$Y > lower & dat_te$Y <= upper)
    
    if(length(indx)==0){
      next
    }else{
    surv_d0_yi[indx,]<-.predict_cox_surv_v2(tmp_haz, newX = newX_f1_obs[indx,], eval.times = eval.times)
    }
    
    rm(tmp_haz) 
  }
  
  # list of length n_te, each element K (fit.times) x l (eval.times)
  SN0_mat <- do.call(rbind, S0N1)
  
  # ---------------------------------------------------------------
  # Convert survival functions to density increments on eval.times grid
  # dF^1[row, j] = S^1_N(n_{j-1}) - S^1_N(n_j)   (probability mass)
  # with S^1_N(n_0) := 1 at n_0 = 0
  # ---------------------------------------------------------------
  dF1 <- cbind(1, SN1_mat)[, seq_len(l), drop=FALSE] - SN1_mat
  # Shape: (n_te * K) x l; each row sums to 1 - S^1_N(last eval.time)
  dF1 <- pmax(dF1, 0)
  
  dF0 <- cbind(1, SN0_mat)[, seq_len(l), drop=FALSE] - SN0_mat
  dF0 <- pmax(dF0, 0)
  
  # ---------------------------------------------------------------
  # Compute r(t_k, n_j, A_i, X_i) = dF^1[i,k,j] / dF^0[i,k,j]
  # Array indexing: rows of SN1_mat are (i, k) pairs, k cycling fastest
  # i.e. row (i-1)*K + k corresponds to participant i at time fit.times[k]
  # ---------------------------------------------------------------
  R_arr <- array(NA_real_, dim = c(n_te, l, K))
  
  for (k in seq_len(K)) {
    # print(k)
    # rows in dF1/dF0 for time index k: k, k+K, k+2K, ..., k+(n_te-1)*K
    # because newX_f1 was built by lapply over fit.times -> k block is rows
    # ((k-1)*n_te + 1) : (k * n_te)
    row_idx <- ((k - 1L) * n_te + 1L) : (k * n_te)   # n_te rows
    
    num <- dF1[row_idx, , drop = FALSE]   # n_te x l
    den <- dF0[row_idx, , drop = FALSE]   # n_te x l
    
    # Ratio: add small epsilon to denominator to avoid 0/0
    r_mat <- num/den
    r_mat[is.nan(r_mat)] <- 1     #enforce 0/0 = 1  # n_te x l
    r_mat[is.infinite(r_mat)] <- 0
    
    # ---------------------------------------------------------------
    # SMOOTHNESS: optional kernel smoothing of r in n (across columns)
    #
    # For each participant i and time t_k, r(t_k, ., A_i, X_i) is a
    # step function in n on the eval.times grid. Kernel smoothing with
    # bandwidth r.smooth.bw replaces each value r[i,j] with a weighted
    # average of nearby r[i,j'] values, preventing abrupt jumps.
    # ---------------------------------------------------------------
    if (!is.null(r.smooth.bw)) {
      r_mat <- .smooth_r_in_n(r_mat, eval.times, bw = r.smooth.bw)
    }
    
    # ---------------------------------------------------------------
    # SMOOTHNESS + STABILITY: hard clip r to [r.clip[1], r.clip[2]]
    #
    # Prevents extreme values of S_{0Y} arising from a single unusual
    # N observation. A value of r.clip = c(0.05, 20) means N can
    # modify the conditional hazard by at most a factor of 20 in either
    # direction relative to the marginal hazard.
    # ---------------------------------------------------------------
    r_mat <- pmin(pmax(r_mat, r.clip[1]), r.clip[2])
    
    R_arr[, , k] <- r_mat
  }
  
  # ---------------------------------------------------------------
  # Compute f^0 and f^1 values over eval.times grid
  # ---------------------------------------------------------------
  f1 <- cbind(1, surv_d1_yi)[,seq_len(l), drop=FALSE] - surv_d1_yi
  f0 <- cbind(1, surv_d0_yi)[,seq_len(l), drop=FALSE] - surv_d0_yi
  
  return(
    list(
      R=R_arr,
      f1 = f1,
      f0 = f0
    )
  )
}

.compute_r <- function(
    dat_tr, dat_te, W_names,
    fit.times, eval.times,
    r.clip, r.smooth.bw, nuis,
    verbose = FALSE
) {
  
  n_te <- nrow(dat_te)
  K    <- length(fit.times)
  l    <- length(eval.times)
  W_te <- dat_te[, W_names, drop = FALSE]
  
  # ---------------------------------------------------------------
  # Fit f^1: survival model for N in the Delta_Y = 1 stratum
  # Covariates: (A, X, tilde_Y) — tilde_Y enters as a covariate
  # because we condition on Y = tilde_Y.
  # ---------------------------------------------------------------
  tr1 <- dat_tr[dat_tr$DeltaY == 1L, ]
  
  if (nrow(tr1) < 10L) {
    warning("Fewer than 10 Delta_Y=1 observations in training fold; r set to 1.")
    return(array(1, dim = c(n_te, l, K)))
  }
  
  # newX for f^1: for each test participant i and each fit.time t_k,
  # we want f^1_N(n | A_i, X_i, Y = t_k). So we expand the test set
  # over fit.times: n_te * K rows with Y = fit.times[k] per block.
  newX_f1 <- do.call(rbind, lapply(fit.times, function(tk) {
    cbind(A = dat_te$A, Y_val = tk, W_te)
  }))   # (n_te * K) x (2 + p)
  
  fit1 <- survival::coxph(
    survival::Surv(N, DeltaN) ~ A + splines::ns(Y, df = nuis$Y.df) + .,
    data    = cbind(tr1[, c("N","DeltaN","A","Y"), drop=FALSE],
                    tr1[, W_names, drop=FALSE]),
    x       = TRUE
  )
  # Predict survival of N on eval.times grid for each row of newX_f1
  # Result: (n_te * K) x l matrix — row indexed by (i, k), col by j
  SN1_mat <- .predict_cox_surv(fit1, newX_f1, eval.times,
                               Y_col = "Y_val", W_names = W_names)
  
  # ---------------------------------------------------------------
  # Fit f^0: survival model for N in the Delta_Y = 0 stratum
  # Covariates: (A, X, tilde_Y) — tilde_Y enters as a covariate
  # because we condition on Y > tilde_Y (i.e. Y was not yet observed).
  # ---------------------------------------------------------------
  tr0 <- dat_tr[dat_tr$DeltaY == 0L, ]
  
  if (nrow(tr0) < 10L) {
    warning("Fewer than 10 Delta_Y=0 observations in training fold; r set to 1.")
    return(array(1, dim = c(n_te, l, K)))
  }
  
  newX_f0 <- newX_f1   # same grid; Y_val now plays role of tilde_Y for Y>tilde_Y
  
  fit0 <- survival::coxph(
    survival::Surv(N, DeltaN) ~ A + splines::ns(Y, df = nuis$Y.df) + .,
    data    = cbind(tr0[, c("N","DeltaN","A","Y"), drop=FALSE],
                    tr0[, W_names, drop=FALSE]),
    x       = TRUE
  )
  SN0_mat <- .predict_cox_surv(fit0, newX_f0, eval.times,
                               Y_col = "Y_val", W_names = W_names)
  
  # ---------------------------------------------------------------
  # Convert survival functions to density increments on eval.times grid
  # dF^1[row, j] = S^1_N(n_{j-1}) - S^1_N(n_j)   (probability mass)
  # with S^1_N(n_0) := 1 at n_0 = 0
  # ---------------------------------------------------------------
  dF1 <- cbind(1, SN1_mat)[, -ncol(SN1_mat) - 1 + seq_len(l), drop=FALSE] - SN1_mat
  # Simpler: dF1[, j] = S^1[, j-1] - S^1[, j], with S^1[, 0] = 1
  dF1 <- cbind(1 - SN1_mat[, 1], -diff(t(SN1_mat)))
  # Shape: (n_te * K) x l; each row sums to 1 - S^1_N(last eval.time)
  dF1 <- pmax(dF1, 0)
  
  dF0 <- cbind(1 - SN0_mat[, 1], -diff(t(SN0_mat)))
  dF0 <- pmax(dF0, 0)
  
  # ---------------------------------------------------------------
  # Compute r(t_k, n_j, A_i, X_i) = dF^1[i,k,j] / dF^0[i,k,j]
  # Array indexing: rows of SN1_mat are (i, k) pairs, k cycling fastest
  # i.e. row (i-1)*K + k corresponds to participant i at time fit.times[k]
  # ---------------------------------------------------------------
  R_arr <- array(NA_real_, dim = c(n_te, l, K))
  
  for (k in seq_len(K)) {
    # rows in dF1/dF0 for time index k: k, k+K, k+2K, ..., k+(n_te-1)*K
    # because newX_f1 was built by lapply over fit.times -> k block is rows
    # ((k-1)*n_te + 1) : (k * n_te)
    row_idx <- ((k - 1L) * n_te + 1L) : (k * n_te)   # n_te rows
    
    num <- dF1[row_idx, , drop = FALSE]   # n_te x l
    den <- dF0[row_idx, , drop = FALSE]   # n_te x l
    
    # Ratio: add small epsilon to denominator to avoid 0/0
    r_mat <- num / pmax(den, 1e-10)       # n_te x l
    
    # ---------------------------------------------------------------
    # SMOOTHNESS: optional kernel smoothing of r in n (across columns)
    #
    # For each participant i and time t_k, r(t_k, ., A_i, X_i) is a
    # step function in n on the eval.times grid. Kernel smoothing with
    # bandwidth r.smooth.bw replaces each value r[i,j] with a weighted
    # average of nearby r[i,j'] values, preventing abrupt jumps.
    # ---------------------------------------------------------------
    if (!is.null(r.smooth.bw)) {
      r_mat <- .smooth_r_in_n(r_mat, eval.times, bw = r.smooth.bw)
    }
    
    # ---------------------------------------------------------------
    # SMOOTHNESS + STABILITY: hard clip r to [r.clip[1], r.clip[2]]
    #
    # Prevents extreme values of S_{0Y} arising from a single unusual
    # N observation. A value of r.clip = c(0.05, 20) means N can
    # modify the conditional hazard by at most a factor of 20 in either
    # direction relative to the marginal hazard.
    # ---------------------------------------------------------------
    r_mat <- pmin(pmax(r_mat, r.clip[1]), r.clip[2])
    
    R_arr[, , k] <- r_mat
  }
  
  R_arr
}


# ===========================================================================
# HELPER: .build_S0Y()
#
# Constructs S_{0Y}(t_k | a0, X_i, n_j) for all (i, j, k) using the
# product-limit Bayes inversion.
#
# MONOTONICITY IN y (BY CONSTRUCTION):
#   S_{0Y}(t_k | a0, X_i, n_j)
#     = prod_{m=1}^{k} max(0, 1 - dLambda_{Y|A,X}(t_m | a0, X_i) * r(t_m, n_j, a0, X_i))
#
#   Each factor is in [0,1]:
#     - dLambda_{Y|A,X}(t_m) = 1 - S_{Y|A,X}(t_m)/S_{Y|A,X}(t_{m-1}) >= 0
#       (guaranteed by survival model)
#     - r >= r.clip[1] > 0  (guaranteed by clipping in .compute_r)
#     - max(0, ...) handles the rare case where r is large and the
#       conditional hazard would exceed 1
#   A cumulative product of factors in [0,1] is non-increasing in k.
#   => MONOTONICITY IS GUARANTEED WITHOUT ISOTONIC REGRESSION.
#
# SMOOTHNESS IN n:
#   r(t_m, n_j, ...) is smooth in n_j (by model choice + clipping).
#   Each factor 1 - h_m * r(t_m, n_j, ...) is therefore smooth in n_j.
#   The product of smooth functions is smooth.
#   => S_{0Y}(t_k | a0, X_i, n_j) is smooth in n_j BY CONSTRUCTION.
#
# INPUTS
#   SY     : n x K matrix, S_{Y|A=a0,X}(t_k | X_i)
#   R      : n x l x K array, r(t_k, n_j, a0, X_i)
#   fit.times, eval.times: grids
#
# OUTPUTS
#   S0Y          : n x l x K array, S_{0Y}(t_k | a0, X_i, n_j)
#   dLambda_cond : n x l x K array, discrete hazard dLambda_{0Y}(t_k|a0,X_i,n_j)
#                  used later inside phi^F computation
# ===========================================================================

.build_S0Y <- function(SY, R, fit.times, eval.times) {
  
  n <- nrow(SY)
  K <- length(fit.times)
  l <- length(eval.times)
  
  # Marginal discrete hazard increments: n x K
  # dLambda_Y[i, k] = 1 - SY[i,k] / SY[i,k-1], with SY[i,0] := 1
  SY_lag      <- cbind(1, SY[, -K, drop = FALSE])    # n x K: SY at t_{k-1}
  dLambda_Y   <- 1 - SY / pmax(SY_lag, 1e-10)        # n x K, in [0,1]
  dLambda_Y   <- pmax(dLambda_Y, 0)
  
  # Conditional discrete hazard: n x l x K
  # dLambda_{0Y}[i, j, k] = dLambda_Y[i,k] * R[i,j,k]
  # (R already clipped to r.clip in .compute_r)
  dLambda_cond <- array(NA_real_, dim = c(n, l, K))
  for (k in seq_len(K)) {
    # dLambda_Y[,k] is n x 1; R[,,k] is n x l
    # outer product via sweep: n x l
    dLambda_cond[, , k] <- sweep(R[, , k], 1, dLambda_Y[, k], FUN = "*")
  }
  # Clip conditional hazard to [0, 1] to ensure valid survival factor
  dLambda_cond <- pmin(pmax(dLambda_cond, 0), 1)
  
  # Product-limit S_{0Y}: n x l x K
  # S_{0Y}[i, j, k] = prod_{m=1}^{k} (1 - dLambda_cond[i, j, m])
  S0Y <- array(NA_real_, dim = c(n, l, K))
  log_S0Y <- array(0, dim = c(n, l, K))
  
  for (k in seq_len(K)) {
    # log(1 - dLambda_cond[,,k]), accumulated
    log_factor <- log(pmax(1 - dLambda_cond[, , k], 1e-15))   # n x l
    if (k == 1L) {
      log_S0Y[, , k] <- log_factor
    } else {
      log_S0Y[, , k] <- log_S0Y[, , k - 1L] + log_factor
    }
  }
  S0Y <- exp(log_S0Y)   # guaranteed in (0, 1], monotone non-increasing in k
  
  list(S0Y = S0Y, dLambda_cond = dLambda_cond)
}


# ===========================================================================
# HELPER: .compute_AIPWCC_IF()
#
# Assembles the three-term AIPWCC influence function for each participant i
# at each time t_k.
#
# phi^obs_i(t_k) = term1_i(t_k) + term2_i(t_k) - term3_i(t_k)
#
# where:
#   term1 = [DeltaN_i / G_C^N(tN_i)] * phi^F_i(tN_i, t_k)
#   term2 = [(1-DeltaN_i) / G_C^N(tN_i)] * Q_0(tN_i | tY_i, dY_i, A_i, X_i, t_k)
#   term3 = int_0^{tN_i} Q_0(u | ..., t_k) / G_C^N(u) * dLambda_C^N(u) du
#
# Q_0(r | ..., t_k) = E[phi^F_i(N, t_k) | N >= r, tY_i, dY_i, A_i, X_i]
#   = sum_{j: n_j >= r} phi^F_i(n_j, t_k) * fN_mix(n_j) / sum_{j: n_j >= r} fN_mix(n_j)
#
# fN_mix(n_j): the posterior mixing density = Delta_Y_i * fN1[i,j] + (1-Delta_Y_i) * fN0[i,j]
#   where fN1[i,j] = f^1(n_j | A_i, X_i, Y=tY_i) and
#         fN0[i,j] = f^0(n_j | A_i, X_i, Y>tY_i)
#   These are passed in as fN1_pred and fN0_pred (n x l matrices).
#
# NOTE ON fN1_pred and fN0_pred:
#   These are the probability mass functions (NOT survival functions) of N
#   on the eval.times grid under f^1 and f^0 respectively. They are computed
#   inside .compute_r as dF1 and dF0 and should be stored separately from R.
#   The current scaffold passes them via cf$fN1_pred and cf$fN0_pred — these
#   need to be populated in .compute_r (TODO marked below).
# ===========================================================================

.compute_AIPWCC_IF <- function(
    time, event, time_N, DeltaN, A, a0, admin_t,
    fit.times, eval.times,
    SY,            # n x K
    S0Y,           # n x l x K
    dLambda0Y,     # n x l x K (conditional hazard increments)
    GY,            # n x K
    GN,            # n x K
    pi_hat,        # n
    fN1_pred,      # n x l: f^1(n_j | A_i, X_i, Y = tY_i)
    fN0_pred,      # n x l: f^0(n_j | A_i, X_i, Y > tY_i)
    verbose = FALSE
) {
  
  n <- length(time)
  K <- length(fit.times)
  l <- length(eval.times)
  
  IF_mat  <- matrix(NA_real_, n, K)
  N_ind   <- findInterval(time_N, eval.times, rightmost.closed = TRUE)
  N_ind   <- pmax(1L, pmin(N_ind, l))
  Y_ind   <- findInterval(time, fit.times,  rightmost.closed = TRUE)
  Y_ind   <- pmax(1L, pmin(Y_ind, K))
  
  # G_C^N step function for each participant (used for term3 integral)
  # dLambda_C^N(t_j) = 1 - GN(t_j)/GN(t_{j-1})
  GN_lag      <- cbind(1, GN[, -K, drop=FALSE])
  dLambda_CN  <- pmax(1 - GN / pmax(GN_lag, 1e-10), 0)   # n x K
  
  if(verbose){
    pb <- progress_bar$new(total = n)
  }
  
  for (i in seq_len(n)) {
    
    Ni_idx  <- N_ind[i]
    Yi_idx  <- Y_ind[i]
    dY_i    <- event[i]
    dN_i    <- DeltaN[i]
    A_i     <- A[i]
    pi_i    <- pi_hat[i]
    IPW     <- as.numeric(A_i == a0) / pi_i
    
    # Posterior mixing weights for participant i on eval.times grid (length l)
    fN_mix <- if (dY_i == 1L) fN1_pred[i, ] else fN0_pred[i, ]
    fN_mix <- pmax(fN_mix, 0)
    fN_sum <- sum(fN_mix)
    if (fN_sum > 0) fN_mix <- fN_mix / fN_sum  else fN_mix <- rep(1/l, l)
    
    # Tail probabilities: P(N >= n_j | ...) = sum_{m >= j} fN_mix[m]
    fN_tail <- rev(cumsum(rev(fN_mix)))   # length l
    
    # GN and dLambda_CN at eval.times for participant i (length l)
    GN_eval   <- c(1, GN[i, findInterval(eval.times, fit.times,
                                         rightmost.closed=TRUE)])
    GN_eval_l <- GN_eval[seq_len(l) + 1L]           # G_C^N at eval.times[j]
    dLCN_eval <- pmax(1 - GN_eval_l / GN_eval[seq_len(l)], 0)
    
    GN_at_Ni <- GN_eval_l[Ni_idx]
    
    for (k in seq_len(K)) {
      
      # S_{0Y}(t_k | a0, X_i, n_j) for all j: S0Y[i, , k] (length l)
      S0Y_ik <- S0Y[i, , k]   # length l
      
      # dLambda_{0Y}(t_k | a0, X_i, n_j) for all j
      dL0Y_ik <- dLambda0Y[i, , k]   # length l
      
      # -------------------------------------------------------
      # phi^F_i(n_j, t_k) for each j: the full-data IF at N=n_j
      # Westling et al. Theorem 1:
      #   phi^F = S_{0Y}(t_k|a0,X,n)
      #           - IPW * S_{0Y}(t_k|a0,X,n)
      #             * [I(tY<=t_k, dY=1) / (S_{0Y}(tY|a0,X,n)*G_C^Y(tY))
      #                - int_0^{min(t_k,tY)} dLambda_{0Y}(u|a0,X,n) /
      #                                      (S_{0Y}(u|a0,X,n)*G_C^Y(u))]
      # -------------------------------------------------------
      
      # inner1[j]: I(tY<=t_k, dY=1) / (S_{0Y}(tY|a0,X,n_j) * G_C^Y(tY))
      GY_at_Yi <- GY[i, Yi_idx]
      S0Y_at_Yi <- S0Y[i, , Yi_idx]   # length l: S_{0Y}(tY | a0, X_i, n_j)
      inner1 <- if (time[i] <= fit.times[k] && dY_i == 1L)
        1 / pmax(S0Y_at_Yi * GY_at_Yi, 1e-10) else rep(0, l)
      
      # inner2[j]: cumulative integral up to min(t_k, tY_i)
      # int_0^{min(t_k,tY_i)} dLambda_{0Y}(u|a0,X,n_j) / (S_{0Y}(u|a0,X,n_j)*G_C^Y(u))
      # Discretized over fit.times[1..min(k, Yi_idx)]
      u_max <- min(k, Yi_idx)
      if (u_max >= 1L) {
        # For each m in 1..u_max: dL0Y[i,,m] / (S0Y[i,,m] * GY[i,m])
        # Result is l x u_max, sum over m -> length l
        integrand <- sapply(seq_len(u_max), function(m) {
          dLambda0Y[i, , m] / pmax(S0Y[i, , m] * GY[i, m], 1e-10)
        })   # l x u_max (or l x 1 if u_max=1)
        inner2 <- if (is.matrix(integrand)) rowSums(integrand) else integrand
      } else {
        inner2 <- rep(0, l)
      }
      
      # phi^F(n_j, t_k) for all j: length l
      phiF_k <- S0Y_ik - IPW * S0Y_ik * (inner1 - inner2)
      
      # -------------------------------------------------------
      # Q_0(n_j | ..., t_k) for j = 1..Ni_idx
      # = E[phi^F(N, t_k) | N >= n_j, ...]
      # = sum_{m>=j} phiF_k[m] * fN_mix[m] / sum_{m>=j} fN_mix[m]
      # -------------------------------------------------------
      num_tail  <- rev(cumsum(rev(phiF_k * fN_mix)))   # length l
      Q_seq     <- num_tail / pmax(fN_tail, 1e-12)     # length l
      
      # -------------------------------------------------------
      # Three AIPWCC terms
      # -------------------------------------------------------
      term1 <- (     dN_i  / pmax(GN_at_Ni, 1e-8)) * phiF_k[Ni_idx]
      term2 <- ((1 - dN_i) / pmax(GN_at_Ni, 1e-8)) * Q_seq[Ni_idx]
      
      # term3: discrete integral int_0^{tN_i} Q_0(u) / G_C^N(u) * dLambda_C^N(u) du
      # sum over j = 1..Ni_idx of Q_seq[j] * dLCN_eval[j] / GN_eval_l[j]
      # Note: term3 uses Q_0 evaluated at each eval.time up to tN_i
      idx3  <- seq_len(Ni_idx)
      term3 <- sum(Q_seq[idx3] / pmax(GN_eval_l[idx3], 1e-8) * dLCN_eval[idx3])
      
      IF_mat[i, k] <- term1 + term2 - term3
    }
    
    if(verbose){
      pb$tick()
    }
  }
  
  IF_mat
}


# ===========================================================================
# HELPER: .smooth_r_in_n()
#
# Applies Gaussian kernel smoothing to r(t_k, ., a0, X_i) across the
# eval.times grid, for each row of r_mat.
#
# INPUT:  r_mat  n_te x l matrix of raw likelihood ratios
# OUTPUT: n_te x l matrix of smoothed likelihood ratios
#
# SMOOTHNESS NOTE:
#   This is an optional second layer of smoothing on top of the model's
#   own smoothness. Use when the Cox model is fit with few observations
#   or when eval.times has very fine spacing. The bandwidth r.smooth.bw
#   should be chosen on the scale of eval.times (e.g. if eval.times is
#   in days, r.smooth.bw = 30 smooths over a 30-day window).
# ===========================================================================

.smooth_r_in_n <- function(r_mat, eval.times, bw) {
  l <- length(eval.times)
  # Gaussian kernel weights: K[j, j'] = dnorm((n_j - n_{j'}) / bw)
  K_mat <- outer(eval.times, eval.times, function(a, b) dnorm((a - b) / bw))
  K_mat <- K_mat / rowSums(K_mat)   # row-normalize
  # Smooth: r_smooth[i, ] = r_mat[i, ] %*% t(K_mat) = (K_mat %*% r_mat[i,]^T)^T
  r_smooth <- tcrossprod(r_mat, K_mat)   # n_te x l
  r_smooth
}


# ===========================================================================
# HELPER: .predict_cox_surv()
#
# Predicts the survival function of N on a grid (eval.times) for new
# covariate values (newX), using a fitted coxph object.
#
# OUTPUT: nrow(newX) x length(eval.times) matrix of survival probabilities
# ===========================================================================

.predict_cox_surv_v2 <- function(fit, newX, eval.times) {
  
  # fit = fit1
  # newX = newX_f1_grid
  
  newX <- as.data.frame(newX)
  
  # Drop columns not in the formula to avoid spurious warnings
  formula_vars <- all.vars(fit$formula)
  keep <- intersect(colnames(newX), formula_vars)
  newX <- newX[, keep, drop = FALSE]
  
  sf <- survival::survfit(fit, newdata = newX, se.fit = FALSE)
  
  # sf$surv dimensions depend on nrow(newX):
  #   nrow=1: vector of length n_times
  #   nrow>1: matrix n_times x nrow(newX)
  n_new    <- nrow(newX)
  n_times  <- length(sf$time)
  
  surv_mat <- if (n_new == 1L) {
    matrix(sf$surv, nrow = 1L)            # 1 x n_times
  } else {
    t(matrix(sf$surv, nrow = n_times))    # n_new x n_times
  }
  
  # Interpolate to eval.times: step function, right-continuous (left-limit at jumps)
  surv_on_grid <- apply(surv_mat, 1, function(s) {
    stepfun(sf$time, c(1, s), right = FALSE)(eval.times)
  })  # l x n_new
  
  # Transpose and ensure non-increasing
  out <- t(surv_on_grid)   # n_new x l
  
  # Enforce non-increasing column-wise (can be violated by interpolation)
  # Apply cumulative minimum across columns for each row
  t(apply(out, 1, cummin))
}

# ===========================================================================
# HELPER: .reverse.nuisance.options()
# ===========================================================================

.reverse.nuisance.options <- function(
    cross.fit       = TRUE,
    V               = 5L,
    folds           = NULL,
    Y.df            = 3L,           # spline df for tilde_Y in f^1 and f^0
    event.SL.library = c("survSL.km","survSL.coxph","survSL.expreg"),
    cens.SL.library  = c("survSL.km","survSL.coxph","survSL.expreg"),
    prop.SL.library  = NULL,
    prop.trunc       = 0.05,
    prop.pred        = NULL,
    survSL.control   = list(initWeightAlg = "survSL.km", verbose = FALSE),
    survSL.cvControl = list(V = 5L)
) {
  list(
    cross.fit        = cross.fit,
    V                = V,
    folds            = folds,
    Y.df             = Y.df,
    event.SL.library = event.SL.library,
    cens.SL.library  = cens.SL.library,
    prop.SL.library  = prop.SL.library,
    prop.trunc       = prop.trunc,
    prop.pred        = prop.pred,
    survSL.control   = survSL.control,
    survSL.cvControl = survSL.cvControl
  )
}


# ===========================================================================
# HELPER: .make_df()
# ===========================================================================

.make_df <- function(time, event, time_N, DeltaN_eff, treat, admin_t, W, idx) {
  df <- data.frame(
    Y       = time,
    DeltaY  = event,
    N       = time_N,
    DeltaN  = DeltaN_eff,
    A       = treat,
    admin_t = admin_t
  )
  df <- cbind(df, W)
  df[idx, ]
}


# ===========================================================================
# TODO LIST — items needed to make this scaffold fully runnable
#
# [TODO-1] fN1_pred and fN0_pred storage:
#   Inside .compute_r, dF1 and dF0 are computed as (n_te*K) x l matrices.
#   The posterior mixing weights in .compute_AIPWCC_IF need f^1(n_j | A_i, X_i, Y=tY_i)
#   and f^0(n_j | A_i, X_i, Y>tY_i) evaluated at tY_i = OBSERVED tilde_Y_i,
#   NOT at the fit.times grid. So .compute_r should be called TWICE:
#     (a) Once on the fit.times grid (for building S_{0Y}) -> R array
#     (b) Once evaluated at each participant's own tY_i (for Q_0 weights)
#         -> fN1_pred (n x l) and fN0_pred (n x l)
#   These should be returned from .compute_r and stored in cf$fN1_pred, cf$fN0_pred.
#
# [TODO-2] Arm-specific SY storage:
#   The cross-fitting loop stores SY_a0 as the A=1 block and SY_alt as A=0.
#   This is correct only when the target arm is a0=1. For a0=0, the
#   swap `SY_use <- if (a0==1L) cf$SY_a0 else cf$SY_alt` handles it,
#   but R array cf$R is computed using A=dat_te$A (observed arm), which
#   should be re-computed for each target arm a0 using A=a0 in newX_f1/f0.
#   Currently .compute_r uses dat_te$A. You should pass a0 into .compute_r
#   and set A = a0 in all newX constructions so r reflects the target arm.
#
# [TODO-3] .predict_cox_surv edge cases:
#   If eval.times contains values beyond the last observed N time in the
#   training fold, survfit extrapolates at the last observed value (flat).
#   Add a check and warning. For parametric models (exponential, Weibull)
#   this extrapolation is better-behaved.
#
# ===========================================================================