# =============================================================================
# AIPWCC_CFsurvival.r 
#
# =============================================================================

AIPWCC_CFsurvival_new <- function(
    time,
    event,
    time_N,
    event_N,
    admin_t,
    treat,
    W,
    id,
    fit.times        = sort(unique(time[time > 0])),
    fit.treat        = c(0, 1),
    nuisance.options = list(),
    conf.band        = FALSE,
    conf.level       = 0.95,
    contrasts        = c("risk.ratio"),
    verbose          = TRUE
) {
  
  nuis <- do.call("CFsurvival.nuisance.options", nuisance.options)
  n    <- length(time)
  W    <- as.data.frame(W)
  
  # ------------------------------------------------------------------
  # Validate / expand admin_t
  # ------------------------------------------------------------------
  if (length(admin_t) == 1L)  admin_t <- rep(admin_t, n)
  stopifnot(length(admin_t) == n)
  
  # ------------------------------------------------------------------
  # KEY: recode DeltaN for the AIPWCC estimator.
  # Administrative censoring of N (time_N == admin_t, event_N == 0)
  # carries no missing-data uncertainty — N is known to be >= admin_t.
  # Treat these as observed events (DeltaN_eff = 1) so they enter
  # term1 rather than the augmentation path.
  # Random dropout (event_N == 0, time_N < admin_t) stays DeltaN_eff = 0.
  # ------------------------------------------------------------------
  DeltaN_eff <- as.integer(event_N == 1L | time_N == admin_t)
  
  # ------------------------------------------------------------------
  # Folds
  # ------------------------------------------------------------------
  if (!nuis$cross.fit) { nuis$V <- 1L; nuis$folds <- rep(1L, n) }
  if (!is.null(nuis$folds)) {
    nuis$folds <- as.integer(factor(nuis$folds))
    nuis$V     <- length(unique(nuis$folds))
  }
  if (nuis$cross.fit && is.null(nuis$folds)) {
    nuis$folds <- rep(NA_integer_, n)
    for (ev in c(0L, 1L)) {
      idx <- which(event == ev)
      nuis$folds[idx] <- sample(rep(seq_len(nuis$V), length.out = length(idx)))
    }
  }
  
  fit.times   <- sort(fit.times[fit.times > 0])
  k           <- length(fit.times)
  if (is.null(nuis$eval.times)) nuis$eval.times <- fit.times
  eval.times  <- nuis$eval.times
  sub.eval.times <- round(seq(min(eval.times), max(eval.times), length.out=100))
  l           <- length(eval.times)
  l.sub       <- length(sub.eval.times)
  W_names     <- colnames(W)
  
  contrasts <- unique(tolower(contrasts))
  contrasts <- intersect(contrasts, c("surv.diff","surv.ratio","risk.ratio","nnt"))
  
  result <- list(fit.times = fit.times, fit.treat = fit.treat,
                 surv.df = data.frame())
  
  # ==================================================================
  # PROPENSITY SCORE
  # ==================================================================
  if (is.null(nuis$prop.pred)) {
    if (verbose) message("Estimating propensities...")
    nuis$prop.pred <- rep(NA_real_, n)
    for (v in seq_len(nuis$V)) {
      tr <- nuis$folds != v;  te <- nuis$folds == v
      pf <- .estimate.propensity(
        A = treat[tr], W = W[tr,, drop=FALSE], newW = W[te,, drop=FALSE],
        SL.library = nuis$prop.SL.library, fit.treat = fit.treat,
        prop.trunc = nuis$prop.trunc, save.fit = FALSE, verbose = FALSE)
      nuis$prop.pred[te] <- pf$prop.pred
    }
  }
  
  # ==================================================================
  # CROSS-FITTED NUISANCE ESTIMATION
  # ==================================================================
  if (verbose) message("Estimating conditional survivals...")
  stopifnot(nuis$V > 1L)
  
  # Accumulators (fold order; reindexed after loop)
  acc <- list(
    SY1       = matrix(NA_real_, 0L, k),   # S_Y(t|A=1,X)
    SY0       = matrix(NA_real_, 0L, k),   # S_Y(t|A=0,X)
    GY1       = matrix(NA_real_, 0L, k),   # G_C_Y(t|A=1,X)
    GY0       = matrix(NA_real_, 0L, k),   # G_C_Y(t|A=0,X)
    GN1       = matrix(NA_real_, 0L, l),   # G_C_N(t|A=1,X)
    GN0       = matrix(NA_real_, 0L, l),   # G_C_N(t|A=0,X)
    SYN1      = matrix(NA_real_, 0L, k),   # S_Y(t|A=1,X,N)
    SYN0      = matrix(NA_real_, 0L, k),   # S_Y(t|A=0,X,N)
    fN1       = matrix(NA_real_, 0L, l.sub),   # S_Y(t|A=1,X,N)
    fN0       = matrix(NA_real_, 0L, l.sub),   # S_Y(t|A=0,X,N)
    test_ids = integer(0L)                # id per dat_test row
  )
  
  for (v in seq_len(nuis$V)) {
    
    if (verbose) message(sprintf("  Fold %d / %d", v, nuis$V))
    
    tr_idx <- which(nuis$folds != v)
    te_idx <- which(nuis$folds == v)
    
    # Convenience data frames
    make_df <- function(idx) as.data.frame(cbind(
      id=id, Y=time, DeltaY=event,
      N=time_N, DeltaN=event_N, DeltaN_eff=DeltaN_eff,
      A=treat, admin_t=admin_t, W
    ))[idx, ]
    
    dat_tr <- make_df(tr_idx)
    dat_te <- make_df(te_idx)
    n_te   <- nrow(dat_te)
    W_tr   <- dat_tr[, W_names, drop=FALSE]
    W_te   <- dat_te[, W_names, drop=FALSE]
    
    newAW_te <- rbind(cbind(A=0L, W_te), 
                      cbind(A=1L, W_te))
    
    # ----------------------------------------------------------------
    # Step 1: S_Y | A, X  (pooled arms, A as covariate)
    # ----------------------------------------------------------------
    if (verbose) message("    [1] S_N | A, X")
    SY_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr$Y,
      event            = dat_tr$DeltaY,
      X                = cbind(A = dat_tr$A, W_tr),
      newX             = as.data.frame(newAW_te),
      new.times        = fit.times,
      event.SL.library = nuis$event.SL.library,
      cens.SL.library  = nuis$event.SL.library,
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl
    )
    SY_pred <- SY_fit$event.SL.predict          # (2*n_te x k) where k is fit.times
    SY0 <- pmin(pmax(SY_pred[seq_len(n_te),         , drop=FALSE], 0), 1) # (n_te x k) where k is fit.times
    SY1 <- pmin(pmax(SY_pred[seq_len(n_te) + n_te,  , drop=FALSE], 0), 1) # (n_te x k) where k is fit.times
    
    # ----------------------------------------------------------------
    # Step 2: S_N | A, X  (pooled arms, A as covariate)
    # ----------------------------------------------------------------
    if (verbose) message("    [2] S_N | A, X")
    SN_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr$N,
      event            = dat_tr$DeltaN,
      X                = cbind(A = dat_tr$A, W_tr),
      newX             = as.data.frame(newAW_te),
      new.times        = eval.times,
      event.SL.library = nuis$event.SL.library,
      cens.SL.library  = nuis$event.SL.library,
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl
    )
    SN_pred <- SN_fit$event.SL.predict          # (2*n_te x k) where k is fit.times
    SN0 <- pmin(pmax(SN_pred[seq_len(n_te),         , drop=FALSE], 0), 1) # (n_te x k) where k is fit.times
    SN1 <- pmin(pmax(SN_pred[seq_len(n_te) + n_te,  , drop=FALSE], 0), 1) # (n_te x k) where k is fit.times
    
    # ----------------------------------------------------------------
    # Step 3: G_C^Y — random-dropout censoring of Y process
    # random_cens_Y = 1 iff Y is censored AND time_Y < admin_t
    # ----------------------------------------------------------------
    if (verbose) message("    [3] G_C^Y")
    rcY  <- as.integer(dat_tr$DeltaY == 0L & dat_tr$Y < dat_tr$admin_t) #dropout indicator
    
    GY_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr$Y,
      event            = rcY,
      X                = cbind(A = dat_tr$A, W_tr),
      newX             = as.data.frame(newAW_te),
      new.times        = fit.times,
      event.SL.library = nuis$cens.SL.library,
      cens.SL.library  = nuis$cens.SL.library,
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl
    )
    GY_pred <- GY_fit$event.SL.predict          # (2*n_te x k) where k is fit.times
    GY0 <- pmax(GY_pred[seq_len(n_te),         , drop=FALSE], nuis$cens.trunc) # (n_te x k) where k is fit.times
    GY1 <- pmax(GY_pred[seq_len(n_te) + n_te,  , drop=FALSE], nuis$cens.trunc) # (n_te x k) where k is fit.times
    GY0[GY0 == 0] <- min(GY0[GY0 > 0])
    GY1[GY1 == 0] <- min(GY1[GY1 > 0])
    
    # ----------------------------------------------------------------
    # Step 4: G_C^N — random-dropout censoring of N process
    # random_cens_N = 1 iff N is censored AND time_N < admin_t
    # ----------------------------------------------------------------
    if (verbose) message("    [4] G_C^N")
    rcN  <- as.integer(dat_tr$DeltaN == 0L & dat_tr$N < dat_tr$admin_t)
    
    GN_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr$N,
      event            = rcN,
      X                = cbind(A = dat_tr$A, W_tr),
      newX             = as.data.frame(newAW_te),
      new.times        = eval.times,
      event.SL.library = nuis$cens.SL.library,
      cens.SL.library  = nuis$cens.SL.library,
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl
    )
    GN_pred <- GN_fit$event.SL.predict
    GN0 <- pmax(GN_pred[seq_len(n_te),        , drop=FALSE], nuis$cens.trunc)
    GN1 <- pmax(GN_pred[seq_len(n_te) + n_te, , drop=FALSE], nuis$cens.trunc)
    GN0[GN0 == 0] <- min(GN0[GN0 > 0])
    GN1[GN1 == 0] <- min(GN1[GN1 > 0])
    
    # ----------------------------------------------------------------
    # Step 5: S_{0Y}(y \mid a, x, n) using IPCW regression
    # ----------------------------------------------------------------
    if (verbose) message("    [5] S_{0Y} | A, X, N")
    if(verbose) message("    [5] NOTE: This may take a minute")
    #identify complete cases in training data
    dat_tr_cc_indx <- which(dat_tr$DeltaN_eff==1)
    dat_tr_cc <- dat_tr[dat_tr_cc_indx, ]
  
    #predict censoring survival probabilities (weights for SL)
    new.times <- sort(unique(dat_tr_cc$N))
    l<-predict.survSuperLearner(GN_fit, 
                                newdata = dat_tr_cc, 
                                new.times = new.times, onlySL=TRUE)
    match_indx <- match(dat_tr_cc$N, new.times)
    G0N_cc <- sapply(1:length(match_indx), FUN = function(x){
      l$event.SL.predict[x,match_indx[x]]
    })
    
    #Construct evaluation data
    #newAW_te_expand <- rbind(cbind(A=0L, W_te), 
    #                  cbind(A=1L, W_te))
    newAW_te_expand <- tidyr::crossing(newAW_te, data.frame(N=sub.eval.times))
    n_te_expand <- sum(newAW_te_expand$A==0)
    
    SYN_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr_cc$Y,
      event            = dat_tr_cc$DeltaY,
      X                = bind_cols(cbind(A = dat_tr_cc$A, 
                                   W_tr %>% slice(dat_tr_cc_indx)), 
                                   N=dat_tr$N[dat_tr_cc_indx]),
      newX             = as.data.frame(newAW_te_expand),
      new.times        = fit.times,
      event.SL.library = nuis$event.SL.library,
      cens.SL.library  = c('survSL.km'),
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl, 
      obsWeights       = 1/G0N_cc
    )
    gc(verbose=FALSE)
    SYN_pred <- SYN_fit$event.SL.predict          # (2*n_te_expand x k) where k is fit.times
    SYN0 <- pmin(pmax(SYN_pred[seq_len(n_te_expand),         , drop=FALSE], 0), 1) # (n_te_expand x k) where k is fit.times
    SYN1 <- pmin(pmax(SYN_pred[seq_len(n_te_expand) + n_te_expand,  , drop=FALSE], 0), 1) # (n_te_expand x k) where k is fit.times
    
    # ----------------------------------------------------------------
    # Step 6: f_{0N}(n | A, X, tilde{Y}, Delta_Y) (By Bayes Inversion)
    # NOTE on grid:
    #   The posterior is computed on sub.eval.times (50 points) for memory
    #   efficiency. The Q_0 integral in the EIF only needs the posterior
    #   evaluated at each participant's own tN_i — which is looked up by
    #   findInterval(tN_i, sub.eval.times). No interpolation needed.
    # ----------------------------------------------------------------
    if (verbose) message("    [6] f_{0N}(n | A, X, tilde{Y}, Delta_Y) (By Bayes Inversion)")
   
    l_sub <- length(sub.eval.times)   # 50
    sub_idx <- findInterval(sub.eval.times, eval.times, rightmost.closed = TRUE)
    sub_idx <- pmax(1L, pmin(sub_idx, length(eval.times)))
    
    # n_te x l_sub matrices
    SN0_sub <- SN0[, sub_idx, drop = FALSE]
    SN1_sub <- SN1[, sub_idx, drop = FALSE]
    dFN0 <- cbind(1, SN0_sub[, -l_sub, drop = FALSE]) - SN0_sub  # n_te x l_sub
    dFN1 <- cbind(1, SN1_sub[, -l_sub, drop = FALSE]) - SN1_sub
    dFN0 <- pmax(dFN0, 0)
    dFN1 <- pmax(dFN1, 0)
    
    # -----------------------------------------------------------------------
    # 6b. Index of each participant's observed Y time in fit.times grid
    # -----------------------------------------------------------------------
    Y_ind <- findInterval(dat_te$Y, fit.times, rightmost.closed = TRUE)
    Y_ind <- pmax(1L, pmin(Y_ind, length(fit.times)))
    
    .compute_posterior_pmf <- function(
    SYN_arm,      # (n_te * l_sub) x K — S0Y for this arm (A=0 or A=1 block)
    SY_arm,       # n_te x K — marginal S_Y for this arm
    dFN_arm,      # n_te x l_sub — marginal N PMF for this arm
    Y_ind,        # n_te integer vector — index of tY_i in fit.times
    dY,           # n_te integer vector — Delta_Y
    n_te, l_sub, K
    ) {
      # Precompute marginal Y density increment at each participant's Y_ind[i]
      # dSY_marg[i] = SY(t_{Y_ind[i]-1}) - SY(t_{Y_ind[i]})
      # with SY(t_0) := 1
      SY_lag <- cbind(1, SY_arm[, -K, drop = FALSE])   # n_te x K: SY at t_{k-1}
      # dSY_marg[i] = SY_lag[i, Y_ind[i]] - SY_arm[i, Y_ind[i]]
      dSY_marg <- SY_lag[cbind(seq_len(n_te), Y_ind)] -
        SY_arm[cbind(seq_len(n_te), Y_ind)]  # n_te vector
      dSY_marg <- pmax(dSY_marg, 1e-15)
      
      # Marginal S_Y at tY_i (for dY=0 case)
      SY_at_Yi <- SY_arm[cbind(seq_len(n_te), Y_ind)]   # n_te vector
      SY_at_Yi <- pmax(SY_at_Yi, 1e-15)
      
      # Posterior weight matrix: n_te x l_sub (unnormalized)
      L_mat <- matrix(NA_real_, nrow = n_te, ncol = l_sub)
      
      # --- Vectorized extraction of S0Y and dS0Y at each participant's Y_ind ---
      # For participant i, the l_sub rows in SYN_arm are:
      #   rows ((i-1)*l_sub + 1) : (i*l_sub)
      # We need column Y_ind[i] (for survival) and columns Y_ind[i]-1, Y_ind[i]
      # (for density increment).
      #
      # Build a row index vector of length (n_te * l_sub) that maps each
      # SYN_arm row to its participant index:
      part_idx <- rep(seq_len(n_te), each = l_sub)  # length n_te*l_sub
      
      # S0Y at tY_i for each (participant, N-grid) pair: length n_te*l_sub
      S0Y_at_Yi_vec <- SYN_arm[cbind(seq_len(n_te * l_sub),
                                     Y_ind[part_idx])]
      
      # S0Y at t_{Y_ind[i]-1} for density increment (use lag column)
      Y_ind_lag <- pmax(Y_ind - 1L, 1L)   # clamp to 1 (at k=1, lag = k = 1 so diff = 0)
      S0Y_lag_vec  <- SYN_arm[cbind(seq_len(n_te * l_sub), Y_ind_lag[part_idx])]
      # When Y_ind[i] == 1, the lag equals current => density increment = 0,
      # which correctly handles the left boundary.
      dS0Y_vec <- pmax(S0Y_lag_vec - S0Y_at_Yi_vec, 0)  # n_te*l_sub vector
      
      # Reshape to n_te x l_sub (participant-major: row i = participant i)
      S0Y_at_Yi_mat <- matrix(S0Y_at_Yi_vec, nrow = n_te, ncol = l_sub, byrow = TRUE)
      dS0Y_mat      <- matrix(dS0Y_vec,      nrow = n_te, ncol = l_sub, byrow = TRUE)
      
      # Likelihood for dY=1: dS0Y / dSY_marg (both at tY_i)
      L_dy1 <- sweep(dS0Y_mat, 1, dSY_marg, "/")  # n_te x l_sub
      
      # Likelihood for dY=0: S0Y(tY_i | n_j) / SY(tY_i)
      L_dy0 <- sweep(S0Y_at_Yi_mat, 1, SY_at_Yi, "/")   # n_te x l_sub
      
      # Select by dY indicator
      dY_mat <- matrix(dY, nrow = n_te, ncol = l_sub)   # broadcast dY across columns
      L_mat  <- dY_mat * L_dy1 + (1 - dY_mat) * L_dy0   # n_te x l_sub
      
      # Unnormalized posterior: marginal N PMF * likelihood
      W_mat <- dFN_arm * pmax(L_mat, 0)   # n_te x l_sub, non-negative
      
      # Normalize row-wise to get posterior PMF
      W_row_sum <- rowSums(W_mat)
      # Participants with zero total weight get uniform posterior (no information)
      zero_rows  <- W_row_sum < 1e-15
      W_mat[!zero_rows, ] <- W_mat[!zero_rows, ] /
        W_row_sum[!zero_rows]
      W_mat[zero_rows, ]  <- 1 / l_sub
      
      W_mat   # n_te x l_sub posterior PMF
    }
    
    # Compute posteriors for both arms
    fN0_post <- .compute_posterior_pmf(
      SYN_arm = SYN0,          # n_te * l_sub x K (A=0 block from SYN_pred split)
      SY_arm  = SY0,
      dFN_arm = dFN0,
      Y_ind   = Y_ind,
      dY      = dat_te$DeltaY,
      n_te    = n_te,
      l_sub   = l_sub,
      K       = length(fit.times)
    )   # n_te x l_sub

    fN1_post <- .compute_posterior_pmf(
      SYN_arm = SYN1,
      SY_arm  = SY1,
      dFN_arm = dFN1,
      Y_ind   = Y_ind,
      dY      = dat_te$DeltaY,
      n_te    = n_te,
      l_sub   = l_sub,
      K       = length(fit.times)
    )   # n_te x l_sub
    
    # ----------------------------------------------------------------
    # Accumulate (fold order)
    # ----------------------------------------------------------------
    acc$SY0     <- rbind(acc$SY0,     SY0)          # n_te x K per fold -> n x K total
    acc$SY1     <- rbind(acc$SY1,     SY1)
    acc$GY1     <- rbind(acc$GY1,     GY1)
    acc$GY0     <- rbind(acc$GY0,     GY0)
    acc$GN1     <- rbind(acc$GN1,     GN1)
    acc$GN0     <- rbind(acc$GN0,     GN0)
    
    acc$SYN0    <- rbind(acc$SYN0,    SYN0)          # (n_te*l_sub) x K
    acc$SYN1    <- rbind(acc$SYN1,    SYN1)
    
    acc$fN0     <- rbind(acc$fN0,     fN0_post)      # n_te x l_sub
    acc$fN1     <- rbind(acc$fN1,     fN1_post)

    acc$test_ids <- c(acc$test_ids, dat_te$id)
    
    rm(dat_tr, dat_te, dat_tr_cc, dFN0, dFN1, fN0_post, fN1_post, GN_pred, GN0, GN1, GY_pred, GY0, GY1, l, newAW_te, newAW_te_expand, SN_fit, SN_pred, SN0, SN1_sub, SY_pred, SY0, SY1, SYN_fit, SYN_pred, SYN1, W_te, W_tr)
    gc()
    if (verbose) message(sprintf("  Fold %d complete", v))
  }
  
  # ------------------------------------------------------------------
  # Reindex ALL fold-accumulated matrices to ORIGINAL DATA ORDER
  # ------------------------------------------------------------------
  orig_ids     <- id
  orig_order <- match(orig_ids, acc$test_ids)   # per-participant index
  
  nuis$SY0     <- acc$SY0[orig_order,  , drop = FALSE]   # n x K
  nuis$SY1     <- acc$SY1[orig_order,  , drop = FALSE]
  nuis$GY0     <- acc$GY0[orig_order,  , drop = FALSE]   # n x K
  nuis$GY1     <- acc$GY1[orig_order,  , drop = FALSE]   # n x K
  nuis$GN0     <- acc$GN0[orig_order,  , drop = FALSE]   # n x K
  nuis$GN1     <- acc$GN1[orig_order,  , drop = FALSE]   # n x K
  nuis$fN0     <- acc$fN0[orig_order,  , drop = FALSE]   # n x l_sub
  nuis$fN1     <- acc$fN1[orig_order,  , drop = FALSE]
  #nuis$dFN0    <- acc$dFN0[orig_order, , drop = FALSE]   # n x l_sub
  #nuis$dFN1    <- acc$dFN1[orig_order, , drop = FALSE]
  expanded_orig_order <- as.vector(
    t(matrix((orig_order - 1L) * l_sub, nrow = n, ncol = l_sub) +
        matrix(seq_len(l_sub), nrow = n, ncol = l_sub, byrow = TRUE))
  )
  nuis$SYN0 <- acc$SYN0[expanded_orig_order, , drop = FALSE]
  nuis$SYN1 <- acc$SYN1[expanded_orig_order, , drop = FALSE]
  
  stopifnot(all(unique(nuis$id) == orig_ids))
  
  # ==================================================================
  # COMPUTE COUNTERFACTUAL SURVIVALS
  # ==================================================================
  if (verbose) message("Computing counterfactual survivals...")
  
  .compute_IF_a0 <- function(
    a0,            # target arm (0 or 1)
    time, event,   # tY, dY
    time_N, DeltaN_eff,  # tN, dN (DeltaN_eff: admin cens treated as observed)
    treat,         # observed A
    pi_hat,        # P(A=1|X), length n
    SY,            # n x K: S_Y(t|a0,X)
    SYN,           # (n*l_sub) x K: S0Y(t|a0,X,n_j), participant-major
    GY,            # n x K: G_C^Y(t|a0,X)
    GN,            # n x l: G_C^N(n|A_i,X_i) — on eval.times, NOT sub.eval.times
    fN_post,       # n x l_sub: posterior PMF of N
    fit.times,     # K-vector
    sub.eval.times,# l_sub-vector
    eval.times,    # l-vector (for GN lookup)
    n, K, l_sub
  ) {
    
    # a0             = a0;
    # time           = time;
    # event          = event;
    # time_N         = time_N;
    # DeltaN_eff     = DeltaN_eff;
    # treat          = treat;
    # pi_hat         = nuis$prop.pred;
    # SY             = SY_use;
    # SYN            = SYN_use;
    # GY             = GY_use;
    # GN             = GN_use;
    # fN_post        = fN_use;
    # fit.times      = fit.times;
    # sub.eval.times = sub.eval.times;
    # eval.times     = eval.times;
    # n = n; 
    # K = k; 
    # l_sub = l_sub
    
    pi_use  <- if (a0 == 1L) pi_hat else 1 - pi_hat
    A_match <- as.numeric(treat == a0)
    
    # Map each participant's tN_i and tY_i to grid indices
    N_ind <- findInterval(time_N, sub.eval.times, rightmost.closed = TRUE)
    N_ind <- pmax(1L, pmin(N_ind, l_sub))
    Y_ind <- findInterval(time, fit.times, rightmost.closed = TRUE)
    Y_ind <- pmax(1L, pmin(Y_ind, K))
    
    # G_C^N at sub.eval.times for each participant
    # GN is n x l on eval.times; interpolate to sub.eval.times
    sub_in_eval <- findInterval(sub.eval.times, eval.times, rightmost.closed = TRUE)
    sub_in_eval <- pmax(1L, pmin(sub_in_eval, ncol(GN)))
    #print(sub_in_eval)
    GN_sub <- GN[, sub_in_eval, drop = FALSE]   # n x l_sub
    
    # dLambda_C^N at sub.eval.times: 1 - GN(n_j)/GN(n_{j-1}), with GN(n_0)=1
    GN_sub_lag   <- cbind(1, GN_sub[, -l_sub, drop = FALSE])
    dLCN_sub     <- pmax(1 - GN_sub / pmax(GN_sub_lag, 1e-10), 0)   # n x l_sub
    
    # Discrete hazard of Y given (a0, X) — marginal: n x K
    SY_lag       <- cbind(1, SY[, -K, drop = FALSE])
    dLambda_Y    <- pmax(1 - SY / pmax(SY_lag, 1e-10), 0)            # n x K
    
    # Discrete hazard of Y given (a0, X, N=n_j):
    # dLambda_0Y[i,j,k] extracted on-the-fly inside the participant loop below
    # (avoids storing n x l_sub x K array)
    
    IF_mat <- matrix(NA_real_, nrow = n, ncol = K)
    
    if(verbose){
      pb <- progress_bar$new(n, format="  EIF [:bar] :percent ETA: :eta")
    }
    
    for (i in seq_len(n)) {
      
      # Rows in SYN for participant i: ((i-1)*l_sub+1) : (i*l_sub)
      i_rows  <- ((i - 1L) * l_sub + 1L) : (i * l_sub)
      S0Y_i   <- SYN[i_rows, , drop = FALSE]   # l_sub x K
      
      # G_C^Y at tY_i for this participant
      GY_at_Yi <- GY[i, Y_ind[i]]
      
      # S0Y at tY_i for all N-grid points: l_sub vector
      S0Y_at_Yi <- S0Y_i[, Y_ind[i]]
      
      # Posterior PMF and tail sums
      fN_i   <- fN_post[i, ]                          # l_sub
      fN_tail <- rev(cumsum(rev(fN_i)))               # P(N >= n_j | ...)
      
      # GN at tN_i for this participant (on sub.eval.times grid)
      GN_at_Ni <- GN_sub[i, N_ind[i]]
      
      dN_i <- DeltaN_eff[i]
      IPW  <- A_match[i] / pmax(pi_use[i], 1e-8)
      
      for (k in seq_len(K)) {
        
        t0 <- fit.times[k]
        
        # S0Y(t0 | a0, x_i, n_j) for all j: l_sub vector
        S0Y_ik <- S0Y_i[, k]   # l_sub
        
        # Conditional discrete hazard at t0 for all j
        # dLambda_0Y[j, k] = 1 - S0Y[j,k] / S0Y[j,k-1]
        if (k == 1L) {
          dL0Y_ik <- pmax(1 - S0Y_ik, 0)           # S0Y(t_{k-1}) = 1
        } else {
          dL0Y_ik <- pmax(1 - S0Y_ik / pmax(S0Y_i[, k - 1L], 1e-15), 0)
        }
        
        # inner1: event indicator term (length l_sub)
        inner1 <- if (time[i] <= t0 && event[i] == 1L)
          1 / pmax(S0Y_at_Yi * GY_at_Yi, 1e-10) else rep(0, l_sub)
        
        # inner2: compensator integral up to min(t0, tY_i)
        # sum_{m=1}^{min(k, Y_ind[i])} dLambda_0Y[j,m] / (S0Y[j,m] * GY[i,m])
        u_max <- min(k, Y_ind[i])
        if (u_max >= 1L) {
          # Build l_sub x u_max matrix of integrand values
          integrand_mat <- sapply(seq_len(u_max), function(m) {
            dL0Y_m <- if (m == 1L) pmax(1 - S0Y_i[, 1L], 0) else
              pmax(1 - S0Y_i[, m] / pmax(S0Y_i[, m - 1L], 1e-15), 0)
            dL0Y_m / pmax(S0Y_i[, m] * GY[i, m], 1e-10)
          })  # l_sub x u_max (or l_sub vector if u_max=1)
          inner2 <- if (is.matrix(integrand_mat)) rowSums(integrand_mat) else integrand_mat
        } else {
          inner2 <- rep(0, l_sub)
        }
        
        # phi^F(n_j, t0) for all j: l_sub vector
        phiF <- S0Y_ik - IPW * S0Y_ik * (inner1 - inner2)
        
        # Q_0(n_j) = E[phi^F(N) | N >= n_j, ...]
        num_cs  <- rev(cumsum(rev(phiF * fN_i)))           # l_sub
        Q_seq   <- num_cs / pmax(fN_tail, 1e-12)           # l_sub
        
        # Three AIPWCC terms
        term1 <- (     dN_i  / pmax(GN_at_Ni, 1e-8)) * phiF[N_ind[i]]
        term2 <- ((1 - dN_i) / pmax(GN_at_Ni, 1e-8)) * Q_seq[N_ind[i]]
        
        # term3: int_0^{tN_i} Q_0(u)/G_C^N(u) dLambda_C^N(u)
        idx3  <- seq_len(N_ind[i])
        term3 <- sum(Q_seq[idx3] / pmax(GN_sub[i, idx3], 1e-8) * dLCN_sub[i, idx3])
        
        IF_mat[i, k] <- term1 + term2 - term3
      }
      
      pb$tick()
    }
    
    IF_mat
  }
  
  if (0L %in% fit.treat) {
    
    a0=0
    
    SY_use  <- nuis$SY0
    SYN_use <- nuis$SYN0
    GY_use  <- nuis$GY0
    GN_use  <- nuis$GN0
    fN_use  <- nuis$fN0
    
    IF_mat <- .compute_IF_a0(
      a0             = a0,
      time           = time,
      event          = event,
      time_N         = time_N,
      DeltaN_eff     = DeltaN_eff,
      treat          = treat,
      pi_hat         = nuis$prop.pred,
      SY             = SY_use,
      SYN            = SYN_use,
      GY             = GY_use,
      GN             = GN_use,
      fN_post        = fN_use,
      fit.times      = fit.times,
      sub.eval.times = sub.eval.times,
      eval.times     = eval.times,
      n = n, 
      K = k, 
      l_sub = l_sub
    )
    
    surv_hat  <- colMeans(IF_mat, na.rm = TRUE)
    IF_center <- sweep(IF_mat, 2, surv_hat, "-")
    
    z   <- qnorm((1 + conf.level) / 2)
    se  <- sqrt(apply(IF_center, 2, var, na.rm = TRUE) / n)
    
    # Logit-scale CI
    logit_s   <- log(pmax(surv_hat, 1e-8) / pmax(1 - surv_hat, 1e-8))
    se_logit  <- se / pmax(surv_hat * (1 - surv_hat), 1e-8)
    
    result$surv.df <- rbind(result$surv.df, data.frame(
      time               = c(0, fit.times),
      trt                = a0,
      surv               = c(1, pmin(1, pmax(0, surv_hat))),
      se                 = c(0, se),
      ptwise.lower       = c(1, plogis(logit_s - z * se_logit)),
      ptwise.upper       = c(1, plogis(logit_s + z * se_logit))
    ))
    
  }
  
  if (1L %in% fit.treat) {
    
    a0=1
    
    SY_use  <- nuis$SY1
    SYN_use <- nuis$SYN1
    GY_use  <- nuis$GY1
    GN_use  <- nuis$GN1
    fN_use  <- nuis$fN1
    
    IF_mat <- .compute_IF_a0(
      a0             = a0,
      time           = time,
      event          = event,
      time_N         = time_N,
      DeltaN_eff     = DeltaN_eff,
      treat          = treat,
      pi_hat         = nuis$prop.pred,
      SY             = SY_use,
      SYN            = SYN_use,
      GY             = GY_use,
      GN             = GN_use,
      fN_post        = fN_use,
      fit.times      = fit.times,
      sub.eval.times = sub.eval.times,
      eval.times     = eval.times,
      n = n, 
      K = k, 
      l_sub = l_sub
    )
    
    surv_hat  <- colMeans(IF_mat, na.rm = TRUE)
    IF_center <- sweep(IF_mat, 2, surv_hat, "-")
    
    z   <- qnorm((1 + conf.level) / 2)
    se  <- sqrt(apply(IF_center, 2, var, na.rm = TRUE) / n)
    
    # Logit-scale CI
    logit_s   <- log(pmax(surv_hat, 1e-8) / pmax(1 - surv_hat, 1e-8))
    se_logit  <- se / pmax(surv_hat * (1 - surv_hat), 1e-8)
    
    result$surv.df <- rbind(result$surv.df, data.frame(
      time               = c(0, fit.times),
      trt                = a0,
      surv               = c(1, pmin(1, pmax(0, surv_hat))),
      se                 = c(0, se),
      ptwise.lower       = c(1, plogis(logit_s - z * se_logit)),
      ptwise.upper       = c(1, plogis(logit_s + z * se_logit))
    ))
    
  }
  
  # if (identical(sort(fit.treat), c(0,1))) {
  #   q01 <- quantile(time[event==1], .01)
  #   q99 <- max(fit.times[!is.na(surv.0$surv.iso) & surv.0$surv.iso>=.01 |
  #                          !is.na(surv.1$surv.iso) & surv.1$surv.iso>=.01])
  #   
  #   if ("risk.ratio" %in% contrasts){
  #     
  #     .risk.ratio(fit.times=fit.times, 
  #                 surv.0=surv.0$surv, 
  #                 IF.vals.0 = surv.0$IF.vals.center, 
  #                 surv.1=surv.1$surv, 
  #                 IF.vals.1 = surv.1$IF.vals.center, 
  #                 conf.band=conf.band, 
  #                 band.end.pts=c(q01, q99), 
  #                 conf.level=conf.level)$risk.ratio.df
  #     
  #     result$risk.ratio <- .risk.ratio(fit.times=fit.times, 
  #                                      surv.0=surv.0$surv, 
  #                                      IF.vals.0 = surv.0$IF.vals.center, 
  #                                      surv.1=surv.1$surv, 
  #                                      IF.vals.1 = surv.1$IF.vals.center, 
  #                                      conf.band=conf.band, 
  #                                      band.end.pts=c(q01, q99), 
  #                                      conf.level=conf.level)$risk.ratio.df
  #   }
  # }
  
  result$data      <- data.frame(time, event, treat)
  result$nuisance  <- nuis
  result
}

