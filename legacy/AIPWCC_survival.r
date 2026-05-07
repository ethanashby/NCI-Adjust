# =============================================================================
# AIPWCC_CFsurvival.r  –  complete rewrite
#
# Summary of all changes incorporated:
#
#  (1) DeltaN recoding: administrative censoring of N is treated as an
#      observed event (DeltaN_eff = 1) throughout the AIPWCC EIF computation.
#      Random dropout of N (DeltaN=0, time_N < admin_t) remains DeltaN_eff=0.
#
#  (2) Separate G models for Y-censoring (GY) and N-censoring (GN):
#      - GY: G_C^Y(t|A,X) fitted using random-dropout-only indicator for Y
#      - GN: G_C^N(t|A,X) fitted using random-dropout-only indicator for N
#      Both fitted via dedicated survSuperLearner calls on un-augmented data.
#
#  (3) Row-ordering fix: all fold-accumulated nuisance matrices are reindexed
#      to original data order after the cross-fitting loop.
#
#  (4) .get.survival.AIPWCC_fast uses GY for phi^F (vals.i, G_y) and GN
#      for the AIPWCC terms (terms 1/2/3).
#
#  (5) Efficiency: id_to_rows lookup replaces O(n_aug) which() scans;
#      denom_cs hoisted outside k loop; dLambda computed vectorially.
#
#  (6) S_N model: fitted on pooled arms (A included as covariate) to borrow
#      strength and avoid degenerate KM priors in small subgroups.
# =============================================================================

AIPWCC_CFsurvival <- function(
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
    verbose          = TRUE,
    FN_fit           = NULL
) {
  
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
  # FN_fit=FN_fit;
  # nuisance.options = CFsurvival.nuisance.options(
  #   cross.fit = TRUE,
  #   V = 10,
  #   folds = NULL,
  #   eval.times = sort(unique(data$N1[data$N1>0 & data$N1 < c(round(L$admin_date*365))])),
  #   event.SL.library = c("survSL.km",
  #                        "survSL.coxph",
  #                        "survSL.expreg"),
  #   event.pred.0 = NULL,
  #   event.pred.1 = NULL,
  #   cens.SL.library = c("survSL.km",
  #                       "survSL.coxph",
  #                       "survSL.expreg"),
  #   cens.trunc = 0,
  #   cens.pred.0 = NULL,
  #   cens.pred.1 = NULL,
  #   survSL.control = list(initWeightAlg = "survSL.km", verbose = FALSE),
  #   survSL.cvControl = list(V = 5),
  #   save.nuis.fits = FALSE,
  #   prop.SL.library = NULL,
  #   prop.trunc = 0,
  #   prop.pred = prop.fit,
  #   verbose = TRUE
  # )
  
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
  l           <- length(eval.times)
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
    event0   = matrix(NA_real_, 0L, k),   # S_Y(t|A=0,X,N): aug rows
    event1   = matrix(NA_real_, 0L, k),   # S_Y(t|A=1,X,N): aug rows
    censY0   = matrix(NA_real_, 0L, k),   # G_C^Y(t|A=0,X): per participant
    censY1   = matrix(NA_real_, 0L, k),
    censN0   = matrix(NA_real_, 0L, k),   # G_C^N(t|A=0,X): per participant
    censN1   = matrix(NA_real_, 0L, k),   
    fN0      = matrix(NA_real_, 0L, l),   # posterior f_N | A=0
    fN1      = matrix(NA_real_, 0L, l),
    aug_ids  = integer(0L),               # id per aug_test row
    test_ids = integer(0L)                # id per dat_test row
  )
  
  #fit marginal model for $f_{N \mid A, X}$ OUTSIDE this function
  stopifnot(!is.null(FN_fit))
  
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
    
    # ----------------------------------------------------------------
    # Step 1: S_N | A, X  (pooled arms, A as covariate)
    # Smooth library only — no KM — to ensure positive prior density.
    # ----------------------------------------------------------------
    if (verbose) message("    [1] S_N | A, X")

    # ----------------------------------------------------------------
    # Step 2: G_C^Y — random-dropout censoring of Y process
    # random_cens_Y = 1 iff Y is censored AND time_Y < admin_t
    # ----------------------------------------------------------------
    if (verbose) message("    [2] G_C^Y")
    rcY  <- as.integer(dat_tr$DeltaY == 0L & dat_tr$Y < dat_tr$admin_t) #dropout indicator
    newAW_te <- rbind(cbind(A=0L, W_te), cbind(A=1L, W_te))
    
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
    # Step 3: G_C^N — random-dropout censoring of N process
    # random_cens_N = 1 iff N is censored AND time_N < admin_t
    # ----------------------------------------------------------------
    if (verbose) message("    [3] G_C^N")
    rcN  <- as.integer(dat_tr$DeltaN == 0L & dat_tr$N < dat_tr$admin_t)
    
    GN_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr$N,
      event            = rcN,
      X                = cbind(A = dat_tr$A, W_tr),
      newX             = as.data.frame(newAW_te),
      new.times        = fit.times,
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
    # Step 4: Augment training fold (redistribute-to-right weights)
    # For each participant, expand censored N obs into pseudo-rows
    # weighted by P(N = n_j | N > tilde_N, A, X) using S_N model.
    # ----------------------------------------------------------------
    if (verbose) message("    [4] Augmenting training fold")
    
    aug_tr <- .augment_fold(
      dat      = dat_tr,
      W_names  = W_names,
      eval.times = eval.times,
      survest_N  = FN_fit$survest_N
    )
    
    # ----------------------------------------------------------------
    # Step 5: Expand test fold over eval.times grid for S_Y prediction
    # ----------------------------------------------------------------
    aug_te <- .expand_test(dat_te, eval.times, W_names)
    
    # ----------------------------------------------------------------
    # Step 6: S_Y | A, X, N on augmented training data
    # ----------------------------------------------------------------
    if (verbose) message("    [5] S_Y | A, X, N")
    
    AW_ev_tr  <- cbind(A = aug_tr$A, aug_tr[, c("N","DeltaN", W_names), drop=FALSE])
    AW_ev_te  <- cbind(A = aug_te$A, aug_te[, c("N","DeltaN", W_names), drop=FALSE])
    AW_ce_tr  <- cbind(A = aug_tr$A, aug_tr[, W_names, drop=FALSE])
    # newX_cens: one row per unique test participant (distinct W only)
    AW_ce_te_0 <- cbind(A=0L, W_te)
    AW_ce_te_1 <- cbind(A=1L, W_te)
    AW_ce_te   <- rbind(AW_ce_te_0, AW_ce_te_1)
    
    SY_fit <- survSuperLearner2(
      time             = aug_tr$Y,
      event            = aug_tr$DeltaY,
      X_event          = AW_ev_tr,
      X_cens           = AW_ce_tr,
      newX_event       = rbind(cbind(A=0L, aug_te[,c("N","DeltaN",W_names),drop=FALSE]),
                               cbind(A=1L, aug_te[,c("N","DeltaN",W_names),drop=FALSE])),
      newX_cens        = AW_ce_te,
      new.times        = fit.times,
      event.SL.library = nuis$event.SL.library,
      cens.SL.library  = nuis$cens.SL.library,
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl,
      obsWeights       = aug_tr$w,
      id               = aug_tr$id
    )
    
    n_aug_te <- nrow(aug_te)
    SY0 <- SY_fit$event.SL.predict[seq_len(n_aug_te),          , drop=FALSE] #n_aug_te x k
    SY1 <- SY_fit$event.SL.predict[seq_len(n_aug_te) + n_aug_te,, drop=FALSE] #n_aug_te x k
    SY0[SY0 == 0] <- min(SY0[SY0 > 0])
    SY1[SY1 == 0] <- min(SY1[SY1 > 0])
    
    # ----------------------------------------------------------------
    # Step 7: Posterior f_N via Bayes rule on test fold
    # ----------------------------------------------------------------
    if (verbose) message("    [6] Posterior f_N")
    
    post <- .compute_posterior_fN(
      dat_te   = dat_te,
      aug_te   = aug_te,
      SY0      = SY0,
      SY1      = SY1,
      survest_N = FN_fit$survest_N,
      fit.times  = fit.times,
      eval.times = eval.times,
      W_names    = W_names
    )
    
    # ----------------------------------------------------------------
    # Accumulate (fold order)
    # ----------------------------------------------------------------
    acc$event0   <- rbind(acc$event0,   SY0)
    acc$event1   <- rbind(acc$event1,   SY1)
    acc$censY0   <- rbind(acc$censY0,   GY0)
    acc$censY1   <- rbind(acc$censY1,   GY1)
    acc$censN0   <- rbind(acc$censN0,   GN0)
    acc$censN1   <- rbind(acc$censN1,   GN1)
    acc$fN0      <- rbind(acc$fN0,      post$fN0)
    acc$fN1      <- rbind(acc$fN1,      post$fN1)
    acc$aug_ids  <- c(acc$aug_ids,      aug_te$id)
    acc$test_ids <- c(acc$test_ids,     dat_te$id)
    
    rm(aug_tr, aug_te, SY_fit, GY_fit, GN_fit)
    gc()
    if (verbose) message(sprintf("  Fold %d complete", v))
  }
  
  # ------------------------------------------------------------------
  # Reindex ALL fold-accumulated matrices to ORIGINAL DATA ORDER
  # ------------------------------------------------------------------
  orig_ids     <- id
  fold_to_orig <- match(orig_ids, acc$test_ids)   # per-participant index
  
  nuis$censY.pred.0 <- acc$censY0[fold_to_orig, , drop=FALSE]
  nuis$censY.pred.1 <- acc$censY1[fold_to_orig, , drop=FALSE]
  nuis$censN.pred.0 <- acc$censN0[fold_to_orig, , drop=FALSE]
  nuis$censN.pred.1 <- acc$censN1[fold_to_orig, , drop=FALSE]
  nuis$fN_0         <- acc$fN0[fold_to_orig,    , drop=FALSE]
  nuis$fN_1         <- acc$fN1[fold_to_orig,    , drop=FALSE]
  
  # aug-row-aligned matrices: sort so participant blocks match orig_ids order
  aug_ord           <- order(match(acc$aug_ids, orig_ids))
  nuis$id           <- acc$aug_ids[aug_ord]
  nuis$event.pred.0 <- acc$event0[aug_ord, , drop=FALSE]
  nuis$event.pred.1 <- acc$event1[aug_ord, , drop=FALSE]
  
  stopifnot(all(unique(nuis$id) == orig_ids))
  
  # ==================================================================
  # COMPUTE COUNTERFACTUAL SURVIVALS
  # ==================================================================
  if (verbose) message("Computing counterfactual survivals...")
  
  if (0L %in% fit.treat) {
    surv.0 <- .get.survival.AIPWCC(
      Y=time, 
      DeltaY=event, 
      N=time_N, 
      DeltaN=DeltaN_eff,
      A=1L-treat, 
      admin_t=admin_t,
      fit.times=fit.times, 
      eval.times=eval.times, 
      ids=nuis$id,
      S.hats=nuis$event.pred.0,
      GY.hats_a0=nuis$censY.pred.0, 
      GY.hats_alt=nuis$censY.pred.1,
      GN.hats_a0=nuis$censN.pred.0, 
      GN.hats_alt=nuis$censN.pred.1,
      g.hats=1-nuis$prop.pred,
      fN_a0=nuis$fN_0, 
      fN_alt=nuis$fN_1,
      isotonize=TRUE, 
      verbose=verbose
    )
  }
  
  if (1L %in% fit.treat) {
    surv.1 <- .get.survival.AIPWCC(
      Y=time, 
      DeltaY=event, 
      N=time_N, 
      DeltaN=DeltaN_eff,
      A=treat, 
      admin_t=admin_t,
      fit.times=fit.times, 
      eval.times=eval.times, 
      ids=nuis$id,
      S.hats=nuis$event.pred.1,
      GY.hats_a0=nuis$censY.pred.1, 
      GY.hats_alt=nuis$censY.pred.0,
      GN.hats_a0=nuis$censN.pred.1, 
      GN.hats_alt=nuis$censN.pred.0,
      g.hats=nuis$prop.pred,
      fN_a0=nuis$fN_1, 
      fN_alt=nuis$fN_0,
      isotonize=TRUE, 
      verbose=verbose
    )
  }
  
  # ==================================================================
  # CONFIDENCE INTERVALS AND CONTRASTS  (unchanged from original)
  # ==================================================================
  .add_curve <- function(surv_obj, trt, event_sub, treat_sub) {
    ft  <- surv_obj$times
    q01 <- quantile(time[event==1 & treat==treat_sub], .01)
    iso <- surv_obj$surv.iso
    q99 <- max(ft[!is.na(iso) & iso >= .01])
    ci  <- .surv.confints(ft, surv_obj$surv, surv_obj$IF.vals.center,
                          conf.band=conf.band, band.end.pts=c(q01,q99),
                          conf.level=conf.level)
    df  <- data.frame(time=c(0, ft), trt=trt, surv=c(1, iso),
                      se=c(0, ci$res$se), se.logit=c(0, ci$res$se.logit),
                      ptwise.lower=c(1, ci$res$ptwise.lower),
                      ptwise.upper=c(1, ci$res$ptwise.upper),
                      ptwise.logit.lower=c(1, ci$res$ptwise.logit.lower),
                      ptwise.logit.upper=c(1, ci$res$ptwise.logit.upper))
    list(df=df, ci=ci)
  }
  
  if (0L %in% fit.treat) {
    r0 <- .add_curve(surv.0, 0L, event[treat==0], 0L)
    result$surv.df    <- rbind(result$surv.df, r0$df)
    result$IF.vals.0  <- surv.0$IF.vals.center
  }
  if (1L %in% fit.treat) {
    r1 <- .add_curve(surv.1, 1L, event[treat==1], 1L)
    result$surv.df    <- rbind(result$surv.df, r1$df)
    result$IF.vals.1  <- surv.1$IF.vals.center
  }
  
  if (identical(sort(fit.treat), c(0,1))) {
    q01 <- quantile(time[event==1], .01)
    q99 <- max(fit.times[!is.na(surv.0$surv.iso) & surv.0$surv.iso>=.01 |
                           !is.na(surv.1$surv.iso) & surv.1$surv.iso>=.01])
    
    if ("risk.ratio" %in% contrasts){
      
      .risk.ratio(fit.times=fit.times, 
                  surv.0=surv.0$surv, 
                  IF.vals.0 = surv.0$IF.vals.center, 
                  surv.1=surv.1$surv, 
                  IF.vals.1 = surv.1$IF.vals.center, 
                  conf.band=conf.band, 
                  band.end.pts=c(q01, q99), 
                  conf.level=conf.level)$risk.ratio.df
      
      result$risk.ratio <- .risk.ratio(fit.times=fit.times, 
                                       surv.0=surv.0$surv, 
                                       IF.vals.0 = surv.0$IF.vals.center, 
                                       surv.1=surv.1$surv, 
                                       IF.vals.1 = surv.1$IF.vals.center, 
                                       conf.band=conf.band, 
                                       band.end.pts=c(q01, q99), 
                                       conf.level=conf.level)$risk.ratio.df
    }
  }
  
  result$data      <- data.frame(time, event, treat)
  result$nuisance  <- nuis
  result
}


# ==================================================================
# HELPER: augment training fold with redistribute-to-right weights
# ==================================================================
.augment_fold <- function(dat, W_names, eval.times, survest_N) {
  
  rows <- vector("list", nrow(dat))
  
  for (i in seq_len(nrow(dat))) {
    ri          <- dat[i, ]
    Ai          <- ri$A
    N0          <- ri$N
    admin_t_i   <- ri$admin_t
    xi          <- ri[, W_names, drop=FALSE]
    
    # Unique event times in same arm + admin_t as support
    N_support <- c(sort(unique(dat$N[dat$DeltaN==1L & dat$A==Ai])), admin_t_i)
    
    base <- data.frame(id=ri$id, Y=ri$Y, DeltaY=ri$DeltaY,
                       A=Ai, admin_t=admin_t_i)
    base <- cbind(base, xi)
    
    # Admin-censored or observed N → single row, weight 1
    if (ri$DeltaN == 1L || N0 == admin_t_i) {
      rows[[i]] <- cbind(base, N=N0, DeltaN=ri$DeltaN, w=1)
      next
    }
    
    # Randomly censored N → redistribute to right
    N_grid <- N_support[N_support > N0]
    if (length(N_grid) == 0L) {
      rows[[i]] <- cbind(base, N=admin_t_i, DeltaN=0L, w=1)
      next
    }
    
    eps   <- 1e-6 * max(1, N0)
    S_all <- survest_N(t = c(N0 - eps, N_grid), x = xi, a = Ai)
    #S_lag <- c(S_all[1], S_all[-length(S_all)])[-1]  # left limits at N_grid
    # equivalently: S_lag = c(S_N0_minus, S_grid[-last])
    S_N0_minus <- S_all[1]
    S_grid     <- S_all[-1] #right hand limits
    S_lag2     <- c(S_N0_minus, S_grid[-length(S_grid)]) #left hand limits
    w_raw <- pmax(0, S_lag2 - S_grid)
    w_sum <- sum(w_raw)
    w     <- if (w_sum > 0) w_raw / w_sum else rep(1/length(N_grid), length(N_grid))
    
    r <- cbind(base[rep(1L, length(N_grid)), ],
               N      = N_grid,
               DeltaN = as.integer(N_grid < admin_t_i),
               w      = w)
    rows[[i]] <- r[w > 0, ]
  }
  
  do.call(rbind, rows)
}


# ==================================================================
# HELPER: expand test fold over eval.times grid
# ==================================================================
.expand_test <- function(dat_te, eval.times, W_names) {
  rows <- vector("list", nrow(dat_te))
  for (i in seq_len(nrow(dat_te))) {
    ri    <- dat_te[i, ]
    times_i <- eval.times[eval.times <= ri$admin_t]
    l_i   <- length(times_i)
    base  <- ri[rep(1L, l_i), c("id","Y","DeltaY","A","admin_t", W_names),
                drop=FALSE]
    base$N      <- times_i
    base$DeltaN <- c(rep(1L, l_i - 1L), 0L)
    rows[[i]]   <- base
  }
  do.call(rbind, rows)
}


# ==================================================================
# HELPER: Bayes-rule posterior f_N on test fold
# ==================================================================
.compute_posterior_fN <- function(dat_te, aug_te, SY0, SY1,
                                  survest_N, fit.times, eval.times, W_names) {
  
  n_te <- nrow(dat_te)
  l    <- length(eval.times)
  fN0  <- matrix(NA_real_, n_te, l)
  fN1  <- matrix(NA_real_, n_te, l)
  
  # Pre-build id → aug_te row index lookup
  id_to_aug <- split(seq_len(nrow(aug_te)), aug_te$id)
  
  # Helper: find last fit.time <= t_val (right-continuous step function index)
  find_idx <- function(t_val, grid) {
    idx <- findInterval(t_val, grid, rightmost.closed=TRUE)
    pmax(1L, pmin(idx, length(grid)))
  }
  
  for (i in seq_len(n_te)) {
    Yi      <- dat_te$Y[i]
    dYi     <- dat_te$DeltaY[i]
    idi     <- dat_te$id[i]
    xi      <- dat_te[i, W_names, drop=FALSE]
    aug_idx <- id_to_aug[[as.character(idi)]]
    
    # Prior: discrete PMF of N on eval.times under each arm
    SN1 <- survest_N(eval.times, xi, 1L)
    SN0 <- survest_N(eval.times, xi, 0L)
    pr1 <- pmax(0, c(1, SN1[-l]) - SN1)
    pr0 <- pmax(0, c(1, SN0[-l]) - SN0)
    
    # Likelihood: L(Yi, dYi | A=a, X, N=n_j) for each aug_te row (= N grid pt)
    if (dYi == 0L) {
      # Censored Y: likelihood = S_Y(Yi | A, X, N)
      y_idx <- find_idx(Yi, fit.times)
      lk1   <- SY1[aug_idx, y_idx]
      lk0   <- SY0[aug_idx, y_idx]
    } else {
      # Observed event: likelihood = S_Y(t_{j-1}) - S_Y(t_j)
      y_lo  <- find_idx(Yi, fit.times)
      S_lo1 <- if (y_lo == 1L) rep(1, length(aug_idx)) else SY1[aug_idx, y_lo-1L]
      S_lo0 <- if (y_lo == 1L) rep(1, length(aug_idx)) else SY0[aug_idx, y_lo-1L]
      lk1   <- pmax(0, S_lo1 - SY1[aug_idx, y_lo])
      lk0   <- pmax(0, S_lo0 - SY0[aug_idx, y_lo])
      if (all(lk1 == 0)) lk1 <- SY1[aug_idx, y_lo]
      if (all(lk0 == 0)) lk0 <- SY0[aug_idx, y_lo]
    }
    
    # Posterior ∝ likelihood × prior
    u1 <- lk1 * pr1;  d1 <- sum(u1) + 1e-12;  fN1[i,] <- u1 / d1
    u0 <- lk0 * pr0;  d0 <- sum(u0) + 1e-12;  fN0[i,] <- u0 / d0
  }
  
  list(fN0=fN0, fN1=fN1)
}


# ==================================================================
# MAIN EIF COMPUTATION
# ==================================================================
.get.survival.AIPWCC <- function(
    Y, DeltaY, N, DeltaN,   # DeltaN = DeltaN_eff (admin cens recoded to 1)
    A,                       # A[i]=1 iff participant in target arm a0
    admin_t,
    fit.times, eval.times, ids,
    S.hats,                  # S_Y(t|a0,X,N): aug-aligned  (n_aug x K)
    GY.hats_a0,              # G_C^Y(t|a0,X):  per-part    (n x K)
    GY.hats_alt,             # G_C^Y(t|1-a0,X):per-part    (n x K)
    GN.hats_a0,              # G_C^N(t|a0,X):  per-part    (n x K)
    GN.hats_alt,             # G_C^N(t|1-a0,X):per-part    (n x K)
    g.hats,                  # pi(a0|X)                     (n)
    fN_a0,                   # posterior P(N|Y,dY,X,A=a0)   (n x l)
    fN_alt,                  # posterior P(N|Y,dY,X,A=1-a0) (n x l)
    isotonize = TRUE,
    verbose   = TRUE
) {
  
  fit.times <- sort(fit.times[fit.times > 0])
  n  <- length(Y)
  K  <- length(fit.times)
  l  <- length(eval.times)
  
  # Sort K-column matrices
  ord        <- order(fit.times)
  fit.times  <- fit.times[ord]
  S.hats     <- S.hats[,       ord, drop=FALSE]
  GY.hats_a0 <- GY.hats_a0[,  ord, drop=FALSE]
  GY.hats_alt<- GY.hats_alt[, ord, drop=FALSE]
  GN.hats_a0 <- GN.hats_a0[,  ord, drop=FALSE]
  GN.hats_alt<- GN.hats_alt[, ord, drop=FALSE]
  
  # Grid index lookups (pre-computed, avoid repeated which() inside loop)
  Y_ind <- findInterval(Y, fit.times, rightmost.closed=TRUE)
  Y_ind <- pmax(1L, pmin(Y_ind, K))
  N_ind <- findInterval(N, eval.times, rightmost.closed=TRUE)
  N_ind <- pmax(1L, pmin(N_ind, l))
  
  # Fast aug-row lookup
  id_to_rows <- split(seq_along(ids), ids)
  unique_ids <- unique(ids)
  
  IF.vals <- matrix(NA_real_, n, K)
  
  if (verbose) {
    pb <- progress_bar$new(n, format="  EIF [:bar] :percent ETA: :eta")
  }
  
  for (i in seq_len(n)) {
    
    N_indx <- N_ind[i]
    fN.i   <- if (A[i] == 1L) fN_a0[i,] else fN_alt[i,]
    
    # S_Y rows for this participant
    rows_i   <- id_to_rows[[as.character(unique_ids[i])]]
    Shi      <- S.hats[rows_i,, drop=FALSE]   # (l x K)
    
    # G_C^Y step function (for phi^F: vals.i and G_y)
    GYi_vec  <- GY.hats_a0[i,]
    GYi_fun  <- stepfun(fit.times, c(1, GYi_vec), right=TRUE)
    
    # phi^F compensator: vals[j,k] = int_0^{min(t_k,Y_i)} d(1/S_Y) / G_C^Y
    S_inv  <- 1 / Shi
    dS_inv <- S_inv[,-1,drop=FALSE] - S_inv[,-K,drop=FALSE]   # (l x K-1)
    # divide column m by G_C^Y(fit.times[m]), i.e. left endpoint of interval
    vals   <- sweep(dS_inv, 2, 1/GYi_vec[-K], FUN="*")
    vals[, fit.times[-1] > Y[i]] <- 0
    vals   <- cbind(0, t(apply(vals, 1, cumsum)))              # (l x K)
    
    S_y  <- Shi[, Y_ind[i]]       # S_Y(Y_i | a0, X, N_j) for each j
    G_y  <- GYi_fun(Y[i])         # G_C^Y(Y_i | a0, X)
    
    # G_C^N on eval.times grid (for AIPWCC terms 1/2/3)
    if (A[i] == 1L) {
      GNi_fun <- stepfun(fit.times, c(1, GN.hats_a0[i,]),  right=TRUE)
    } else {
      GNi_fun <- stepfun(fit.times, c(1, GN.hats_alt[i,]), right=TRUE)
    }
    GN_at_Ni <- GNi_fun(N[i])
    GN_eval  <- c(1, GNi_fun(eval.times))          # length l+1
    GN_inv   <- 1 / pmax(GN_eval, 1e-8)            # length l+1
    dLambdaN <- 1 - GN_eval[-1] / GN_eval[-(l+1)] # length l
    
    # fN tail-sum denominator (k-independent, hoist outside k loop)
    denom_cs <- rev(cumsum(rev(fN.i)))   # P(N >= eval.times[j] | ...)
    
    IPW_factor <- as.numeric(A[i] == 1L) / g.hats[i]
    
    for (k in seq_len(K)) {
      
      S_t0 <- Shi[, k]   # S_Y(t0 | a0, X, N_j), j=1..l
      
      inner1 <- if (Y[i] <= fit.times[k] && DeltaY[i] == 1L)
        1 / (S_y * G_y) else rep(0, l)
      
      # phi^F as a function of N grid point j
      IF_k <- S_t0 - IPW_factor * S_t0 * (inner1 - vals[, k])
      
      # Q_seq[j] = E[phi^F(N) | N >= eval.times[j], obs, A=A_i]
      num_cs <- rev(cumsum(rev(IF_k * fN.i)))
      Q_seq  <- num_cs / pmax(denom_cs, 1e-12)
      
      # AIPWCC terms
      # DeltaN[i] = DeltaN_eff: admin-censored N has DeltaN_eff=1 → term1
      # randomly censored N has DeltaN_eff=0 → term2
      term1 <- (     DeltaN[i]  / GN_at_Ni) * IF_k[N_indx]
      term2 <- ((1 - DeltaN[i]) / GN_at_Ni) * Q_seq[N_indx]
      
      idx3  <- seq_len(N_indx)
      term3 <- sum(Q_seq[idx3] * GN_inv[idx3] * dLambdaN[idx3])
      
      IF.vals[i, k] <- term1 + term2 - term3
    }
    
    if (verbose) pb$tick()
  }
  
  surv           <- colMeans(IF.vals, na.rm=TRUE)
  IF.vals.center <- sweep(IF.vals, 2, surv, "-")
  
  res <- list(
    times          = fit.times,
    surv           = pmin(1, pmax(0, surv)),
    IF.vals        = IF.vals,
    IF.vals.center = IF.vals.center
  )
  
  if (isotonize) {
    idx              <- !is.na(surv)
    res$surv.iso     <- surv
    res$surv.iso[idx] <- 1 - isoreg(fit.times[idx], 1 - surv[idx])$yf
  }
  
  res
}
