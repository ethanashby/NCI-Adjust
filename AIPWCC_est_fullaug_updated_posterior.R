# =============================================================================
# AIPWCC_CFsurvival.r 
# More memory efficient -- evaluates EIF at each inner CV loop to avoid memory exhaustion limits
# ALSO implements new EIF calculation
# ALSO offers EM-style update to S_{0Y}
# NOTE that the projection is given by 
# Delta_N * Psi^F / K(N | V) + int_0^{infty} Q_0(u)/K(u |V) dM_C(u | V)
# Where V is all observed data A,X,TildeY, DeltaY
# We CANNOT treat K(u|V)=K(u|A,X) because there is some partial information about C
# =============================================================================

AIPWCC_newEIF_em <- function(
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
    single_cens_mechanism = TRUE,
    em               = TRUE,
    em_fast          = TRUE,
    max_iter         = 10,
    tol              = 1E-4
) {
  
  # time = data$Y;
  # event = data$DeltaY;
  # time_N = data$N;
  # event_N = data$DeltaN;
  # treat = data$A;
  # W = data %>%  dplyr::select("X");
  # fit.times = sort(unique(data$Y));
  # fit.treat = c(0,1);
  # admin_t = round(365*L$admin_date);
  # id = 1:nrow(data);
  # nuisance.options = .CFsurvival.nuis.opt(
  #   cross.fit = TRUE,
  #   V = 5,
  #   folds = NULL,
  #   eval.times = sort(unique(c(0, data$N))),
  #   event.SL.library = c("survSL.coxph",
  #                        "survSL.expreg",
  #                        "survSL.km",
  #                        "survSL.coxph.ps"),
  #   event.SL.library_SY = c("survSL.coxph",
  #                           "survSL.coxph.ps",
  #                           "survSL.rfsrc"),
  #   event.pred.0 = NULL,
  #   event.pred.1 = NULL,
  #   cens.SL.library = c("survSL.coxph",
  #                       "survSL.expreg",
  #                       "survSL.km",
  #                       "survSL.coxph.ps"),
  #   cens.trunc = 0,
  #   cens.pred.0 = NULL,
  #   cens.pred.1 = NULL,
  #   survSL.control = list(initWeightAlg = "survSL.coxph", verbose = FALSE),
  #   survSL.cvControl = list(V = 5),
  #   save.nuis.fits = FALSE,
  #   prop.SL.library = NULL,
  #   prop.trunc = 0,
  #   prop.pred = prop.fit,
  #   verbose = TRUE
  # );
  # conf.band = FALSE;
  # conf.level = 0.95;
  # contrasts = c("risk.ratio");
  # verbose = TRUE;
  # single_cens_mechanism = TRUE;
  # em = TRUE;
  # em_fast = TRUE
  # max_iter         = 3
  # tol              = 1E-4
  
  nuis <- do.call(".CFsurvival.nuis.opt", nuisance.options)
  n    <- length(time)
  W    <- as.data.frame(W)
  
  # ------------------------------------------------------------------
  # Validate / expand admin_t
  # ------------------------------------------------------------------
  if (length(admin_t) == 1L)  admin_t <- rep(admin_t, n)
  stopifnot(length(admin_t) == n)
  
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
    for(a in c(0L, 1L)){
      for (ev in c(0L, 1L)) {
        idx <- which(event == ev & treat == a)
        nuis$folds[idx] <- sample(rep(seq_len(nuis$V), length.out = length(idx)))
      }
    }
  }
  
  fit.times   <- sort(fit.times[fit.times > 0])
  k           <- length(fit.times)
  if(is.null(nuis$eval.times)){
    nuis$eval.times <- fit.times
  }
  admin_times <- unique(admin_t)
  nuis$eval.times <- sort(unique(c(nuis$eval.times, admin_times)))
  eval.times <- nuis$eval.times
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
    IF0      = matrix(NA_real_, 0L, k),   # EIF rows for arm 0, per participant
    IF1      = matrix(NA_real_, 0L, k),   # EIF rows for arm 1, per participant
    test_ids = integer(0L)
  )
  
  for (v in seq_len(nuis$V)) {
    #v=1
    if (verbose) message(sprintf("  Fold %d / %d", v, nuis$V))
    tr_idx <- which(nuis$folds != v)
    te_idx <- which(nuis$folds == v)
    
    # Convenience data frames
    make_df <- function(idx) as.data.frame(cbind(
      id = id, 
      Y = time, 
      DeltaY = event,
      N = time_N, # FIX: assigned N to dataset where admin_censored participants are replaced by sentinel value 
      DeltaN = event_N, 
      A = treat, 
      admin_t = admin_t, 
      W
    ))[idx, ]
    
    dat_tr <- make_df(tr_idx)
    dat_te <- make_df(te_idx)
    n_te   <- nrow(dat_te)
    W_tr   <- dat_tr[, W_names, drop=FALSE]
    W_te   <- dat_te[, W_names, drop=FALSE]
    
    # ----------------------------------------------------------------
    # Step 0: generate test datasets
    # ----------------------------------------------------------------
    
    # Global augmentation over entire eval.times grid
    newAW_te <- rbind(cbind(A=0L, W_te), cbind(A=1L, W_te))
    aug_te <- .expand_test(dat_te, eval.times, W_names)
    n_aug_te <- nrow(aug_te)
    
    # ----------------------------------------------------------------
    # Step 1: S_N | A, X  (pooled arms, A as covariate)
    # Smooth library only — no KM — to ensure positive prior density.
    # ----------------------------------------------------------------
    if (verbose) message("    [1] S_N | A, X")
    SN_fit <- survSuperLearner::survSuperLearner(
      time             = dat_tr$N,
      event            = dat_tr$DeltaN,
      X                = W_tr,
      newX             = as.data.frame(W_te),
      new.times        = eval.times,
      event.SL.library = nuis$event.SL.library[nuis$event.SL.library!="survSL.km"], #ENFORCING SMOOTHNESSS
      cens.SL.library  = nuis$event.SL.library, #[nuis$event.SL.library!="survSL.km"], #No interest in censoring, keep simple
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl
    )
    
    #survest function
    survest_N <- function(fit=SN_fit, t, x) {
      x <- as.data.frame(x)
      colnames(x) <- colnames(W)
      if (nrow(x) == 1L) { x <- rbind(x, x); drop1 <- TRUE } else drop1 <- FALSE
      newdat <- bind_cols(x)
      out <- .pSSL_event_only(object = fit, 
                              newdata = newdat, 
                              new.times = t, 
                              onlySL = TRUE)$event.SL.predict
      if (drop1) out <- out[1,, drop=FALSE]
      as.vector(out)
    }
    
    # ----------------------------------------------------------------
    # Step 2: Augment training fold (redistribute-to-right weights)
    # For each participant, expand censored N obs into pseudo-rows
    # weighted by P(N = n_j | N > tilde_N, A, X) using S_N model.
    # ----------------------------------------------------------------
    if (verbose) message("    [2] Augmenting training fold via redistribute-to-the-right")
    
    aug_tr <- .augment_fold(
      dat      = dat_tr,
      W_names  = W_names,
      eval.times = eval.times,
      survest_N  = survest_N
    )
    
    # ----------------------------------------------------------------
    # Step 6: S_Y | A, X, N on augmented training data
    # ----------------------------------------------------------------
    # TRIED: `stackG` in `survML` but just takes FAR too long computationally
    # You could try SuperSurv in future
    
    if (verbose) message("    [3] S_Y | A, X, N")
    
    if(em==FALSE){
      
      if(verbose) message("        No EM-style updates (S_{0Y} may be misspecified)")
      
      AW_ev_tr  <- cbind(A = aug_tr$A, aug_tr[, c("N", W_names), drop=FALSE])
      AW_ce_tr  <- cbind(A = aug_tr$A, aug_tr[, W_names, drop=FALSE])
      AW_ev_te  <- cbind(A = aug_te$A, aug_te[, c("N", W_names), drop=FALSE])
      # newX_cens: one row per unique test participant (distinct W only)
      AW_ce_te_0 <- cbind(A=0L, W_te)
      AW_ce_te_1 <- cbind(A=1L, W_te)
      AW_ce_te   <- rbind(AW_ce_te_0, AW_ce_te_1)
      
      SY_fit <- survSuperLearner2(
        time             = aug_tr$Y,
        event            = aug_tr$DeltaY,
        X_event          = AW_ev_tr,
        X_cens           = AW_ce_tr,
        newX_event       = cbind(A=0L, aug_te[,c("N",W_names),drop=FALSE] %>% slice(1:2)),
        newX_cens        = AW_ce_te[1:2,],
        new.times        = fit.times,
        event.SL.library = nuis$event.SL.library_SY[nuis$event.SL.library_SY!="survSL.km"], #ENFORCING SMOOTHNESSS
        cens.SL.library  = nuis$cens.SL.library[nuis$cens.SL.library!="survSL.km"], #ENFORCING SMOOTHNESSS
        verbose          = FALSE,
        control          = nuis$survSL.control,
        cvControl        = nuis$survSL.cvControl,
        obsWeights       = aug_tr$w,
        id               = aug_tr$id
      )
      
      SY0 <- .pSSL_event_only(SY_fit, 
                              newdata = cbind(A=0L, aug_te[,c("N",W_names),drop=FALSE]), 
                              new.times = fit.times, 
                              onlySL = TRUE)$event.SL.predict %>%  bigmemory::as.big.matrix()
      
      gc()
      
      SY1 <- .pSSL_event_only(SY_fit, 
                              newdata = cbind(A=1L, aug_te[,c("N",W_names),drop=FALSE]), 
                              new.times = fit.times, 
                              onlySL = TRUE)$event.SL.predict %>%  bigmemory::as.big.matrix()
      
      gc()
      
      rm(SY_fit)
      gc()
      
    }else{
      
      if(verbose) message("        Iterative EM-style updates to S_{0Y}")  
      
      gc()
      
      # This code throws an error: `Error in eval(extras, data, env) : argument "weights" is missing, with no default`
      em_out <- .fit_SY_em(
        aug_tr            = aug_tr,
        aug_te            = aug_te,
        dat_te            = dat_te,
        fit.times         = fit.times,
        W_names           = W_names,
        nuis              = nuis,
        max_iter          = max_iter,             # 1 disables EM; 2 = one update
        tol               = tol,
        fast_iter_library = "survSL.coxph",           # e.g. c("survSL.coxph") for speed
        verbose           = verbose, 
        weight_norm = "L2")
      
      SY0    <- bigmemory::as.big.matrix(em_out$SY0)
      SY1    <- bigmemory::as.big.matrix(em_out$SY1)
      aug_tr <- em_out$aug_tr_final         # downstream code uses aug_tr$id only
      rm(em_out)
      gc()
      
    }
    # ----------------------------------------------------------------
    # Step 2: G_C^Y — random-dropout censoring of Y process
    # random_cens_Y = 1 iff Y is censored AND time_Y < admin_t
    # ----------------------------------------------------------------
    
    if(single_cens_mechanism==TRUE){
      
      if (verbose) message("SINGLE censoring mechanism specified")
      if (verbose) message("    [4] G_C")
      rc <- pmax(as.integer(dat_tr$DeltaY == 0L & dat_tr$Y < dat_tr$admin_t),
                 as.integer(dat_tr$DeltaN == 0L & dat_tr$N < dat_tr$admin_t))
      
      C <- case_when(
        dat_tr$DeltaY == 0L & dat_tr$DeltaN==0L ~ pmax(dat_tr$Y, dat_tr$N),
        dat_tr$DeltaY == 1L & dat_tr$DeltaN==0L ~ dat_tr$N,
        dat_tr$DeltaY == 0L & dat_tr$DeltaN==1L ~ dat_tr$Y,
        dat_tr$DeltaY == 1L & dat_tr$DeltaN==1L ~ pmax(dat_tr$Y, dat_tr$N)
      )
      
      # USE CENS TO GET LEFT CONT
      GC_fit <- survSuperLearner::survSuperLearner(
        time             = C,
        event            = 1-rc,
        X                = cbind(A = dat_tr$A, W_tr),
        newX             = as.data.frame(newAW_te),
        new.times        = sort(unique(c(eval.times, fit.times))),
        event.SL.library = nuis$cens.SL.library,
        cens.SL.library  = nuis$cens.SL.library,
        verbose          = FALSE,
        control          = nuis$survSL.control,
        cvControl        = nuis$survSL.cvControl
      )
      GC_pred <- GC_fit$cens.SL.predict          # (2*n_te x k) where k is ALL times
      if(0 %in% eval.times){
        GC_pred[, 1] <- 1  # S(0) = 1 by construction
      }
      GC0 <- pmax(GC_pred[seq_len(n_te),         , drop=FALSE], nuis$cens.trunc) # (n_te x k) where k is ALL times
      GC1 <- pmax(GC_pred[seq_len(n_te) + n_te,  , drop=FALSE], nuis$cens.trunc) # (n_te x k) where k is ALL times
      GC0[GC0 == 0] <- min(GC0[GC0 > 0])
      GC1[GC1 == 0] <- min(GC1[GC1 > 0])
      
    }else{
      
      if (verbose) message("TWO censoring mechanism specified")
      if (verbose) message("    [4] G_{CY}")
      
      # Dropout censoring Y
      rc <- as.integer(dat_tr$DeltaY == 0L & dat_tr$Y < dat_tr$admin_t)
      
      C <- dat_tr$Y
      
      #USE CENS TO GET LEFT CONTINUOUS
      GC_fit <- survSuperLearner::survSuperLearner(
        time             = C,
        event            = 1-rc,
        X                = cbind(A = dat_tr$A, W_tr),
        newX             = as.data.frame(newAW_te),
        new.times        = sort(unique(c(eval.times, fit.times))),
        event.SL.library = nuis$cens.SL.library,
        cens.SL.library  = nuis$cens.SL.library,
        verbose          = FALSE,
        control          = nuis$survSL.control,
        cvControl        = nuis$survSL.cvControl
      )
      GC_pred <- GC_fit$cens.SL.predict          # (2*n_te x k) where k is ALL times
      if(0 %in% eval.times){
        GC_pred[, 1] <- 1  # S(0) = 1 by construction
      }
      GCY0 <- pmax(GC_pred[seq_len(n_te),         , drop=FALSE], nuis$cens.trunc) # (n_te x k) where k is ALL times
      GCY1 <- pmax(GC_pred[seq_len(n_te) + n_te,  , drop=FALSE], nuis$cens.trunc) # (n_te x k) where k is ALL times
      GCY0[GCY0 == 0] <- min(GCY0[GCY0 > 0])
      GCY1[GCY1 == 0] <- min(GCY1[GCY1 > 0])
      
      # Dropout censoring N
      # Direct approach as opposed to data augmentation approach
      if (verbose) message("    [5] G_{CN}")
      rc <- as.integer(dat_tr$DeltaN == 0L & dat_tr$N < dat_tr$admin_t)
      C <- dat_tr$N
      # Prediction grid: both arms × test participants
      newAW_GN <- as.data.frame(cbind(A = dat_te$A, W_te, Y = dat_te$Y, DeltaY = dat_te$DeltaY))
      
      #USE CENS TO GET LEFT CONTINUOUS
      GN_fit <- survSuperLearner::survSuperLearner(
        time             = C,
        event            = 1-rc,
        X                = cbind(A = dat_tr$A, W_tr, Y=dat_tr$Y, DeltaY=dat_tr$DeltaY),
        newX             = newAW_GN,
        new.times        = eval.times,
        event.SL.library = unique(c(nuis$cens.SL.library, "survSL.rfsrc")),
        cens.SL.library  = unique(c(nuis$cens.SL.library, "survSL.rfsrc")),
        verbose          = FALSE,
        control          = nuis$survSL.control,
        cvControl        = nuis$survSL.cvControl
      )
      GN_pred <- GN_fit$cens.SL.predict          # (2*n_te x length(eval.times))
      if(0 %in% eval.times){
        GN_pred[, 1] <- 1  # S(0) = 1 by construction
      }
    }
    
    # ----------------------------------------------------------------
    # Step 7: Posterior f_N via Bayes rule on test fold
    # ----------------------------------------------------------------
    if (verbose) message("    [6] Posterior f_{0N}(n | a, x, delta_y, tilde{y}) via bayes inversion")
    post <- .compute_posterior_fN_2(
      dat_te   = dat_te,
      aug_te   = aug_te,
      SY0      = SY0,
      SY1      = SY1,
      survest_N = survest_N,
      fit.times  = fit.times,
      eval.times = eval.times,
      W_names    = W_names
    )
    
    # ----------------------------------------------------------------
    # Accumulate (fold order)
    # ----------------------------------------------------------------
    if (0L %in% fit.treat) {
      if(single_cens_mechanism==TRUE){
        eif0_fold<-.get.survival.AIPWCC(
          Y       = dat_te$Y,
          DeltaY  = dat_te$DeltaY,
          N       = dat_te$N,
          DeltaN  = dat_te$DeltaN,
          A       = 1L - dat_te$A,
          admin_t = dat_te$admin_t,
          fit.times  = fit.times,
          eval.times = eval.times,
          ids        = aug_te$id,
          S.hats     = SY0,                    # (n_aug_te x K) — fold-local only
          G.hats  = GC0,
          GN.hats = NULL,
          g.hats  = 1 - nuis$prop.pred[te_idx],
          fN_a0   = post$fN0,
          fN_alt  = post$fN1,
          isotonize = TRUE,
          verbose   = TRUE,
          single_cens_mechanism = TRUE)
      }else{
        eif0_fold<-.get.survival.AIPWCC(
          Y       = dat_te$Y,
          DeltaY  = dat_te$DeltaY,
          N       = dat_te$N,
          DeltaN  = dat_te$DeltaN,
          A       = 1L - dat_te$A,
          admin_t = dat_te$admin_t,
          fit.times  = fit.times,
          eval.times = eval.times,
          ids        = aug_te$id,
          S.hats     = SY0,                    # (n_aug_te x K) — fold-local only
          G.hats  = GCY0,
          GN.hats = GN_pred,
          g.hats  = 1 - nuis$prop.pred[te_idx],
          fN_a0   = post$fN0,
          fN_alt  = post$fN1,
          isotonize = TRUE,
          verbose   = TRUE,
          single_cens_mechanism = FALSE)  
      }
      if(any(is.na(eif0_fold$IF.vals))){
        print("NA values found")
        print(paste0(dat_te$id[which(is.na(eif0_fold$IF.vals[,ncol(eif0_fold$IF.vals)]))]))
        stop("Missing IF values")
      }
      
      acc$IF0 <- rbind(acc$IF0, eif0_fold$IF.vals)
    }
    
    if (1L %in% fit.treat) {
      if(single_cens_mechanism==TRUE){
        eif1_fold <- .get.survival.AIPWCC(
          Y       = dat_te$Y,
          DeltaY  = dat_te$DeltaY,
          N       = dat_te$N,
          DeltaN  = dat_te$DeltaN,
          A       = dat_te$A,
          admin_t = dat_te$admin_t,
          fit.times  = fit.times,
          eval.times = eval.times,
          ids        = aug_te$id,
          S.hats     = SY1,                    # (n_aug_te x K) — fold-local only
          G.hats  = GC1,
          GN.hats = NULL,
          g.hats  = nuis$prop.pred[te_idx],
          fN_a0   = post$fN1,
          fN_alt  = post$fN0,
          isotonize = TRUE,
          verbose   = TRUE,
          single_cens_mechanism = TRUE
        )
      }else{
        eif1_fold<-.get.survival.AIPWCC(
          Y       = dat_te$Y,
          DeltaY  = dat_te$DeltaY,
          N       = dat_te$N,
          DeltaN  = dat_te$DeltaN,
          A       = dat_te$A,
          admin_t = dat_te$admin_t,
          fit.times  = fit.times,
          eval.times = eval.times,
          ids        = aug_te$id,
          S.hats     = SY1,                    # (n_aug_te x K) — fold-local only
          G.hats  = GCY1,
          GN.hats = GN_pred,
          g.hats  = nuis$prop.pred[te_idx],
          fN_a0   = post$fN1,
          fN_alt  = post$fN0,
          isotonize = TRUE,
          verbose   = TRUE,
          single_cens_mechanism = FALSE
        )  
      }
      if(any(is.na(eif1_fold$IF.vals))){
        print("NA values found")
        print(paste0(dat_te$id[which(is.na(eif1_fold$IF.vals[,ncol(eif1_fold$IF.vals)]))]))
        stop("Missing IF values")
      }
      acc$IF1 <- rbind(acc$IF1, eif1_fold$IF.vals)
    }
    
    acc$test_ids <- c(acc$test_ids,     dat_te$id)
    
    rm(SY0, SY1, aug_tr, aug_te, post, SN_pred,
       eif0_fold, eif1_fold, GC_fit, GC_pred, GC0, GC1, GCY1, CGY0, GN_pred, SN_fit, SN0, SN1)
    gc()
    if (verbose) message(sprintf("  Fold %d complete", v))
  }
  
  # ------------------------------------------------------------------
  # Reindex ALL fold-accumulated matrices to ORIGINAL DATA ORDER
  # ------------------------------------------------------------------
  orig_ids     <- id
  fold_to_orig <- match(orig_ids, acc$test_ids)   # per-participant index
  
  make_surv_from_IF <- function(IF_mat, isotonize) {
    IF_mat    <- IF_mat[fold_to_orig, , drop=FALSE]
    surv      <- colMeans(IF_mat, na.rm=TRUE)
    IF_center <- sweep(IF_mat, 2, surv, "-")
    res <- list(
      times          = fit.times,
      surv           = pmin(1, pmax(0, surv)),
      IF.vals        = IF_mat,
      IF.vals.center = IF_center
    )
    if (isotonize) {
      idx <- !is.na(res$surv)
      res$surv.iso     <- res$surv
      res$surv.iso[idx] <- 1 - isoreg(fit.times[idx], 1 - res$surv[idx])$yf
    }
    res
  }
  
  if (verbose) message("Computing counterfactual survivals...")
  if (0L %in% fit.treat) surv.0 <- make_surv_from_IF(acc$IF0, isotonize=TRUE)
  if (1L %in% fit.treat) surv.1 <- make_surv_from_IF(acc$IF1, isotonize=TRUE)
  
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
      
      # .risk.ratio(fit.times=fit.times, 
      #             surv.0=surv.0$surv, 
      #             IF.vals.0 = surv.0$IF.vals.center, 
      #             surv.1=surv.1$surv, 
      #             IF.vals.1 = surv.1$IF.vals.center, 
      #             conf.band=conf.band, 
      #             band.end.pts=c(q01, q99), 
      #             conf.level=conf.level)$risk.ratio.df
      
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

survSL.rfsrc.fast <- function (time, event, X, newX, new.times, obsWeights, id, ...)
{
  data <- data.frame(time, event)
  data <- cbind(data, X)
  fit.rfsrc <- randomForestSRC::rfsrc.fast(Surv(time, event) ~ .,
                                           data = data, importance = FALSE, case.wt = obsWeights, forest=TRUE, save.memory = TRUE,
                                           ...)
  survs <- predict(fit.rfsrc, newdata = newX, importance = "none")$survival
  pred <- t(sapply(1:nrow(survs), function(i) {
    stats::approx(c(0, fit.rfsrc$time.interest), c(1, survs[i,
    ]), method = "constant", xout = new.times, rule = 2)$y
  }))
  fit <- list(object = fit.rfsrc)
  class(fit) <- c("survSL.rfsrc")
  out <- list(pred = pred, fit = fit)
  return(out)
}

# Flexible Cox

survSL.coxph.ps <- function(time, event, X, newX, new.times, cts.num = 10, obsWeights, ...) {
  
  bad <- !is.finite(time) |
    !is.finite(event) |
    !is.finite(obsWeights)
  
  if (any(bad)) {
    stop("Non-finite values detected")
  }
  
  if (any(obsWeights <= 0)) {
    stop("Non-positive weights detected")
  }
  
  # time=t
  # event=event
  # X=data.frame("X"=X)
  # obsWeights = rep(1,1000)
  
  # Identify continuous variables for spline smoothing
  cts.x <- which(apply(X, 2, function(x) length(unique(x)) > cts.num))
  
  # Build formula using pspline() instead of s()
  if (length(cts.x) > 0 && length(cts.x) < ncol(X)) {
    gam.model <- as.formula(paste(
      "Surv(time, event) ~",
      paste(paste0("pspline(", colnames(X[, cts.x, drop = FALSE]), ", df=3)"),
            collapse = " + "),
      "+",
      paste(colnames(X[, -cts.x, drop = FALSE]),
            collapse = " + ")
    ))
  } else if (length(cts.x) > 0) {
    
    gam.model <- as.formula(paste(
      "Surv(time, event) ~",
      paste(paste0("pspline(", colnames(X[, cts.x, drop = FALSE]), ", df=3)"),
            collapse = " + ")
    ))
    
  } else {
    
    gam.model <- as.formula(
      paste("Surv(time, event) ~",
            paste(colnames(X), collapse = " + "))
    )
  }
  
  X_fit       <- as.data.frame(X)
  X_fit$time  <- time
  X_fit$event <- event
  
  # Fit using standard coxph, which natively supports case weights
  fit.cox <- tryCatch(
    survival::coxph(gam.model, data = X_fit, weights = obsWeights, x = TRUE, y=TRUE),
    error = function(e) NULL,
    warning = function(w){
      if(grepl("singular|did not converge|NA/NaN", w$message, ignore.case=TRUE)){
        NULL
      } else {
        suppressWarnings(survival::coxph(gam.model, data = X_fit, weights = obsWeights, x = TRUE, y=TRUE))
      }
    }
  )
  
  fit_usable <- !is.null(fit.cox) && 
    !any(is.na(coef(fit.cox))) && 
    !any(is.infinite(coef(fit.cox)))
  
  if (!fit_usable) {
    # Fall back to linear Cox in N
    fit.cox <- coxph(as.formula(paste("Surv(time, event) ~", paste(colnames(X), collapse = " + "))), 
                     data = X_fit, weights = obsWeights, x = TRUE, y = TRUE)
  }
  
  # Predict survival curves for newX
  #sfit <- survival::survfit(fit.cox, newdata = as.data.frame(newX))
  
  # survfit returns a matrix of dimensions (times x newX rows)
  # Transpose to get (newX rows x times) and interpolate to new.times
  # pred <- t(apply(t(sfit$surv), 1, function(s) {
  #   stats::approx(
  #     x = sfit$time,
  #     y = s,
  #     xout = new.times,
  #     method = "linear",
  #     yleft  = 1
  #   )$y
  # }))
  pred <- t(summary(survival::survfit(fit.cox, newdata = newX, 
                                      se.fit = FALSE, conf.int = FALSE), times = new.times)$surv)
  
  if (ncol(pred) < length(new.times)) {
    pred <- cbind(pred, matrix(pred[, ncol(pred)], nrow = nrow(pred), 
                               ncol = length(new.times) - ncol(pred)))
  }
  
  fit <- list(object = fit.cox)
  class(fit) <- c("survSL.coxph.ps")
  # fit_obj <- list(
  #   fit       = fit.cox,
  #   new.times = new.times
  # )
  # class(fit_obj) <- "survSL.coxph.ps"
  
  return(list(pred = pred, fit = fit))
}

predict.survSL.coxph.ps <- function(object, newX, new.times, ...) {
  
  pred <- t(summary(survival::survfit(formula = object$object, 
                                      newdata = newX, se.fit = FALSE, conf.int = FALSE), 
                    times = new.times)$surv)
  
  if (ncol(pred) < length(new.times)) {
    pred <- cbind(pred, matrix(pred[, ncol(pred)], nrow = nrow(pred), 
                               ncol = length(new.times) - ncol(pred)))
  }
  # sfit <- survival::survfit(object$fit, newdata = as.data.frame(newX), se.fit=FALSE, conf.int=FALSE)
  # 
  # pred <- t(apply(t(sfit$surv), 1, function(s) {
  #   stats::approx(
  #     x = sfit$time,
  #     y = s,
  #     xout = new.times,
  #     method = "linear",
  #     yleft  = 1
  #   )$y
  # }))
  return(pred)
}

# DEFUNCT: I don't think this incorporates weights correctly.
survSL.gam.w <- function(time, event, X, newX, new.times, cts.num = 5,
                         obsWeights, ...) {
  
  # Build formula
  cts.x <- apply(X, 2, function(x) length(unique(x)) > cts.num)
  if (sum(cts.x) > 0 && sum(!cts.x) > 0) {
    gam.model <- as.formula(paste(
      "time ~",
      paste(paste0("s(", colnames(X[, cts.x, drop = FALSE]), ")"), collapse = "+"),
      "+",
      paste(colnames(X[, !cts.x, drop = FALSE]), collapse = "+")
    ))
  } else if (sum(cts.x) > 0) {
    gam.model <- as.formula(paste(
      "time ~",
      paste(paste0("s(", colnames(X[, cts.x, drop = FALSE]), ")"), collapse = "+")
    ))
  } else {
    gam.model <- as.formula(paste("time ~", paste(colnames(X), collapse = "+")))
  }
  
  X_fit       <- as.data.frame(X)
  X_fit$time  <- time
  X_fit$event <- event
  
  fit.gam <- mgcv::gam(gam.model,
                       family  = mgcv::cox.ph(),
                       data    = X_fit,
                       weights = event * obsWeights)
  
  # ----- Build a baseline survival curve via a single reference participant -----
  # Pick the column-wise mean of training X as the reference. Doesn't matter
  # which point we pick as long as we record the linear predictor at it.
  ref_X <- as.data.frame(lapply(X_fit[, setdiff(names(X_fit), c("time", "event")),
                                      drop = FALSE], function(col) {
                                        if (is.numeric(col)) mean(col) else as(names(sort(table(col), decreasing = TRUE))[1],
                                                                               class(col))
                                      }))
  
  # Evaluate baseline survival on new.times for the reference participant
  ref_grid       <- ref_X[rep(1L, length(new.times)), , drop = FALSE]
  ref_grid$time  <- new.times
  S_ref_on_grid  <- as.numeric(predict(fit.gam, newdata = ref_grid, type = "response"))
  
  # Compute the eta(x) - eta(ref) shift for each newX row:
  # Use lpmatrix predictions at any common time (S depends on time only via H0).
  # Trick: S(t | x) = S_ref(t)^{exp(eta(x) - eta(ref))}.
  # We need (eta(x) - eta(ref)) for each newX row. Get it via type="link" predictions
  # at a *single* time point (any time -- the link includes the same eta term).
  one_t <- new.times[1]
  newX_fixed_time          <- as.data.frame(newX)
  newX_fixed_time$time     <- one_t
  ref_fixed_time           <- ref_X
  ref_fixed_time$time      <- one_t
  
  link_newX <- as.numeric(predict(fit.gam, newdata = newX_fixed_time, type = "link"))
  link_ref  <- as.numeric(predict(fit.gam, newdata = ref_fixed_time,  type = "link"))
  
  eta_shift <- link_newX - link_ref   # length n_newX, the (eta(x) - eta(ref)) per row
  
  # Combine: S(t | x) = S_ref(t)^{exp(eta_shift)}.
  # Avoid 0^anything = 0 issue; clamp S_ref away from {0, 1}.
  S_ref_safe <- pmax(pmin(S_ref_on_grid, 1 - 1e-12), 1e-12)
  pred       <- t(outer(log(S_ref_safe), exp(eta_shift), "*"))
  pred       <- exp(pred)
  # pred[i, j] = S(new.times[j] | newX[i, ])
  
  fit_obj <- list(
    fit       = fit.gam,
    ref_X     = ref_X,
    new.times = new.times
  )
  class(fit_obj) <- "survSL.gam.w"
  
  list(pred = pred, fit = fit_obj)
}

predict.survSL.gam.w <- function(object, newX, new.times, ...) {
  fit.gam <- object$fit
  ref_X   <- object$ref_X
  
  # Baseline survival at new.times for the reference participant
  ref_grid       <- ref_X[rep(1L, length(new.times)), , drop = FALSE]
  ref_grid$time  <- new.times
  S_ref_on_grid  <- as.numeric(predict(fit.gam, newdata = ref_grid, type = "response"))
  
  # Linear predictor shift for newX
  one_t                    <- new.times[1]
  newX_fixed_time          <- as.data.frame(newX)
  newX_fixed_time$time     <- one_t
  ref_fixed_time           <- ref_X
  ref_fixed_time$time      <- one_t
  
  link_newX <- as.numeric(predict(fit.gam, newdata = newX_fixed_time, type = "link"))
  link_ref  <- as.numeric(predict(fit.gam, newdata = ref_fixed_time,  type = "link"))
  
  eta_shift <- link_newX - link_ref
  
  # Combine
  S_ref_safe <- pmax(pmin(S_ref_on_grid, 1 - 1e-12), 1e-12)
  pred <- t(outer(log(S_ref_safe), exp(eta_shift), "*"))
  exp(pred)
}

# Idea: used smoothed random forest
survSL.rfsrc.smooth <- function(time, event, X, newX, new.times, obsWeights, id,
                                ntree = 500,
                                nodesize = NULL,
                                mtry = NULL,
                                smooth_h = 30,
                                ...) {
  
  # time=aug_tr$Y
  # event=aug_tr$DeltaY
  # X = aug_tr %>% dplyr::select(X, N, DeltaN)
  # ntree = 500
  # nodesize = NULL
  # mtry = NULL
  # newX = aug_te[1:20,]
  # new.times=fit.times
  # id = aug_tr$id
  # obsWeights = aug_tr$w
  # smooth_h=14
  
  # ----- Tune nodesize and mtry based on training set size -----
  n_train <- length(time)
  
  if (is.null(nodesize)) {
    # Force much larger leaves than the default 15
    # Rule of thumb: nodesize ~ sqrt(n) gives good smoothing for survival
    # Floor at 50 to avoid tiny-leaf instability
    nodesize <- max(50L, floor(sqrt(n_train)))
  }
  
  if (is.null(mtry)) {
    # Use a smaller fraction of variables per split for stability
    # Default is sqrt(p); we use p/3 for survival to encourage smoother trees
    mtry <- max(1L, floor(ncol(X) / 3))
  }
  
  # ----- Fit RF -----
  train_data <- data.frame(time = time, event = event, X)
  
  fit <- randomForestSRC::rfsrc(
    Surv(time, event) ~ .,
    data           = train_data,
    case.wt        = obsWeights,
    ntree          = ntree,
    nodesize       = nodesize,
    mtry           = mtry,
    samptype       = "swor",       # sampling without replacement is more stable
    sampsize       = function(x) min(x, max(100, floor(0.632 * x))),
    splitrule      = "logrank",
    importance     = "none",       # skip importance for speed
    save.memory    = TRUE#,
    #...
  )
  
  # ----- Predict on newX at the original RF event time grid -----
  #pred_obj  <- predict(fit, newdata = as.data.frame(newX), na.action = "na.impute", importance="none")$survival
  rf_times  <- fit$time.interest
  rf_surv   <- predict(fit, newdata = as.data.frame(newX), na.action = "na.impute", importance="none")$survival
  
  # ----- Smooth each row's survival curve using Gaussian kernel on the time axis -----
  # The convolution preserves S(0) ≈ 1 and decays smoothly
  #smooth_h = 30
  if(is.null(smooth_h)){
    mu = sum(X$N * obsWeights)/sum(obsWeights)
    smooth_h = 1.06 * sqrt(sum(obsWeights * (X$N - mu)^2)/sum(obsWeights)) * (sum(obsWeights)^2/sum(obsWeights^2))^(-1/5)
  }
  
  smooth_row <- function(s) {
    if (length(rf_times) < 2L) return(s)
    smoothed <- sapply(seq_along(rf_times), function(j) {
      w <- dnorm(rf_times - rf_times[j], sd = smooth_h)
      sum(w * s) / sum(w)
    })
    # Enforce monotonicity (right-continuous, non-increasing)
    rev(cummax(rev(smoothed)))
  }
  
  rf_surv_smooth <- t(apply(rf_surv, 1, smooth_row))
  
  # ----- Interpolate smoothed predictions onto new.times grid -----
  pred <- t(apply(rf_surv_smooth, 1, function(s) {
    approx(x      = rf_times,
           y      = s,
           xout   = new.times,
           method = "linear",
           yleft  = 1,
           yright = min(s),
           ties   = "ordered")$y
  }))
  
  # Ensure predictions are in [0, 1] and monotone in t
  pred <- pmin(pmax(pred, 0), 1)
  pred <- t(apply(pred, 1, function(s) rev(cummax(rev(s)))))
  
  fit_obj <- list(
    fit      = fit,
    smooth_h = smooth_h,
    rf_times = rf_times
  )
  class(fit_obj) <- "survSL.rfsrc.smooth"
  
  list(pred = pred, fit = fit_obj)
}

predict.survSL.rfsrc.smooth <- function(object, newX, new.times, ...) {
  fit       <- object$fit
  smooth_h  <- object$smooth_h
  
  pred_obj  <- predict(fit, newdata = as.data.frame(newX), na.action = "na.impute")
  rf_times  <- pred_obj$time.interest
  rf_surv   <- pred_obj$survival
  
  smooth_row <- function(s) {
    if (length(rf_times) < 2L) return(s)
    smoothed <- sapply(seq_along(rf_times), function(j) {
      w <- dnorm(rf_times - rf_times[j], sd = smooth_h)
      sum(w * s) / sum(w)
    })
    rev(cummin(rev(smoothed)))
  }
  
  rf_surv_smooth <- t(apply(rf_surv, 1, smooth_row))
  
  pred <- t(apply(rf_surv_smooth, 1, function(s) {
    approx(x      = rf_times,
           y      = s,
           xout   = new.times,
           method = "linear",
           yleft  = 1,
           yright = min(s),
           ties   = "ordered")$y
  }))
  
  pred <- pmin(pmax(pred, 0), 1)
  t(apply(pred, 1, function(s) rev(cummin(rev(s)))))
}


.CFsurvival.nuis.opt <- function(cross.fit = TRUE, V = ifelse(cross.fit, 10, 1), folds = NULL, eval.times = NULL,
                                 event.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) c(alg, "survscreen.glmnet", "survscreen.marg", "All") ),  event.pred.0 = NULL, event.pred.1 = NULL,
                                 event.SL.library_SY = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) c(alg, "survscreen.glmnet", "survscreen.marg", "All") ),
                                 cens.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) c(alg, "survscreen.glmnet", "survscreen.marg", "All") ),  cens.trunc=0, cens.pred.0 = NULL, cens.pred.1 = NULL,
                                 survSL.control = list(initWeightAlg = "survSL.rfsrc", verbose=FALSE), survSL.cvControl = list(V = 10), save.nuis.fits = FALSE,
                                 prop.SL.library = lapply(c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger", "SL.xgboost"), function(alg) c(alg, "screen.glmnet", "screen.corRank", "All") ), prop.trunc=0, prop.pred = NULL,
                                 verbose=FALSE) {
  list(cross.fit = cross.fit, V = V, folds = folds, eval.times = eval.times,
       event.SL.library = event.SL.library, event.pred.0 = event.pred.0, event.pred.1 = event.pred.1,
       cens.SL.library = cens.SL.library, 
       event.SL.library_SY = event.SL.library_SY,
       cens.trunc=cens.trunc, 
       cens.pred.0 = cens.pred.0, 
       cens.pred.1 = cens.pred.1,
       survSL.control = survSL.control,
       survSL.cvControl = survSL.cvControl,
       save.nuis.fits = save.nuis.fits,
       prop.SL.library = prop.SL.library,  prop.trunc=prop.trunc, prop.pred = prop.pred, verbose=verbose)
}

# This function only ouputs the *event* probabilities (ignore censoring)
.pSSL_event_only <- function (object, newdata, new.times, X = NULL, time = NULL, 
                              event = NULL, onlySL = FALSE, threshold = 1e-04, ...){
  
  if (missing(newdata)) {
    out <- list(survpred = object$survSL.predict, 
                survlibrary.predict = object$survlibrary.predict)
    return(out)
  }
  if (!object$control$saveFitLibrary) {
    stop("This SuperLearner fit was created using control$saveFitLibrary = FALSE, so new predictions cannot be made.")
  }
  
  event.k <- nrow(object$event.libraryNames) # 1
  event.pred <- array(0, dim = c(nrow(newdata), length(new.times), event.k))
  dimnames(event.pred)[[3]] <- apply(object$event.libraryNames, 1, paste, collapse = "_")
  dimnames(event.pred)[[2]] <- new.times
  
  if (onlySL) {
    event.whichLibrary <- which(object$event.coef > threshold)
    event.coef <- object$event.coef
    event.coef[-event.whichLibrary] <- 0
    event.coef <- event.coef/sum(event.coef)
  } else {
    event.whichLibrary <- seq(event.k)
    event.coef <- object$event.coef
  }
  for (mm in event.whichLibrary) {
    newdataMM <- subset(newdata, select = object$event.whichScreen[object$event.SL.library$library[mm, 2], ])
    XMM <- if (is.null(X)) {
      NULL
    }
    else {
      subset(X, select = object$event.whichScreen[object$event.SL.library$library[mm, 2], ])
    }
    event.pred[, , mm] <- do.call(paste0("predict", ".", class(object$event.fitLibrary[[mm]])),
                                  list(object = object$event.fitLibrary[[mm]], 
                                       newX = newdataMM, new.times = new.times, X = XMM, 
                                       time = time, event = event, ...))
  }
  event.SL.predict <- matrix(NA, nrow = nrow(newdata), ncol = length(new.times))
  if(event.k==1){
    for (j in seq(length(new.times))) {
      event.SL.predict[, j] <- event.pred[, j, ] * event.coef
    }
  }else{
    for (j in seq(length(new.times))) {
      event.SL.predict[, j] <- event.pred[, j, ] %*% event.coef
    }
  }
  
  out <- list(event.SL.predict = event.SL.predict, event.library.predict = event.pred)
  return(out)
}

# ------------------------------------------------------------------
# .compute_aug_likelihood
#
#   Compute L(Y_tilde, Delta_Y | A, X, N) for each row of aug_tr using a
#   fitted survSuperLearner2 object. For uncensored Y we approximate
#   f_Y(y | A, X, N) by the discrete mass S_Y(t_{k-1}) - S_Y(t_k) on the
#   fit.times grid; for censored Y we use S_Y(Y_tilde | A, X, N) directly.
#
# Returns: numeric vector of length nrow(aug_tr).
# ------------------------------------------------------------------
.compute_aug_likelihood <- function(SY_fit, aug_tr, fit.times, W_names) {
  
  newdat <- cbind(A = aug_tr$A,
                  aug_tr[, c("N", W_names), drop = FALSE])
  
  SY_pred <- .pSSL_event_only(
    object    = SY_fit,
    newdata   = newdat,
    new.times = fit.times,
    onlySL    = TRUE
  )$event.SL.predict
  # SY_pred is (nrow(aug_tr) x K)
  
  K <- length(fit.times)
  #Y_idx <- findInterval(aug_tr$Y, c(0, fit.times), rightmost.closed = FALSE) #FIX
  Y_idx <- findInterval(aug_tr$Y, fit.times, left.open = TRUE)+1
  Y_idx <- pmax(1L, pmin(Y_idx, K))
  
  n_aug <- nrow(aug_tr)
  L_vec <- numeric(n_aug)
  
  # Vectorize the two cases for speed
  is_event <- aug_tr$DeltaY == 1L
  
  # Censored Y: L = S_Y(Y_tilde | A, X, N)
  if (any(!is_event)) {
    idx_c   <- which(!is_event)
    L_vec[idx_c] <- SY_pred[cbind(idx_c, Y_idx[idx_c])]
  }
  
  # Observed event: L = S_Y(t_{yi-1}) - S_Y(t_yi)
  if (any(is_event)) {
    idx_e <- which(is_event)
    yi_e  <- Y_idx[idx_e]
    S_hi  <- SY_pred[cbind(idx_e, yi_e)]
    S_lo  <- ifelse(yi_e == 1L, 1, SY_pred[cbind(idx_e, pmax(yi_e - 1L, 1L))])
    L_vec[idx_e] <- pmax(0, S_lo - S_hi)
  }
  
  L_vec
}


# ------------------------------------------------------------------
# .update_aug_weights
#
#   Given a previous S_Y fit, update aug_tr$w to be proportional to
#   L(Y_tilde, Delta_Y | A, X, N) * w0, renormalized per participant.
#   Stores the original prior weights in aug_tr$w0 on first call so
#   subsequent iterations multiply against the *same* prior, not the
#   running posterior.
#
#   For single-row participants (DeltaN_eff = 1, no augmentation), w
#   stays at 1.
# ------------------------------------------------------------------
.update_aug_weights <- function(aug_tr, SY_fit_prev, fit.times, W_names,
                                floor_w = 1e-12) {
  
  L_vec <- .compute_aug_likelihood(SY_fit_prev, aug_tr, fit.times, W_names)
  
  # Save original prior weights only on first call
  if (!"w0" %in% colnames(aug_tr)) aug_tr$w0 <- aug_tr$w
  
  # Per-participant renormalization
  ids        <- as.character(aug_tr$id)
  rows_by_id <- split(seq_along(ids), ids)
  
  w_new <- aug_tr$w
  
  for (uid in names(rows_by_id)) {
    rows <- rows_by_id[[uid]]
    if (length(rows) == 1L) {
      w_new[rows] <- 1
    } else {
      w_raw <- L_vec[rows] * aug_tr$w0[rows]
      w_sum <- sum(w_raw)
      if (w_sum > floor_w) {
        w_new[rows] <- w_raw / w_sum
      } else {
        # Likelihood collapsed to ~0 on every grid point: fall back to prior
        prior_sum <- sum(aug_tr$w0[rows])
        w_new[rows] <- if (prior_sum > 0) aug_tr$w0[rows] / prior_sum
        else rep(1 / length(rows), length(rows))
      }
    }
  }
  
  aug_tr$w <- w_new
  aug_tr
}

# ------------------------------------------------------------------
# .fit_SY_em
#
#   EM-style fit of S_Y on augmented training data, with iterative
#   refinement of redistribute-to-right weights.
#
# Arguments
#   aug_tr, aug_te, dat_te : per-fold dataframes, as produced by
#       .augment_fold and .expand_test in AIPWCC_survival_2.r
#   fit.times              : event times grid (length K)
#   W_names                : names of baseline covariates
#   nuis                   : nuisance options object (CFsurvival.nuisance.options)
#   max_iter               : maximum number of S_Y fits. max_iter=1 disables
#                            EM (matches original behavior). max_iter=2 does
#                            exactly one weight update and one refit.
#   tol                    : convergence threshold on max |w^(t+1) - w^(t)|.
#   fast_iter_library      : optional simpler SL library used for non-final
#                            iterations. If NULL, use full library throughout.
#                            Final fit always uses full library.
#                            E.g., c("survSL.coxph") or c("survSL.expreg").
#   verbose                : print progress
#
# Returns: list with
#   SY_fit         : final survSuperLearner2 fit
#   SY0, SY1       : predictions on aug_te under arms 0 and 1
#                    (matrices matching the originals)
#   aug_tr_final   : augmented training data with final weights and w0
#   history        : list of per-iteration diagnostics
#                    (iter, max_w_diff, n_aug, library)
#   converged      : TRUE if max_w_diff < tol within max_iter iterations
#   n_iter         : total number of S_Y fits performed
# ------------------------------------------------------------------
.fit_SY_em <- function(
    aug_tr, 
    aug_te, 
    dat_te,
    fit.times, 
    W_names, 
    nuis,
    max_iter          = 2L,
    tol               = 0.05,
    fast_iter_library = NULL,
    verbose           = TRUE,
    logN              = FALSE,
    weight_norm       = c("max", "L2"),
    ...
) {
  
  stopifnot(max_iter >= 1L)
  
  AW_ev_tr <- cbind(A = aug_tr$A, aug_tr[, c("N", W_names), drop = FALSE])
  if(logN==TRUE){
    AW_ev_tr$N <- log(AW_ev_tr$N)
  }
  AW_ce_tr <- cbind(A = aug_tr$A, aug_tr[, W_names, drop = FALSE])
  
  W_te     <- dat_te[, W_names, drop = FALSE]
  AW_ce_te <- rbind(cbind(A = 0L, W_te), cbind(A = 1L, W_te))
  
  full_event_lib <- nuis$event.SL.library_SY[nuis$event.SL.library_SY != "survSL.km"]
  full_cens_lib  <- nuis$cens.SL.library
  
  has_fast <- !is.null(fast_iter_library) && length(fast_iter_library) >= 1L
  
  # Helper: one survSuperLearner2 fit using current aug_tr
  fit_once <- function(event_lib, cens_lib) {
    survSuperLearner2(
      time             = aug_tr$Y,
      event            = aug_tr$DeltaY,
      X_event          = AW_ev_tr,
      X_cens           = AW_ce_tr,
      newX_event       = cbind(A = 0L,
                               aug_te[, c("N",W_names), drop = FALSE] %>%
                                 dplyr::slice(1:2)),
      newX_cens        = AW_ce_te[1:2, ],
      new.times        = fit.times,
      event.SL.library = event_lib,
      cens.SL.library  = cens_lib,
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl,
      obsWeights       = aug_tr$w,
      id               = aug_tr$id
    )
  }
  
  history    <- list()
  SY_fit     <- NULL
  converged  <- FALSE
  
  update_indices <- which(duplicated(aug_tr$id))
  
  for (iter in seq_len(max_iter)) {
    
    is_final  <- (iter == max_iter)
    use_fast  <- has_fast && !is_final
    event_lib <- if (use_fast) fast_iter_library else full_event_lib
    cens_lib  <- if (use_fast) fast_iter_library else full_cens_lib
    
    if (verbose) message(sprintf(
      "    [EM iter %d/%d] Fitting S_Y (n_aug=%d, library=%s)",
      iter, max_iter, nrow(aug_tr), if (use_fast) "fast" else "full"
    ))
    
    w_in   <- aug_tr$w   # weights in this fit
    SY_fit <- fit_once(event_lib, cens_lib)
    
    if (is_final) {
      history[[iter]] <- list(
        iter       = iter,
        max_w_diff = NA_real_,
        n_aug      = nrow(aug_tr),
        library    = if (use_fast) "fast" else "full",
        note       = "final iter; no further weight update"
      )
      break
    }
    
    # E-step: update weights using current SY_fit
    aug_tr_new <- .update_aug_weights(aug_tr, SY_fit, fit.times, W_names)
    if(weight_norm=="max"){
      
      max_w_diff <- max(abs(aug_tr_new$w - w_in), na.rm = TRUE)
      
      if (verbose) message(sprintf(
        "      max |w^(t+1) - w^(t)| = %.4g  (tol = %.4g)", max_w_diff, tol
      ))
      
    }else if(weight_norm=="L2"){
      
      max_w_diff <- sqrt(mean((aug_tr_new$w[update_indices] - w_in[update_indices])^2)/mean(w_in[update_indices]^2))
      
      if (verbose) message(sprintf(
        "      sqrt(SUM((w^(t+1) - w^(t))^2)/ SUM(w^(t)^2)) = %.4g  (tol = %.4g)", max_w_diff, tol
      ))
      
    }
    
    history[[iter]] <- list(
      iter       = iter,
      max_w_diff = max_w_diff,
      n_aug      = nrow(aug_tr),
      library    = if (use_fast) "fast" else "full"
    )
    
    if (max_w_diff < tol) {
      converged <- TRUE
      
      if (use_fast) {
        # Refit once more with the full library on the converged weights,
        # so the SY_fit returned is library-consistent with the final SE.
        if (verbose) message(sprintf(
          "    [EM iter %d/%d] Converged. Refitting on full library...",
          iter, max_iter
        ))
        aug_tr <- aug_tr_new
        SY_fit <- fit_once(full_event_lib, full_cens_lib)
      } else if (verbose) {
        message(sprintf("    [EM iter %d/%d] Converged.", iter, max_iter))
      }
      break
    }
    
    aug_tr <- aug_tr_new
  }
  
  # Final predictions on aug_te under both arms
  SY0 <- .pSSL_event_only(
    object    = SY_fit,
    newdata   = cbind(A = 0L, aug_te[, c("N",W_names), drop = FALSE]),
    new.times = fit.times,
    onlySL    = TRUE
  )$event.SL.predict
  
  SY1 <- .pSSL_event_only(
    object    = SY_fit,
    newdata   = cbind(A = 1L, aug_te[, c("N",W_names), drop = FALSE]),
    new.times = fit.times,
    onlySL    = TRUE
  )$event.SL.predict
  
  list(
    #SY_fit       = SY_fit,
    SY0          = SY0,
    SY1          = SY1,
    aug_tr_final = aug_tr,
    history      = history,
    converged    = converged,
    n_iter       = length(history)
  )
}

# ==================================================================
# HELPER: augment training fold with redistribute-to-right weights
# NOTE: this is global augmentation that expands all censored observations (dropout & administrative) to values to the right
# ==================================================================

.augment_fold <- function(dat, W_names, eval.times, survest_N) {
  
  rows <- vector("list", nrow(dat))
  
  for (i in seq_len(nrow(dat))) {
    ri          <- dat[i, ]
    Ai          <- ri$A
    N0          <- ri$N
    #admin_t_i   <- ri$admin_t
    xi          <- ri[, W_names, drop=FALSE]
    
    # Unique event times in same arm + admin_t as support
    # N_support <- sort(unique(c(eval.times, admin_t_i))) # FIXED
    
    base <- data.frame(id=ri$id, Y=ri$Y, DeltaY=ri$DeltaY, A=Ai)
    base <- cbind(base, xi)
    
    # Observed N → single row, weight 1
    if (ri$DeltaN == 1L) {
      rows[[i]] <- cbind(base, N=N0, w=1)
      next
    }
    
    # ELSE IF N censored → redistribute to right
    N_grid <- eval.times[eval.times > N0]
    if (length(N_grid) == 0L) {
      rows[[i]] <- cbind(base, N=N0, w=1)
      next
    }
    
    S_all <- survest_N(t = c(N0, N_grid), x = xi) #S(lhs, Ngrid)
    S_N0_minus <- S_all[1]
    S_grid     <- S_all[-1] #right hand limits
    S_lag2     <- c(S_N0_minus, S_grid[-length(S_grid)]) #left hand limits
    w_raw <- pmax(0, S_lag2 - S_grid)
    w_sum <- sum(w_raw)
    w     <- if (w_sum > 0) w_raw / w_sum else rep(1/length(N_grid), length(N_grid))
    
    r <- cbind(base[rep(1L, length(N_grid)), ],
               N      = N_grid,
               w      = w)
    rows[[i]] <- r[w > 0, ]
  }
  
  do.call(rbind, rows)
}

# ==================================================================
# HELPER: expand test fold over eval.times grid
# ==================================================================

.expand_test <- function(dat_te, eval.times, W_names) {
  l <- length(eval.times)
  rows <- vector("list", nrow(dat_te))
  for (i in seq_len(nrow(dat_te))) {
    ri    <- dat_te[i, ]
    base  <- ri[rep(1L, l), c("id","Y","DeltaY","A", W_names), drop = FALSE]
    base$N      <- eval.times
    rows[[i]]   <- base
  }
  do.call(rbind, rows)
}

# ==================================================================
# HELPER: Bayes-rule posterior f_N on test fold
# ==================================================================
.compute_posterior_fN_2 <- function(dat_te, aug_te, SY0, SY1,
                                    survest_N, fit.times, eval.times, W_names) {
  
  n_te <- nrow(dat_te)
  l    <- length(eval.times) #NOTE THAT l includes the 0 appended to the left hand side of eval.times, before which there are no events
  fN0  <- matrix(NA_real_, n_te, l-1)
  fN1  <- matrix(NA_real_, n_te, l-1)
  times_mass <- eval.times[eval.times>0] #Times to assign N mass
  stopifnot(length(times_mass)==(l-1))
  
  # Pre-build id → aug_te row index lookup
  id_to_aug <- split(seq_len(nrow(aug_te)), aug_te$id)
  
  # Helper: find last fit.time <= t_val (right-continuous step function index)
  find_idx <- function(t_val, grid) {
    idx <- findInterval(t_val, grid, left.open=TRUE)+1 #FIXED
    pmax(1L, pmin(idx, length(grid)))
  }
  
  for (i in seq_len(n_te)) {
    Yi      <- dat_te$Y[i]
    dYi     <- dat_te$DeltaY[i]
    idi     <- dat_te$id[i]
    xi      <- dat_te[i, W_names, drop=FALSE]
    admin_t_i <- dat_te$admin_t[i]
    aug_idx <- id_to_aug[[as.character(idi)]]
    
    # Prior: discrete PMF of N on eval.times under each arm
    SN1 <- survest_N(t=eval.times, x=xi) # predict survival over eval.times grid (l x 1)
    SN0 <- survest_N(t=eval.times, x=xi) # predict survival over eval.times grid (l x 1)
    pr1_i <- pmax(0, (SN1[-l] - SN1[-1])) # length (l-1)
    pr0_i <- pmax(0, (SN0[-l] - SN0[-1])) # length (l-1)
    # Append remaining tail mass to final observed time
    pr1_i[l-1] <- pr1_i[l-1] + SN1[l]
    pr0_i[l-1] <- pr0_i[l-1] + SN0[l]
    
    # Likelihood: L(Yi, dYi | A=a, X, N=n_j) for each aug_te row (= N grid pt)
    if (dYi == 0L) {
      # Censored: numerator = S_{0Y}(Yi | a,x,n)
      y_idx <- find_idx(Yi, fit.times)
      lk1 <- SY1[aug_idx, y_idx][-1]          # length l-1
      lk0 <- SY0[aug_idx, y_idx][-1]
    } else {
      # Event: numerator = f_{0Y}(Yi | a,x,n) ≈ S(t_{j-1}) - S(t_j)
      y_lo <- pmin(find_idx(Yi, fit.times) + 1L, length(fit.times))
      f_1 <- t(sapply(seq_along(aug_idx), function(x) -diff(c(1, SY1[aug_idx[x], ]))))
      f_0 <- t(sapply(seq_along(aug_idx), function(x) -diff(c(1, SY0[aug_idx[x], ]))))
      
      pos_indx <- which(colSums(f_1)>0)
      pos_times <- fit.times[pos_indx]
      closest_gridpoint_indx <- which.min(abs(pos_times - Yi))
      lk1 <- f_1[,pos_indx[closest_gridpoint_indx]][-1]
      
      pos_indx <- which(colSums(f_0)>0)
      pos_times <- fit.times[pos_indx]
      closest_gridpoint_indx <- which.min(abs(pos_times - Yi))
      lk0 <- f_0[,pos_indx[closest_gridpoint_indx]][-1]
    }
    
    finish <- function(lk, pr) {
      u <- lk * pr
      u[!is.finite(u) | u < 0] <- 0
      d <- sum(u)
      if (d > 0) u / d
      else if (sum(pr) > 0) pr / sum(pr)
      else rep(1 / length(pr), length(pr))
    }
    fN1[i, ] <- finish(lk1, pr1_i)
    fN0[i, ] <- finish(lk0, pr0_i)
  }
  
  list(fN0=fN0, fN1=fN1)
}


.compute_posterior_fN <- function(dat_te, aug_te, SY0, SY1,
                                  survest_N, fit.times, eval.times, W_names) {
  
  n_te <- nrow(dat_te)
  l    <- length(eval.times)
  fN0  <- matrix(NA_real_, n_te, l-1)
  fN1  <- matrix(NA_real_, n_te, l-1)
  
  # Pre-build id → aug_te row index lookup
  id_to_aug <- split(seq_len(nrow(aug_te)), aug_te$id)
  
  # Helper: find last fit.time <= t_val (right-continuous step function index)
  find_idx <- function(t_val, grid) {
    #old idx <- findInterval(t_val, grid, rightmost.closed=FALSE)
    idx <- findInterval(t_val, grid, left.open=TRUE)+1 #FIXED
    pmax(1L, pmin(idx, length(grid)))
  }
  
  for (i in seq_len(n_te)) {
    Yi      <- dat_te$Y[i]
    dYi     <- dat_te$DeltaY[i]
    idi     <- dat_te$id[i]
    xi      <- dat_te[i, W_names, drop=FALSE]
    admin_t_i <- dat_te$admin_t[i]
    aug_idx <- id_to_aug[[as.character(idi)]]
    
    # Prior: discrete PMF of N on eval.times under each arm
    SN1 <- c(survest_N(t=eval.times, x=xi))
    SN0 <- c(survest_N(t=eval.times, x=xi))
    pr1_i <- pmax(0, (SN1[-length(SN1)] - SN1[-1]))
    pr0_i <- pmax(0, (SN0[-length(SN0)] - SN0[-1]))
    
    # This removes positive posterior mass from values greater than the administrative censoring time
    #pr1_i[eval.times > admin_t_i] <- 0
    #pr0_i[eval.times > admin_t_i] <- 0
    # Renormalize
    #pr1_i <- pr1_i / sum(pr1_i)
    #pr0_i <- pr0_i / sum(pr0_i)
    
    # Likelihood: L(Yi, dYi | A=a, X, N=n_j) for each aug_te row (= N grid pt)
    if (dYi == 0L) {
      # Censored Y: likelihood = S_Y(Yi | A, X, N)
      y_idx <- find_idx(Yi, fit.times)
      lk1   <- SY1[aug_idx, y_idx]
      lk0   <- SY0[aug_idx, y_idx]
    } else {
      # Observed event: likelihood = S_Y(t_{j-1}) - S_Y(t_j)
      y_lo  <- find_idx(Yi, fit.times)
      #S_lo1 <- if (y_lo == 1L) rep(1, length(aug_idx)) else SY1[aug_idx, y_lo-1L]
      #S_lo0 <- if (y_lo == 1L) rep(1, length(aug_idx)) else SY0[aug_idx, y_lo-1L]
      f_1 <- sapply(1:length(aug_idx), FUN=function(x){-diff(c(1, SY1[aug_idx[x], ]))})
      f_0 <- sapply(1:length(aug_idx), FUN=function(x){-diff(c(1, SY0[aug_idx[x], ]))})
      
      lk1   <- f_1[y_lo,]
      lk0   <- f_0[y_lo,]
      #if (all(lk1 == 0)) lk1 <- pr1_i
      #if (all(lk0 == 0)) lk0 <- pr0_i
    }
    
    # Posterior ∝ likelihood × prior
    u1 <- lk1[-1] * pr1_i
    d1 <- sum(u1)
    if (d1 > 0) {
      fN1[i, ] <- u1 / d1
    } else if (sum(pr1_i) > 0) { #fallback to prior if divide by 0 error
      fN1[i, ] <- pr1_i / sum(pr1_i)
    } else {
      fN1[i, ] <- rep(1 / (l - 1), l - 1)
    }
    
    u0 <- lk0[-1] * pr0_i
    d0 <- sum(u0)
    if (d0 > 0) {
      fN0[i, ] <- u0 / d0
    } else if (sum(pr0_i) > 0) { #fallback to prior if divide by 0 error
      fN0[i, ] <- pr0_i / sum(pr0_i)
    } else {
      fN0[i, ] <- rep(1 / (l - 1), l - 1)
    }
  }
  
  list(fN0=fN0, fN1=fN1)
}


# ==================================================================
# MAIN EIF COMPUTATION
# ==================================================================
.get.survival.AIPWCC <- function(
    Y, DeltaY, N, DeltaN,
    A,                       # A[i] = 1 iff participant in target arm a0
    admin_t,
    fit.times, eval.times,
    ids,
    S.hats,
    G.hats, # G_Y predictions for compensator and inner1
    GN.hats, #G_N predictions : (n x l) for the participant's target arm.
    g.hats,                  # length n -- pi(a0 | x)
    fN_a0, fN_alt,           # n x l posterior matrices
    isotonize = TRUE,
    verbose   = TRUE,
    single_cens_mechanism = TRUE) {
  
  # Y       = dat_te$Y;
  # DeltaY  = dat_te$DeltaY;
  # N       = dat_te$N;
  # DeltaN  = dat_te$DeltaN;
  # A       = 1L - dat_te$A;
  # admin_t = dat_te$admin_t;
  # fit.times  = fit.times;
  # eval.times = eval.times;
  # ids        = aug_te$id;
  # S.hats     = SY0;                # (n_aug_te x K) — fold-local only
  # G.hats  = GCY0;
  # GN.hats = GN_pred;
  # g.hats  = 1 - nuis$prop.pred[te_idx];
  # fN_a0   = post$fN0;
  # fN_alt  = post$fN1;
  # isotonize = TRUE;
  # verbose   = TRUE;
  # single_cens_mechanism = FALSE
  
  fit.times <- sort(fit.times[fit.times > 0])
  n  <- length(Y)
  K  <- length(fit.times)
  l  <- length(eval.times)
  #Y_ind <- pmax(1L, pmin(findInterval(Y, fit.times), K))
  N_ind <- pmax(1L, pmin(findInterval(N, eval.times), l))
  
  IF.vals <- matrix(NA_real_, n, K)
  
  # Fast aug-row lookup
  id_to_rows <- split(seq_along(ids), ids)
  unique_ids <- unique(ids)
  
  IF.vals <- matrix(NA_real_, n, K)
  MARTS <- matrix(NA_real_, n, K)
  
  if (verbose) {
    pb <- progress_bar$new(n, format="  EIF [:bar] :percent ETA: :eta")
  }
  
  for (i in seq_len(n)) {
    N_indx = N_ind[i] # ONLY USED FOR INDEXING INTO IF & Q_SEQ
    fN.i  <- if (A[i] == 1L) fN_a0[i,] else fN_alt[i,] # l x 1
    
    # S_Y rows for this participant
    rows_i   <- id_to_rows[[as.character(unique_ids[i])]]
    Shi      <- S.hats[rows_i,, drop=FALSE]   # (l x K)
    
    # G_C^Y step function (for phi^F: vals.i and G_y)
    GYi_vec  <- G.hats[i,]
    # CONFIRMED: this code gives the left hand limit at time t
    #GYi_fun_lhl  <- stepfun(sort(unique(c(eval.times, fit.times))), c(1, GYi_vec), right=TRUE)
    # CONFIRMED: this code gives GY[i]
    GYi_fun_lhl  <- stepfun(sort(unique(c(eval.times, fit.times))), c(1, GYi_vec), right=FALSE)
    # This iwll be used inside the integral calculation and nowhere else
    # To be consistent with Westling
    fit.times.indx <- sort(which(sort(unique(c(eval.times, fit.times))) %in% fit.times))
    GYi_int <- GYi_vec[fit.times.indx] # length K actual censoring survival times
    
    # phi^F compensator
    # Based on cancellation algebra and discrete hazard approximation
    # FOLLOWING WESTLING, we should use the right hand limits in the integral
    S_inv  <- 1 / Shi
    S_inv[is.infinite(S_inv)] <- 0
    # #Matches Ted's diff, skips the [0, t_1] baseline hazard
    dS_inv <- S_inv[, -1, drop=FALSE] - S_inv[, -K, drop=FALSE] #l x (K-1)
    G_left <- GYi_fun_lhl(fit.times[-K])
    #G_left = GYi_int[-K]
    vals_inner <- sweep(dS_inv, 2, 1/G_left, FUN="*")
    vals_inner[,fit.times[-1] > Y[i]] <- 0
    vals <- cbind(0, t(apply(vals_inner, 1, cumsum)))
    
    # SOMEHOW MY CODE IS WRONG
    #S_inv_padded <- cbind(1, S_inv) # pad by 1 lhs, K+1
    #dS_inv <- S_inv_padded[,-1,drop=FALSE] - S_inv_padded[,-(K+1),drop=FALSE] #K
    #G_right <- GYi_fun_rhl(fit.times) #length K avoid specifying a right hand limit function
    #vals_inner <- sweep(dS_inv, 2, 1 / GYi_int, FUN = "*")
    #vals_inner[, fit.times > Y[i]] <- 0
    #Take cumulative sum to express integral, NO NEED to add zero now to first increment
    # if(K==1){
    # vals <- vals_inner
    # }else{
    # vals <- t(apply(vals_inner, 1, cumsum))  # l x (K)
    # }
    
    # PRIOR MISTAKE WAS INDEXING WITH "i". Use "j" instead
    # matches Shi[,which(fit.times==Y[i])]
    S_y <- sapply(
      1:nrow(Shi),
      FUN=function(j){
        pmax(stepfun(fit.times, c(1, Shi[j,]), right=FALSE)(Y[i]), 1E-8)
      }
    )
    # DOES MATCH GYi_int[which(fit.times==Y[i])]
    G_y <- pmax(GYi_fun_lhl(Y[i]), 1E-8)
    
    # K_N on eval.times grid
    # NOTE: we use right=FALSE because we want to ensure that the survival function of censoring is RIGHT CONTINUOUS NOW
    # P(C > t) rather than P(C \geq t) to match Tsiatis' text, which is K(r) = \prod_{u \leq r} (1-\lambda_C(r)) = P(C > r)
    if (single_cens_mechanism) {
      # Original: K_N derived from G_Y under common censoring
      K_N_right <- sapply(eval.times,
                          FUN = function(t) {
                            DeltaY[i] * GYi_fun_lhl(pmax(t, Y[i])) / GYi_fun_lhl(Y[i]) +
                              (1 - DeltaY[i]) * as.numeric(t <= Y[i])
                          })
      # DELTAY=0 convention correct, left-continuous
      # DELTAY=1 conventionc correct, left-continuous
      
      KN_fun_lhl <- stepfun(eval.times, c(1, K_N_right), right=TRUE) #left-continuous step function
      KNi <- KN_fun_lhl(N[i])
    }else{
      # New: K_N evaluated from estimated G_N at participant i's (A_i, X_i, Y_i, DeltaY_i)
      KN_fun_lhl <- stepfun(eval.times, c(1,GN.hats[i, ]), right=TRUE) # left-continuous step function
      KN_fun_rhl <- stepfun(eval.times, c(1,GN.hats[i, ]), right=FALSE) # left-continuous step function
      K_N_right    <- GN.hats[i, ]                    # = K_N(t_j),   length l, the post-jump value
      KNi <- KN_fun_lhl(N[i])
    }
    
    # Matching Westling args
    # int_0^{Ntilde} Q_0(u | ...) * (1/G_0(u) - 1/G_0(u-)) du
    #KN_vec      <- c(1, K_N_right)             # = K_N(t_j^-), length l
    #dLambdaN <- 1 - KN_vec[-1]/KN_vec[-length(KN_vec)]
    #dLambdaN[is.na(dLambdaN) | dLambdaN < 0 | is.nan(dLambdaN)] <- 0
    #KN_inv  <- 1 / KN_vec
    #KN_inv[is.infinite(KN_inv)] <- 0
    #KN_inv_extended <- c(1, KN_inv)
    # SLIGHT DIFF FROM WESTLING -- includes first compensator jump rather than just assuming first jump is 0
    #dKN_inv <- diff(KN_inv_extended) # EQUIVALENT TO KN_inv_extended[-1] -  KN_inv_extended[-l] OR dKN_inv <- c(KN_inv[1] - 1, KN_inv[-1] - KN_inv[-l]) # (l)
    # If matching Ted's code
    #dKN_inv <- KN_inv[-1] - KN_inv[-l] #l-1
    
    denom_cs   <- rev(cumsum(rev(fN.i)))
    IPW_factor <- as.numeric(A[i] == 1L) / g.hats[i]
    
    for (k in seq_len(K)) {
      S_t0 <- Shi[, k] # l x 1
      
      inner1 <- if (Y[i] <= fit.times[k] && DeltaY[i] == 1L) {
        1 / (S_y * G_y)
      } else {
        rep(0, l)
      }
      
      IF_k <- S_t0 - IPW_factor * S_t0 * (inner1 - vals[, k]) # l x 1
      
      num_cs <- rev(cumsum(rev(IF_k[-1] * fN.i))) # l x 1
      Q_seq  <- num_cs / pmax(denom_cs, 1e-8)
      
      term1 <- (    DeltaN[i]  / KNi) * IF_k[N_indx]
      
      int_N <- Q_seq * (1/KN_fun_lhl(eval.times[-1])) #length l
      int_N[is.infinite(int_N) | is.nan(int_N) | is.na(int_N)] <- 0
      # P(C_N >= 2) - P(C_N >= 1)
      dCP <- as.numeric(N[i] == eval.times[-1] & DeltaN[i]==0)
      k_left <- KN_fun_lhl(eval.times)
      #k_right <- KN_fun_rhl(eval.times[-1])
      #dComp <- c(0, 1- (kl_left[-1]/kl_left[-length(kl)]))
      dComp <- 1 - k_left[-1]/k_left[-l]
      dComp[eval.times[-1] > N[i]] <- 0
      dMart <- dCP - dComp
      
      aug_term <- sum(int_N * dMart)
      
      IF.vals[i, k] <- term1 + aug_term
      MARTS[i, k] <- mean(dMart)
    }
    
    if (verbose) pb$tick()
  }
  
  surv <- colMeans(IF.vals, na.rm = TRUE)
  IF.vals.center <- sweep(IF.vals, 2, surv, "-")
  
  res <- list(
    times = fit.times,
    surv = pmin(1, pmax(0, surv)),
    IF.vals = IF.vals,
    IF.vals.center = IF.vals.center
  )
  if (isotonize) {
    idx <- !is.na(surv)
    res$surv.iso <- surv
    res$surv.iso[idx] <- 1 - isoreg(fit.times[idx], 1 - surv[idx])$yf
  }
  res
}