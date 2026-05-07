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
    logN             = TRUE
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
  # nuisance.options = CFsurvival.nuisance.options(
  #   cross.fit = TRUE,
  #   V = 5,
  #   folds = NULL,
  #   eval.times = sort(unique(c(data$N1, 730))),
  #   event.SL.library = c("survSL.rfsrc"),
  #   event.pred.0 = NULL,
  #   event.pred.1 = NULL,
  #   cens.SL.library = c("survSL.coxph", "survSL.km", "survSL.expreg"),
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
  # single_cens_mechanism = TRUE;
  # em = TRUE;
  # logN = FALSE
  
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
    IF0      = matrix(NA_real_, 0L, k),   # EIF rows for arm 0, per participant
    IF1      = matrix(NA_real_, 0L, k),   # EIF rows for arm 1, per participant
    #censY0   = matrix(NA_real_, 0L, k),
    #censY1   = matrix(NA_real_, 0L, k),
    #censN0   = matrix(NA_real_, 0L, l),
    #censN1   = matrix(NA_real_, 0L, l),
    #fN0      = matrix(NA_real_, 0L, l),
    #fN1      = matrix(NA_real_, 0L, l),
    test_ids = integer(0L)
  )
  
  #fit marginal model for $f_{N \mid A, X}$ OUTSIDE this function
  # stopifnot(!is.null(FN_fit))
  
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
    # Step 0: generate test datasets
    # ----------------------------------------------------------------
    
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
      X                = cbind(A = dat_tr$A, W_tr),
      newX             = as.data.frame(newAW_te),
      new.times        = eval.times,
      event.SL.library = nuis$event.SL.library[nuis$event.SL.library!="survSL.km"], #ENFORCING SMOOTHNESSS
      cens.SL.library  = nuis$cens.SL.library[nuis$cens.SL.library!="survSL.km"], #ENFORCING SMOOTHNESSS
      verbose          = FALSE,
      control          = nuis$survSL.control,
      cvControl        = nuis$survSL.cvControl
    )
    SN_pred <- SN_fit$event.SL.predict          # (2*n_te x k) where k is fit.times
    SN0 <- SN_pred[seq_len(n_te),         , drop=FALSE] # (n_te x k) where k is fit.times
    SN1 <- SN_pred[seq_len(n_te) + n_te,  , drop=FALSE] # (n_te x k) where k is fit.times
    SN0[SN0 == 0] <- min(SN0[SN0 > 0])
    SN1[SN1 == 0] <- min(SN1[SN1 > 0])
    
    pSSL <- function (object, newdata, new.times, X = NULL, time = NULL, 
              event = NULL, onlySL = FALSE, threshold = 1e-04, ...){
      
      if (missing(newdata)) {
        out <- list(survpred = object$survSL.predict, cens.pred = object$cens.SL.predict, 
                    survlibrary.predict = object$survlibrary.predict, 
                    cens.library.predict = object$cens.library.predict)
        return(out)
      }
      if (!object$control$saveFitLibrary) {
        stop("This SuperLearner fit was created using control$saveFitLibrary = FALSE, so new predictions cannot be made.")
      }
      event.k <- nrow(object$event.libraryNames) # 1
      event.pred <- array(0, dim = c(nrow(newdata), length(new.times), 
                                     event.k))
      dimnames(event.pred)[[3]] <- apply(object$event.libraryNames, 
                                         1, paste, collapse = "_")
      dimnames(event.pred)[[2]] <- new.times
      cens.k <- nrow(object$cens.libraryNames)
      cens.pred <- array(0, dim = c(nrow(newdata), length(new.times), 
                                    cens.k))
      dimnames(cens.pred)[[3]] <- apply(object$cens.libraryNames, 
                                        1, paste, collapse = "_")
      dimnames(cens.pred)[[2]] <- new.times
      if (onlySL) {
        event.whichLibrary <- which(object$event.coef > threshold)
        event.coef <- object$event.coef
        event.coef[-event.whichLibrary] <- 0
        event.coef <- event.coef/sum(event.coef)
        cens.whichLibrary <- which(object$cens.coef > threshold)
        cens.coef <- object$cens.coef
        cens.coef[-cens.whichLibrary] <- 0
        cens.coef <- cens.coef/sum(cens.coef)
      }
      else {
        event.whichLibrary <- seq(event.k)
        event.coef <- object$event.coef
        cens.whichLibrary <- seq(cens.k)
        cens.coef <- object$cens.coef
      }
      for (mm in event.whichLibrary) {
        newdataMM <- subset(newdata, select = object$event.whichScreen[object$event.SL.library$library[mm, 
                                                                                                       2], ])
        XMM <- if (is.null(X)) {
          NULL
        }
        else {
          subset(X, select = object$event.whichScreen[object$event.SL.library$library[mm, 
                                                                                      2], ])
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
      for (mm in cens.whichLibrary) {
        newdataMM <- subset(newdata, select = object$cens.whichScreen[object$cens.SL.library$library[mm, 
                                                                                                     2], ])
        XMM <- if (is.null(X)) {
          NULL
        }
        else {
          subset(X, select = object$cens.whichScreen[object$cens.SL.library$library[mm, 
                                                                                    2], ])
        }
        cens.pred[, , mm] <- do.call(paste0("predict", ".", class(object$cens.fitLibrary[[mm]])),
                                     list(object = object$cens.fitLibrary[[mm]], 
                                                     newX = newdataMM, new.times = new.times, X = XMM, 
                                                     time = time, event = event, ...))
      }
      cens.SL.predict <- matrix(NA, nrow = nrow(newdata), ncol = length(new.times))
      for (j in seq(length(new.times))) {
        cens.SL.predict[, j] <- cens.pred[, j, ] %*% cens.coef
      }
      out <- list(event.SL.predict = event.SL.predict, event.library.predict = event.pred, 
                  cens.SL.predict = cens.SL.predict, cens.library.predict = cens.pred)
      return(out)
    }
    
    #survest function
    survest_N <- function(fit=SN_fit, t, x, a) {
      # x <- newAW_te[1,2]
      # t=c(100, 200)
      # a=1
      x <- as.data.frame(x)
      colnames(x) <- colnames(W)
      if (nrow(x) == 1L) { x <- rbind(x, x); drop1 <- TRUE } else drop1 <- FALSE
      newdat <- bind_cols(A = a, x)
      out <- pSSL(object = fit, 
                  newdata = newdat, 
                  new.times = t, 
                  onlySL = TRUE)$event.SL.predict
      if (drop1) out <- out[1,, drop=FALSE]
      as.vector(out)
    }
    
    # ----------------------------------------------------------------
    # Step 4: Augment training fold (redistribute-to-right weights)
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
    if (verbose) message("    [3] S_Y | A, X, N")
    
    if(em==FALSE){
      
      if(verbose) message("        No EM-style updates (S_{0Y} may be misspecified)")
      
      AW_ev_tr  <- cbind(A = aug_tr$A, aug_tr[, c("N","DeltaN", W_names), drop=FALSE])
      AW_ce_tr  <- cbind(A = aug_tr$A, aug_tr[, W_names, drop=FALSE])
      AW_ev_te  <- cbind(A = aug_te$A, aug_te[, c("N","DeltaN", W_names), drop=FALSE])
      # newX_cens: one row per unique test participant (distinct W only)
      AW_ce_te_0 <- cbind(A=0L, W_te)
      AW_ce_te_1 <- cbind(A=1L, W_te)
      AW_ce_te   <- rbind(AW_ce_te_0, AW_ce_te_1)
      
      SY_fit <- survSuperLearner2(
        time             = aug_tr$Y,
        event            = aug_tr$DeltaY,
        X_event          = AW_ev_tr,
        X_cens           = AW_ce_tr,
        newX_event       = cbind(A=0L, aug_te[,c("N","DeltaN",W_names),drop=FALSE] %>% slice(1:2)),
        newX_cens        = AW_ce_te[1:2,],
        new.times        = fit.times,
        event.SL.library = nuis$event.SL.library[nuis$event.SL.library!="survSL.km"], #ENFORCING SMOOTHNESSS
        cens.SL.library  = nuis$cens.SL.library[nuis$cens.SL.library!="survSL.km"], #ENFORCING SMOOTHNESSS
        verbose          = FALSE,
        control          = nuis$survSL.control,
        cvControl        = nuis$survSL.cvControl,
        obsWeights       = aug_tr$w,
        id               = aug_tr$id
      )
      
        SY0 <- predict.survSuperLearner.2(SY_fit, 
                                          newdata = cbind(A=0L, aug_te[,c("N","DeltaN",W_names),drop=FALSE]), 
                                          new.times = fit.times, onlySL = TRUE)$event.SL.predict %>% 
          bigmemory::as.big.matrix()
        
        gc()
        
        SY1 <- predict.survSuperLearner.2(SY_fit, 
                                          newdata = cbind(A=1L, aug_te[,c("N","DeltaN",W_names),drop=FALSE]), 
                                          new.times = fit.times, onlySL = TRUE)$event.SL.predict %>% 
          bigmemory::as.big.matrix()                  
        
        rm(SY_fit)
        gc()
      
    }else{
      
      if(verbose) message("        Iterative EM-style updates to S_{0Y}")  
      
      if(logN==FALSE){
        em_out <- .fit_SY_em(
        aug_tr            = aug_tr,
        aug_te            = aug_te,
        dat_te            = dat_te,
        fit.times         = fit.times,
        W_names           = W_names,
        nuis              = nuis,
        max_iter          = 10,             # 1 disables EM; 2 = one update
        tol               = 1/length(fit.times),
        fast_iter_library = NULL,           # e.g. c("survSL.coxph") for speed
        verbose           = verbose)
      }else{
        em_out <- .fit_SY_em(
          aug_tr            = aug_tr %>% mutate(N=log(N)),
          aug_te            = aug_te %>% mutate(N=log(N)),
          dat_te            = dat_te,
          fit.times         = fit.times,
          W_names           = W_names,
          nuis              = nuis,
          max_iter          = 10,             # 1 disables EM; 2 = one update
          tol               = 1/length(fit.times),
          fast_iter_library = NULL,           # e.g. c("survSL.coxph") for speed
          verbose           = verbose)
      }
      
      SY0    <- bigmemory::as.big.matrix(em_out$SY0)
      SY1    <- bigmemory::as.big.matrix(em_out$SY1)
      aug_tr <- em_out$aug_tr_final         # downstream code uses aug_tr$id only
      gc()
      
    }
    # ----------------------------------------------------------------
    # Step 2: G_C^Y — random-dropout censoring of Y process
    # random_cens_Y = 1 iff Y is censored AND time_Y < admin_t
    # ----------------------------------------------------------------
    
    
    if(single_cens_mechanism){
      
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
      
      GC_fit <- survSuperLearner::survSuperLearner(
        time             = C,
        event            = rc,
        X                = cbind(A = dat_tr$A, W_tr),
        newX             = as.data.frame(newAW_te),
        new.times        = sort(unique(c(eval.times, fit.times))),
        event.SL.library = nuis$cens.SL.library,
        cens.SL.library  = nuis$cens.SL.library,
        verbose          = FALSE,
        control          = nuis$survSL.control,
        cvControl        = nuis$survSL.cvControl
      )
      GC_pred <- GC_fit$event.SL.predict          # (2*n_te x k) where k is ALL times
      GC0 <- pmax(GC_pred[seq_len(n_te),         , drop=FALSE], nuis$cens.trunc) # (n_te x k) where k is ALL times
      GC1 <- pmax(GC_pred[seq_len(n_te) + n_te,  , drop=FALSE], nuis$cens.trunc) # (n_te x k) where k is ALL times
      GC0[GC0 == 0] <- min(GC0[GC0 > 0])
      GC1[GC1 == 0] <- min(GC1[GC1 > 0])
      
      
      
    }else{
      
      ### TO FILL OUT LATER
      
      
    }
    
    # ----------------------------------------------------------------
    # Step 7: Posterior f_N via Bayes rule on test fold
    # ----------------------------------------------------------------
    if (verbose) message("    [6] Posterior f_{0N}(n | a, x, delta_y, tilde{y}) via bayes inversion")
    
    post <- .compute_posterior_fN(
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
      eif0_fold <- .get.survival.AIPWCC(
        Y       = dat_te$Y,
        DeltaY  = dat_te$DeltaY,
        N       = dat_te$N,
        DeltaN  = dat_te$DeltaN_eff,
        A       = 1L - dat_te$A,
        admin_t = dat_te$admin_t,
        fit.times  = fit.times,
        eval.times = eval.times,
        ids        = aug_te$id,
        S.hats     = SY0,                    # (n_aug_te x K) — fold-local only
        GY.hats_a0  = GC0,
        GY.hats_alt = GC1,
        GN.hats_a0  = GC0,
        GN.hats_alt = GC1,
        g.hats  = 1 - nuis$prop.pred[te_idx],
        fN_a0   = post$fN0,
        fN_alt  = post$fN1,
        isotonize = FALSE,
        verbose   = TRUE,
        single_cens_mechanism = TRUE
      )
      acc$IF0 <- rbind(acc$IF0, eif0_fold$IF.vals)
    }
    
    if (1L %in% fit.treat) {
      eif1_fold <- .get.survival.AIPWCC(
        Y       = dat_te$Y,
        DeltaY  = dat_te$DeltaY,
        N       = dat_te$N,
        DeltaN  = dat_te$DeltaN_eff,
        A       = dat_te$A,
        admin_t = dat_te$admin_t,
        fit.times  = fit.times,
        eval.times = eval.times,
        ids        = aug_te$id,
        S.hats     = SY1,
        GY.hats_a0  = GC1,
        GY.hats_alt = GC0,
        GN.hats_a0  = GC1,
        GN.hats_alt = GC0,
        g.hats  = nuis$prop.pred[te_idx],
        fN_a0   = post$fN1,
        fN_alt  = post$fN0,
        isotonize = FALSE,
        verbose   = TRUE, 
        single_cens_mechanism = TRUE
      )
      acc$IF1 <- rbind(acc$IF1, eif1_fold$IF.vals)
    }
    
    #acc$event0   <- rbind(acc$event0,   SY0)
    #acc$event1   <- rbind(acc$event1,   SY1)
    #acc$censY0   <- rbind(acc$censY0,   GY0)
    #acc$censY1   <- rbind(acc$censY1,   GY1)
    #acc$censN0   <- rbind(acc$censN0,   GN0)
    #acc$censN1   <- rbind(acc$censN1,   GN1)
    #acc$fN0      <- rbind(acc$fN0,      post$fN0)
    #acc$fN1      <- rbind(acc$fN1,      post$fN1)
    #acc$aug_ids  <- c(acc$aug_ids,      aug_te$id)
    acc$test_ids <- c(acc$test_ids,     dat_te$id)
    
    rm(SY0, SY1, aug_tr, aug_te, post,
       eif0_fold, eif1_fold, GC_fit, GC0, GC1)
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
      idx <- !is.na(surv)
      res$surv.iso     <- surv
      res$surv.iso[idx] <- 1 - isoreg(fit.times[idx], 1 - surv[idx])$yf
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
                  aug_tr[, c("N", "DeltaN", W_names), drop = FALSE])
  
  SY_pred <- predict.survSuperLearner.2(
    object    = SY_fit,
    newdata   = newdat,
    new.times = fit.times,
    onlySL    = TRUE
  )$event.SL.predict
  # SY_pred is (nrow(aug_tr) x K)
  
  K <- length(fit.times)
  Y_idx <- findInterval(aug_tr$Y, c(0, fit.times), rightmost.closed = FALSE) #FIXED
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
    logN              = FALSE
) {
  
  stopifnot(max_iter >= 1L)
  
  AW_ev_tr <- cbind(A = aug_tr$A, aug_tr[, c("N", "DeltaN", W_names), drop = FALSE])
  if(logN==TRUE){
  AW_ev_tr$N <- log(AW_ev_tr$N)
  }
  AW_ce_tr <- cbind(A = aug_tr$A, aug_tr[, W_names, drop = FALSE])
  
  W_te     <- dat_te[, W_names, drop = FALSE]
  AW_ce_te <- rbind(cbind(A = 0L, W_te), cbind(A = 1L, W_te))
  
  full_event_lib <- nuis$event.SL.library[nuis$event.SL.library != "survSL.km"]
  full_cens_lib  <- nuis$cens.SL.library [nuis$cens.SL.library  != "survSL.km"]
  
  has_fast <- !is.null(fast_iter_library) && length(fast_iter_library) >= 1L
  
  # Helper: one survSuperLearner2 fit using current aug_tr
  fit_once <- function(event_lib, cens_lib) {
    survSuperLearner2(
      time             = aug_tr$Y,
      event            = aug_tr$DeltaY,
      X_event          = AW_ev_tr,
      X_cens           = AW_ce_tr,
      newX_event       = cbind(A = 0L,
                               aug_te[, c("N","DeltaN",W_names), drop = FALSE] %>%
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
    max_w_diff <- max(abs(aug_tr_new$w - w_in), na.rm = TRUE)
    
    history[[iter]] <- list(
      iter       = iter,
      max_w_diff = max_w_diff,
      n_aug      = nrow(aug_tr),
      library    = if (use_fast) "fast" else "full"
    )
    
    if (verbose) message(sprintf(
      "      max |w^(t+1) - w^(t)| = %.4g  (tol = %.4g)", max_w_diff, tol
    ))
    
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
  SY0 <- predict.survSuperLearner.2(
    object    = SY_fit,
    newdata   = cbind(A = 0L, aug_te[, c("N","DeltaN",W_names), drop = FALSE]),
    new.times = fit.times,
    onlySL    = TRUE
  )$event.SL.predict
  
  SY1 <- predict.survSuperLearner.2(
    object    = SY_fit,
    newdata   = cbind(A = 1L, aug_te[, c("N","DeltaN",W_names), drop = FALSE]),
    new.times = fit.times,
    onlySL    = TRUE
  )$event.SL.predict
  
  list(
    SY_fit       = SY_fit,
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
    idx <- findInterval(t_val, grid, rightmost.closed=FALSE) #FIXED
    pmax(1L, pmin(idx, length(grid)))
  }
  
  for (i in seq_len(n_te)) {
    Yi      <- dat_te$Y[i]
    dYi     <- dat_te$DeltaY[i]
    idi     <- dat_te$id[i]
    xi      <- dat_te[i, W_names, drop=FALSE]
    aug_idx <- id_to_aug[[as.character(idi)]]
    
    # Prior: discrete PMF of N on eval.times under each arm
    SN1 <- survest_N(t=eval.times, x=xi, a=1)
    SN0 <- survest_N(t=eval.times, x=xi, a=0)
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
      if (all(lk1 == 0)) lk1 <- pr1
      if (all(lk0 == 0)) lk0 <- pr0
    }
    
    # Posterior ∝ likelihood × prior
    
    #ENFORCE 0/0 = 1 and add 1E-12 to any zero denominator
    '|' <- function(a,b) ifelse(a==0 & b==0, 1, ifelse(a!=0 & b==0, a/(b+1e-12), a/b))
    
    u1 <- lk1 * pr1;  
    d1 <- sum(u1) 
    fN1[i,] <- u1 | d1
    u0 <- lk0 * pr0;  
    d0 <- sum(u0);  
    fN0[i,] <- u0 | d0
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
    verbose   = TRUE,
    single_cens_mechanism = TRUE
) {
  
  # Y       = dat_te$Y;
  # DeltaY  = dat_te$DeltaY;
  # N       = dat_te$N;
  # DeltaN  = dat_te$DeltaN_eff;
  # A       = dat_te$A;
  # admin_t = dat_te$admin_t;
  # fit.times  = fit.times;
  # eval.times = eval.times;
  # ids        = aug_te$id;
  # S.hats     = SY1;
  # GY.hats_a0  = GC1;
  # GY.hats_alt = GC0;
  # GN.hats_a0  = GC1;
  # GN.hats_alt = GC0;
  # g.hats  = nuis$prop.pred[te_idx];
  # fN_a0   = post$fN1;
  # fN_alt  = post$fN0;
  # isotonize = FALSE;
  # verbose   =TRUE;
  # single_cens_mechanism = TRUE
  
  fit.times <- sort(fit.times[fit.times > 0])
  n  <- length(Y)
  K  <- length(fit.times)
  l  <- length(eval.times)
  
  # Sort K-column matrices
  ord        <- order(fit.times)
  fit.times  <- fit.times[ord]
  S.hats     <- S.hats[,       ord, drop=FALSE]
  # GY.hats_a0 <- GY.hats_a0[,  ord, drop=FALSE]
  # GY.hats_alt<- GY.hats_alt[, ord, drop=FALSE]
  # GN.hats_a0 <- GN.hats_a0[,  ord, drop=FALSE]
  # GN.hats_alt<- GN.hats_alt[, ord, drop=FALSE]
  
  # Grid index lookups (pre-computed, avoid repeated which() inside loop)
  Y_ind <- findInterval(Y, fit.times, rightmost.closed=FALSE) # FIXED
  Y_ind <- pmax(1L, pmin(Y_ind, K))
  N_ind <- findInterval(N, eval.times, rightmost.closed=FALSE) # FIXED
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
    fN.i   <- if (A[i] == 1L) fN_a0[i,] else fN_alt[i,] # l x 1
    
    # S_Y rows for this participant
    rows_i   <- id_to_rows[[as.character(unique_ids[i])]]
    Shi      <- S.hats[rows_i,, drop=FALSE]   # (l x K)
    
    # G_C^Y step function (for phi^F: vals.i and G_y)
    GYi_vec  <- GY.hats_a0[i,]
    if(single_cens_mechanism==TRUE){
      GYi_fun  <- stepfun(sort(unique(c(eval.times, fit.times))), c(1, GYi_vec), right=TRUE)
    }
    
    # ---- phi^F compensator (Ted's pattern, adapted for (l x K) Shi) ----
    S_inv  <- 1 / Shi # (l x K)
    dS_inv <- S_inv[,-1,drop=FALSE] - S_inv[,-K,drop=FALSE]   # (l x K-1)
    # Divisor: G_Y at left limit of each interval (t_{m-1}, t_m]
    G_left = GYi_fun(fit.times[-K]) #length (K-1)
    vals_inner   <- sweep(dS_inv, 2, 1/G_left, FUN="*") # (l x (K-1))
    # Zero out increments at event times t_m > Y_i
    vals_inner[, fit.times[-1] > Y[i]] <- 0 # length K-1 mask, matches columns
    # Prepend 0 and take cumsum to get length-K cumulative compensator
    vals <- cbind(0, t(apply(vals_inner, 1, cumsum)))              # (l x K)
    
    S_y  <- Shi[, Y_ind[i]]       # S_Y(Y_i | a0, X, N_j) for each j
    G_y  <- GYi_fun(Y[i])         # G_C^Y(Y_i | a0, X)
    
    # ---- K_N on eval.times grid (no 0 prepended) ---
    if (single_cens_mechanism) {
      K_N <- sapply(eval.times,
                    FUN = function(t) {
                      DeltaY[i] * GYi_fun(pmax(t, Y[i])) / GYi_fun(Y[i]) +
                        (1 - DeltaY[i]) * as.numeric(t <= Y[i])              # CHANGED to <=
                    })
    }
    
    KN_fun <- stepfun(eval.times, c(1, K_N), right=TRUE)
    
    # ---- Hazard increments and at-risk weights ----
    KN_inv_left   <- 1 / KN_fun(eval.times[-l])                               # length (l-1)
    KN_inv_left = ifelse(is.infinite(KN_inv_left), 0, KN_inv_left)
    dLambdaN <- 1 - KN_fun(eval.times[-1]) / KN_fun(eval.times[-l])           # length (l-1)
    dLambdaN <- ifelse(is.nan(dLambdaN), 0, dLambdaN)
    

    # fN tail-sum denominator (k-independent, hoist outside k loop)
    denom_cs <- rev(cumsum(rev(fN.i)))   # P(N >= eval.times[j] | ...)
    
    IPW_factor <- as.numeric(A[i] == 1L) / g.hats[i]
    
    for (k in seq_len(K)) {
      
      S_t0 <- Shi[, k]   # S_Y(t0 | a0, X, N_j), j=1..l
      
      inner1 <- if (Y[i] <= fit.times[k] && DeltaY[i] == 1L){
        1 / (S_y * G_y) 
      }else{ 
        rep(0, l)
      }
      
      # phi^F as a function of N grid point j
      IF_k <- S_t0 - IPW_factor * S_t0 * (inner1 - vals[, k])
      
      # Q_seq[j] = E[phi^F(N) | N >= eval.times[j], obs, A=A_i]
      num_cs <- rev(cumsum(rev(IF_k * fN.i)))
      Q_seq  <- num_cs / pmax(denom_cs, 1e-12) #length l
      
      # AIPWCC terms
      # DeltaN[i] = DeltaN_eff: admin-censored N has DeltaN_eff=1 → term1
      # randomly censored N has DeltaN_eff=0 → term2
      term1 <- (     DeltaN[i]  / K_N[N_indx]) * IF_k[N_indx]
      term2 <- ((1 - DeltaN[i]) / K_N[N_indx]) * Q_seq[N_indx]
      
      # Term 3: martingale convention, integrate over intervals 1..N_indx-1
      # (the jump at tilde_N_i is captured by Term 2, not Term 3)
      if (N_indx >= 2L) {
        idx3  <- seq_len(N_indx - 1L)
        term3 <- sum(Q_seq[idx3] * KN_inv_left[idx3] * dLambdaN[idx3])
      } else {
        term3 <- 0
      }
      
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
