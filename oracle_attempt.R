# =============================================================================
# AIPWCC_newEIF_em_oracle.r 
# More memory efficient -- evaluates EIF at each inner CV loop to avoid memory exhaustion limits
# ALSO implements new EIF calculation
# USES ORACLE NUISANCE PARAMS
# =============================================================================

AIPWCC_newEIF_oracle <- function(
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
    #nuisance.options = list(),
    conf.band        = FALSE,
    conf.level       = 0.95,
    contrasts        = c("risk.ratio"),
    verbose          = TRUE,
    dgm_params       = list(n, trt_prop, RR, RR_C, per_year_incidence_Y, per_year_incidence_N, per_year_incidence_C, rho, admin_date, rho_X)
) {
  # 
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
  # dgm_params<-list(
  #   n=100000,
  #   trt_prop = 0.5,
  #   RR = 0.75,
  #   RR_C = 1,
  #   per_year_incidence_Y = 0.06,
  #   per_year_incidence_N = 0.12,
  #   per_year_incidence_C = 0.10,
  #   rho = 0,
  #   admin_date = 2,
  #   rho_X=0
  # )
  
  #nuis <- do.call("CFsurvival.nuisance.options", nuisance.options)
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
  N_for_lookup <- ifelse(event_N == 0L & time_N == admin_t, 731L, time_N)
  
  W_names     <- colnames(W)
  # Define grids once
  admin_t    <- dgm_params$admin_date * 365
  #fit.times  <- seq(1L, admin_t, by=1L)   # times where we evaluate S_Y
  eval.times <- seq(1L, admin_t+1, by=1L)   # values of N where we stratify
  K <- k <- length(fit.times)
  L <- l <- length(eval.times)
  
  contrasts <- unique(tolower(contrasts))
  contrasts <- intersect(contrasts, c("surv.diff","surv.ratio","risk.ratio","nnt"))
  
  result <- list(fit.times = fit.times, fit.treat = fit.treat,
                 surv.df = data.frame())
  
  make_df <- function(idx) as.data.frame(cbind(
    id=id, Y=time, DeltaY=event,
    N=time_N, DeltaN=event_N, DeltaN_eff=DeltaN_eff,
    A=treat, admin_t=admin_t, W
  ))[idx, ]
  
  data <- make_df(1:length(time))
  
  # ==================================================================
  # CROSS-FITTED NUISANCE ESTIMATION
  # ==================================================================
  
  big_data_event <- create_multivar_survival(n=dgm_params$n, 
                                   trt_prob=dgm_params$trt_prop, 
                                   RR = dgm_params$RR, 
                                   RR_C = dgm_params$RR_C, 
                                   per_year_incidence_Y = dgm_params$per_year_incidence_Y, 
                                   per_year_incidence_N = dgm_params$per_year_incidence_N, 
                                   per_year_incidence_C = 1E-16,
                                   J = 1, 
                                   rho = dgm_params$rho, 
                                   admin_date = dgm_params$admin_date,
                                   X=TRUE,
                                   rho_X = dgm_params$rho_X, 
                                   convert_days = TRUE)
  
  ##############
  # Fit Nuisances to Huge Dataset (Oracle)
  
  big_data_event$N1 <- ifelse(big_data_event$DeltaN1 == 0L & big_data_event$N1 == admin_t, 731L, big_data_event$N1)
  big_data_event$DeltaN1 <- ifelse(big_data_event$N1==731, 1, big_data_event$DeltaN1)
  
  fit_S_N <- survfit(Surv(N1, DeltaN1) ~ 1, data=big_data_event)
  
  # Evaluate S_N at the grid
  S_N_grid <- summary(fit_S_N, times=eval.times, extend=TRUE)$surv  # length L
  
  f_N_grid <- pmax(0, c(1, S_N_grid[-l]) - S_N_grid)  # length L
  f_N_grid <- f_N_grid / sum(f_N_grid) 
  
  # ---- Step 4: oracle S_Y(t | a, n) -- this is the core nuisance ----
  # Stratify by (A, N_bin) where N_bin is a small window around eval.times[s].
  # Store as a 3D array: [arm, n_index, fit_time_index]
  # This avoids the closure bug entirely.
  
  bin_half_width <- 0L
  
  S_Y_array <- array(NA_real_, dim=c(2, l, k),
                     dimnames=list(arm=c("0","1"), n_idx=NULL, t_idx=NULL))
  
  if(verbose){
  pb <- progress_bar$new(length(eval.times)*2, format="  S_{0Y}(y | a,x,n) [:bar] :percent ETA: :eta")
  }
  for (a in c(0L, 1L)) {
    for (s in seq_along(eval.times)) {
      n_target <- eval.times[s]
      
      if (n_target <= 730L) {
        # Conditional on a real N event at this time
        sub <- big_data_event %>%
          filter(A == a,
                 abs(N1 - n_target) <= bin_half_width,
                 DeltaN1 == 1L)
        if (nrow(sub) < 50L) {
          sub <- big_data_event %>%
            filter(A == a, abs(N1 - n_target) <= 1L, DeltaN1 == 1L)
        }
      } else {
        # n_target = 731: post-admin sentinel.
        # Conditional on N > admin_t = 730. In the aux data with no random
        # censoring of N, this is exactly the admin-censored cohort.
        sub <- big_data_event %>%
          filter(A == a, N1 == 731 & DeltaN1==1L)
      }
      
      if (nrow(sub) < 1L) {
        # Pathologically sparse: store NA, will need to handle downstream
        next
      }
      
      fit_SY <- survfit(Surv(Y, DeltaY) ~ 1, data=sub)
      S_Y_array[a + 1L, s, ] <- summary(fit_SY, times=fit.times, extend=TRUE)$surv
      if(verbose){pb$tick()}
    }
    #message(sprintf("Completed arm %d", a))
  }
  
  # ---- Step 5: lookup function for S_Y(t | a, n) ----
  # Vectorized: takes vectors t and n (and scalar or vector a), returns matrix of survivals.
  S_Y_oracle <- function(t, n, a) {
    # Find indices in fit.times and eval.times grids
    t_idx <- pmax(1L, pmin(findInterval(t, fit.times), k))
    n_idx <- pmax(1L, pmin(findInterval(n, eval.times), l))
    a_idx <- as.integer(a) + 1L
    # Vectorized lookup
    S_Y_array[cbind(a_idx, n_idx, t_idx)]
  }
  
  # Test it
  #S_Y_oracle(t=1:730, n=c(730), a=c(0))
  #S_Y_oracle(t=1:730, n=c(731), a=c(0))
  
  # ---- Step 6: lookup for f_N(n) ----
  # Since N _|_ (A, X), this is a marginal -- doesn't depend on a or x
  f_N_oracle <- function(n) {
    n_idx <- pmax(1L, pmin(findInterval(n, eval.times), L))
    f_N_grid[n_idx]
  }
  
  
  # ---- Step 7: posterior f_{0N}(n | a, x, y, dy) for a single participant ----
  # Returns a vector of length L (one entry per eval.times point)
  
  f_N_posterior_oracle <- function(y_tilde, delta_y, a, x = NULL) {
    # x is unused under our DGM (rho_X = 0 makes N _|_ X). Kept in signature
    # for interface consistency.
    
    # Get S_Y at all (n_idx, y_tilde) pairs for this arm
    y_idx <- pmax(1L, pmin(findInterval(y_tilde, fit.times), k))
    S_Y_at_y      <- S_Y_array[a + 1L, , y_idx]              # length L
    S_Y_at_y_prev <- if (y_idx > 1L) S_Y_array[a + 1L, , y_idx - 1L] else rep(1, l)
    f_Y_at_y      <- pmax(0, S_Y_at_y_prev - S_Y_at_y)        # length L
    
    # Marginal at this y_tilde: S_Y(y | a) = sum_n S_0Y(y | a, n) f_N(n)
    S_Y_marg <- sum(S_Y_at_y * f_N_grid)
    f_Y_marg <- sum(f_Y_at_y * f_N_grid)
    
    # Likelihood ratio L(y, dy | a, n) / L(y, dy | a)
    if (delta_y == 1L) {
      LR <- if (f_Y_marg > 1e-12) f_Y_at_y / f_Y_marg else rep(0, L)
    } else {
      LR <- if (S_Y_marg > 1e-12) S_Y_at_y / S_Y_marg else rep(0, L)
    }
    
    # Posterior: prior x LR, then renormalize (in case of small numeric error)
    posterior <- f_N_grid * LR
    s <- sum(posterior)
    if (s > 1e-12) posterior <- posterior / s else posterior <- f_N_grid
    
    posterior}
  
  # Test it
  # post <- f_N_posterior_oracle(y_tilde=150, delta_y=1L, a=1L)
  # plot(eval.times, post, type='l',
  #      main="Oracle posterior f_N(n | y_tilde=300, delta_y=1, a=1)")
  # plot(eval.times, f_N_grid, type='l', col='red', add=TRUE)
  
  S_Y_predict_oracle <- function(newdata, fit.times) {
    # newdata is a data frame with columns A, X, N (and possibly DeltaN)
    out <- matrix(NA_real_, nrow=nrow(newdata), ncol=length(fit.times))
    for (i in seq_len(nrow(newdata))) {
      out[i, ] <- S_Y_oracle(t=fit.times,
                             n=rep(newdata$N[i], length(fit.times)),
                             a=rep(newdata$A[i], length(fit.times)))
    }
    out
  }
  
  f_N_post_predict_oracle <- function(newdata) {
    out <- matrix(NA_real_, nrow=nrow(newdata), ncol=length(eval.times))
    for (i in seq_len(nrow(newdata))) {
      out[i, ] <- f_N_posterior_oracle(
        y_tilde = newdata$Y[i],
        delta_y = newdata$DeltaY[i],
        a       = newdata$A[i]
      )
    }
    out
  }
  
  lambda_C <- -log(1 - dgm_params$per_year_incidence_C) / 365   # daily hazard for 10% annual incidence
  G_0Y_oracle <- function(t, a, x = NULL) {
    exp(-lambda_C * t)   # arm-independent in your DGM
  }

  SY_block_a0 <- S_Y_array[1,,]
  SY_block_a1 <- S_Y_array[2,,]
   
  fN_1 <- f_N_post_predict_oracle(data %>% mutate(A=1)) # n x l
  
  fN_0 <- f_N_post_predict_oracle(data %>% mutate(A=0)) # n xl
  
  preds <- predict(glm(A ~ X, data, family="binomial"), type="response")
  
  GY_1 <- G_0Y_oracle(fit.times, a=1) # k x 1
  
  GY_0 <- G_0Y_oracle(fit.times, a=0) # k x 1
  
  if(1 %in% fit.treat){
    IF1 <- .get.survival.AIPWCC.oracle(
      Y = data$Y, 
      DeltaY = data$DeltaY, 
      N = N_for_lookup, 
      DeltaN = data$DeltaN,
      A = data$A,                       # A[i]=1 iff participant in target arm a0
      admin_t = rep((data$admin_t+1), nrow(data)),
      fit.times = fit.times, 
      eval.times = eval.times, 
      SY_block_target = SY_block_a1,                  
      SY_block_alt = SY_block_a0,                 
      GY_vec = GY_1,
      g.hats = preds,                  # pi(a0|X)                     (n)
      fN_a0 = fN_1,                   # posterior P(N|Y,dY,X,A=a0)   (n x l)
      fN_alt = fN_0,                  # posterior P(N|Y,dY,X,A=1-a0) (n x l)
      isotonize = TRUE,
      verbose   = TRUE,
      single_cens_mechanism = TRUE
    )  
  }
  
  if(0 %in% fit.treat){
    IF0 <- .get.survival.AIPWCC.oracle(
      Y = data$Y, 
      DeltaY = data$DeltaY, 
      N = N_for_lookup, 
      DeltaN = data$DeltaN,
      A = 1-data$A,                       # A[i]=1 iff participant in target arm a0
      admin_t = rep((data$admin_t+1), nrow(data)),
      fit.times = fit.times, 
      eval.times = eval.times, 
      SY_block_target = SY_block_a0,                  
      SY_block_alt = SY_block_a1,                 
      GY_vec = GY_0,
      g.hats = 1-preds,                  # pi(a0|X)                     (n)
      fN_a0 = fN_0,                   # posterior P(N|Y,dY,X,A=a0)   (n x l)
      fN_alt = fN_1,                  # posterior P(N|Y,dY,X,A=1-a0) (n x l)
      isotonize = TRUE,
      verbose   = TRUE,
      single_cens_mechanism = TRUE
    )  
  }
  
  # ------------------------------------------------------------------
  # Reindex ALL fold-accumulated matrices to ORIGINAL DATA ORDER
  # ------------------------------------------------------------------
  
  make_surv_from_IF <- function(IF_mat, isotonize) {
    #IF_mat    <- IF_mat[fold_to_orig, , drop=FALSE]
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
  if (0L %in% fit.treat) surv.0 <- make_surv_from_IF(IF0$IF.vals, isotonize=TRUE)
  if (1L %in% fit.treat) surv.1 <- make_surv_from_IF(IF1$IF.vals, isotonize=TRUE)
  
  # ==================================================================
  # CONFIDENCE INTERVALS AND CONTRASTS  (unchanged from original)
  # ==================================================================
  .add_curve <- function(surv_obj, trt, event_sub, treat_sub) {
    ft  <- surv_obj$times
    q01 <- quantile(time[event==1 & treat==treat_sub], .01)
    iso <- surv_obj$surv.iso
    q99 <- max(ft[!is.na(iso) & iso >= .01])
    ci  <- .surv.confints(ft, 
                          surv_obj$surv, 
                          surv_obj$IF.vals.center,
                          conf.band=conf.band, 
                          band.end.pts=c(q01,q99),
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
  #result$nuisance  <- nuis
  result
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
.get.survival.AIPWCC.oracle <- function(
    Y, DeltaY, N, DeltaN,
    A,                       # A[i] = 1 iff participant in target arm a0
    admin_t,
    fit.times, eval.times,
    SY_block_target,         # l x K -- target arm's S_Y(t | n) block
    SY_block_alt,            # l x K -- opposite arm (unused in current IF, kept for symmetry)
    GY_vec,                  # length K -- oracle G_Y(t)
    g.hats,                  # length n -- pi(a0 | x)
    fN_a0, fN_alt,           # n x l posterior matrices
    isotonize = TRUE,
    verbose   = TRUE,
    single_cens_mechanism = TRUE
) {
  
  # Y = data$Y;
  # DeltaY = data$DeltaY;
  # N = N_for_lookup;
  # DeltaN = data$DeltaN;
  # A = 1-data$A;                    # A[i]=1 iff participant in target arm a0
  # admin_t = rep((data$admin_t+1), nrow(data));
  # fit.times = fit.times;
  # eval.times = eval.times;
  # SY_block_target = SY_block_a0;
  # SY_block_alt = SY_block_a1;
  # GY_vec = GY_0;
  # g.hats = 1-preds;                 # pi(a0|X)                     (n)
  # fN_a0 = fN_0;                  # posterior P(N|Y,dY,X,A=a0)   (n x l)
  # fN_alt = fN_1;                  # posterior P(N|Y,dY,X,A=1-a0) (n x l)
  # isotonize = TRUE;
  # verbose   = TRUE;
  # single_cens_mechanism = TRUE
  
  fit.times <- sort(fit.times[fit.times > 0])
  n  <- length(Y)
  K  <- length(fit.times)
  l  <- length(eval.times)
  
  Y_ind <- pmax(1L, pmin(findInterval(Y, fit.times), K))
  N_ind <- pmax(1L, pmin(findInterval(N, eval.times), l))
  
  IF.vals <- matrix(NA_real_, n, K)
  
  # Hoist out everything that doesn't depend on participant i
  # GY_vec is length K, same for all participants
  GYi_fun <- stepfun(fit.times, c(1, GY_vec), right = TRUE)
  G_left  <- GYi_fun(fit.times[-K])  # length K-1
  
  if (verbose) pb <- progress_bar$new(total = n, format = "  EIF [:bar] :percent")
  
  for (i in seq_len(n)) {
    N_indx <- N_ind[i]
    fN.i   <- if (A[i] == 1L) fN_a0[i, ] else fN_alt[i, ]
    
    # *** KEY CHANGE: use the shared block, not a per-participant slice ***
    # Under oracle nuisances with rho_X = 0, every participant in the target arm
    # has the same S_Y(t | n) block. So Shi is just SY_block_target.
    Shi <- SY_block_target  # l x K, shared by reference (no copy)
    
    # phi^F compensator
    S_inv  <- 1 / Shi
    S_inv[is.infinite(S_inv)]<-0
    dS_inv <- S_inv[, -1, drop = FALSE] - S_inv[, -K, drop = FALSE]   # l x (K-1)
    vals_inner <- sweep(dS_inv, 2, 1 / G_left, FUN = "*")             # l x (K-1)
    vals_inner[, fit.times[-1] > Y[i]] <- 0
    vals <- cbind(0, t(apply(vals_inner, 1, cumsum)))                 # l x K
    
    S_y <- Shi[, Y_ind[i]]
    G_y <- GYi_fun(Y[i])
    
    # K_N on eval.times grid
    if (single_cens_mechanism) {
      K_N <- sapply(eval.times,
                    FUN = function(t) {
                      DeltaY[i] * GYi_fun(pmax(t, Y[i])) / GYi_fun(Y[i]) +
                        (1 - DeltaY[i]) * as.numeric(t <= Y[i])
                    })
    }
    
    KN_fun       <- stepfun(eval.times, c(1, K_N), right = TRUE)
    KN_inv_left  <- 1 / KN_fun(eval.times[-l])
    KN_inv_left  <- ifelse(is.infinite(KN_inv_left), 0, KN_inv_left)
    dLambdaN     <- 1 - KN_fun(eval.times[-1]) / KN_fun(eval.times[-l])
    dLambdaN     <- ifelse(is.nan(dLambdaN), 0, dLambdaN)
    
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
      #IF_k <- ifelse(is.na(IF_k, na.rm = TRUE), na.rm=TRUE)
      
      num_cs <- rev(cumsum(rev(IF_k * fN.i))) # l x 1
      Q_seq  <- num_cs / pmax(denom_cs, 1e-12)
      
      if(K_N[(N_indx)]==0){
        KNi <- K_N[(N_indx-1)]
      }else{
        KNi <- K_N[N_indx]
      }
      
      term1 <- (    DeltaN[i]  / KNi) * IF_k[N_indx]
      term2 <- ((1 - DeltaN[i]) /KNi) * Q_seq[N_indx]
      
      if (N_indx >= 2L) {
        idx3  <- seq_len(N_indx-1)
        term3 <- sum(Q_seq[idx3] * KN_inv_left[idx3] * dLambdaN[idx3])
      } else {
        term3 <- 0
      }
      
      IF.vals[i, k] <- term1 + term2 - term3
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
