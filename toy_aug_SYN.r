
source("./create_multivar_survival_rand_admin.r")
source("./survSuperLearner_funcs.r")
source("./CFsurvival_funcs.r")
library(CFsurvival)
library(survSuperLearner)
library(survival)
library(tidyverse)
library(randomForestSRC)

set.seed(47)
data <- create_multivar_survival_var_admin(n=10000000, 
                                 trt_prob=0.5, 
                                 RR = 0.75, 
                                 RR_C = 1, 
                                 per_year_incidence_Y = 0.10, 
                                 per_year_incidence_N = 0.12, 
                                 per_year_incidence_C = 1E-16,
                                 J = 1, 
                                 rho = 0.9, 
                                 admin_date = 2,
                                 X=TRUE,
                                 rho_X = 0, 
                                 convert_days = TRUE
)

set.seed(14)
data_small <- create_multivar_survival_var_admin(n=3000, 
                                 trt_prob=0.5, 
                                 RR = 0.75, 
                                 RR_C = 1, 
                                 per_year_incidence_Y = 0.10, 
                                 per_year_incidence_N = 0.12, 
                                 per_year_incidence_C = 0.10,
                                 J = 1, 
                                 rho = 0.9, 
                                 admin_date = 2,
                                 X=TRUE,
                                 rho_X = 0, 
                                 convert_days = TRUE
)

make_df <- function(idx) as.data.frame(cbind(
  id=id, Y=time, DeltaY=event,
  N=time_N, DeltaN=event_N, DeltaN_eff=DeltaN_eff,
  A=treat, admin_t=admin_t, W
))[idx, ]

id=1:nrow(data_small)
time=data_small$Y
event = data_small$DeltaY
time_N = data_small$N1
event_N = data_small$DeltaN1
DeltaN_eff = ifelse(data_small$DeltaN1==1 | data_small$N1==730, 1, 0)
treat = data_small$A
admin_t = data_small$admin_t
W = data_small %>% dplyr::select("X")
W_names = "X"

# ==================================================================
# HELPER: augment training fold with redistribute-to-right weights
# ==================================================================
.augment_fold <- function(dat, W_names, eval.times, survest_N) {
  
  rows <- vector("list", nrow(dat))
  
  for (i in seq_len(nrow(dat))) {
    print(i)
    ri          <- dat[i, ]
    Ai          <- ri$A
    N0          <- ri$N
    admin_t_i   <- ri$admin_t
    xi          <- ri[, W_names, drop=FALSE]
    
    # Unique event times in same arm + admin_t as support
    N_support <- sort(unique(c(eval.times, admin_t_i))) #c(sort(unique(dat$N[dat$DeltaN==1L & dat$A==Ai])), admin_t_i)
    
    base <- data.frame(id=ri$id, Y=ri$Y, DeltaY=ri$DeltaY,
                       A=Ai, admin_t=admin_t_i)
    base <- cbind(base, xi)
    
    # Admin-censored or observed N → single row, weight 1
    if (ri$DeltaN == 1L) {
      rows[[i]] <- cbind(base, N=N0, w=1)
      next
    }
    
    # Randomly censored N → redistribute to right
    N_grid <- N_support[N_support > N0]
    if (length(N_grid) == 0L) {
      rows[[i]] <- cbind(base, N=max(N_support), w=1)
      next
    }
    
    eps   <- 1e-6 * max(1, N0)
    S_all <- survest_N(t = c(N0 - eps, N_grid), x = xi)
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
               #DeltaN = as.integer(N_grid < admin_t_i),
               w      = w)
    rows[[i]] <- r[w > 0, ]
  }
  
  do.call(rbind, rows)
}

.fit_SY_em_old <- function(
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

.compute_aug_likelihood <- function(SY_fit, aug_tr, fit.times, W_names) {
  
  newdat <- cbind(A = aug_tr$A,
                  aug_tr[, c("N", W_names), drop = FALSE])
  
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

dat_tr <- make_df(idx=id)
W_tr   <- dat_tr[, W_names, drop=FALSE]
eval.times <- sort(unique(c(dat_tr$N[dat_tr$DeltaN==1])))
#eval.times <- round(seq(min(dat_tr$N), max(dat_tr$N), length.out=65))

# ----------------------------------------------------------------
# Step 1: S_N | A, X  (pooled arms, A as covariate)
# Smooth library only — no KM — to ensure positive prior density.
# ----------------------------------------------------------------
message("    [1] S_N | A, X")
SN_fit <- survSuperLearner::survSuperLearner(
  time             = dat_tr$N,
  event            = dat_tr$DeltaN,
  X                = cbind(W_tr),
  newX             = dat_tr[1:2,] %>% dplyr::select(X),
  new.times        = eval.times,
  event.SL.library = c("survSL.expreg", "survSL.coxph"),
  cens.SL.library  = c("survSL.expreg", "survSL.coxph"),
  verbose          = FALSE,
  control          = list(initWeightAlg = "survSL.km", verbose = FALSE),
  cvControl        = list(V = 5)
)

survest_N <- function(fit=SN_fit, t, x) {
  # x <- newAW_te[1,2]
  # t=c(100, 200)
  # a=1
  x <- as.data.frame(x)
  colnames(x) <- colnames(W)
  if (nrow(x) == 1L) { x <- rbind(x, x); drop1 <- TRUE } else drop1 <- FALSE
  newdat <- bind_cols(x)
  out <- predict.survSuperLearner(object = fit, 
                                  newdata = newdat, 
                                  new.times = t, 
                                  onlySL = TRUE)$event.SL.predict
  if (drop1) out <- out[1,, drop=FALSE]
  as.vector(out)
}

message("    [2] Augmenting training fold via redistribute-to-the-right")

aug_tr <- .augment_fold(
  dat      = dat_tr,
  W_names  = W_names,
  eval.times = eval.times,
  survest_N  = survest_N
)

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
  fit.cox <- survival::coxph(gam.model, data = X_fit, weights = obsWeights, x = TRUE)
  
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

# ----------------------------------------------------------------
# Step 6: S_Y | A, X, N on augmented training data
# ----------------------------------------------------------------
message("    [3] S_Y | A, X, N")
message("        Iterative EM-style updates to S_{0Y}")  
  
nuis<-list()
nuis$event.SL.library <- c("survSL.coxph.ps")
nuis$cens.SL.library <- c("survSL.km", "survSL.coxph", "survSL.expreg")
nuis$survSL.control <- list(initWeightAlg = "survSL.km", verbose = FALSE)
nuis$survSL.cvControl <- list(V=5)

fit.times <- sort(unique(c(dat_tr$Y, 730)))

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

nuis$event.SL.library_SY <- c("survSL.coxph.ps", "survSL.rfsrc.smooth")
nuis$cens.SL.library <- c("survSL.coxph.ps", "survSL.rfsrc.smooth")

em_out <- .fit_SY_em_old(
    aug_tr            = aug_tr %>% filter(id<=500),
    aug_te            = aug_tr %>% filter(id %in% seq(501, 1000)),
    dat_te            = aug_tr %>% filter(id %in% seq(501, 1000)),
    fit.times         = fit.times,
    W_names           = W_names,
    nuis              = nuis,
    max_iter          = 2,             # 1 disables EM; 2 = one update
    tol               = 1E-3,
    fast_iter_library = c('survSL.coxph.ps'),           # e.g. c("survSL.coxph") for speed
    verbose           = TRUE,
    logN              = FALSE,
    weight_norm       = "L2"
)
  
SY0<-em_out$SY0 # l x K

SY1 <- em_out$SY1

list_p <- vector(mode="list", length=1)
counter=1
eval.times.new <- eval.times
max_indx <- which.max(eval.times.new)
for(s in c(1, 15, 40, 80, 100, 150, 200, 300, 350, 400, 417)){
  
t = eval.times.new[s]
print(counter)

km<-survfit(Surv(Y, DeltaY) ~ 1, data=
              data %>% filter(abs(N1-t)<=0 & DeltaN1==1 & A==0))

km_est <- summary(km, times = unique(sort(fit.times)), extend=TRUE)$surv

list_p[[counter]]<-bind_cols('time'=fit.times, 
                             "est"=SY0[s,], "km"=km_est) %>%
  pivot_longer(cols=2:3) %>%
  mutate(name=factor(name, levels=c("km", "est"), labels=c("Oracle", "RF"))) %>%
  ggplot(aes(x=time, y=value, color=name))+
  geom_step()+
  ylim(c(0, 1))+
  labs(title=paste0("S_{0Y}(y | a, x, n=", t, ")"))
  
counter<-counter+1
}
library(patchwork)
wrap_plots(list_p)

# New approach
#library(devtools)
#devtools::install_github("cwolock/CFsurvival")
#library(CFsurvival)
# library(survML)
# 
# f_y_stack_SuperLearner <- function(time,
#                                    event,
#                                    X,
#                                    censored,
#                                    time_grid,
#                                    bin_size,
#                                    isotonize = TRUE,
#                                    SL_control,
#                                    time_basis){
#   
#   if (!is.null(censored)){
#     if (censored == TRUE){
#       time <- time[!as.logical(event)]
#       X <- X[!as.logical(event),]
#       obsWeights <- SL_control$obsWeights[!as.logical(event)]
#     } else if (censored == FALSE){
#       time <- time[as.logical(event)]
#       X <- X[as.logical(event),]
#       obsWeights <- SL_control$obsWeights[as.logical(event)]
#     }
#   } else{
#     time <- time
#     X <- X
#     obsWeights <- SL_control$obsWeights
#   }
#   
#   if (is.null(time_grid)){
#     bin_variable <- time
#     
#     # if user gives bin size, set time grid based on quantiles. otherwise, every observed time
#     if (!is.null(bin_size)){
#       time_grid <- sort(unique(stats::quantile(bin_variable, type = 1, probs = seq(0, 1, by = bin_size))))
#       time_grid <- c(0, time_grid) # 013123 changed this to try to get better predictions at time 0
#       #time_grid[1] <- 0 # manually set first point to 0, instead of first observed time
#     } else{
#       time_grid <- sort(unique(bin_variable))
#       time_grid <- c(0, time_grid)
#     }
#     
#   }
#   
#   cv_folds <- split(sample(1:length(time)), rep(1:SL_control$V, length = length(time)))
#   
#   X <- as.matrix(X)
#   time <- as.matrix(time)
#   dat <- data.frame(X, time)
#   
#   ids <- seq(1:length(time))
#   
#   if (!is.null(obsWeights)){
#     stackX <- as.matrix(data.frame(X,
#                                    obsWeights = obsWeights,
#                                    ids = ids))
#   } else{
#     stackX <- as.matrix(data.frame(X,
#                                    ids = ids))
#   }
#   
#   stacked <- stack_cdf(time = time,
#                        X = stackX,
#                        time_grid = time_grid)
#   
#   # change t to dummy variable
#   if (time_basis == "dummy"){
#     stacked$t <- factor(stacked$t)
#     dummy_mat <- stats::model.matrix(~-1 + t, data=stacked)
#     risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
#     colnames(dummy_mat) <- risk_set_names
#     stacked$t <- NULL
#     stacked <- cbind(dummy_mat, stacked)
#   }
#   
#   long_obsWeights <- stacked$obsWeights
#   stacked_ids <- stacked$ids
#   stacked$obsWeights <- NULL
#   stacked$ids <- NULL
#   .Y <- stacked[,ncol(stacked)]
#   .X <- stacked[,-ncol(stacked)]
#   
#   get_validRows <- function(fold_sample_ids){
#     validRows <- which(stacked_ids %in% fold_sample_ids)
#     return(validRows)
#   }
#   
#   validRows <- lapply(cv_folds, get_validRows)
#   
#   if (is.null(SL_control$method)){
#     SL_control$method <- "method.NNLS"
#   }
#   if (is.null(SL_control$V)){
#     SL_control$V <- 10
#   }
#   if (is.null(SL_control$SL.library)){
#     SL_control$SL.library <- c("SL.mean")
#   }
#   if (is.null(SL_control$stratifyCV)){
#     SL_control$stratifyCV <- FALSE
#   }
#   
#   fit <- SuperLearner::SuperLearner(Y = .Y,
#                                     X = .X,
#                                     SL.library = SL_control$SL.library,
#                                     family = stats::binomial(),
#                                     method = SL_control$method,
#                                     verbose = FALSE,
#                                     obsWeights = long_obsWeights,
#                                     cvControl = list(V = SL_control$V,
#                                                      validRows = validRows,
#                                                      stratifyCV = SL_control$stratifyCV))
#   
#   fit <- list(reg.object = fit,
#               time_grid = time_grid,
#               isotonize = isotonize,
#               time_basis = time_basis)
#   class(fit) <- c("f_y_stack_SuperLearner")
#   return(fit)
# }
# 
# predict.f_y_stack_SuperLearner <- function(fit, newX, newtimes){
#   
#   time_grid <- fit$time_grid
#   
#   if (fit$time_basis == "continuous"){
#     get_stacked_pred <- function(t){
#       new_stacked <- data.frame(t = t, newX)
#       preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
#       return(preds)
#     }
#     
#     predictions <- matrix(apply(X = matrix(newtimes), FUN = get_stacked_pred, MARGIN = 1),
#                           nrow = nrow(newX))
#   } else{
#     get_preds <- function(t){
#       dummies <- matrix(0, ncol = length(time_grid), nrow = nrow(newX))
#       index <- max(which(time_grid <= t))
#       dummies[,index] <- 1
#       new_stacked <- cbind(dummies, newX)
#       risk_set_names <- paste0("risk_set_", seq(1, (length(time_grid))))
#       colnames(new_stacked)[1:length(time_grid)] <- risk_set_names
#       new_stacked <- data.frame(new_stacked)
#       preds <- stats::predict(fit$reg.object, newdata=new_stacked)$pred
#       return(preds)
#     }
#     
#     predictions <- matrix(apply(X = matrix(newtimes), FUN = get_preds, MARGIN = 1),
#                           nrow = nrow(newX))
#   }
#   
#   if (fit$isotonize){
#     iso.cdf.ests <- t(apply(predictions, MARGIN = 1, FUN = Iso::pava))
#   } else{
#     iso.cdf.ests <- predictions
#   }
#   
#   return(iso.cdf.ests)
# }
# 
# stack_cdf <- function(time, X, time_grid){
#   # we will treat time as continuous
#   ncol_stacked <- ncol(X) + 2 # covariates, time, binary outcome
#   stacked <- matrix(NA, ncol = ncol_stacked, nrow = 1)
#   for (i in 1:(length(time_grid))){ # can change this to not do anything in last time bin
#     event_indicators <- matrix(ifelse(time <= time_grid[i], 1, 0))
#     t <- time_grid[i]
#     newdata <- as.matrix(cbind(t, X, event_indicators))
#     stacked <- rbind(stacked, newdata)
#   }
#   stacked <- stacked[-1,]
#   colnames(stacked)[ncol(stacked)] <- "event_indicators"
#   
#   stacked <- data.frame(stacked)
#   return(stacked)
# }
# 
# vrows = rep(1:5, 300)
# 
# valid_rows <- sapply(1:length(vrows), FUN=function(x){
#   
#   slice = which(aug_tr$id==x)
#   return(slice)
# })
# 
# valid_rows <- sapply(1:5, FUN=function(x){
#   
#   unlist(valid_rows[which(vrows==x)])
#   
# })
# 
# tmp<-f_y_stack_SuperLearner(time = aug_tr$Y,
#                        event = aug_tr$DeltaY,
#                        X     = aug_tr %>% dplyr::select("A", "X", "N", "DeltaN"),
#                        censored = FALSE,
#                        time_grid = sort(unique(aug_tr$Y)),
#                        bin_size = 1,
#                        time_basis = "continuous",
#                        isotonize=TRUE,
#                        SL_control = list(
#                          obsWeights = aug_tr$w,
#                          V         = 5,
#                          validRows = valid_rows,
#                          SL.library = c("SL.mean", "SL.earth", "SL.xgboost")
#                        )
#                        )
# 
# fit.times = sort(unique(data_small$Y))
# 
# l<-predict.f_y_stack_SuperLearner(fit=tmp, 
#                                   newX = data.frame(A=1, X=0, N=eval.times, DeltaN=as.numeric(eval.times<730)), 
#                                   newtimes=fit.times)
# 
# surv <- 1-l
# 
# km<-survfit(Surv(Y, DeltaY) ~ 1, data=data %>% filter(abs(N1-eval.times[5])<=0 & DeltaN1==1 & A==1))
# 
# km_est <- summary(km, times = unique(sort(c(fit.times, 730))), extend=TRUE)$surv
# 
# bind_cols('time'=unique(sort(c(fit.times, 730))), "est"=surv[5,], "km"=km_est) %>%
#   pivot_longer(cols=2:3) %>%
#   ggplot(aes(x=time, y=value, color=name))+
#   geom_step()+
#   ylim(c(0, 1))