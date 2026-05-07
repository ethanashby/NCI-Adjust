# =============================================================================
# AIPWCC_oracle_v2.R
#
# Memory-efficient oracle nuisance computation for AIPWCC, designed for the
# normal-copula DGM in create_multivar_survival.r.
#
# KEY INSIGHTS THAT ENABLE MEMORY SAVINGS
# ----------------------------------------
# 1. In the DGM with rho_X = 0 and randomized A, the oracle nuisances are:
#       pi(a | x)        = 0.5
#       G_C(t | a, x)    = function of t only (C is independent of A, X)
#       f_N(n | a, x)    = function of n only (N is independent of A, X)
#       S_Y(t | a, x, n) = function of (a, t, n) only -- the only non-trivial one
#    So oracle nuisances live on a small grid: K x l x 2 doubles for S_Y plus
#    a few K-vectors. ~1-5 MB total. We compute them ONCE and reuse.
#
# 2. .get.survival.AIPWCC expects S.hats as an (n_aug x K) matrix indexed by
#    (participant, fit.time). But for *oracle* nuisances under our DGM,
#    every participant in arm a gets the same S_Y(t | a, n_j) at grid point j.
#    So we DON'T need an n_aug x K matrix in memory -- we can pass a
#    closure / lookup that indexes (a, n_grid_idx, t_grid_idx) into a small
#    (l x K x 2) tensor. We achieve this by passing a "virtual" big.matrix-
#    like object, or more pragmatically, by building a SMALL (l x K) matrix
#    and writing a thin wrapper around .get.survival.AIPWCC that uses
#    the FOLD-AGNOSTIC indexing.
#
# 3. The auxiliary "near-oracle" KM fits need ~1M rows max, not 25M. With
#    continuous N, we *bin* N onto eval.times and stratify within bins,
#    rather than matching N1 == t exactly.
#
# STRUCTURE
# ---------
#   Step A: simulate auxiliary dataset (1M rows -- adjust by memory)
#   Step B: estimate oracle marginal nuisances
#             G_Y(t), G_N(t), f_N on eval.times, S_Y(t | a, N_bin) on (eval x fit) grid
#   Step C: free auxiliary dataset
#   Step D: build the per-participant nuisance objects WITHOUT replicating
#             huge matrices: each is just a small lookup into the grids
#   Step E: call .get.survival.AIPWCC fold-by-fold with these compact inputs
#
# DEPENDS ON
# ----------
#   AIPWCC_survival_2.r  (for .get.survival.AIPWCC, .expand_test, helpers)
#   create_multivar_survival.r
#   survival, dplyr
#
# USAGE
# -----
#   res <- AIPWCC_oracle(data,
#                        fit.times  = sort(unique(data$Y[data$Y > 0])),
#                        eval.times = sort(unique(c(data$N1, 730))),
#                        admin_t    = 730,
#                        n_aux      = 1e6,
#                        seed_aux   = 42)
#   res$risk.ratio    # for log RR coverage diagnostics
# =============================================================================


# ------------------------------------------------------------------
# Step B helper: build oracle nuisance grids from an aux dataset
# Returns a list of small objects (no big matrices).
# ------------------------------------------------------------------
.build_oracle_nuisances <- function(aux, fit.times, eval.times, admin_t) {
  
  K <- length(fit.times)
  l <- length(eval.times)
  
  # admin censoring recoding -- match what the main estimator does
  aux <- transform(aux,
                   DeltaY_eff = ifelse(Y < admin_t,  DeltaY,  1L),
                   DeltaN_eff = ifelse(N1 < admin_t, DeltaN1, 1L))
  
  # ---- censoring of Y: G_Y(t) ---------------------------------------
  fit_GY <- survival::survfit(survival::Surv(Y, 1L - DeltaY_eff) ~ 1, data = aux)
  GY_fun <- stats::stepfun(fit_GY$time, c(1, fit_GY$surv))
  GY_grid <- pmax(GY_fun(fit.times), 1e-8)   # length K
  
  # ---- censoring of N: G_N(t) on eval.times grid --------------------
  fit_GN <- survival::survfit(survival::Surv(N1, 1L - DeltaN_eff) ~ 1, data = aux)
  GN_fun <- stats::stepfun(fit_GN$time, c(1, fit_GN$surv))
  GN_grid <- pmax(GN_fun(eval.times), 1e-8)   # length l
  
  # ---- marginal S_N(t), f_N(t) on eval.times grid -------------------
  fit_SN <- survival::survfit(survival::Surv(N1, DeltaN1) ~ 1, data = aux)
  SN_fun <- stats::stepfun(fit_SN$time, c(1, fit_SN$surv))
  SN_grid <- SN_fun(eval.times)               # length l
  fN_grid <- pmax(0, c(1, SN_grid[-l]) - SN_grid)
  fN_grid <- fN_grid / sum(fN_grid)
  
  # ---- S_Y(t | a, N=n_j) on (l x K x 2) grid ------------------------
  # Bin N onto eval.times rather than match exactly.
  aux$N_bin <- pmin(findInterval(aux$N1, eval.times,
                                 rightmost.closed = TRUE), l)
  aux$N_bin <- pmax(aux$N_bin, 1L)
  
  SY_array <- array(NA_real_, dim = c(l, K, 2),
                    dimnames = list(NULL, NULL, c("a0","a1")))
  
  for (a in c(0L, 1L)) {
    for (j in seq_len(l)) {
      sub <- aux[aux$A == a & aux$N_bin == j, , drop = FALSE]
      if (nrow(sub) < 2L) {
        # Fall back to wider neighborhood if a bin is sparse
        widen <- 1L
        repeat {
          sub <- aux[aux$A == a &
                       aux$N_bin >= max(1L, j - widen) &
                       aux$N_bin <= min(l, j + widen), , drop = FALSE]
          if (nrow(sub) >= 50L || widen >= 5L) break
          widen <- widen + 1L
        }
      }
      if (nrow(sub) < 2L) {
        SY_array[j, , a + 1L] <- 1
      } else {
        fit_SY <- tryCatch(
          survival::survfit(survival::Surv(Y, DeltaY) ~ 1, data = sub),
          error = function(e) NULL)
        if (is.null(fit_SY)) {
          SY_array[j, , a + 1L] <- 1
        } else {
          SY_fun <- stats::stepfun(fit_SY$time, c(1, fit_SY$surv))
          SY_array[j, , a + 1L] <- SY_fun(fit.times)
        }
      }
    }
  }
  
  list(
    fit.times = fit.times,
    eval.times = eval.times,
    K = K, l = l,
    GY_grid = GY_grid,         # length K
    GN_grid = GN_grid,         # length l
    fN_grid = fN_grid,         # length l
    SN_grid = SN_grid,         # length l
    SY_array = SY_array        # l x K x 2
  )
}


# ------------------------------------------------------------------
# Step D helper: build per-participant compact inputs to
# .get.survival.AIPWCC.
#
# .get.survival.AIPWCC expects:
#   S.hats     : (n_aug x K)  -- aug-aligned via id_to_rows
#   GY.hats_*  : (n   x K)    -- per-participant
#   GN.hats_*  : (n   x l)    -- per-participant
#   fN_*       : (n   x l)    -- per-participant posterior
#
# Under oracle DGM with rho_X=0:
#   GY_a0[i,] = GY_alt[i,] = GY_grid     (same K-vector for everyone)
#   GN_a0[i,] = GN_alt[i,] = GN_grid     (same l-vector for everyone)
#   S.hats[ rows_for_id_i , ] depends only on a (=A_i) and bin(N_j):
#       SY_array[ N_bin , , A_i + 1 ]   -- l x K matrix per participant
#   fN_a0[i,] = posterior of N given (A_i, Y_i, dY_i, X_i) computed from
#               the small grids above -- length l per participant
# ------------------------------------------------------------------
.build_compact_inputs <- function(data, oracle, ex_data) {
  
  fit.times  <- oracle$fit.times
  eval.times <- oracle$eval.times
  K <- oracle$K
  l <- oracle$l
  n <- nrow(data)
  
  # -------- per-participant censoring nuisances (replicated rows) -----
  # These are TINY in absolute terms (n*K and n*l doubles) -- 2400*730 ~ 14 MB
  # so replication is fine. The big-matrix-blowup only happens for S.hats.
  GY_mat <- matrix(oracle$GY_grid, n, K, byrow = TRUE)  # n x K
  GN_mat <- matrix(oracle$GN_grid, n, l, byrow = TRUE)  # n x l
  
  # -------- S.hats: aug-aligned --------------------------------------
  # ex_data has nrow ~ n * l rows; one row per (participant, eval.times[j])
  # S.hats[r, k] = SY_array[ ex_data$N_bin[r] , k , ex_data$A[r] + 1 ]
  #
  # We build this directly without going through any l x K replicate trick.
  # Memory: n_aug x K doubles. With n=2400, l=731, K=730, this is
  #   2400 * 731 * 730 * 8 bytes ~ 10.2 GB.  STILL TOO BIG.
  #
  # WHY THIS MATTERS:
  #   .get.survival.AIPWCC pulls S.hats[rows_for_id_i, ] which is l x K,
  #   meaning for each participant we need an l x K block. Replicating
  #   gives you n*l*K total. You can't avoid l*K per participant, but you
  #   CAN avoid replicating across participants by sharing memory:
  #   under oracle nuisances every participant in arm a uses the SAME
  #   l x K block.
  #
  # SOLUTION: pass two l x K blocks (one per arm) and use a thin wrapper
  # that reroutes id_to_rows lookup to these shared blocks. We implement
  # this via an environment-based "virtual" matrix.
  
  list(
    GY_mat = GY_mat,
    GN_mat = GN_mat,
    SY_block_a0 = oracle$SY_array[, , "a0"],   # l x K, shared
    SY_block_a1 = oracle$SY_array[, , "a1"]    # l x K, shared
  )
}


# ------------------------------------------------------------------
# Posterior f_N | (A, X, Y, dY) per participant, using oracle quantities.
# Memory: n x l doubles ~ 2400*731*8 ~ 14 MB. Fine.
# ------------------------------------------------------------------
.compute_posterior_oracle <- function(data, oracle) {
  
  fit.times  <- oracle$fit.times
  eval.times <- oracle$eval.times
  K <- oracle$K
  l <- oracle$l
  n <- nrow(data)
  
  fN0 <- matrix(NA_real_, n, l)
  fN1 <- matrix(NA_real_, n, l)
  
  prior <- oracle$fN_grid                  # length l, same for all (i, a)
  SY0   <- oracle$SY_array[, , "a0"]       # l x K
  SY1   <- oracle$SY_array[, , "a1"]       # l x K
  
  Y_idx <- pmax(1L, pmin(findInterval(data$Y, fit.times,
                                      rightmost.closed = TRUE), K))
  
  for (i in seq_len(n)) {
    Yi  <- data$Y[i]
    dYi <- data$DeltaY[i]
    yk  <- Y_idx[i]
    
    if (dYi == 0L) {
      lk0 <- SY0[, yk]
      lk1 <- SY1[, yk]
    } else {
      if (yk == 1L) {
        S_lo0 <- rep(1, l); S_lo1 <- rep(1, l)
      } else {
        S_lo0 <- SY0[, yk - 1L]; S_lo1 <- SY1[, yk - 1L]
      }
      lk0 <- pmax(0, S_lo0 - SY0[, yk])
      lk1 <- pmax(0, S_lo1 - SY1[, yk])
    }
    
    u0 <- lk0 * prior; u1 <- lk1 * prior
    s0 <- sum(u0);     s1 <- sum(u1)
    fN0[i, ] <- if (s0 > 1e-12) u0 / s0 else prior
    fN1[i, ] <- if (s1 > 1e-12) u1 / s1 else prior
  }
  
  list(fN0 = fN0, fN1 = fN1)
}


# ------------------------------------------------------------------
# Memory-efficient version of .get.survival.AIPWCC for the oracle case
# where S_Y(t | a, X, N) depends only on (a, N, t) -- so all participants
# in arm a share a single (l x K) S_Y block. We avoid the n_aug x K
# replication entirely.
#
# This is a near-verbatim copy of the original loop body, with the
# only change being the line:
#       Shi <- S.hats[rows_i, , drop=FALSE]
# replaced by:
#       Shi <- if (A[i] == 1L) SY_block_a0 else SY_block_alt
# (note: a0 = the target arm passed in; we follow the same A[i]==1
#  iff in target arm convention as the original)
# ------------------------------------------------------------------
.get.survival.AIPWCC.oracle <- function(
    Y, DeltaY, N, DeltaN, A, admin_t,
    fit.times, eval.times,
    SY_block_a0,        # l x K, target arm
    SY_block_alt,       # l x K, opposite arm (unused here, kept for symmetry)
    GY_vec,             # length K, oracle G_Y
    GN_vec,             # length l, oracle G_N
    g_hats,             # length n, oracle pi(a0 | x) = 0.5
    fN_a0, fN_alt,      # n x l posterior matrices
    isotonize = TRUE,
    verbose = TRUE
) {
  fit.times <- sort(fit.times[fit.times > 0])
  n <- length(Y)
  K <- length(fit.times)
  l <- length(eval.times)
  
  Y_ind <- pmax(1L, pmin(findInterval(Y, fit.times, rightmost.closed = TRUE), K))
  N_ind <- pmax(1L, pmin(findInterval(N, eval.times, rightmost.closed = TRUE), l))
  
  IF.vals <- matrix(NA_real_, n, K)
  
  # GN-related quantities -- same for everyone under oracle, hoist out
  GN_eval_full <- c(1, GN_vec)                  # length l+1
  GN_inv       <- 1 / pmax(GN_eval_full, 1e-8)  # length l+1
  dLambdaN     <- 1 - GN_eval_full[-1] / GN_eval_full[-(l + 1)]   # length l
  GNi_fun      <- stats::stepfun(eval.times, GN_eval_full, right = TRUE)
  
  # GY step function -- same for everyone under oracle
  GYi_fun  <- stats::stepfun(fit.times, c(1, GY_vec), right = FALSE)
  GYi_vec  <- GY_vec
  
  if (verbose) {
    pb <- progress::progress_bar$new(total = n,
                                     format = "  EIF [:bar] :percent ETA: :eta")
  }
  
  for (i in seq_len(n)) {
    
    Shi <- if (A[i] == 1L) SY_block_a0 else SY_block_alt   # l x K
    
    fN.i <- if (A[i] == 1L) fN_a0[i, ] else fN_alt[i, ]
    
    # phi^F compensator: vals[j,k] = int_0^{min(t_k, Y_i)} d(1/S_Y) / G_C^Y
    S_inv  <- 1 / pmax(Shi, 1e-12)
    dS_inv <- S_inv[, -1, drop = FALSE] - S_inv[, -K, drop = FALSE]
    vals   <- sweep(dS_inv, 2, 1 / pmax(GYi_vec[-K], 1e-8), FUN = "*")
    vals[, fit.times[-1] > Y[i]] <- 0
    vals   <- cbind(0, t(apply(vals, 1, cumsum)))             # l x K
    
    S_y    <- Shi[, Y_ind[i]]                                  # length l
    G_y    <- GYi_fun(Y[i])
    
    GN_at_Ni <- GNi_fun(N[i])
    
    denom_cs   <- rev(cumsum(rev(fN.i)))
    IPW_factor <- as.numeric(A[i] == 1L) / g_hats[i]
    
    for (k in seq_len(K)) {
      
      S_t0 <- Shi[, k]
      
      inner1 <- if (Y[i] <= fit.times[k] && DeltaY[i] == 1L)
        1 / pmax(S_y * G_y, 1e-12) else rep(0, l)
      
      IF_k  <- S_t0 - IPW_factor * S_t0 * (inner1 - vals[, k])
      
      num_cs <- rev(cumsum(rev(IF_k * fN.i)))
      Q_seq  <- num_cs / pmax(denom_cs, 1e-12)
      
      term1 <- (    DeltaN[i]   / GN_at_Ni) * IF_k[N_ind[i]]
      term2 <- ((1 - DeltaN[i]) / GN_at_Ni) * Q_seq[N_ind[i]]
      
      idx3  <- seq_len(N_ind[i])
      term3 <- sum(Q_seq[idx3] * GN_inv[idx3] * dLambdaN[idx3])
      
      IF.vals[i, k] <- term1 + term2 - term3
    }
    
    if (verbose) pb$tick()
  }
  
  surv           <- colMeans(IF.vals, na.rm = TRUE)
  IF.vals.center <- sweep(IF.vals, 2, surv, "-")
  
  res <- list(
    times          = fit.times,
    surv           = pmin(1, pmax(0, surv)),
    IF.vals        = IF.vals,
    IF.vals.center = IF.vals.center
  )
  if (isotonize) {
    idx <- !is.na(surv)
    res$surv.iso <- surv
    res$surv.iso[idx] <- 1 - stats::isoreg(fit.times[idx], 1 - surv[idx])$yf
  }
  res
}


# ==================================================================
# MAIN ENTRY POINT
# ==================================================================
AIPWCC_oracle <- function(
    data,
    fit.times,
    eval.times,
    admin_t       = 730,
    n_aux         = 1e6,
    seed_aux      = 42L,
    aux_dgm_args  = list(trt_prob = 0.5, RR = 0.75, RR_C = 1,
                         per_year_incidence_Y = 0.06,
                         per_year_incidence_N = 0.12,
                         per_year_incidence_C = 0.10,
                         J = 1, rho = 0.8, admin_date = 2,
                         X = TRUE, rho_X = 0, convert_days = TRUE),
    conf.level    = 0.95,
    verbose       = TRUE
) {
  
  stopifnot(all(c("Y","DeltaY","N1","DeltaN1","A") %in% names(data)))
  
  n <- nrow(data)
  fit.times  <- sort(unique(fit.times[fit.times > 0]))
  eval.times <- sort(unique(eval.times))
  K <- length(fit.times); l <- length(eval.times)
  
  if (verbose) message(sprintf(
    "AIPWCC_oracle: n=%d, K=%d (fit.times), l=%d (eval.times), n_aux=%d",
    n, K, l, n_aux))
  
  # ----- Step A: simulate aux data ------------------------------------
  if (verbose) message("  [A] Simulating auxiliary dataset...")
  set.seed(seed_aux)
  aux <- do.call(create_multivar_survival, c(list(n = n_aux), aux_dgm_args))
  
  # ----- Step B: oracle nuisances on grids ----------------------------
  if (verbose) message("  [B] Building oracle nuisance grids...")
  oracle <- .build_oracle_nuisances(aux, fit.times, eval.times, admin_t)
  
  # ----- Step C: free aux ---------------------------------------------
  rm(aux); invisible(gc(verbose = FALSE))
  if (verbose) message(sprintf(
    "  [C] Aux freed. Oracle grids: ~%.2f MB",
    (object.size(oracle) / 1024^2)))
  
  # ----- Step D: per-participant nuisances (compact) ------------------
  if (verbose) message("  [D] Computing posterior f_N for each participant...")
  post <- .compute_posterior_oracle(data, oracle)
  
  # DeltaN_eff: admin censoring → 1 (matches main estimator)
  DeltaN_eff <- as.integer(data$DeltaN1 == 1L | data$N1 == admin_t)
  
  g_hats <- rep(0.5, n)   # randomized
  
  # ----- Step E: arm-specific IF using compact get.survival.AIPWCC ----
  if (verbose) message("  [E] Arm 1 IF computation...")
  eif1 <- .get.survival.AIPWCC.oracle(
    Y = data$Y, DeltaY = data$DeltaY,
    N = data$N1, DeltaN = DeltaN_eff,
    A = data$A,                       # A==1 iff in target arm 1
    admin_t = rep(admin_t, n),
    fit.times = fit.times, eval.times = eval.times,
    SY_block_a0  = oracle$SY_array[, , "a1"],   # target arm = 1
    SY_block_alt = oracle$SY_array[, , "a0"],
    GY_vec  = oracle$GY_grid,
    GN_vec  = oracle$GN_grid,
    g_hats  = g_hats,
    fN_a0   = post$fN1,
    fN_alt  = post$fN0,
    isotonize = TRUE,
    verbose = verbose
  )
  
  if (verbose) message("  [E] Arm 0 IF computation...")
  eif0 <- .get.survival.AIPWCC.oracle(
    Y = data$Y, DeltaY = data$DeltaY,
    N = data$N1, DeltaN = DeltaN_eff,
    A = 1L - data$A,                  # A==1 iff in target arm 0
    admin_t = rep(admin_t, n),
    fit.times = fit.times, eval.times = eval.times,
    SY_block_a0  = oracle$SY_array[, , "a0"],   # target arm = 0
    SY_block_alt = oracle$SY_array[, , "a1"],
    GY_vec  = oracle$GY_grid,
    GN_vec  = oracle$GN_grid,
    g_hats  = 1 - g_hats,
    fN_a0   = post$fN0,
    fN_alt  = post$fN1,
    isotonize = TRUE,
    verbose = verbose
  )
  
  # ----- Step F: log RR + IF-based SE ---------------------------------
  if (verbose) message("  [F] Building risk-ratio output...")
  
  # Use the existing risk.ratio helper if loaded. Otherwise compute inline.
  if (exists(".risk.ratio", mode = "function")) {
    rr <- .risk.ratio(
      fit.times    = fit.times,
      surv.0       = eif0$surv,
      IF.vals.0    = eif0$IF.vals.center,
      surv.1       = eif1$surv,
      IF.vals.1    = eif1$IF.vals.center,
      conf.band    = FALSE,
      band.end.pts = c(min(fit.times), max(fit.times)),
      conf.level   = conf.level
    )$risk.ratio.df
  } else {
    risk.0 <- 1 - eif0$surv
    risk.1 <- 1 - eif1$surv
    IF_logRR <- (-eif1$IF.vals.center) / pmax(risk.1, 1e-8) -
      (-eif0$IF.vals.center) / pmax(risk.0, 1e-8)
    se_logRR <- apply(IF_logRR, 2, function(x) sqrt(mean(x^2) / n))
    log_rr   <- log(pmax(risk.1, 1e-8) / pmax(risk.0, 1e-8))
    qz       <- qnorm(1 - (1 - conf.level)/2)
    rr <- data.frame(
      time             = fit.times,
      log.risk.ratio   = log_rr,
      se.log.ratio     = se_logRR,
      ptwise.lower     = exp(log_rr - qz * se_logRR),
      ptwise.upper     = exp(log_rr + qz * se_logRR)
    )
  }
  
  list(
    risk.ratio = rr,
    surv.0     = eif0,
    surv.1     = eif1,
    oracle     = oracle,
    posterior  = post
  )
}


# =============================================================================
# DIAGNOSTICS
# =============================================================================
#
# At target landmark (e.g. t = 730):
#   row <- subset(res$risk.ratio, time == max(time))
#   row$log.risk.ratio                # point estimate
#   row$se.log.ratio                  # IF-based SE
#
# Run this across many simulated datasets:
#   logRR_vec <- replicate(400, {
#     d <- create_multivar_survival(n = 2400, ...)
#     r <- AIPWCC_oracle(d, fit.times = ..., eval.times = ...,
#                        n_aux = 1e6)$risk.ratio
#     row <- subset(r, time == max(time))
#     c(logRR = row$log.risk.ratio, SE = row$se.log.ratio)
#   })
#
#   emp_SD     <- sd(logRR_vec["logRR", ])
#   mean_SE    <- mean(logRR_vec["SE", ])
#   SE_ratio   <- mean_SE / emp_SD
#   coverage   <- mean(abs(logRR_vec["logRR", ] - log(0.75)) <
#                      qnorm(0.975) * logRR_vec["SE", ])
#
# CASE A: SE_ratio ~ 1.00 and coverage ~ 95%
#   ==> EIF derivation correct; original undercoverage is nuisance-driven.
#       Recommend cross-fit aggregation (AIPWCC_aggregate.r).
#
# CASE B: SE_ratio < 1 (e.g. 0.92) even with oracle nuisances
#   ==> EIF derivation has a structural issue. Most likely the
#       Tsiatis 10.4 conditioning -- the censoring filtration in the IF
#       should condition on (A, X, Y_tilde, Delta_Y), not just (A, X).
#       Path A in earlier discussion: rewrite term3/IPCW with V-conditional
#       compensator.
#
# AUX SAMPLE SIZE GUIDANCE
# ------------------------
# n_aux = 1e6 is enough for stable KM estimates within each (a, N_bin) cell
# given l ~ 730 bins -- ~700 obs per cell, halved by arm. If you see noisy
# S_Y blocks (jagged surfaces in N), bump n_aux to 5e6 or coarsen
# eval.times. Memory cost of 1e6 aux: ~80 MB (transient, freed before IF).
# =============================================================================