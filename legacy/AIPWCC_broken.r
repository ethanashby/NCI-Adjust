# AIPWCC version of CFSurvival

#' Estimate counterfactual survival functions while adjusting for a right-censored predictor (negative control infection)
#'
#' This function estimates counterfactual survival functions and contrasts from right-censored data subject to potential confounding, including adjustment for a negative control time-to-event predictor.
#'
#' For a binary treatment measured at baseline, the counterfactual survival through time \code{t} for treatment level \code{a} is defined as the probability of survival up to time \code{t} were all units in the population assigned to treatment \code{a}. This quantity is estimated using an augmented inverse probability of treatment/censoring weighted (AIPTW) estimator. Our estimator of the counterfactual survival function requires further estimation of FOUR nuisance parameters: the conditional probability of receiving treatment given a set of potential confounders (i.e. the treatment propensity), (NEW) the conditional survival function of the negative control event given a set of potential confounders, the conditional survival function of the event given confounders and the negative control event time, and the conditional survival function of censoring given confounders. At the moment, we provide two ways of estimating each of these nuisance parameters. The treatment propensity can be estimated either using a logistic regression or SuperLearner (van der Laan and Polley, 2006). The conditional survival functions of negative control event, primary event, and censoring can be estimated using single estimators or a SuperLearner ensemble. Additionally, the estimates of any of these parameters can be passed in directly if another method is desired.
#'
#'  \code{AIPWCC_CFsurvival} by default returns 95\% pointwise confidence intervals and 95\% uniform confidence band, both of which are based on the influence function of the estimator. Contrasts of the counterfactual survival under treatment and control can also be obtained, including the survival difference, survival ratio, risk ratio, and number needed to treat. Pointwise confidence intervals and uniform confidence bands for these contrasts will also be computed.
#'
#' @param time \code{n x 1} numeric vector of observed right-censored follow-up times for the primary event; i.e. the minimum of the event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of whether a primary event was observed.
#' @param time_N \code{n x 1} numeric vector of observed right-censored follow-up times for the negative control event; i.e. the minimum of the event and censoring times.
#' @param event_N \code{n x 1} numeric vector of status indicators of whether a negative control event was observed.
#' @param admin_t \code{n x 1} numeric vector of administrative censoring times for the negative control event.
#' @param treat \code{n x 1} numeric vector of observed binary treatment/exposure group; either 0 or 1.
#' @param W \code{n x p} numeric matrix of completely-observed potential confounders used to model conditional survival, censoring, and treatment probabilities. Missing values are not allowed.
#' @param id \code{n x 1} a vector of unique participant identifiers. Defaults to 1:n.
#' @param fit.times \code{k x 1} numeric vector of time grid at which the counterfactual survival function should be computed. Only values > 0 are allowed. Note that even if the survival is only desired at one timepoint (e.g. a final study time), fit.times should be a relatively fine grid between 0 and this final time in order for the computation to work. Defauls to all unique positive values in \code{time}.
#' @param fit.treat Optional subset of \code{c(0,1)} for which the counterfactual survival curves are desired. For example, if \code{fit.treat=c(0,1)} (default behavior) then both curves will be estimated, and if \code{fit.treat=0} then only the curve under non-treatment will be estimated.
#' @param nuisance.options List of options for nuisance parameter specification. See \code{\link{CFsurvival.nuisance.options}} for details. Defaults to SuperLearners with rich libraries.
#' @param conf.band Logical indicating whether to compute simultaneous confidence bands.
#' @param conf.level Desired coverage of confidence intervals/bands.
#' @param contrasts Character vector indicating which (if any) contrasts of the survival functions are desired. Can include \code{"surv.diff"} (difference of survival functions), \code{"surv.ratio"} (ratio of survival functions), \code{"risk.ratio"} (ratio of the risk functions), or \code{"nnt"} (number needed to treat; i.e. recipricol of survival difference). If not \code{NULL}, then \code{fit.treat} must be \code{c(0,1)}. Defaults to \code{c("surv.diff", "surv.ratio")}.
#' @param verbose Logical indicating whether progress should be printed to the command line.
#' @return \code{CFsurvival} returns a named list with the following elements:
#' \item{fit.times}{The time points at which the counterfactual survival curves (and contrasts) were fit.}
#' \item{fit.treat}{The treatment values for which the counterfactual survival curves were fit.}
#' \item{surv.df}{A data frame with the estimated counterfactual survival functions, estimated standard errors, confidence interval (ptwise.lower/ptwise.upper is constructed on the survival scale, and ptwise.logit.lower/ptwise.logit.upper is constructed on the logit scale), confidence bands (unif.ew.lower/unif.ew.upper is an equal-width confidence band on the survival scale, and unif.logit.lower/unif.logit.upper is an equi-precision confidence band on the logit scale).}
#' \item{IF.vals.0, IF.vals.1}{n x k matrices with the influence values for the counterfactual survival functions. Rows index observations, columns index time points.}
#' \item{surv.diff.df}{A data frame with the estimated counterfactual survival difference (treatment survival minus control survival), as well as confidence intervals and tests of the null hypothesis that the difference equals zero. Only returned if \code{surv.diffs=TRUE}.}
#' \item{surv.ratio.df}{A data frame with the estimated counterfactual survival ratio (treatment survival divided by control survival), as well as confidence intervals and tests of the null hypothesis that the ratio equals one. Note that standard errors and confidence intervals are first computed on the log scale, then exponentiated. Only returned if \code{surv.ratios=TRUE}.}
#' \item{risk.ratio.df}{A data frame with the estimated counterfactual risk ratio (treatment risk divided by control risk), as well as confidence intervals and tests of the null hypothesis that the ratio equals one. Note that standard errors and confidence intervals are first computed on the log scale, then exponentiated. Only returned if \code{risk.ratios=TRUE}.}
#' \item{nnt.df}{A data frame with the estimated counterfactual number needed to treat (one divided by the survival difference), as well as confidence intervals and tests of the null hypothesis that the ratio equals one. Note that standard errors and confidence intervals are first computed on the log scale, then exponentiated. Only returned if \code{nnt=TRUE}.}
#' \item{surv.0.unif.ew.sim.maxes, surv.0.logit.sim.maxes, surv.1.unif.ew.sim.maxes, surv.1.logit.sim.maxes}{Samples from the approximate distribution of the maximum of the limiting Gaussian processes of the counterfactual survival functions.}
#' \item{surv.0.unif.ew.quant, surv.0.logit.quant, surv.1.unif.ew.quant, surv.1.logit.quant, surv.diff.unif.quant, log.risk.ratio.unif.quant, log.surv.ratio.unif.quant, nnt.unif.quant}{Estimated quantile of the approximate distribution of the maximum of the limiting Gaussian processes of the counterfactual survival functions or contrast function.}
#' \item{prop.fits}{The estimated treatment propensities (if \code{surv.nuis.fits = TRUE}).}
#' \item{surv.fits}{The estimated conditional survival functions (if \code{surv.nuis.fits = TRUE}).}
#' \item{nuisance}{Predicted nuisance function values and details of the nuisance function fits.}
#' \item{data}{The original time, event, and treatment data supplied to the function.}
#' @details \code{surv.df$ptwise.lower} and \code{surv.df$ptwise.upper} are lower and upper endpoints of back-transformed pointwise confidence intervals on the logit scale. Pointwise confidence intervals are only provided for values of \code{fit.times} that are at least as large as the smallest observed event time in the corresponding treatment cohort. \code{surv.df$unif.ew.lower} and \code{surv.df$unif.ew.upper} are lower and upper endpoints of an equal-width confidence band. \code{surv.df$unif.logit.lower} and \code{surv.df$unif.logit.upper} are lower and upper endpoints of an equi-precision confidence band on the logit scale, and trasformed back to the survival scale. This confidence band only covers values of \code{fit.times} such that the corresponding estimated survival function is \code{<= .99} and \code{>= .01}.


AIPWCC_CFsurvival <- function(time, 
                              event,
                              time_N,
                              event_N,
                              admin_t,
                              treat, 
                              W,
                              id,
                              fit.times=sort(unique(time[time > 0])), 
                              fit.treat=c(0,1), 
                              nuisance.options = list(), 
                              conf.band=FALSE, 
                              conf.level=.95, 
                              contrasts = c("risk.ratio"), 
                              verbose=TRUE) {
  
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
  #   V = 8,
  #   folds = NULL,
  #   eval.times = sort(unique(c(data$N1[data$N1>0]))),
  #   #eval.times = seq(min(data$N1[data$N1>0]), max(data$N1[data$N1>0]), length.out=100),
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
  
  .args <- mget(names(formals()), sys.frame(sys.nframe()))
  
  # if(is.null(contrasts)) contrasts <- c("surv.diff", "surv.ratio")
  
  nuis <- do.call("CFsurvival.nuisance.options", nuisance.options)
  
  n <- length(time)
  
  if(length(unique(admin_t))==1){
    admin_t = rep(admin_t, n)
  }else if(length(admin_t)==n){
    admin_t=admin_t
  }else{
    stop("admin_t must be either a n x 1 vector or a scalar")
  }
  
  if(sum(1-event) == 0) {
    message("No censored events; unanticipated errors may occur.")
    nuis$cens.pred.0 <- nuis$cens.pred.1 <- matrix(1, nrow=n, ncol=k)
  }
  if(sum(event) == 0) {
    stop("No uncensored events; cannot perform estimation.")
  }
  
  if(!nuis$cross.fit) {
    if(!is.null(nuis$folds)) {
      message("nuisance.options$cross.fit set to FALSE, but V > 1 or folds specified. Cross-fitting will not be performed.")
    }
    nuis$V <- 1
    nuis$folds <- rep(1, n)
  }
  
  if(!is.null(nuis$folds)) {
    nuis$folds <- as.numeric(factor(nuis$folds))
    nuis$V <- length(unique(nuis$folds))
  }
  if(nuis$cross.fit & is.null(nuis$folds)) {
    
    event.0 <- which(event == 0)
    event.1 <- which(event == 1)
    folds.0 <- sample(rep(1:nuis$V, length = length(event.0)))
    folds.1 <- sample(rep(1:nuis$V, length = length(event.1)))
    nuis$folds <- rep(NA, n)
    nuis$folds[event.0] <- folds.0
    nuis$folds[event.1] <- folds.1
  }
  
  if(is.null(nuis$verbose)) nuis$verbose <- FALSE
  
  contrasts <- unique(tolower(contrasts))
  contrasts <- contrasts[contrasts %in% c("surv.diff", "surv.ratio", "risk.ratio", "nnt")]
  
  .check.input(time=time, 
               event=event, 
               treat=treat, 
               W_event = W,
               W_cens = W,
               fit.times=fit.times, 
               fit.treat=fit.treat, 
               nuisance.options=nuis, 
               conf.band=conf.band, 
               conf.level=conf.level, 
               contrasts=contrasts, 
               verbose=verbose)
  
  if(any(fit.times < 0)) {
    fit.times <- fit.times[fit.times > 0]
    message("fit.times < 0 removed.")
  }
  
  if(is.null(nuis$eval.times)){
    nuis$eval.times <- fit.times #sort(unique(c(0,time[time > 0 & time <= max(fit.times)], max(fit.times))))
  }
  
  k <- length(fit.times)
  
  W <- as.data.frame(W)
  
  surv.df <- data.frame()
  result <- list(fit.times=fit.times, fit.treat=fit.treat, surv.df=surv.df)
  
  #### ESTIMATE PROPENSITY ####
  
  if(is.null(nuis$prop.pred)) {
    if(verbose) message("Estimating propensities...")
    nuis$prop.pred <- rep(NA, n)
    if(nuis$V > 1) {
      if(nuis$save.nuis.fits) result$prop.fits <- vector(mode='list',length=nuis$V)
      for(v in 1:nuis$V) {
        if(verbose) message(paste("Fold ", v, "..."))
        train <- nuis$folds != v
        test <- nuis$folds == v
        
        prop.fit <- .estimate.propensity(A=treat[train], 
                                         W=W[train,,drop=FALSE], 
                                         newW=W[test,, drop=FALSE], 
                                         SL.library=nuis$prop.SL.library, 
                                         fit.treat=fit.treat, 
                                         prop.trunc=nuis$prop.trunc, 
                                         save.fit = nuis$save.nuis.fits,
                                         verbose = nuis$verbose)
        
        nuis$prop.pred[test] <- prop.fit$prop.pred
        if(nuis$save.nuis.fits) result$prop.fits[[v]] <- prop.fit$prop.fit
      }
    } else {
      prop.fit <- .estimate.propensity(A=treat, 
                                       W=W, 
                                       newW=W, 
                                       SL.library=nuis$prop.SL.library, 
                                       fit.treat=fit.treat, 
                                       prop.trunc=nuis$prop.trunc, 
                                       save.fit = nuis$save.nuis.fits, 
                                       verbose = FALSE)
      nuis$prop.pred <- prop.fit$prop.pred
      if(nuis$save.nuis.fits) {
        result$prop.fit <- prop.fit$prop.fit
      }
    }
  }
  
  #### ESTIMATE CONDITIONAL SURVIVALS ####
  
  if((1 %in% fit.treat & (is.null(nuis$event.pred.1) | is.null(nuis$cens.pred.1))) | (0 %in% fit.treat & (is.null(nuis$event.pred.0) | is.null(nuis$cens.pred.0)))) {
    if(verbose) message("Estimating conditional survivals...")
    
    # =============================================================================
    # AIPWCC_censoring_fix_v2.r
    #
    # Corrected estimate_censoring_survival() with process-specific dropout logic.
    #
    # Two separate functions:
    #   estimate_censoring_survival_Y() – fits G_C^rand for the Y process
    #   estimate_censoring_survival_N() – fits G_C^rand for the N process
    #
    # And a wrapper estimate_censoring_survivals() that fits both and returns
    # cens.pred.0 / cens.pred.1 for each, ready to plug into the fold loop.
    #
    # KEY CORRECTION from v1:
    #   The v1 logic said: if Delta_Y=1 OR Delta_N=1, treat as non-dropout.
    #   This is wrong. Example: Y occurs at day 5, N is censored at day 10
    #   (well before admin_t). That participant dropped out of the N-observation
    #   process at day 10, and this should be a random censoring event in
    #   G_C^rand for N, regardless of what happened to Y.
    #
    #   The correct logic is process-specific:
    #
    #   For G_C^Y (censoring of the Y process):
    #     random_cens_Y = (Delta_Y == 0) & (time_Y < admin_t)
    #     g_time_Y      = time_Y
    #     Interpretation: participant dropped out of Y-observation before admin_t.
    #     Admin-censored Y participants (time_Y == admin_t, Delta_Y == 0) are
    #     NOT random dropouts — they simply reached end of follow-up.
    #
    #   For G_C^N (censoring of the N process):
    #     random_cens_N = (Delta_N == 0) & (time_N < admin_t)
    #     g_time_N      = time_N
    #     Interpretation: participant dropped out of N-observation before admin_t.
    #     This correctly captures the case where Y=5, N-censored=10 < admin_t:
    #     random_cens_N = 1, g_time_N = 10. ✓
    #
    #   Note: if the study design ensures a SINGLE dropout process that
    #   simultaneously ends observation of both Y and N (e.g. participant
    #   physically leaves the study), then time_Y == time_N for censored
    #   participants and both indicators agree. In that case combining
    #   information across processes (v1 approach, corrected) would be:
    #     dropout_time = min(time_Y, time_N)
    #     random_cens  = (Delta_Y == 0) & (Delta_N == 0) &
    #                    (time_Y < admin_t) & (time_N < admin_t)
    #   But this is a special case. The process-specific approach below
    #   is always valid and is the safe default.
    # =============================================================================
    
    
    # =============================================================================
    # estimate_censoring_survival_Y()
    #
    # Fits G_C^rand(t | A, X) for the Y process.
    #
    #   random_cens_Y[i] = 1  iff  Delta_Y[i] == 0  AND  time_Y[i] < admin_t[i]
    #   g_time_Y[i]      = time_Y[i]
    #
    # Administrative censorings (time_Y == admin_t, Delta_Y == 0) contribute
    # to the risk set up to admin_t but are NOT censoring events.
    # G_C^rand(admin_t | A, X) is therefore bounded away from 0.
    # =============================================================================
    
    estimate_censoring_survival_Y <- function(
    time_Y,           # observed tilde_Y (length n_train)
    event_Y,          # Delta_Y          (length n_train)
    admin_t,          # admin censoring times (length n_train)
    A,                # treatment indicator (length n_train)
    W,                # covariate data frame (n_train x p)
    newW,             # covariates for prediction (n_pred x p); usually dat_test W
    newA_vals = c(0, 1),  # treatment values to predict at
    new.times,        # time grid for prediction (length k)
    cens.SL.library,
    survSL.control   = NULL,
    survSL.cvControl = list(V = 5),
    cens.trunc       = 0,
    verbose          = FALSE
    ) {
      # ------------------------------------------------------------------
      # Construct dropout indicator for Y process
      # ------------------------------------------------------------------
      random_cens_Y <- as.integer(event_Y == 0 & time_Y < admin_t)
      g_time_Y      <- time_Y   # follow-up time in G model = observed Y time
      
      # ------------------------------------------------------------------
      # Fit G_C^rand on training data
      # ------------------------------------------------------------------
      AW_train <- as.data.frame(cbind(A = A, W))
      
      # Build prediction matrix: stack A=0 and A=1 blocks for newW
      newAW <- do.call(rbind, lapply(newA_vals, function(a) cbind(A = a, newW)))
      
      if (is.null(survSL.control)) survSL.control <- list()
      
      fit_G <- survSuperLearner::survSuperLearner(
        time             = g_time_Y,
        event            = random_cens_Y,
        X                = AW_train,
        newX             = as.data.frame(newAW),
        new.times        = new.times,
        event.SL.library = cens.SL.library,
        cens.SL.library  = cens.SL.library,
        verbose          = verbose,
        control          = survSL.control,
        cvControl        = survSL.cvControl
      )
      
      # ------------------------------------------------------------------
      # Extract and clean predictions
      # ------------------------------------------------------------------
      n_pred <- nrow(newW)
      preds  <- fit_G$event.SL.predict   # (n_pred * length(newA_vals)) x k
      
      trunc_clean <- function(m) {
        m <- pmax(m, cens.trunc)
        if (any(m == 0)) m[m == 0] <- min(m[m > 0])
        m
      }
      
      result <- list(fit = fit_G, g_time = g_time_Y, random_cens = random_cens_Y)
      
      for (j in seq_along(newA_vals)) {
        rows <- ((j - 1) * n_pred + 1):(j * n_pred)
        result[[paste0("cens.pred.", newA_vals[j])]] <- trunc_clean(preds[rows, , drop = FALSE])
      }
      
      result
    }
    
    
    # =============================================================================
    # estimate_censoring_survival_N()
    #
    # Fits G_C^rand(t | A, X) for the N process.
    #
    #   random_cens_N[i] = 1  iff  Delta_N[i] == 0  AND  time_N[i] < admin_t[i]
    #   g_time_N[i]      = time_N[i]
    #
    # This correctly handles the case: Y=5, N-censored=10 < admin_t → dropout
    # at day 10 is a genuine random censoring event in the N process. ✓
    #
    # Administrative censorings of N (time_N == admin_t, Delta_N == 0) are
    # NOT censoring events. G_C^rand(admin_t | A, X) is bounded away from 0.
    # =============================================================================
    
    estimate_censoring_survival_N <- function(
    time_N,           # observed tilde_N (length n_train)
    event_N,          # Delta_N          (length n_train)
    admin_t,          # admin censoring times (length n_train)
    A,                # treatment indicator (length n_train)
    W,                # covariate data frame (n_train x p)
    newW,             # covariates for prediction (n_pred x p); usually dat_test W
    newA_vals = c(0, 1),
    new.times,        # time grid for prediction (length k)
    cens.SL.library,
    survSL.control   = NULL,
    survSL.cvControl = list(V = 5),
    cens.trunc       = 0,
    verbose          = FALSE
    ) {
      # ------------------------------------------------------------------
      # Construct dropout indicator for N process
      # ------------------------------------------------------------------
      random_cens_N <- as.integer(event_N == 0 & time_N < admin_t)
      g_time_N      <- time_N
      
      # ------------------------------------------------------------------
      # Fit G_C^rand on training data
      # ------------------------------------------------------------------
      AW_train <- as.data.frame(cbind(A = A, W))
      newAW    <- do.call(rbind, lapply(newA_vals, function(a) cbind(A = a, newW)))
      
      if (is.null(survSL.control)) survSL.control <- list()
      
      fit_G <- survSuperLearner::survSuperLearner(
        time             = g_time_N,
        event            = random_cens_N,
        X                = AW_train,
        newX             = as.data.frame(newAW),
        new.times        = new.times,
        event.SL.library = cens.SL.library,
        cens.SL.library  = cens.SL.library,
        verbose          = verbose,
        control          = survSL.control,
        cvControl        = survSL.cvControl
      )
      
      # ------------------------------------------------------------------
      # Extract and clean predictions
      # ------------------------------------------------------------------
      n_pred <- nrow(newW)
      preds  <- fit_G$event.SL.predict
      
      trunc_clean <- function(m) {
        m <- pmax(m, cens.trunc)
        if (any(m == 0)) m[m == 0] <- min(m[m > 0])
        m
      }
      
      result <- list(fit = fit_G, g_time = g_time_N, random_cens = random_cens_N)
      
      for (j in seq_along(newA_vals)) {
        rows <- ((j - 1) * n_pred + 1):(j * n_pred)
        result[[paste0("cens.pred.", newA_vals[j])]] <- trunc_clean(preds[rows, , drop = FALSE])
      }
      
      result
    }
    
    estimate_conditional_survival_N <- function(dat, 
                                                N_name, 
                                                DeltaN_name, 
                                                A_name, 
                                                X_names, 
                                                admin_t, 
                                                fit.treat=c(0,1), 
                                                event.SL.library, 
                                                cens.SL.library,
                                                survSL.control, 
                                                survSL.cvControl, 
                                                verbose){
      
      req <- require(survSuperLearner)
      
      if(is.null(X_names)){
        
        X_names = "X"
        dat <- dat %>%
          mutate(
            X=1,
            DeltaN = case_when(
              get(N_name) >= admin_t ~ 1,
              .default = get(DeltaN_name)),
            N = case_when(
              get(N_name) > admin_t ~ admin_t,
              .default = get(N_name))
          )
        
      }else{
        dat=dat %>%
          mutate(
            DeltaN = case_when(
              get(N_name) >= admin_t ~ 1,
              .default = get(DeltaN_name)),
            N = case_when(
              get(N_name) > admin_t ~ admin_t,
              .default = get(N_name))
          ) %>%
          dplyr::select(X_names, N, DeltaN, A_name)
      }
      
      if(all(c(0,1) %in% fit.treat)){
        
        FN_1 <- survSuperLearner::survSuperLearner(time = dat[['N']][dat[[A_name]]==1], 
                                                   event = dat[['DeltaN']][dat[[A_name]]==1],
                                                   X = as.data.frame(dat[which(dat[[A_name]]==1),which(colnames(dat) %in% X_names),drop=FALSE]),
                                                   newX = NULL,
                                                   new.times = sort(unique(dat[['N']])),
                                                   event.SL.library = event.SL.library,
                                                   cens.SL.library = cens.SL.library,
                                                   verbose=verbose, 
                                                   control = survSL.control, 
                                                   cvControl = survSL.cvControl)
        
        FN_0 <- survSuperLearner(time = dat[['N']][dat[[A_name]]==0], 
                                 event = dat[['DeltaN']][dat[[A_name]]==0],
                                 X = as.data.frame(dat[which(dat[[A_name]]==0),which(colnames(dat) %in% X_names),drop=FALSE]),
                                 newX = NULL,
                                 new.times = sort(unique(dat[['N']])),
                                 event.SL.library = event.SL.library,
                                 cens.SL.library = cens.SL.library,
                                 control = survSL.control, 
                                 cvControl = survSL.cvControl,
                                 verbose=verbose)
        
        ### These functions provide a means of prediction for general model for F_{N|X, A}
        survest_1 <- function(t, x){
          
          x <- as.data.frame(x)
          
          if (nrow(x) == 1) {
            x <- bind_rows(x, x)
            drop_first <- TRUE
          } else {
            drop_first <- FALSE
          }
          
          out<-predict.survSuperLearner.2(object = FN_1, 
                                          newdata = x,
                                          new.times = t,
                                          X = as.data.frame(dat[which(dat[[A_name]]==1),which(colnames(dat) %in% X_names),drop=FALSE]), 
                                          time = dat[['N']][dat[[A_name]]==1],
                                          event = dat[['DeltaN']][dat[[A_name]]==1])$event.SL.predict
          
          if(drop_first){
            out <- out[1, , drop = FALSE]
          }
          
          as.vector(out)
        }
        
        survest_0 <- function(t, x){
          
          x <- as.data.frame(x)
          
          if (nrow(x) == 1) {
            x <- bind_rows(x, x)
            drop_first <- TRUE
          } else {
            drop_first <- FALSE
          }
          
          out<-predict.survSuperLearner.2(object = FN_0, 
                                          newdata = x,
                                          new.times = t,
                                          X = as.data.frame(dat[which(dat[[A_name]]==0),which(colnames(dat) %in% X_names),drop=FALSE]), 
                                          time = dat[['N']][dat[[A_name]]==0],
                                          event = dat[['DeltaN']][dat[[A_name]]==0])$event.SL.predict
          
          if(drop_first){
            out <- out[1, , drop = FALSE]
          }
          
          as.vector(out)
          
          
        }
      }
      
      if(all(fit.treat==0)){
        
        FN_0 <- survSuperLearner(time = dat[['N']][dat[[A_name]]==0], 
                                 event = dat[['DeltaN']][dat[[A_name]]==0],
                                 X = as.data.frame(dat[which(dat[[A_name]]==0),which(colnames(dat) %in% X_names),drop=FALSE]),
                                 newX = NULL,
                                 new.times = sort(unique(dat[['N']])),
                                 event.SL.library = event.SL.library,
                                 cens.SL.library = cens.SL.library,
                                 control = survSL.control, 
                                 cvControl = survSL.cvControl,
                                 verbose=verbose)
        
        ### These functions provide a means of prediction for general model for F_{N|X, A}
        survest_1 <- NULL
        
        survest_0 <- function(t, x){
          
          x <- as.data.frame(x)
          
          if (nrow(x) == 1) {
            x <- bind_rows(x, x)
            drop_first <- TRUE
          } else {
            drop_first <- FALSE
          }
          
          out<-predict.survSuperLearner.2(object = FN_0, 
                                          newdata = x,
                                          new.times = t,
                                          X = as.data.frame(dat[which(dat[[A_name]]==0),which(colnames(dat) %in% X_names),drop=FALSE]), 
                                          time = dat[['N']][dat[[A_name]]==0],
                                          event = dat[['DeltaN']][dat[[A_name]]==0])$event.SL.predict
          
          if(drop_first){
            out <- out[1, , drop = FALSE]
          }
          
          as.vector(out)
          
          
        }
      }
      
      if(all(fit.treat==1)){
        
        FN_1 <- survSuperLearner(time = dat[['N']][dat[[A_name]]==1], 
                                 event = dat[['DeltaN']][dat[[A_name]]==1],
                                 X = as.data.frame(dat[which(dat[[A_name]]==1),which(colnames(dat) %in% X_names),drop=FALSE]),
                                 newX = NULL,
                                 new.times = sort(unique(dat[['N']])),
                                 event.SL.library = event.SL.library,
                                 cens.SL.library = cens.SL.library,
                                 control = survSL.control, 
                                 cvControl = survSL.cvControl,
                                 verbose=verbose)
        
        ### These functions provide a means of prediction for general model for F_{N|X, A}
        survest_1 <- function(t, x){
          
          x <- as.data.frame(x)
          
          if (nrow(x) == 1) {
            x <- bind_rows(x, x)
            drop_first <- TRUE
          } else {
            drop_first <- FALSE
          }
          
          out<-predict.survSuperLearner.2(object = FN_1, 
                                          newdata = x,
                                          new.times = t,
                                          X = as.data.frame(dat[which(dat[[A_name]]==1),which(colnames(dat) %in% X_names),drop=FALSE]), 
                                          time = dat[['N']][dat[[A_name]]==1],
                                          event = dat[['DeltaN']][dat[[A_name]]==1])$event.SL.predict
          
          if(drop_first){
            out <- out[1, , drop = FALSE]
          }
          
          as.vector(out)
        }
        
        survest_0 <- NULL
      }
      
      return(
        list(
          #FN_1 = FN_1,
          #FN_0 = FN_0,
          survest_1 = survest_1,
          survest_0 = survest_0
        )
      ) 
    }
    
    augment_dataset <- function(i, df, id, A, Y, DeltaY, W, N, DeltaN, admin_t, SN){
      
      survest_1 <- SN$survest_1
      survest_0 <- SN$survest_0
      
      row_i <- df[i, ]
      
      Ai = row_i[[A]]
      N0 <- row_i[[N]]  # censoring time
      admin_t_tmp <- row_i[[admin_t]]
      
      N_support <- c(sort(unique(df[[N]][df[[DeltaN]]==1 & df[[A]]==Ai])), admin_t_tmp)
      
      # Case 1: N observed
      if (row_i[[DeltaN]] == 1 | row_i[[N]]==admin_t_tmp) {
        if(all(is.na(W) | is.null(W))){
          return(
            data.frame(
              id     = row_i[[id]],
              N      = row_i[[N]],
              DeltaN = row_i[[DeltaN]],
              Y      = row_i[[Y]],
              DeltaY = row_i[[DeltaY]],
              A      = row_i[[A]],
              X      = NA,
              admin_t = admin_t_tmp,
              w      = 1
            ))
        }else{
          return(
            suppressMessages(bind_cols(data.frame(
              id     = row_i[[id]],
              N      = row_i[[N]],
              DeltaN = row_i[[DeltaN]],
              Y      = row_i[[Y]],
              DeltaY = row_i[[DeltaY]],
              A      = row_i[[A]],
              admin_t = admin_t_tmp,
              w      = 1
            ), row_i %>% select(all_of(W)), .name_repair="minimal")))
        }
      }
      
      # Case 2: N censored → redistribution
      
      # plausible N values
      N_grid <- N_support[N_support > N0]
      
      if (length(N_grid) == 0){
        return(suppressMessages(bind_cols(data.frame(
          id     = row_i[[id]],
          N      = admin_t_tmp,
          DeltaN = 0,
          Y      = row_i[[Y]],
          DeltaY = row_i[[DeltaY]],
          A      = row_i[[A]],
          admin_t = admin_t_tmp,
          w      = 1
        ), row_i %>% select(all_of(W)), .name_repair="minimal")))
      }
      
      # conditional survival from fitted model
      # FN_survfun returns S_N(t | a, x)
      
      if(Ai==0){
        S_vals <- as.vector(survest_0(
          t=c(N0, N_grid), 
          x=row_i %>% dplyr::select(all_of(W))
        ))
      }else{
        S_vals <- as.vector(survest_1(
          t=c(N0, N_grid), 
          x=row_i %>% dplyr::select(all_of(W))
        ))
      }
      
      # compute conditional probabilities
      S_N0 <- S_vals[1]
      S_grid <- S_vals
      
      # P(N = n | N > N0)
      w_raw <- -diff(c(S_grid)) / S_N0
      w_raw[w_raw < 0] <- 0
      w <- w_raw / sum(w_raw)
      
      if(all(is.na(W) | is.null(W))){
        return(data.frame(
          id     = row_i[[id]],
          N      = N_grid,
          DeltaN = as.numeric(N_grid < admin_t_tmp),
          Y      = row_i[[Y]],
          DeltaY = row_i[[DeltaY]],
          A      = row_i[[A]],
          X      = NA,
          admin_t = admin_t_tmp,
          w      = w
        ) %>% filter(w>0))
      }else{
        return(suppressMessages(bind_cols(data.frame(
          id     = row_i[[id]],
          N      = N_grid,
          DeltaN = as.numeric(N_grid < admin_t_tmp),
          Y      = row_i[[Y]],
          DeltaY = row_i[[DeltaY]],
          A      = row_i[[A]],
          admin_t = admin_t_tmp,
          w      = w
        ) %>% filter(w>0), row_i %>% select(all_of(W)), .name_repair="minimal")))
      }
    }
    
    # --- helper: find the discrete interval index for a scalar time value
    .find_interval_index <- function(t_val, time_grid) {
      # Returns j such that time_grid[j] is the largest grid value <= t_val.
      # If t_val < time_grid[1] returns 1 (first interval).
      idx <- findInterval(t_val, time_grid, rightmost.closed = TRUE)
      pmax(1L, pmin(idx, length(time_grid)))
    }
    
    compute_posterior_fN <- function(dat_test,
                                     aug_test,
                                     surv_fit_event_pred_0,
                                     surv_fit_event_pred_1,
                                     SN,
                                     fit.times,
                                     eval.times,
                                     W_names,
                                     verbose = TRUE) {
      
      # Returns fN_0 and fN_1: matrices of dimension (nrow(dat_test) x length(eval.times))
      # Row i, column j = P(N = eval.times[j] | Yi, DeltaYi, Xi, A=0or1)
      #
      # Arguments:
      #   dat_test               – original (un-augmented) test fold data frame
      #   aug_test               – test fold expanded over eval.times grid
      #   surv_fit_event_pred_0  – S_Y(t | A=0, X, N): matrix (nrow(aug_test) x length(fit.times))
      #   surv_fit_event_pred_1  – S_Y(t | A=1, X, N): matrix (nrow(aug_test) x length(fit.times))
      #   SN                     – list with survest_0, survest_1 from estimate_conditional_survival_N
      #   fit.times              – grid used for Y survival prediction
      #   eval.times             – grid used for N (rows of aug_test cycle through this)
      #   W_names                – covariate column names
      
      n_test   <- nrow(dat_test)
      l        <- length(eval.times)
      fN_1     <- matrix(NA_real_, nrow = n_test, ncol = l)
      fN_0     <- matrix(NA_real_, nrow = n_test, ncol = l)
      
      if (verbose) {
        pb <- progress_bar$new(
          total  = n_test,
          format = "Posterior fN [:bar] :percent ETA: :eta"
        )
      }
      
      for (i in seq_len(n_test)) {
        
        Yi      <- dat_test[i, "Y"]
        DeltaYi <- dat_test[i, "DeltaY"]
        idi     <- dat_test[i, "id"]
        xi      <- dat_test[i, W_names, drop = FALSE]
        
        # Rows of aug_test that belong to this individual (one per eval.time grid point)
        indx_aug <- which(aug_test[["id"]] == idi)
        
        # ------------------------------------------------------------------
        # Prior: marginal PMF of N on eval.times grid under each arm
        # f_N(n_j | A=a, X=xi)  ∝  S_N(n_{j-1}^-) – S_N(n_j)
        # We use the left-limit trick: prepend 1 (= S at time 0) to S_grid
        # ------------------------------------------------------------------
        SN1_grid <- as.vector(SN$survest_1(t = eval.times, x = xi))  # S_N(t|A=1,X)
        SN0_grid <- as.vector(SN$survest_0(t = eval.times, x = xi))  # S_N(t|A=0,X)
        
        # Discrete PMF: P(N = eval.times[j] | A=a, X)
        prior_1  <- pmax(0, c(1, SN1_grid[-length(SN1_grid)]) - SN1_grid)
        prior_0  <- pmax(0, c(1, SN0_grid[-length(SN0_grid)]) - SN0_grid)
        
        # ------------------------------------------------------------------
        # Likelihood: L(Yi, DeltaYi | A=a, X, N = eval.times[j])
        # One value per row of aug_test (= per grid point j of N)
        # ------------------------------------------------------------------
        # FIX (b): for censored Y, use S_Y evaluated at the last fit.time <= Yi
        # FIX (a): for observed Y,  use the discrete hazard mass at the interval
        #          containing Yi, found robustly via findInterval
        
        if (DeltaYi == 0) {
          # Censored: likelihood = S_Y(Yi | A, X, N)
          # Use last fit.time <= Yi (right-continuous step function)
          y_idx   <- .find_interval_index(Yi, fit.times)
          # S_Y at fit.times[y_idx] per aug_test row (one row per N grid point)
          lik_1   <- surv_fit_event_pred_1[indx_aug, y_idx]  # length l
          lik_0   <- surv_fit_event_pred_0[indx_aug, y_idx]
          
        } else {
          # Observed event: likelihood = f_Y(Yi | A, X, N) ≈ S_Y(t_{j-1}) - S_Y(t_j)
          # where [t_{j-1}, t_j] is the fit.times interval containing Yi.
          # FIX (a): use findInterval to locate Yi robustly.
          y_lo_idx <- .find_interval_index(Yi, fit.times)   # last fit.time <= Yi
          y_hi_idx <- min(y_lo_idx + 1L, length(fit.times))
          
          # S_Y at the two bracketing time points (matrix: nrow=l, ncol=2)
          S_lo_1 <- if (y_lo_idx == 1) rep(1, length(indx_aug)) else
            surv_fit_event_pred_1[indx_aug, y_lo_idx - 1L]
          S_hi_1 <- surv_fit_event_pred_1[indx_aug, y_lo_idx]
          
          S_lo_0 <- if (y_lo_idx == 1) rep(1, length(indx_aug)) else
            surv_fit_event_pred_0[indx_aug, y_lo_idx - 1L]
          S_hi_0 <- surv_fit_event_pred_0[indx_aug, y_lo_idx]
          
          lik_1  <- pmax(0, S_lo_1 - S_hi_1)   # length l
          lik_0  <- pmax(0, S_lo_0 - S_hi_0)
          
          # If all likelihoods are zero (flat survival at Yi) fall back to S_Y(Yi)
          if (all(lik_1 == 0)) lik_1 <- surv_fit_event_pred_1[indx_aug, y_lo_idx]
          if (all(lik_0 == 0)) lik_0 <- surv_fit_event_pred_0[indx_aug, y_lo_idx]
        }
        
        # ------------------------------------------------------------------
        # Posterior: P(N = n_j | Yi, DeltaYi, A=a, X)
        #            ∝  lik_a(n_j) * prior_a(n_j)
        # FIX (c): use the *same-arm* prior for each arm's posterior
        # ------------------------------------------------------------------
        unnorm_1 <- lik_1 * prior_1
        unnorm_0 <- lik_0 * prior_0
        
        denom_1  <- sum(unnorm_1) + 1e-12
        denom_0  <- sum(unnorm_0) + 1e-12
        
        fN_1[i, ] <- unnorm_1 / denom_1
        fN_0[i, ] <- unnorm_0 / denom_0
        
        if (verbose) pb$tick()
      }
      
      list(fN_0 = fN_0, fN_1 = fN_1)
    }
    
    # =============================================================================
    # Revised per-fold block
    #
    # Now fits FOUR censoring models per fold:
    #   G_C^Y(t | A=0, X),  G_C^Y(t | A=1, X)   – for phi^F IPCW weights
    #   G_C^N(t | A=0, X),  G_C^N(t | A=1, X)   – for AIPWCC terms 1/2/3
    #
    # G_C^Y predictions go into nuis$cens.pred.0 / nuis$cens.pred.1
    #   (used as G.hats_a0 and G.hats_1_minus_a0 in .get.survival.AIPWCC_fast
    #    for the phi^F compensator integral and inner1 terms)
    #
    # G_C^N predictions go into nuis$censN.pred.0 / nuis$censN.pred.1
    #   (used as the N-censoring G denominator in terms 1, 2, 3)
    #
    # NOTE: .get.survival.AIPWCC_fast will need a small update to accept
    #   separate G_C^N arguments for the AIPWCC terms. See the comment at
    #   the bottom of this file.
    # =============================================================================
    
    if (nuis$V > 1) {
      
      if (nuis$save.nuis.fits) result$surv.fits <- vector(mode = "list", length = nuis$V)
      nuis$event.coef <- nuis$cens.coef <- list()
      
      # Fold-order accumulators
      acc_event_pred_0  <- matrix(NA, nrow = 0, ncol = k)
      acc_event_pred_1  <- matrix(NA, nrow = 0, ncol = k)
      acc_censY_pred_0  <- matrix(NA, nrow = 0, ncol = k)   # G_C^Y | A=0
      acc_censY_pred_1  <- matrix(NA, nrow = 0, ncol = k)   # G_C^Y | A=1
      acc_censN_pred_0  <- matrix(NA, nrow = 0, ncol = k)   # G_C^N | A=0
      acc_censN_pred_1  <- matrix(NA, nrow = 0, ncol = k)   # G_C^N | A=1
      acc_fN_0          <- matrix(NA, nrow = 0, ncol = length(nuis$eval.times))
      acc_fN_1          <- matrix(NA, nrow = 0, ncol = length(nuis$eval.times))
      acc_aug_ids       <- vector(length = 0)   # aug_test id per row
      acc_test_ids      <- vector(length = 0)   # dat_test id per row
      
      for (v in seq_len(nuis$V)) {
        
        if (verbose) message(paste("Fold", v, "..."))
        
        train_idx <- which(nuis$folds != v)
        test_idx  <- which(nuis$folds == v)
        
        dat_train <- as.data.frame(cbind(
          id = id, Y = time, DeltaY = event,
          N = time_N, DeltaN = event_N,
          A = treat, admin_t = admin_t, W
        ))[train_idx, ]
        
        dat_test <- as.data.frame(cbind(
          id = id, Y = time, DeltaY = event,
          N = time_N, DeltaN = event_N,
          A = treat, admin_t = admin_t, W
        ))[test_idx, ]
        
        W_train <- dat_train %>% dplyr::select(all_of(colnames(W)))
        W_test  <- dat_test  %>% dplyr::select(all_of(colnames(W)))
        n_test  <- nrow(dat_test)
        
        # ------------------------------------------------------------------
        # Step 1(a): Fit S_N | A, X
        # ------------------------------------------------------------------
        if (verbose) message("\tStep 1(a): Estimate S_N | A, X")
        SN <- estimate_conditional_survival_N(
          dat              = dat_train,
          N_name           = "N",
          DeltaN_name      = "DeltaN",
          A_name           = "A",
          X_names          = colnames(W),
          admin_t          = dat_train$admin_t,
          fit.treat        = c(0, 1),
          event.SL.library = c("survSL.coxph", "survSL.expreg",
                               "survSL.weibreg", "survSL.gam"),
          cens.SL.library  = nuis$cens.SL.library,
          survSL.control   = nuis$survSL.control,
          survSL.cvControl = nuis$survSL.cvControl,
          verbose          = FALSE
        )
        
        # ------------------------------------------------------------------
        # Step 1(b): Augment training fold
        # ------------------------------------------------------------------
        if (verbose) message("\tStep 1(b): Augmenting train fold")
        aug_train_list <- vector(mode = "list", length = nrow(dat_train))
        for (i in seq_len(nrow(dat_train))) {
          aug_train_list[[i]] <- augment_dataset(
            i, df = dat_train, id = "id", Y = "Y", DeltaY = "DeltaY",
            A = "A", N = "N", W = colnames(W), DeltaN = "DeltaN",
            admin_t = "admin_t", SN = SN
          )
        }
        aug_train <- do.call(rbind, aug_train_list)
        rm(aug_train_list)
        
        # ------------------------------------------------------------------
        # Step 1(c): Expand test fold over eval.times grid
        # ------------------------------------------------------------------
        if (verbose) message("\tStep 1(c): Expanding test fold over N grid")
        aug_test_list <- vector(mode = "list", length = n_test)
        for (i in seq_len(n_test)) {
          aug_test_list[[i]] <- bind_cols(
            dat_test[i, ] %>% dplyr::select(-c("N", "DeltaN")),
            data.frame(
              N      = nuis$eval.times[nuis$eval.times <= dat_test[i, ]$admin_t],
              DeltaN = c(
                rep(1L, length(nuis$eval.times[nuis$eval.times <= dat_test[i,]$admin_t]) - 1L),
                0L
              )
            )
          )
        }
        aug_test <- do.call(rbind, aug_test_list)
        rm(aug_test_list)
        
        # ------------------------------------------------------------------
        # Step 1(d): Fit dedicated G_C^Y (censoring of Y process)
        #
        # random_cens_Y = (Delta_Y == 0) & (time_Y < admin_t)
        # Adminstrative censorings of Y are NOT censoring events here.
        # ------------------------------------------------------------------
        if (verbose) message("\tStep 1(d): Fitting G_C^Y (censoring of Y process)")
        GY_fit <- estimate_censoring_survival_Y(
          time_Y           = dat_train$Y,
          event_Y          = dat_train$DeltaY,
          admin_t          = dat_train$admin_t,
          A                = dat_train$A,
          W                = W_train,
          newW             = W_test,
          newA_vals        = c(0, 1),
          new.times        = fit.times,
          cens.SL.library  = nuis$cens.SL.library,
          survSL.control   = nuis$survSL.control,
          survSL.cvControl = nuis$survSL.cvControl,
          cens.trunc       = nuis$cens.trunc,
          verbose          = FALSE
        )
        
        # ------------------------------------------------------------------
        # Step 1(e): Fit dedicated G_C^N (censoring of N process)
        #
        # random_cens_N = (Delta_N == 0) & (time_N < admin_t)
        # Correctly handles: Y=5, N-censored=10 < admin_t → dropout at 10. ✓
        # ------------------------------------------------------------------
        if (verbose) message("\tStep 1(e): Fitting G_C^N (censoring of N process)")
        GN_fit <- estimate_censoring_survival_N(
          time_N           = dat_train$N,
          event_N          = dat_train$DeltaN,
          admin_t          = dat_train$admin_t,
          A                = dat_train$A,
          W                = W_train,
          newW             = W_test,
          newA_vals        = c(0, 1),
          new.times        = fit.times,
          cens.SL.library  = nuis$cens.SL.library,
          survSL.control   = nuis$survSL.control,
          survSL.cvControl = nuis$survSL.cvControl,
          cens.trunc       = nuis$cens.trunc,
          verbose          = FALSE
        )
        
        # ------------------------------------------------------------------
        # Step 1(f): Fit S_0Y | A, X, N on augmented training fold.
        # Pass G_C^Y predictions for cens.pred (used inside phi^F).
        # The internal cens model of survSuperLearner2 is still fit but
        # its predictions are discarded.
        # ------------------------------------------------------------------
        if (verbose) message("\tStep 1(f): Fitting S_0Y | A, X, N on augmented data")
        surv.fit <- .estimate.conditional.survival.2(
          Y                  = aug_train[["Y"]],
          Delta              = aug_train[["DeltaY"]],
          A                  = aug_train[["A"]],
          W_event            = aug_train %>% dplyr::select(c("N", "DeltaN", colnames(W))),
          W_cens             = aug_train %>% dplyr::select(all_of(colnames(W))),
          newW_event         = aug_test  %>% dplyr::select(c("N", "DeltaN", colnames(W))),
          newW_cens          = aug_test  %>% dplyr::select(all_of(colnames(W))) %>% distinct(),
          fit.times          = fit.times,
          fit.treat          = fit.treat,
          event.SL.library   = nuis$event.SL.library,
          cens.SL.library    = nuis$cens.SL.library,
          survSL.control     = nuis$survSL.control,
          survSL.cvControl   = nuis$survSL.cvControl,
          cens.trunc         = nuis$cens.trunc,
          save.fit           = nuis$save.nuis.fits,
          verbose            = nuis$verbose,
          obsWeights         = aug_train$w,
          id                 = aug_train$id#,
          #G_pred_0_external  = GY_fit$cens.pred.0,   # G_C^Y | A=0
          #G_pred_1_external  = GY_fit$cens.pred.1    # G_C^Y | A=1
        )
        
        # ------------------------------------------------------------------
        # Step 2: Compute posterior f_N via Bayes rule
        # ------------------------------------------------------------------
        if (verbose) message("\tStep 2: Computing posterior f_N")
        posteriors <- compute_posterior_fN(
          dat_test              = dat_test,
          aug_test              = aug_test,
          surv_fit_event_pred_0 = surv.fit$event.pred.0,
          surv_fit_event_pred_1 = surv.fit$event.pred.1,
          SN                    = SN,
          fit.times             = fit.times,
          eval.times            = nuis$eval.times,
          W_names               = colnames(W),
          verbose               = TRUE
        )
        
        # ------------------------------------------------------------------
        # Accumulate in fold order
        # ------------------------------------------------------------------
        acc_event_pred_0  <- rbind(acc_event_pred_0,  surv.fit$event.pred.0)
        acc_event_pred_1  <- rbind(acc_event_pred_1,  surv.fit$event.pred.1)
        acc_censY_pred_0  <- rbind(acc_censY_pred_0,  GY_fit$cens.pred.0)
        acc_censY_pred_1  <- rbind(acc_censY_pred_1,  GY_fit$cens.pred.1)
        acc_censN_pred_0  <- rbind(acc_censN_pred_0,  GN_fit$cens.pred.0)
        acc_censN_pred_1  <- rbind(acc_censN_pred_1,  GN_fit$cens.pred.1)
        acc_fN_0          <- rbind(acc_fN_0,          posteriors$fN_0)
        acc_fN_1          <- rbind(acc_fN_1,          posteriors$fN_1)
        acc_aug_ids       <- c(acc_aug_ids,            aug_test$id)
        acc_test_ids      <- c(acc_test_ids,           dat_test$id)
        
        if (nuis$save.nuis.fits) result$surv.fits[[v]] <- surv.fit$surv.fit
        nuis$event.coef[[v]] <- surv.fit$event.coef
        nuis$cens.coef[[v]]  <- surv.fit$cens.coef
        
        rm(surv.fit, GY_fit, GN_fit, aug_train, aug_test)
        gc()
        if (verbose) message(paste("Fold", v, "complete"))
      }
      
      # ------------------------------------------------------------------
      # Reindex all fold-accumulated matrices to ORIGINAL DATA ORDER
      # ------------------------------------------------------------------
      orig_ids     <- id
      fold_to_orig <- match(orig_ids, acc_test_ids)
      
      # Per-participant matrices (n x k)
      nuis$cens.pred.0  <- acc_censY_pred_0[fold_to_orig, , drop = FALSE]
      nuis$cens.pred.1  <- acc_censY_pred_1[fold_to_orig, , drop = FALSE]
      nuis$censN.pred.0 <- acc_censN_pred_0[fold_to_orig, , drop = FALSE]
      nuis$censN.pred.1 <- acc_censN_pred_1[fold_to_orig, , drop = FALSE]
      nuis$fN_0         <- acc_fN_0[fold_to_orig,          , drop = FALSE]
      nuis$fN_1         <- acc_fN_1[fold_to_orig,          , drop = FALSE]
      
      # Aug-test-aligned rows (multiple per participant)
      aug_row_order     <- order(match(acc_aug_ids, orig_ids))
      nuis$id           <- acc_aug_ids[aug_row_order]
      nuis$event.pred.0 <- acc_event_pred_0[aug_row_order, , drop = FALSE]
      nuis$event.pred.1 <- acc_event_pred_1[aug_row_order, , drop = FALSE]
      
      stopifnot(all(unique(nuis$id) == orig_ids))
      
    } else {
      stop("Please specify V > 1 to ensure cross-fitting is used.")
    }
  }
  
  #### ESTIMATE CF SURVIVALS ####
  if(verbose) message("Computing counterfactual survivals...")
  
  .get.survival.AIPWCC_fast <- function(
    Y,                      # observed tilde_Y                   (n)
    DeltaY,                 # Delta_Y                            (n)
    N,                      # observed tilde_N                   (n)
    DeltaN,                 # Delta_N                            (n)
    A,                      # A[i]=1 iff participant in target arm a0  (n)
    admin_t,                # administrative censoring times     (n)
    fit.times,              # evaluation grid for survivor curve (K)
    eval.times,             # evaluation grid for N              (l)
    ids,                    # id vector aligned with S.hats rows (n_aug)
    S.hats,                 # S_Y(t|a0,X,N):  aug-aligned matrix (n_aug x K)
    GY.hats_a0,             # G_C^Y(t|a0,X):  per-participant   (n x K)
    GY.hats_1_minus_a0,     # G_C^Y(t|1-a0,X):per-participant   (n x K)
    GN.hats_a0,             # G_C^N(t|a0,X):  per-participant   (n x K)
    GN.hats_1_minus_a0,     # G_C^N(t|1-a0,X):per-participant   (n x K)
    g.hats,                 # pi(a0|X)                           (n)
    fN_a0,                  # posterior P(N=n_j|Y,DeltaY,X,A=a0)    (n x l)
    fN_1_minus_a0,          # posterior P(N=n_j|Y,DeltaY,X,A=1-a0)  (n x l)
    isotonize = TRUE,
    verbose   = TRUE
  ) {
    # ------------------------------------------------------------------
    # Argument key (target arm = a0):
    #
    #  A[i] = 1  iff  participant i received treatment a0.
    #  Caller passes A = 1-treat for the a0=0 curve,
    #                A = treat   for the a0=1 curve.
    #
    #  GY.*  – censoring of Y; drives phi^F (vals.i, G_y)
    #  GN.*  – censoring of N; drives AIPWCC terms 1/2/3
    #
    #  S.hats rows are in aug_test order: for participant with id=u,
    #  the corresponding rows are ids == u (one per eval.times grid point).
    # ------------------------------------------------------------------
    
    # Y                    = time;
    # DeltaY               = event;
    # N                    = time_N;
    # DeltaN               = event_N;
    # A                    = 1-treat;             # A[i]=1 iff treat[i]==1
    # admin_t              = admin_t;
    # fit.times            = fit.times;
    # eval.times           = nuis$eval.times;
    # ids                  = nuis$id;
    # S.hats               = nuis$event.pred.0;
    # GY.hats_a0           = nuis$cens.pred.0;   # G_C^Y | A=1; for phi^F
    # GY.hats_1_minus_a0   = nuis$cens.pred.1;   # G_C^Y | A=0; for phi^F
    # GN.hats_a0           = nuis$censN.pred.0;  # G_C^N | A=1; for AIPWCC
    # GN.hats_1_minus_a0   = nuis$censN.pred.1;  # G_C^N | A=0; for AIPWCC
    # g.hats               = 1-nuis$prop.pred;
    # fN_a0                = nuis$fN_0;
    # fN_1_minus_a0        = nuis$fN_1;
    # isotonize            = TRUE;
    # verbose              = TRUE
    
    
    fit.times <- fit.times[fit.times > 0]
    n  <- length(Y)
    K  <- length(fit.times)
    l  <- length(eval.times)
    
    if (verbose) {
      pb <- progress_bar$new(
        total  = n,
        format = "Computing survival [:bar] :percent ETA: :eta"
      )
    }
    
    # Sort fit.times; reorder all K-column matrices consistently
    ord                  <- order(fit.times)
    fit.times            <- fit.times[ord]
    S.hats               <- S.hats[,            ord, drop = FALSE]
    GY.hats_a0           <- GY.hats_a0[,        ord, drop = FALSE]
    GY.hats_1_minus_a0   <- GY.hats_1_minus_a0[,ord, drop = FALSE]
    GN.hats_a0           <- GN.hats_a0[,        ord, drop = FALSE]
    GN.hats_1_minus_a0   <- GN.hats_1_minus_a0[,ord, drop = FALSE]
    
    # Pre-compute grid indices (avoid repeated which() inside loop)
    # Y_ind[i]: index of first fit.time >= Y[i]  (for S_y, inner1)
    Y_ind <- vapply(Y, function(y) {
      idx <- which(fit.times >= y)
      if (length(idx) == 0L) K else idx[1L]
    }, integer(1L))
    
    # N_ind[i]: index of first eval.time >= N[i]  (for N_indx in AIPWCC terms)
    N_ind <- vapply(N, function(ni) {
      idx <- which(eval.times >= ni)
      if (length(idx) == 0L) l else idx[1L]
    }, integer(1L))
    
    # Pre-build a lookup from id → row indices in S.hats
    # (avoids O(n_aug) scan per participant)
    unique_ids   <- unique(ids)
    id_to_rows   <- split(seq_along(ids), ids)   # named list: id → row indices
    
    IF.vals <- matrix(NA_real_, n, K)
    
    for (i in seq_len(n)) {
      
      id_i    <- unique_ids[i]
      N_indx  <- N_ind[i]
      
      # ------------------------------------------------------------------
      # Posterior fN under participant's ACTUAL arm
      # ------------------------------------------------------------------
      fN.i <- if (A[i] == 1L) fN_a0[i, ] else fN_1_minus_a0[i, ]
      
      # ------------------------------------------------------------------
      # S_Y(t | a0, X, N_j) for each eval.times grid point j  →  (l x K)
      # ------------------------------------------------------------------
      row_idx  <- id_to_rows[[as.character(id_i)]]
      S.hats.i <- S.hats[row_idx, , drop = FALSE]   # (l x K)
      
      # ------------------------------------------------------------------
      # G_C^Y step functions for this participant
      # Used in: vals.i (compensator) and G_y (inner1 denominator)
      # ------------------------------------------------------------------
      GYfun     <- stepfun(fit.times, c(1, GY.hats_a0[i,       ]), right = TRUE)
      GYfun_alt <- stepfun(fit.times, c(1, GY.hats_1_minus_a0[i,]), right = TRUE)
      
      # G_C^N step functions for this participant
      # Used in: terms 1, 2, 3 of AIPWCC
      GNfun     <- stepfun(fit.times, c(1, GN.hats_a0[i,       ]), right = TRUE)
      GNfun_alt <- stepfun(fit.times, c(1, GN.hats_1_minus_a0[i,]), right = TRUE)
      
      # ------------------------------------------------------------------
      # phi^F compensator:
      #   vals.i[j, k] = int_0^{min(t_k, Y_i)} dLambda_Y(u|a0,X,N_j) / (S_Y * G_C^Y)
      #
      # Approximated as a cumulative sum over fit.times increments using
      # the product-integral identity: dLambda_Y / S_Y = d(1/S_Y).
      # Divided by G_C^Y(t | a0, X) — the Y-process censoring survival.
      #
      # Note: G_C^Y is always evaluated under arm a0 here because phi^F
      # is defined under a0 for all participants (the IPW factor I(A=a0)/pi
      # zeroes out the correction for participants with A != a0 anyway).
      # ------------------------------------------------------------------
      GY_vec <- GY.hats_a0[i, ]                       # (K) under a0
      
      S_inv  <- 1 / S.hats.i                          # (l x K)
      dS_inv <- S_inv[, -1, drop = FALSE] -
        S_inv[, -K,  drop = FALSE]             # (l x K-1)
      
      # Multiply each column m by 1/G_C^Y(fit.times[m]):
      # GY_vec[-K] gives G at fit.times[1..K-1] (left endpoint of each interval)
      vals.i <- sweep(dS_inv, 2, 1 / GY_vec[-K], FUN = "*")   # (l x K-1)
      vals.i[, fit.times[-1] > Y[i]] <- 0                       # zero beyond Y[i]
      vals.i <- cbind(0, t(apply(vals.i, 1, cumsum)))           # (l x K)
      
      # S_Y(Y_i | a0, X, N_j) for each j;  G_C^Y(Y_i | a0, X) scalar
      S_y <- S.hats.i[, Y_ind[i]]   # (l)
      G_y <- GYfun(Y[i])            # scalar
      
      # ------------------------------------------------------------------
      # G_C^N on eval.times grid — for AIPWCC terms 1/2/3
      #
      # Select the correct arm based on participant's actual treatment.
      # GN_eval[1]   = 1           (sentinel at t=0)
      # GN_eval[j+1] = G_C^N(eval.times[j] | actual arm, X)
      #
      # dLambdaN[j]  = 1 - GN_eval[j+1]/GN_eval[j]
      #              = discrete hazard increment of N-censoring at eval.times[j]
      #
      # GN_inv[j]    = 1/GN_eval[j]
      #              = 1/G_C^N(eval.times[j-1])  = left-limit weight at j
      # ------------------------------------------------------------------
      if (A[i] == 1L) {
        GN_at_Ni  <- GNfun(N[i])                          # scalar for terms 1/2
        GN_eval   <- c(1, GNfun(eval.times))              # length l+1
      } else {
        GN_at_Ni  <- GNfun_alt(N[i])
        GN_eval   <- c(1, GNfun_alt(eval.times))
      }
      
      GN_inv    <- 1 / pmax(GN_eval, 1E-8)                            # length l+1
      dLambdaN  <- 1 - GN_eval[-1] / GN_eval[-(l+1)]     # length l
      
      # ------------------------------------------------------------------
      # Pre-compute fN tail-sum denominator (k-independent)
      # denom_cs[j] = sum_{m >= j} fN.i[m]
      # ------------------------------------------------------------------
      denom_cs <- rev(cumsum(rev(fN.i)))   # (l)
      
      # ------------------------------------------------------------------
      # Inner loop over target time points t0 = fit.times[k]
      # ------------------------------------------------------------------
      for (k in seq_len(K)) {
        
        S_t0 <- S.hats.i[, k]   # S_Y(t0 | a0, X, N_j) for each j  (l)
        
        # inner1: 1/(S_Y(Y_i) * G_C^Y(Y_i)) if event observed at/before t0
        if (Y[i] <= fit.times[k] && DeltaY[i] == 1L) {
          inner1 <- 1 / (S_y * G_y)   # (l)
        } else {
          inner1 <- rep(0, l)
        }
        
        # Full-data EIF phi^F as a function of N grid point j  (l)
        IF_k <- S_t0 - (as.numeric(A[i] == 1L) / g.hats[i]) *
          S_t0 * (inner1 - vals.i[, k])
        
        # ----------------------------------------------------------------
        # Q_seq[j] = E[ phi^F(N) | N >= eval.times[j], obs data, A=A_i ]
        #          = sum_{m>=j} IF_k[m] * fN.i[m]  /  sum_{m>=j} fN.i[m]
        # ----------------------------------------------------------------
        num_cs <- rev(cumsum(rev(IF_k * fN.i)))           # (l)
        Q_seq  <- num_cs / pmax(denom_cs, 1e-12)          # (l)
        
        # ----------------------------------------------------------------
        # Term 1: (delta_N / G_C^N(N_i | actual arm, X)) * phi^F(N_i)
        #
        # For administratively censored N: DeltaN[i] = 0 by construction
        # (event_N == 0 & time_N == admin_t).  With G_C^N now fit on
        # random-censoring-only, G_C^N(admin_t) is bounded away from 0,
        # so 0/G_C^N(admin_t) = 0.  Term 2 then carries the full weight.
        # No special-casing needed.
        # ----------------------------------------------------------------
        term1 <- (DeltaN[i] / GN_at_Ni) * IF_k[N_indx]
        
        # ----------------------------------------------------------------
        # Term 2: ((1 - delta_N) / G_C^N(N_i | actual arm, X)) * Q(N_i)
        # ----------------------------------------------------------------
        term2 <- ((1 - DeltaN[i]) / GN_at_Ni) * Q_seq[N_indx]
        
        # ----------------------------------------------------------------
        # Term 3: sum_{j=1}^{N_indx} Q(t_j) * (1/G_C^N(t_{j-1})) * dLambdaN(t_j)
        #
        # GN_inv[j] = 1/G_C^N(eval.times[j-1]) = left-limit weight at j.
        # GN_inv[1] = 1/G_C^N(0) = 1  (sentinel).
        # For j = 1..N_indx:
        #   left-limit weight = GN_inv[j]    (positions 1..N_indx)
        #   hazard increment  = dLambdaN[j]  (positions 1..N_indx)
        # ----------------------------------------------------------------
        idx3  <- seq_len(N_indx)
        term3 <- sum(Q_seq[idx3] * GN_inv[idx3] * dLambdaN[idx3])
        
        IF.vals[i, k] <- term1 + term2 - term3
        
      }  # end k loop
      
      if (verbose) pb$tick()
      
    }  # end i loop
    
    surv           <- colMeans(IF.vals, na.rm = TRUE)
    IF.vals.center <- sweep(IF.vals, 2, surv, "-")
    
    res <- list(
      times          = fit.times,
      surv           = pmin(1, pmax(0, surv)),
      IF.vals        = IF.vals,
      IF.vals.center = IF.vals.center
    )
    
    if (isotonize) {
      idx              <- !is.na(res$surv)
      res$surv.iso     <- res$surv
      res$surv.iso[idx] <- 1 - isoreg(res$times[idx], 1 - res$surv[idx])$yf
    }
    
    res
  }
  
  DeltaN_for_AIPWCC <- as.integer(event_N == 1 | time_N == admin_t)
  
  if(0 %in% fit.treat) {
    
    if(verbose){message("Computating Survival Probability for A=0")}
    
    surv.0 <- .get.survival.AIPWCC_fast(
      Y                    = time,
      DeltaY               = event,
      N                    = time_N,
      DeltaN               = DeltaN_for_AIPWCC,
      A                    = 1-treat,             # A[i]=1 iff treat[i]==1
      admin_t              = admin_t,
      fit.times            = fit.times,
      eval.times           = nuis$eval.times,
      ids                  = nuis$id,
      S.hats               = nuis$event.pred.0,
      GY.hats_a0           = nuis$cens.pred.0,   # G_C^Y | A=1, for phi^F
      GY.hats_1_minus_a0   = nuis$cens.pred.1,   # G_C^Y | A=0, for phi^F
      GN.hats_a0           = nuis$censN.pred.0,  # G_C^N | A=1, for AIPWCC
      GN.hats_1_minus_a0   = nuis$censN.pred.1,  # G_C^N | A=0, for AIPWCC
      g.hats               = 1-nuis$prop.pred,
      fN_a0                = nuis$fN_0,
      fN_1_minus_a0        = nuis$fN_1,
      isotonize            = TRUE,
      verbose              = TRUE
    )
    
    surv.df.0 <- data.frame(time=unique(c(0,fit.times)), 
                            trt=0, 
                            surv=c(1, surv.0$surv.iso))
    
    result$IF.vals.0 <- surv.0$IF.vals.center
    
    q.01 <- quantile(time[event == 1 & treat == 0], .01) #min(fit.times[surv.0$surv <= .99])
    
    q.99.inds <- !is.na(surv.0$surv.iso)
    q.99.inds[q.99.inds &  surv.0$surv.iso < .01] <- FALSE
    q.99 <- max(fit.times[q.99.inds])
    min.obs.time <- min(time[event == 1 & treat == 0])
    
    c.int <- .surv.confints(fit.times, 
                            surv.0$surv, 
                            surv.0$IF.vals.center, 
                            conf.band = conf.band, 
                            band.end.pts = c(q.01, q.99), 
                            conf.level=conf.level)
    
    surv.df.0$se <- c(0,c.int$res$se)
    surv.df.0$se.logit <- c(0,c.int$res$se.logit)
    surv.df.0$ptwise.lower <- c(1,c.int$res$ptwise.lower)
    surv.df.0$ptwise.upper <- c(1,c.int$res$ptwise.upper)
    surv.df.0$ptwise.logit.lower <- c(1,c.int$res$ptwise.logit.lower)
    surv.df.0$ptwise.logit.upper <- c(1,c.int$res$ptwise.logit.upper)
    
    #surv.df.0$se[surv.df.0$time < min.obs.time & surv.df.0$time > 0] <- surv.df.0$se.logit[surv.df.0$time < min.obs.time & surv.df.0$time > 0] <- surv.df.0$ptwise.lower[surv.df.0$time < min.obs.time & surv.df.0$time > 0] <- surv.df.0$ptwise.upper[surv.df.0$time < min.obs.time & surv.df.0$time > 0] <- NA
    
    if(conf.band) {
      surv.df.0$unif.ew.lower <- c(1,c.int$res$unif.ew.lower)
      surv.df.0$unif.ew.upper <- c(1,c.int$res$unif.ew.upper)
      result$surv.0.unif.ew.quant <- c.int$unif.ew.quant
      result$surv.0.ew.sim.maxes <- c.int$ew.sim.maxes
      
      surv.df.0$unif.logit.lower <- c(1,c.int$res$unif.logit.lower)
      surv.df.0$unif.logit.upper <- c(1,c.int$res$unif.logit.upper)
      result$surv.0.unif.logit.quant <- c.int$unif.logit.quant
      result$surv.0.logit.sim.maxes <- c.int$logit.sim.maxes
    }
    result$surv.df <- rbind(result$surv.df, surv.df.0)
  }
  
  if(1 %in% fit.treat) {
    
    if(verbose){message("Computing Survival Probability for A=1")}
    
    surv.1 <- .get.survival.AIPWCC_fast(
      Y                    = time,
      DeltaY               = event,
      N                    = time_N,
      DeltaN               = DeltaN_for_AIPWCC,
      A                    = treat,             # A[i]=1 iff treat[i]==1
      admin_t              = admin_t,
      fit.times            = fit.times,
      eval.times           = nuis$eval.times,
      ids                  = nuis$id,
      S.hats               = nuis$event.pred.1,
      GY.hats_a0           = nuis$cens.pred.1,   # G_C^Y | A=1, for phi^F
      GY.hats_1_minus_a0   = nuis$cens.pred.0,   # G_C^Y | A=0, for phi^F
      GN.hats_a0           = nuis$censN.pred.1,  # G_C^N | A=1, for AIPWCC
      GN.hats_1_minus_a0   = nuis$censN.pred.0,  # G_C^N | A=0, for AIPWCC
      g.hats               = nuis$prop.pred,
      fN_a0                = nuis$fN_1,
      fN_1_minus_a0        = nuis$fN_0,
      isotonize            = TRUE,
      verbose              = TRUE
    )
    
    surv.df.1 <- data.frame(time=unique(c(0,fit.times)), trt=1, surv=c(1, surv.1$surv.iso))
    result$IF.vals.1 <- surv.1$IF.vals.center
    
    q.01 <- quantile(time[event == 1 & treat == 1], .01)#min(fit.times[surv.1$surv <= .99])
    q.99 <- max(fit.times[surv.1$surv.iso >= .01])
    min.obs.time <- min(time[event == 1 & treat == 1])
    
    c.int <- .surv.confints(fit.times, 
                            surv.1$surv, 
                            surv.1$IF.vals.center, 
                            conf.band = conf.band, 
                            band.end.pts = c(q.01, q.99), 
                            conf.level=conf.level)
    
    surv.df.1$se <- c(0,c.int$res$se)
    surv.df.1$se.logit <- c(0,c.int$res$se.logit)
    surv.df.1$ptwise.lower <- c(1,c.int$res$ptwise.lower)
    surv.df.1$ptwise.upper <- c(1,c.int$res$ptwise.upper)
    surv.df.1$ptwise.logit.lower <- c(1,c.int$res$ptwise.logit.lower)
    surv.df.1$ptwise.logit.upper <- c(1,c.int$res$ptwise.logit.upper)
    
    #surv.df.1$se[surv.df.1$time < min.obs.time & surv.df.1$time > 0] <- surv.df.1$se.logit[surv.df.1$time < min.obs.time & surv.df.1$time > 0] <- surv.df.1$ptwise.lower[surv.df.1$time < min.obs.time & surv.df.1$time > 0] <- surv.df.1$ptwise.upper[surv.df.1$time < min.obs.time & surv.df.1$time > 0] <- NA
    
    if(conf.band) {
      surv.df.1$unif.ew.lower <- c(1,c.int$res$unif.ew.lower)
      surv.df.1$unif.ew.upper <- c(1,c.int$res$unif.ew.upper)
      result$surv.1.unif.ew.quant <- c.int$unif.ew.quant
      result$surv.1.ew.sim.maxes <- c.int$ew.sim.maxes
      
      surv.df.1$unif.logit.lower <- c(1,c.int$res$unif.logit.lower)
      surv.df.1$unif.logit.upper <- c(1,c.int$res$unif.logit.upper)
      result$surv.1.unif.logit.quant <- c.int$unif.logit.quant
      result$surv.1.logit.sim.maxes <- c.int$logit.sim.maxes
    }
    result$surv.df <- rbind(result$surv.df, surv.df.1)
  }
  
  #### ESTIMATE CF CONTRASTS ####
  if(identical(sort(fit.treat), c(0,1))) {
    q.01 <- quantile(time[event == 1], .01)#min(fit.times[surv.0$surv <= .99])
    q.99 <- max(fit.times[surv.0$surv.iso >= .01 | surv.1$surv.iso >= .01])
  }
  
  if('surv.diff' %in% contrasts & identical(sort(fit.treat), c(0,1))) {
    if(verbose) message("Computing survival differences...")
    surv.diff <- .surv.difference(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, band.end.pts=c(q.01, q.99), conf.level=conf.level)
    result <- c(result, surv.diff)
  }
  
  if('surv.ratio' %in% contrasts & identical(sort(fit.treat), c(0,1))) {
    if(verbose) message("Computing survival ratios...")
    surv.ratio <- .surv.ratio(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, band.end.pts=c(q.01, q.99),conf.level=conf.level)
    result <- c(result, surv.ratio)
  }
  
  if('risk.ratio' %in% contrasts  & identical(sort(fit.treat), c(0,1))) {
    if(verbose) message("Computing risk ratios...")
    risk.ratio <- .risk.ratio(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, band.end.pts=c(q.01, q.99), conf.level=conf.level)
    result <- c(result, risk.ratio)
  }
  
  if('nnt' %in% contrasts  & identical(sort(fit.treat), c(0,1))) {
    if(verbose) message("Computing number needed to treat...")
    nnt <- .nnt(fit.times=fit.times, surv.0=surv.0$surv, IF.vals.0 = surv.0$IF.vals, surv.1=surv.1$surv, IF.vals.1 = surv.1$IF.vals, conf.band=conf.band, band.end.pts=c(q.01, q.99), conf.level=conf.level)
    result <- c(result, nnt)
  }
  
  result$data <- data.frame(time, event, treat)
  result$nuisance <- nuis
  
  return(result)
}
