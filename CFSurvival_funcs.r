# Key functions adapted from the CFSurvival Package

#' Estimate counterfactual survival functions while adjusting for a negative control right-censored predictor
#'
#' This function estimates counterfactual survival functions and contrasts from right-censored data subject to potential confounding, including adjustment for a negative control time-to-event predictor.
#'
#' For a binary treatment measured at baseline, the counterfactual survival through time \code{t} for treatment level \code{a} is defined as the probability of survival up to time \code{t} were all units in the population assigned to treatment \code{a}. This quantity is estimated using an augmented inverse probability of treatment/censoring weighted (AIPTW) estimator. Our estimator of the counterfactual survival function requires further estimation of FOUR nuisance parameters: the conditional probability of receiving treatment given a set of potential confounders (i.e. the treatment propensity), (NEW) the conditional survival function of the negative control event given a set of potential confouners, the conditional survival function of the event given confounders and the negative control event time, and the conditional survival function of censoring given confounders. At the moment, we provide two ways of estimating each of these nuisance parameters. The treatment propensity can be estimated either using a logistic regression or SuperLearner (van der Laan and Polley, 2006). The conditional survival functions of negative control event, primary event, and censoring can be estimated using single estimators or a SuperLearner ensemble. Additionally, the estimates of any of these parameters can be passed in directly if another method is desired.
#'
#'  \code{CFsurvival} by default returns 95\% pointwise confidence intervals and 95\% uniform confidence band, both of which are based on the influence function of the estimator. Contrasts of the counterfactual survival under treatment and control can also be obtained, including the survival difference, survival ratio, risk ratio, and number needed to treat. Pointwise confidence intervals and uniform confidence bands for these contrasts will also be computed.
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

CFsurvival3 <- function(time, 
                        event,
                        time_N,
                        event_N,
                        information_level,
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
                        verbose=TRUE,
                        maxiter = 3) {
  # time = to_fit_1$Y;
  # event = to_fit_1$DeltaY;
  # treat = to_fit_1$A;
  # W = W;
  # time_N = to_fit_1$N1;
  # event_N = to_fit_1$DeltaN1;
  # information_level = 2;
  # admin_t = to_fit_1$admin_t;
  # id = to_fit_1$ptid;
  # fit.times = sort(unique(to_fit_1$Y[to_fit_1$Y>0])); 
  # fit.treat = c(0,1);
  # nuisance.options = CFsurvival.nuisance.options(
  #   cross.fit = TRUE,
  #   V = 5,
  #   folds = NULL,
  #   eval.times = NULL,
  #   event.SL.library = lapply(c("survSL.km", 
  #                               "survSL.coxph", 
  #                               "survSL.gam"), function(alg)
  #                                 c(alg, "survscreen.glmnet", "All")),
  #   event.pred.0 = NULL,
  #   event.pred.1 = NULL,
  #   cens.SL.library = lapply(c("survSL.km", 
  #                              "survSL.coxph", 
  #                              "survSL.gam"), function(alg)
  #                                c(alg, "survscreen.glmnet", "All")),
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
  # );
  # conf.band = FALSE;
  # conf.level = 0.95;
  # contrasts = c("risk.ratio");
  # verbose = TRUE
  
  if(!(information_level %in% c(0,1,2,3))){
    stop("information_level must be one of 0, 1, 2, or 3")
  }
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
  #if(any(fit.times == 0)) fit.times <- fit.times[fit.times > 0]
  #if(any(fit.times > max(time[event == 1]))) {
  #  fit.times <- fit.times[fit.times <= max(time[event == 1])]
  #  message("fit.times > max(time[event == 1]) removed.")
  #}
  
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
        train <- nuis$fold != v
        test <- nuis$fold == v
        
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
    nuis$w <- vector(length=0)
    nuis$id <- vector(length=0)
    if(0 %in% fit.treat & is.null(nuis$event.pred.0)) {
      do.event.pred.0 <- TRUE
      nuis$event.pred.0 <- matrix(NA, nrow=0, ncol=k)
    } else do.event.pred.0 <- FALSE
    if(1 %in% fit.treat & is.null(nuis$event.pred.1)) {
      do.event.pred.1 <- TRUE
      nuis$event.pred.1 <- matrix(NA, nrow=0, ncol=k)
    } else do.event.pred.1 <- FALSE
    if(0 %in% fit.treat & is.null(nuis$cens.pred.0)) {
      do.cens.pred.0 <- TRUE
      nuis$cens.pred.0 <- matrix(NA, nrow=0, ncol=k)
    } else do.cens.pred.0 <- FALSE
    if(1 %in% fit.treat & is.null(nuis$cens.pred.1)) {
      do.cens.pred.1 <- TRUE
      nuis$cens.pred.1 <- matrix(NA, nrow=0, ncol=k)
    } else{do.cens.pred.1 <- FALSE}
    if(nuis$V > 1) {
      if(nuis$save.nuis.fits) result$surv.fits <- vector(mode='list',length=nuis$V)
      nuis$event.coef <- nuis$cens.coef <- list()
      for(v in 1:nuis$V) {
        if(verbose) message(paste("Fold ", v, "..."))
        train <- which(nuis$fold != v)
        test <- which(nuis$fold == v)
        
        ### All nuisance parameters must be fit to the training fold
        dat_train<- as.data.frame(cbind("id"= id, 
                                        "Y"=time, 
                                        "DeltaY"=event, 
                                        "N"=time_N, 
                                        'DeltaN'=event_N, 
                                        "A"=treat, 
                                        'admin_t' = admin_t, 
                                        W))[train,]
        
        dat_test<- as.data.frame(cbind("id" = id, 
                                       "Y"=time, 
                                       "DeltaY"=event, 
                                       "N"=time_N, 
                                       'DeltaN'=event_N, 
                                       "A"=treat, 
                                       'admin_t' = admin_t, 
                                       W))[test,]
        
        augment_dataset <- function(i, df, id, A, Y, DeltaY, W, N, DeltaN, admin_t, FN, information_level){
          
          survest_1 <- FN$survest_1
          survest_0 <- FN$survest_0
          
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
          if(information_level == 2){
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
          }else{
            if(Ai==0){
              S_vals <- as.vector(survest_0(
                c(N0, N_grid)
              ))
            }else{
              S_vals <- as.vector(survest_1(
                c(N0, N_grid)
              ))
            }
          }
          
          # compute conditional probabilities
          S_N0 <- S_vals[1]
          S_grid <- S_vals
          
          # P(N = n | N > N0)
          w_raw <- -diff(c(S_grid)) / S_N0
          w_raw[w_raw < 0] <- 0
          w <- w_raw / sum(w_raw)
          
          if(all(is.na(W) | is.null(W))){
            data.frame(
              id     = row_i[[id]],
              N      = N_grid,
              DeltaN = as.numeric(N_grid < admin_t_tmp),
              Y      = row_i[[Y]],
              DeltaY = row_i[[DeltaY]],
              A      = row_i[[A]],
              X      = NA,
              admin_t = admin_t_tmp,
              w      = w
            ) %>% filter(w>0)
          }else{
            suppressMessages(bind_cols(data.frame(
              id     = row_i[[id]],
              N      = N_grid,
              DeltaN = as.numeric(N_grid < admin_t_tmp),
              Y      = row_i[[Y]],
              DeltaY = row_i[[DeltaY]],
              A      = row_i[[A]],
              admin_t = admin_t_tmp,
              w      = w
            ) %>% filter(w>0), row_i %>% select(all_of(W)), .name_repair="minimal"))
          }
        }
        
        if(information_level==0){
          if(verbose){
            message("Level 0 -- fitting (unconditional) censored covariate distribution N")
            warning("Global KM estimator relies on completely at random censoring. Level 1 may be preferable from standpoint of assumptions")
          }
          
        dat_train=dat_train %>%
            mutate(
              DeltaN = case_when(
                N >= admin_t ~ 1,
                .default = DeltaN),
              N = case_when(
                N > admin_t ~ admin_t,
                .default = N)
            )
          
        KM<-survfit(Surv(N, DeltaN)~1, data=dat_train)
        
        FN<-list()
        FN$survest_1 <- stepfun(KM$time, c(1, KM$surv))
        FN$survest_0 <- stepfun(KM$time, c(1, KM$surv))
        
        }else if(information_level==1){
          if(verbose){
            message("Level 1 -- fitting (treatment-conditional) censored covariate distribution N | A")
            message("Using treatment-stratified KM estimator")
          }
          
          dat_train=dat_train %>%
            mutate(
              DeltaN = case_when(
                N >= admin_t ~ 1,
                .default = DeltaN),
              N = case_when(
                N > admin_t ~ admin_t,
                .default = N)
            )
          
        KM<-survfit(Surv(N, DeltaN)~A, data=dat_train)
        
        FN<-list()
        FN$survest_0 <- stepfun(KM[1]$time, c(1, KM[1]$surv))
        FN$survest_1 <- stepfun(KM[2]$time, c(1, KM[2]$surv))
        
        }else if(information_level==2){
          if(verbose){
            message("Level 2 -- fitting (treatment and covariate-conditional) censored covariate distribution N | A, X")
            message("Using survSuperLearner ensemble")
          }
          
          FN<-estimate_conditional_survival_N(dat=dat_train, 
                                              N_name="N", 
                                              DeltaN_name = "DeltaN", 
                                              A_name="A", 
                                              X_names = colnames(W), 
                                              admin_t = admin_t, 
                                              event.SL.library = nuis$event.SL.library,
                                              cens.SL.library = nuis$cens.SL.library,
                                              survSL.control = nuis$survSL.control,
                                              survSL.cvControl = nuis$survSL.cvControl,
                                              verbose=FALSE)
        }else{
          if(verbose){
            message("\tLevel 3 -- fitting (full observed-data-conditional) censored covariate distribution N | A, X, Delta_Y, tilde{Y}")
            message("\tNOTE: This approach will require iterative nuisance fitting")
            message("\tStep 1: estimate prior model -- fitting (treatment and covariate-conditional) censored covariate distribution N | A, X")
          }
          
          FN<-estimate_conditional_survival_N(dat=dat_train, 
                                              N_name="N", 
                                              DeltaN_name = "DeltaN", 
                                              A_name="A", 
                                              X_names = colnames(W), 
                                              admin_t = admin_t, 
                                              event.SL.library = nuis$event.SL.library,
                                              cens.SL.library = nuis$cens.SL.library,
                                              survSL.control = nuis$survSL.control,
                                              survSL.cvControl = nuis$survSL.cvControl,
                                              verbose=FALSE)
          
        }
        
        if(information_level %in% c(0,1,2)){
          message(paste0("Information Level=", information_level, ". Conducting dataset augmentation."))
          if(nuis$verbose){
          message("Augmenting training dataset using redistribute-to-the-right via N|X")
          pb <- progress_bar$new(total = nrow(dat_train),
                                 format = "[:bar] :percent ETA: :eta")
        }
          aug_train_list <- vector(mode="list", length=nrow(dat_train))
          for(i in 1:nrow(dat_train)){
          aug_train_list[[i]]<-augment_dataset(i, 
                                               df=dat_train, 
                                               id="id",
                                               Y="Y", 
                                               DeltaY='DeltaY', 
                                               A="A", 
                                               N="N", 
                                               W=colnames(W), 
                                               DeltaN='DeltaN', 
                                               admin_t='admin_t',
                                               FN = FN,
                                               information_level = information_level)
          if(nuis$verbose){pb$tick()}
        }
        
          ### Create augmented training data
          aug_train <- do.call(rbind, aug_train_list)
          rm(aug_train_list)
        
          if(nuis$verbose){
          message("Augmenting test dataset using redistribute-to-the-right via N|X")
          pb <- progress_bar$new(total = nrow(dat_test),
                                 format = "[:bar] :percent ETA: :eta")}
        
          aug_test_list <- vector(mode="list", length=nrow(dat_test))
          for(i in 1:nrow(dat_test)){
          aug_test_list[[i]]<-augment_dataset(i, 
                                              df=dat_test, 
                                              id = 'id',
                                              Y="Y", 
                                              DeltaY='DeltaY', 
                                              A="A", 
                                              N="N", 
                                              W=colnames(W), 
                                              DeltaN='DeltaN', 
                                              admin_t='admin_t',
                                              FN = FN,
                                              information_level = information_level)
          if(nuis$verbose){pb$tick()}
        }
        
          ### Create augmented training data
          aug_test <- do.call(rbind, aug_test_list)
          rm(aug_test_list)
        
          message("Fit conditional survival distributions Y|N,A,X")
          surv.fit <- .estimate.conditional.survival.2(Y=aug_train[["Y"]], 
                                                     Delta=aug_train[["DeltaY"]], 
                                                     A=aug_train[["A"]], 
                                                     W_event=aug_train %>% dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                     W_cens=aug_train %>% dplyr::select(colnames(W)),
                                                     newW_event = aug_test %>% dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                     newW_cens = aug_test %>% dplyr::select(colnames(W)),
                                                     event.SL.library=nuis$event.SL.library, 
                                                     fit.times=fit.times, 
                                                     fit.treat=fit.treat, 
                                                     cens.SL.library=nuis$cens.SL.library, 
                                                     survSL.control=nuis$survSL.control, 
                                                     survSL.cvControl = nuis$survSL.cvControl, 
                                                     cens.trunc=nuis$cens.trunc, 
                                                     save.fit = nuis$save.nuis.fits, 
                                                     verbose = nuis$verbose,
                                                     obsWeights = aug_train$w, 
                                                     id = aug_train$id)
        }
        
        if(information_level == 3){
          
          augment_dataset_unweighted <- function(i, 
                                                 df, 
                                                 N_support, 
                                                 id, 
                                                 A, 
                                                 Y, 
                                                 DeltaY,
                                                 W, 
                                                 N, 
                                                 DeltaN, 
                                                 admin_t){
            
            row_i <- df[i, ]
            
            Ai = row_i[[A]]
            N0 <- row_i[[N]]  # censoring time
            admin_t_tmp <- row_i[[admin_t]]
            
            #N_support <- c(sort(unique(df[[N]][df[[DeltaN]]==1])), admin_t_tmp)
            
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
                    admin_t = admin_t_tmp
                  ))
              }else{
                return(
                  bind_cols(data.frame(
                    id     = row_i[[id]],
                    N      = row_i[[N]],
                    DeltaN = row_i[[DeltaN]],
                    Y      = row_i[[Y]],
                    DeltaY = row_i[[DeltaY]],
                    A      = row_i[[A]],
                    admin_t = admin_t_tmp
                  ), row_i %>% select(all_of(W))))
              }
            }
            
            # Case 2: N censored → redistribution
            
            # plausible N values
            N_grid <- N_support[N_support > N0]
            
            if (length(N_grid) == 0){
              return(bind_cols(data.frame(
                id     = row_i[[id]],
                N      = admin_t_tmp,
                DeltaN = 0,
                Y      = row_i[[Y]],
                DeltaY = row_i[[DeltaY]],
                A      = row_i[[A]],
                admin_t = admin_t_tmp
              ), row_i %>% select(all_of(W))))
            }else{
              return(bind_cols(data.frame(
                id     = row_i[[id]],
                N      = N_grid,
                DeltaN = c(rep(1, length(N_grid)-1), 0),
                Y      = row_i[[Y]],
                DeltaY = row_i[[DeltaY]],
                A      = row_i[[A]],
                admin_t = admin_t_tmp
              ), row_i[[N]]))
            }
          }
          
          # Augment training data without weights
          aug_train<-lapply(1:nrow(dat_train), FUN=function(x){
            augment_dataset_unweighted(i=x, 
                                       df=dat_train, 
                                       N_support = sort(unique(c(time_N[event_N==1 & time_N <= admin_t[train[x]]], admin_t[train[x]]))),
                                       id='id', 
                                       A="A", 
                                       Y="Y", 
                                       DeltaY='DeltaY', 
                                       W=colnames(W), 
                                       N="N", 
                                       DeltaN='DeltaN', 
                                       admin_t='admin_t')
          })
          aug_train <- do.call(rbind, aug_train)
          
          # Augment test data without weights
          aug_test<-lapply(1:nrow(dat_test), FUN=function(x){
            augment_dataset_unweighted(i=x, 
                                       df=dat_test, 
                                       N_support = sort(unique(c(time_N[event_N==1 & time_N <= admin_t[test[x]]], admin_t[test[x]]))),
                                       id='id', 
                                       A="A", 
                                       Y="Y", 
                                       DeltaY='DeltaY', 
                                       W=colnames(W), 
                                       N="N", 
                                       DeltaN='DeltaN', 
                                       admin_t='admin_t')
          })
          aug_test <- do.call(rbind, aug_test)
          
          message(paste0("Aug train dimension (", nrow(aug_train), "x", ncol(aug_train), ")"))
          message(paste0("Aug test dimension (", nrow(aug_test), "x", ncol(aug_test), ")"))
          
          #Augment test data without weights
          
          message("\tFit preliminary conditional survival distributions Y|N,A,X with DeltaN==1")
          
          indx_full <- which(dat_train$DeltaN==1 | dat_train$N==admin_t[train])
          # time grid for weight calculation
          grid <- round(seq(min(unique(dat_train[["Y"]])), 
                            max(unique(dat_train[["Y"]])), 
                            length=200), 
                        0)
          
          surv.fit <- .estimate.conditional.survival.2(Y=dat_train[["Y"]][indx_full], 
                                                       Delta=dat_train[["DeltaY"]][indx_full], 
                                                       A=dat_train[["A"]][indx_full], 
                                                       W_event=dat_train %>% 
                                                         dplyr::select(c("N", "DeltaN", colnames(W))) %>%
                                                         slice(indx_full), 
                                                       W_cens=dat_train %>% 
                                                         dplyr::select(colnames(W)) %>% 
                                                         slice(indx_full),
                                                       newW_event = aug_train %>% 
                                                         dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                       newW_cens = aug_train %>% 
                                                         dplyr::select(colnames(W)),
                                                       event.SL.library=nuis$event.SL.library, 
                                                       fit.times=grid, 
                                                       fit.treat=fit.treat, 
                                                       cens.SL.library=c('survSL.km'), 
                                                       survSL.control=nuis$survSL.control, 
                                                       survSL.cvControl = nuis$survSL.cvControl, 
                                                       cens.trunc=nuis$cens.trunc, 
                                                       save.fit = nuis$save.nuis.fits, 
                                                       verbose = FALSE, 
                                                       obsWeights = rep(1, length(indx_full)),
                                                       id = 1:length(indx_full))
        
        message("\tIterative nuisance updating -- (i) conditional distn Y|N,A,X and (ii) weights")
        calc_weights_old <- function(i, aug_df, df){
            
            indx <- which(aug_df$id==i)
            slice <- aug_df %>% slice(indx)
          
            if(nrow(slice)==0 | nrow(slice)==1){
              
              return(slice %>% mutate(w=1))
              
            }else{
              
              if(all(slice[["A"]]==1)){
                SN <- FN$survest_1(t=c(df[["N"]][df$id==i], slice$N), 
                                   x = slice %>% 
                                     dplyr::select("X") %>% 
                                     distinct()) #Conditional survival function predictions
                PN = -diff(SN) # Failure probabilities at each slice
                
                #what if the value of Y is not in fit.times?
                # Essential that all times >0 are included in fit.times
                
                #surv.fit$event.pred.1[indx,]
                
                #fit <- stepfun(grid, c(1, surv.fit$event.pred.1[1,]))
                Yind <- which(surv.fit$fit.times==unique(slice$Y))
                
                if(is.na(Yind) | length(Yind)==0){
                  w<-1/nrow(slice)
                }else if(unique(slice$DeltaY)==1){
                  SY <- surv.fit$event.pred.1[indx,(Yind-1):Yind] #row corresponds to id, columns correspond to (T-1, T), hence, P(T=t) = P(T > t-1) - P(T > t) = S(t-1) - S(t)
                  LY<-apply(SY, MARGIN=1, FUN=function(x){x[1]-x[2]}) # This is the likelihood contribution when can event occurs
                  w <- LY * PN / sum(LY * PN)
                }else{
                  LY <- surv.fit$event.pred.1[indx,Yind]  #This is the likelihood contribution under censoring
                  w <- LY * PN / sum(LY * PN)
                }
                return(slice %>% mutate(w = w))
              }else{
                SN <- FN$survest_0(t=c(df[['N']][df$id==i], slice$N), 
                                   x = slice %>% 
                                     dplyr::select("X") %>% 
                                     distinct())
                PN = -diff(SN)
                
                Yind <- which(surv.fit$fit.times==unique(slice$Y))
                
                #what if the value of Y is not in fit.times?
                # Essential that all times >0 are included in fit.times
                Yind <- which(surv.fit$fit.times==unique(slice$Y))
                
                if(is.na(Yind) | length(Yind)==0){
                  w<-1/nrow(slice)
                }else if(unique(slice$DeltaY)==1){
                  SY <- surv.fit$event.pred.0[indx,(Yind-1):Yind] #row corresponds to id, columns correspond to (T-1, T), hence, P(T=t) = P(T > t-1) - P(T > t) = S(t-1) - S(t)
                  LY<-apply(SY, MARGIN=1, FUN=function(x){x[1]-x[2]}) # This is the likelihood contribution when can event occurs
                  w <- LY * PN / sum(LY * PN)
                }else{
                  LY <- surv.fit$event.pred.0[indx,Yind]  #This is the likelihood contribution under censoring
                  w <- LY * PN / sum(LY * PN)
                }
                
                w <- LY * PN / sum(LY * PN)
                
                return(slice %>% mutate(w = w))
                
              }
              
            }
            
        }
        
        calc_weights <- function(i, aug_df, df, grid){
          
          # i=26
          # aug_df = aug_test
          # df = dat_test
          # grid=grid
          
          indx <- which(aug_df$id==i)
          slice <- aug_df %>% slice(indx)
          
          if(nrow(slice)==0 | nrow(slice)==1){
            
            return(slice %>% mutate(w=1))
            
          }else{
            
            if(all(slice[["A"]]==1)){
              
              SN <- FN$survest_1(t=c(df[["N"]][df$id==i], slice$N), 
                                 x = slice %>% 
                                   dplyr::select("X") %>% 
                                   distinct()) #Conditional survival function predictions
              PN = -diff(SN) # Failure probabilities at each slice
              
              Y_tmp <- unique(slice[["Y"]])
              Delta_tmp <- unique(slice[["DeltaY"]])
              LY <- rep(NA, length(indx))
              for(i in 1:length(indx)){
                
                fit <- stepfun(grid, c(1, surv.fit$event.pred.1[indx[i],]))
                if(Delta_tmp==1){
                  LY[i] = fit(Y_tmp) - fit(Y_tmp-1)
                }else{
                  LY[i] = fit(Y_tmp)
                }
              }
              
              w <- LY * PN / sum(LY * PN)
              return(slice %>% mutate(w = w))
            }else{
              SN <- FN$survest_0(t=c(df[["N"]][df$id==i], slice$N), 
                                 x = slice %>% 
                                   dplyr::select("X") %>% 
                                   distinct()) #Conditional survival function predictions
              PN = -diff(SN) # Failure probabilities at each slice
              
              Y_tmp <- unique(slice[["Y"]])
              Delta_tmp <- unique(slice[["DeltaY"]])
              LY <- rep(NA, length(indx))
              for(i in 1:length(indx)){
                
                fit <- stepfun(grid, c(1, surv.fit$event.pred.0[indx[i],]))
                if(Delta_tmp==1){
                  LY[i] = fit(Y_tmp) - fit(Y_tmp-1)
                }else{
                  LY[i] = fit(Y_tmp)
                }
              }
              
              w <- LY * PN / sum(LY * PN)
              return(slice %>% mutate(w = w))
              
            }
            
          }
          
        }
        
        for(j in 1:maxiter){
          message(paste0("\t--Iteration ", j))
        
          message(paste0("\tWeight estimation"))
          out<-vector(mode="list", length=length(dat_train$id))
          iter=0
          for(i in dat_train$id){
            iter=iter+1
            out[[iter]]<-calc_weights(i=i, aug_df = aug_train, df = dat_train, grid=grid)
          }
        
          aug_train <- do.call(rbind, out) %>%
            filter(w>0)
          
          rm(out)
          gc()
          rm(surv.fit)
        
          message(paste0("\tUpdating distn Y|N,A,X via weighted estimation"))
          surv.fit <- .estimate.conditional.survival.2(Y=aug_train[["Y"]], 
                                                     Delta=aug_train[["DeltaY"]], 
                                                     A=aug_train[["A"]], 
                                                     W_event=aug_train %>% 
                                                       dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                     W_cens=aug_train %>% 
                                                       dplyr::select(colnames(W)),
                                                     newW_event = aug_train %>% 
                                                       dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                     newW_cens = aug_train %>% 
                                                       dplyr::select(colnames(W)),
                                                     fit.times=grid, 
                                                     fit.treat=fit.treat, 
                                                     event.SL.library=nuis$event.SL.library, 
                                                     cens.SL.library=c("survSL.km"), 
                                                     survSL.control=nuis$survSL.control, 
                                                     survSL.cvControl = nuis$survSL.cvControl, 
                                                     cens.trunc=nuis$cens.trunc, 
                                                     save.fit = nuis$save.nuis.fits, 
                                                     verbose = FALSE, 
                                                     obsWeights = aug_train$w,
                                                     id = aug_train$id)
        gc()
        }
        
        message("\tFinal calculation on test data")
        surv.fit <- .estimate.conditional.survival.2(Y=aug_train[["Y"]], 
                                                     Delta=aug_train[["DeltaY"]], 
                                                     A=aug_train[["A"]], 
                                                     W_event=aug_train %>% 
                                                       dplyr::select(c("N", 
                                                                       "DeltaN", 
                                                                       colnames(W))), 
                                                     W_cens=aug_train %>% 
                                                       dplyr::select(colnames(W)),
                                                     newW_event = aug_test %>% 
                                                       dplyr::select(c("N", 
                                                                       "DeltaN", 
                                                                       colnames(W))), 
                                                     newW_cens = aug_test %>% 
                                                       dplyr::select(colnames(W)),
                                                     event.SL.library=nuis$event.SL.library, 
                                                     fit.times=fit.times, 
                                                     fit.treat=fit.treat, 
                                                     cens.SL.library=nuis$cens.SL.library, 
                                                     survSL.control=nuis$survSL.control, 
                                                     survSL.cvControl = nuis$survSL.cvControl, 
                                                     cens.trunc=nuis$cens.trunc, 
                                                     save.fit = nuis$save.nuis.fits, 
                                                     verbose = nuis$verbose, 
                                                     obsWeights = aug_train$w,
                                                     id = aug_train$id)
        
        out<-vector(mode="list", length=length(dat_test$id))
        iter=0
        for(i in dat_test$id){
          iter=iter+1
          out[[iter]]<-calc_weights(i=i, aug_df = aug_test, df = dat_test, grid=fit.times)
        }
        aug_test <- do.call(rbind, out)
        rm(out)
        gc()
        
        }
          
        if(do.event.pred.0) nuis$event.pred.0 <- rbind(nuis$event.pred.0, surv.fit$event.pred.0)
        if(do.event.pred.1) nuis$event.pred.1 <- rbind(nuis$event.pred.1, surv.fit$event.pred.1)
        if(do.cens.pred.0) nuis$cens.pred.0 <- rbind(nuis$cens.pred.0, surv.fit$cens.pred.0)
        if(do.cens.pred.1) nuis$cens.pred.1 <- rbind(nuis$cens.pred.1, surv.fit$cens.pred.1)
        nuis$w <- c(nuis$w, aug_test$w)
        nuis$id <- c(nuis$id, aug_test$id)
        if(nuis$save.nuis.fits) result$surv.fits[[v]] <- surv.fit$surv.fit
        nuis$event.coef[[v]] <- surv.fit$event.coef
        nuis$cens.coef[[v]] <- surv.fit$cens.coef
        
        rm(surv.fit)
        gc()
        if(verbose){message(paste("Fold ", v, "Complete"))}
      }
    } else {
      
      if(nuis$save.nuis.fits) result$surv.fits <- vector(mode='list',length=nuis$V)
      nuis$event.coef <- nuis$cens.coef <- list()
      
      ### We MUST fit nuisance parameters on the training fold
      ### This code fits N | X to training data and outputs augmented training dataset
      dat<- as.data.frame(cbind("id"= id, 
                                "Y"=time, 
                                "DeltaY"=event, 
                                "N"=time_N, 
                                'DeltaN'=event_N, 
                                "A"=treat, 
                                'admin_t'=admin_t,
                                W))
      
      if(nuis$verbose){message("Fitting censored covariate distribution N|X")}
      
      FN<-estimate_conditional_survival_N(dat=dat, 
                                          N_name="N", 
                                          DeltaN_name = "DeltaN", 
                                          A_name="A", 
                                          X_names = colnames(W), 
                                          admin_t = admin_t, 
                                          event.SL.library = nuis$event.SL.library,
                                          cens.SL.library = nuis$cens.SL.library,
                                          survSL.control = nuis$survSL.control,
                                          survSL.cvControl = nuis$survSL.cvControl,
                                          verbose=FALSE)
      
      augment_dataset <- function(i, df, id, A, Y, DeltaY, W, N, DeltaN, admin_t, FN){
        
        survest_1 <- FN$survest_1
        survest_0 <- FN$survest_0
        
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
              bind_cols(data.frame(
                id     = row_i[[id]],
                N      = row_i[[N]],
                DeltaN = row_i[[DeltaN]],
                Y      = row_i[[Y]],
                DeltaY = row_i[[DeltaY]],
                A      = row_i[[A]],
                admin_t = admin_t_tmp,
                w      = 1
              ), row_i %>% select(all_of(W)))
            )
          }
        }
        
        # Case 2: N censored → redistribution
        
        # plausible N values
        N_grid <- N_support[N_support > N0]
        
        if (length(N_grid) == 0){
          return(bind_cols(data.frame(
            id     = row_i[[id]],
            N      = admin_t_tmp,
            DeltaN = 0,
            Y      = row_i[[Y]],
            DeltaY = row_i[[DeltaY]],
            A      = row_i[[A]],
            admin_t = admin_t_tmp,
            w      = 1
          ), row_i %>% select(all_of(W))))
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
          data.frame(
            id     = row_i[[id]],
            N      = N_grid,
            DeltaN = as.numeric(N_grid < admin_t_tmp),
            Y      = row_i[[Y]],
            DeltaY = row_i[[DeltaY]],
            A      = row_i[[A]],
            X      = NA,
            admin_t = admin_t_tmp,
            w      = w
          ) %>% filter(w>0)
        }else{
          bind_cols(data.frame(
            id     = row_i[[id]],
            N      = N_grid,
            DeltaN = as.numeric(N_grid < admin_t_tmp),
            Y      = row_i[[Y]],
            DeltaY = row_i[[DeltaY]],
            A      = row_i[[A]],
            admin_t = admin_t_tmp,
            w      = w
          ) %>% filter(w>0), row_i %>% select(all_of(W)))
        }
      }
      
      if(nuis$verbose){
        message("Augmenting training dataset using redistribute-to-the-right via N|X")
        pb <- progress_bar$new(total = nrow(dat),
                               format = "[:bar] :percent ETA: :eta")
      }
      aug_list <- vector(mode="list", length=nrow(dat))
      for(i in 1:nrow(dat)){
        aug_train_list[[i]]<-augment_dataset(i, 
                                             df=dat, 
                                             id="id",
                                             Y="Y", 
                                             DeltaY='DeltaY', 
                                             A="A", 
                                             N="N", 
                                             W='X', 
                                             DeltaN='DeltaN', 
                                             admin_t='admin_t',
                                             FN = FN)
        if(nuis$verbose){pb$tick()}
      }
      
      ### Create augmented training data
      aug <- do.call(rbind, aug_list)
      rm(aug_list)
      
      message("Fit conditional survival distributions Y|N,X")
      surv.fit <- .estimate.conditional.survival.2(Y=aug[["Y"]], 
                                                   Delta=aug[["DeltaY"]], 
                                                   A=aug[["A"]], 
                                                   W_event=aug %>% dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                   W_cens=aug %>% dplyr::select(colnames(W)),
                                                   newW_event = aug %>% dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                   newW_cens = aug %>% dplyr::select(colnames(W)),
                                                   event.SL.library=nuis$event.SL.library, 
                                                   fit.times=nuis$eval.times, 
                                                   fit.treat=fit.treat, 
                                                   cens.SL.library=nuis$cens.SL.library, 
                                                   survSL.control=nuis$survSL.control, 
                                                   survSL.cvControl = nuis$survSL.cvControl, 
                                                   cens.trunc=nuis$cens.trunc, 
                                                   save.fit = nuis$save.nuis.fits, 
                                                   verbose = nuis$verbose,
                                                   obsWeights = aug$w, 
                                                   id = aug$id)
      
      if(do.event.pred.0) nuis$event.pred.0 <- surv.fit$event.pred.0
      if(do.event.pred.1) nuis$event.pred.1 <- surv.fit$event.pred.1
      if(do.cens.pred.0) nuis$cens.pred.0 <- surv.fit$cens.pred.0
      if(do.cens.pred.1) nuis$cens.pred.1 <- surv.fit$cens.pred.1
      nuis$w <- c(nuis$w, aug_test$w)
      nuis$id <- c(nuis$id, aug_test$id)
      if(nuis$save.nuis.fits) result$surv.fits <- surv.fit$surv.fit
      nuis$event.coef <- surv.fit$event.coef
      nuis$cens.coef <- surv.fit$cens.coef
    }
  }
  
  #### ESTIMATE CF SURVIVALS ####
  if(verbose) message("Computing counterfactual survivals...")
  if(0 %in% fit.treat) {
    
    if(verbose){message("Computating Survival Probability for A=0")}
    
    surv.0 <- .get.survival.2(Y=time[nuis$id], 
                              Delta=event[nuis$id], 
                              A=1-treat[nuis$id], 
                              fit.times=fit.times, 
                              eval.times=nuis$eval.times, 
                              S.hats=nuis$event.pred.0, 
                              G.hats=nuis$cens.pred.0, 
                              g.hats=1-nuis$prop.pred[nuis$id], 
                              obsWeights = nuis$w, 
                              id = nuis$id, 
                              isotonize = TRUE)
    
    surv.df.0 <- data.frame(time=unique(c(0,fit.times)), trt=0, surv=c(1, surv.0$surv.iso))
    result$IF.vals.0 <- surv.0$IF.vals
    
    q.01 <- quantile(time[event == 1 & treat == 0], .01)#min(fit.times[surv.0$surv <= .99])
    
    q.99.inds <- !is.na(surv.0$surv.iso)
    q.99.inds[q.99.inds &  surv.0$surv.iso < .01] <- FALSE
    q.99 <- max(fit.times[q.99.inds])
    min.obs.time <- min(time[event == 1 & treat == 0])
    
    c.int <- .surv.confints(fit.times, 
                            surv.0$surv, 
                            surv.0$IF.vals, 
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
    
    surv.1 <- .get.survival.2(Y=time[nuis$id], 
                              Delta=event[nuis$id], 
                              A=treat[nuis$id], 
                              fit.times=fit.times, 
                              eval.times=nuis$eval.times, 
                              S.hats=nuis$event.pred.1, 
                              G.hats=nuis$cens.pred.1, 
                              g.hats=nuis$prop.pred[nuis$id], 
                              obsWeights = nuis$w, 
                              id = nuis$id, 
                              isotonize = TRUE)
    
    surv.df.1 <- data.frame(time=unique(c(0,fit.times)), trt=1, surv=c(1, surv.1$surv.iso))
    result$IF.vals.1 <- surv.1$IF.vals
    
    q.01 <- quantile(time[event == 1 & treat == 1], .01)#min(fit.times[surv.1$surv <= .99])
    q.99 <- max(fit.times[surv.1$surv.iso >= .01])
    min.obs.time <- min(time[event == 1 & treat == 1])
    
    c.int <- .surv.confints(fit.times, 
                            surv.1$surv, 
                            surv.1$IF.vals, 
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

CFsurvival2 <- function(time, 
                        event,
                        time_N,
                        event_N,
                        admin_t,
                        treat, 
                        W,
                        id,
                        fit.times=sort(unique(time[time > 0 & time < max(time[event == 1])])), 
                        fit.treat=c(0,1), 
                        nuisance.options = list(), 
                        conf.band=FALSE, 
                        conf.level=.95, 
                        contrasts = c("risk.ratio"), 
                        verbose=TRUE) {
  
  .args <- mget(names(formals()), sys.frame(sys.nframe()))
  
  # if(is.null(contrasts)) contrasts <- c("surv.diff", "surv.ratio")
  
  nuis <- do.call("CFsurvival.nuisance.options", nuisance.options)
  
  n <- length(time)
  
  if(length(unique(admin_t))==1){
    admin_t = rep(admin_t, n)
  }else if(length(unique(admin_t))==n){
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
  if(any(fit.times == 0)) fit.times <- fit.times[fit.times > 0]
  if(any(fit.times > max(time[event == 1]))) {
    fit.times <- fit.times[fit.times <= max(time[event == 1])]
    message("fit.times > max(time[event == 1]) removed.")
  }
  
  
  if(is.null(nuis$eval.times)) nuis$eval.times <- sort(unique(c(0,time[time > 0 & time <= max(fit.times)], max(fit.times))))
  k <- length(nuis$eval.times)
  
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
        train <- nuis$fold != v
        test <- nuis$fold == v
        
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
    nuis$w <- vector(length=0)
    nuis$id <- vector(length=0)
    if(0 %in% fit.treat & is.null(nuis$event.pred.0)) {
      do.event.pred.0 <- TRUE
      nuis$event.pred.0 <- matrix(NA, nrow=0, ncol=k)
    } else do.event.pred.0 <- FALSE
    if(1 %in% fit.treat & is.null(nuis$event.pred.1)) {
      do.event.pred.1 <- TRUE
      nuis$event.pred.1 <- matrix(NA, nrow=0, ncol=k)
    } else do.event.pred.1 <- FALSE
    if(0 %in% fit.treat & is.null(nuis$cens.pred.0)) {
      do.cens.pred.0 <- TRUE
      nuis$cens.pred.0 <- matrix(NA, nrow=0, ncol=k)
    } else do.cens.pred.0 <- FALSE
    if(1 %in% fit.treat & is.null(nuis$cens.pred.1)) {
      do.cens.pred.1 <- TRUE
      nuis$cens.pred.1 <- matrix(NA, nrow=0, ncol=k)
    } else{do.cens.pred.1 <- FALSE}
    if(nuis$V > 1) {
      if(nuis$save.nuis.fits) result$surv.fits <- vector(mode='list',length=nuis$V)
      nuis$event.coef <- nuis$cens.coef <- list()
      for(v in 1:nuis$V) {
        if(verbose) message(paste("Fold ", v, "..."))
        train <- nuis$fold != v
        test <- nuis$fold == v
        
        ### We MUST fit nuisance parameters on the training fold
        ### This code fits N | X to training data and outputs augmented training dataset
        dat_train<- as.data.frame(cbind("id"= id, "Y"=time, "DeltaY"=event, "N"=time_N, 'DeltaN'=event_N, "A"=treat, 'admin_t' = admin_t, W))[train,]
        dat_test<- as.data.frame(cbind("id" = id, "Y"=time, "DeltaY"=event, "N"=time_N, 'DeltaN'=event_N, "A"=treat, 'admin_t' = admin_t, W))[test,]
        
        if(nuis$verbose){message("Fitting censored covariate distribution N|X")}
        
        FN<-estimate_conditional_survival_N(dat=dat_train, 
                                            N_name="N", 
                                            DeltaN_name = "DeltaN", 
                                            A_name="A", 
                                            X_names = colnames(W), 
                                            admin_t = admin_t, 
                                            event.SL.library = nuis$event.SL.library,
                                            cens.SL.library = nuis$cens.SL.library,
                                            survSL.control = nuis$survSL.control,
                                            survSL.cvControl = nuis$survSL.cvControl,
                                            verbose=FALSE)
        
        augment_dataset <- function(i, df, id, A, Y, DeltaY, W, N, DeltaN, admin_t, FN){
          
          survest_1 <- FN$survest_1
          survest_0 <- FN$survest_0
          
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
                bind_cols(data.frame(
                  id     = row_i[[id]],
                  N      = row_i[[N]],
                  DeltaN = row_i[[DeltaN]],
                  Y      = row_i[[Y]],
                  DeltaY = row_i[[DeltaY]],
                  A      = row_i[[A]],
                  admin_t = admin_t_tmp,
                  w      = 1
                ), row_i %>% select(all_of(W))))
            }
          }
          
          # Case 2: N censored → redistribution
          
          # plausible N values
          N_grid <- N_support[N_support > N0]
          
          if (length(N_grid) == 0){
            return(bind_cols(data.frame(
              id     = row_i[[id]],
              N      = admin_t_tmp,
              DeltaN = 0,
              Y      = row_i[[Y]],
              DeltaY = row_i[[DeltaY]],
              A      = row_i[[A]],
              admin_t = admin_t_tmp,
              w      = 1
            ), row_i %>% select(all_of(W))))
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
            data.frame(
              id     = row_i[[id]],
              N      = N_grid,
              DeltaN = as.numeric(N_grid < admin_t_tmp),
              Y      = row_i[[Y]],
              DeltaY = row_i[[DeltaY]],
              A      = row_i[[A]],
              X      = NA,
              admin_t = admin_t_tmp,
              w      = w
            ) %>% filter(w>0)
          }else{
            bind_cols(data.frame(
              id     = row_i[[id]],
              N      = N_grid,
              DeltaN = as.numeric(N_grid < admin_t_tmp),
              Y      = row_i[[Y]],
              DeltaY = row_i[[DeltaY]],
              A      = row_i[[A]],
              admin_t = admin_t_tmp,
              w      = w
            ) %>% filter(w>0), row_i %>% select(all_of(W)))
          }
        }
        
        if(nuis$verbose){
        message("Augmenting training dataset using redistribute-to-the-right via N|X")
        pb <- progress_bar$new(total = nrow(dat_train),
                               format = "[:bar] :percent ETA: :eta")
        }
        aug_train_list <- vector(mode="list", length=nrow(dat_train))
        for(i in 1:nrow(dat_train)){
          aug_train_list[[i]]<-augment_dataset(i, 
                                               df=dat_train, 
                                               id="id",
                                               Y="Y", 
                                               DeltaY='DeltaY', 
                                               A="A", 
                                               N="N", 
                                               W='X', 
                                               DeltaN='DeltaN', 
                                               admin_t='admin_t',
                                               FN = FN)
          if(nuis$verbose){pb$tick()}
        }
        
        ### Create augmented training data
        aug_train <- do.call(rbind, aug_train_list)
        rm(aug_train_list)
        
        if(nuis$verbose){
          message("Augmenting test dataset using redistribute-to-the-right via N|X")
          pb <- progress_bar$new(total = nrow(dat_test),
                               format = "[:bar] :percent ETA: :eta")
        }
        
        aug_test_list <- vector(mode="list", length=nrow(dat_test))
        for(i in 1:nrow(dat_test)){
          aug_test_list[[i]]<-augment_dataset(i, 
                                              df=dat_test, 
                                              id = 'id',
                                              Y="Y", 
                                              DeltaY='DeltaY', 
                                              A="A", 
                                              N="N", 
                                              W='X', 
                                              DeltaN='DeltaN', 
                                              admin_t='admin_t',
                                              FN = FN)
          if(nuis$verbose){pb$tick()}
        }
        
        ### Create augmented training data
        aug_test <- do.call(rbind, aug_test_list)
        rm(aug_test_list)
        
        message("Fit conditional survival distributions Y|N,X")
        surv.fit <- .estimate.conditional.survival.2(Y=aug_train[["Y"]], 
                                                     Delta=aug_train[["DeltaY"]], 
                                                     A=aug_train[["A"]], 
                                                     W_event=aug_train %>% dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                     W_cens=aug_train %>% dplyr::select(colnames(W)),
                                                     newW_event = aug_test %>% dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                     newW_cens = aug_test %>% dplyr::select(colnames(W)),
                                                     event.SL.library=nuis$event.SL.library, 
                                                     fit.times=nuis$eval.times, 
                                                     fit.treat=fit.treat, 
                                                     cens.SL.library=nuis$cens.SL.library, 
                                                     survSL.control=nuis$survSL.control, 
                                                     survSL.cvControl = nuis$survSL.cvControl, 
                                                     cens.trunc=nuis$cens.trunc, 
                                                     save.fit = nuis$save.nuis.fits, 
                                                     verbose = nuis$verbose,
                                                     obsWeights = aug_train$w, 
                                                     id = aug_train$id)
        
        if(do.event.pred.0) nuis$event.pred.0 <- rbind(nuis$event.pred.0, surv.fit$event.pred.0)
        if(do.event.pred.1) nuis$event.pred.1 <- rbind(nuis$event.pred.1, surv.fit$event.pred.1)
        if(do.cens.pred.0) nuis$cens.pred.0 <- rbind(nuis$cens.pred.0, surv.fit$cens.pred.0)
        if(do.cens.pred.1) nuis$cens.pred.1 <- rbind(nuis$cens.pred.1, surv.fit$cens.pred.1)
        nuis$w <- c(nuis$w, aug_test$w)
        nuis$id <- c(nuis$id, aug_test$id)
        if(nuis$save.nuis.fits) result$surv.fits[[v]] <- surv.fit$surv.fit
        nuis$event.coef[[v]] <- surv.fit$event.coef
        nuis$cens.coef[[v]] <- surv.fit$cens.coef
      }
    } else {
      
      if(nuis$save.nuis.fits) result$surv.fits <- vector(mode='list',length=nuis$V)
      nuis$event.coef <- nuis$cens.coef <- list()
      
      ### We MUST fit nuisance parameters on the training fold
      ### This code fits N | X to training data and outputs augmented training dataset
      dat<- as.data.frame(cbind("id"= id, 
                                "Y"=time, 
                                "DeltaY"=event, 
                                "N"=time_N, 
                                'DeltaN'=event_N, 
                                "A"=treat, 
                                'admin_t'=admin_t,
                                W))
      
      if(nuis$verbose){message("Fitting censored covariate distribution N|X")}
      
      FN<-estimate_conditional_survival_N(dat=dat, 
                                          N_name="N", 
                                          DeltaN_name = "DeltaN", 
                                          A_name="A", 
                                          X_names = colnames(W), 
                                          admin_t = admin_t, 
                                          event.SL.library = nuis$event.SL.library,
                                          cens.SL.library = nuis$cens.SL.library,
                                          survSL.control = nuis$survSL.control,
                                          survSL.cvControl = nuis$survSL.cvControl,
                                          verbose=FALSE)
      
      augment_dataset <- function(i, df, id, A, Y, DeltaY, W, N, DeltaN, admin_t, FN){
        
        survest_1 <- FN$survest_1
        survest_0 <- FN$survest_0
        
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
              bind_cols(data.frame(
                id     = row_i[[id]],
                N      = row_i[[N]],
                DeltaN = row_i[[DeltaN]],
                Y      = row_i[[Y]],
                DeltaY = row_i[[DeltaY]],
                A      = row_i[[A]],
                admin_t = admin_t_tmp,
                w      = 1
              ), row_i %>% select(all_of(W)))
            )
          }
        }
        
        # Case 2: N censored → redistribution
        
        # plausible N values
        N_grid <- N_support[N_support > N0]
        
        if (length(N_grid) == 0){
          return(
            bind_cols(data.frame(
            id     = row_i[[id]],
            N      = admin_t_tmp,
            DeltaN = 0,
            Y      = row_i[[Y]],
            DeltaY = row_i[[DeltaY]],
            A      = row_i[[A]],
            admin_t = admin_t_tmp,
            w      = 1
          ), row_i %>% select(all_of(W)))
        )
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
          data.frame(
            id     = row_i[[id]],
            N      = N_grid,
            DeltaN = as.numeric(N_grid < admin_t_tmp),
            Y      = row_i[[Y]],
            DeltaY = row_i[[DeltaY]],
            A      = row_i[[A]],
            X      = NA,
            admin_t = admin_t_tmp,
            w      = w
          ) %>% filter(w>0)
        }else{
          bind_cols(
          data.frame(
            id     = row_i[[id]],
            N      = N_grid,
            DeltaN = as.numeric(N_grid < admin_t_tmp),
            Y      = row_i[[Y]],
            DeltaY = row_i[[DeltaY]],
            A      = row_i[[A]],
            admin_t = admin_t_tmp,
            w      = w
          ) %>% filter(w>0), row_i %>% select(all_of(W)))
        }
      }
      
      if(nuis$verbose){
        message("Augmenting training dataset using redistribute-to-the-right via N|X")
        pb <- progress_bar$new(total = nrow(dat),
                               format = "[:bar] :percent ETA: :eta")
      }
      aug_list <- vector(mode="list", length=nrow(dat))
      for(i in 1:nrow(dat)){
        aug_train_list[[i]]<-augment_dataset(i, 
                                             df=dat, 
                                             id="id",
                                             Y="Y", 
                                             DeltaY='DeltaY', 
                                             A="A", 
                                             N="N", 
                                             W='X', 
                                             DeltaN='DeltaN', 
                                             admin_t='admin_t',
                                             FN = FN)
        if(nuis$verbose){pb$tick()}
      }
      
      ### Create augmented training data
      aug <- do.call(rbind, aug_list)
      rm(aug_list)
      
      message("Fit conditional survival distributions Y|N,X")
      surv.fit <- .estimate.conditional.survival.2(Y=aug[["Y"]], 
                                                   Delta=aug[["DeltaY"]], 
                                                   A=aug[["A"]], 
                                                   W_event=aug %>% dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                   W_cens=aug %>% dplyr::select(colnames(W)),
                                                   newW_event = aug %>% dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                   newW_cens = aug %>% dplyr::select(colnames(W)),
                                                   event.SL.library=nuis$event.SL.library, 
                                                   fit.times=nuis$eval.times, 
                                                   fit.treat=fit.treat, 
                                                   cens.SL.library=nuis$cens.SL.library, 
                                                   survSL.control=nuis$survSL.control, 
                                                   survSL.cvControl = nuis$survSL.cvControl, 
                                                   cens.trunc=nuis$cens.trunc, 
                                                   save.fit = nuis$save.nuis.fits, 
                                                   verbose = nuis$verbose,
                                                   obsWeights = aug$w, 
                                                   id = aug$id)
      
      if(do.event.pred.0) nuis$event.pred.0 <- surv.fit$event.pred.0
      if(do.event.pred.1) nuis$event.pred.1 <- surv.fit$event.pred.1
      if(do.cens.pred.0) nuis$cens.pred.0 <- surv.fit$cens.pred.0
      if(do.cens.pred.1) nuis$cens.pred.1 <- surv.fit$cens.pred.1
      nuis$w <- c(nuis$w, aug_test$w)
      nuis$id <- c(nuis$id, aug_test$id)
      if(nuis$save.nuis.fits) result$surv.fits <- surv.fit$surv.fit
      nuis$event.coef <- surv.fit$event.coef
      nuis$cens.coef <- surv.fit$cens.coef
    }
  }
  
  #### ESTIMATE CF SURVIVALS ####
  if(verbose) message("Computing counterfactual survivals...")
  if(0 %in% fit.treat) {
    
    if(verbose){message("Computating Survival Probability for A=0")}
    
    surv.0 <- .get.survival.2(Y=time[nuis$id], 
                              Delta=event[nuis$id], 
                              A=1-treat[nuis$id], 
                              fit.times=fit.times, 
                              eval.times=nuis$eval.times, 
                              S.hats=nuis$event.pred.0, 
                              G.hats=nuis$cens.pred.0, 
                              g.hats=1-nuis$prop.pred[nuis$id], 
                              obsWeights = nuis$w, 
                              id = nuis$id, 
                              isotonize = TRUE)
    
    surv.df.0 <- data.frame(time=unique(c(0,fit.times)), trt=0, surv=c(1, surv.0$surv.iso))
    result$IF.vals.0 <- surv.0$IF.vals
    
    q.01 <- quantile(time[event == 1 & treat == 0], .01)#min(fit.times[surv.0$surv <= .99])
    
    q.99.inds <- !is.na(surv.0$surv.iso)
    q.99.inds[q.99.inds &  surv.0$surv.iso < .01] <- FALSE
    q.99 <- max(fit.times[q.99.inds])
    min.obs.time <- min(time[event == 1 & treat == 0])
    
    c.int <- .surv.confints(fit.times, 
                            surv.0$surv, 
                            surv.0$IF.vals, 
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
    
    if(verbose){message("Computating Survival Probability for A=1")}

    surv.1 <- .get.survival.2(Y=time[nuis$id], 
                              Delta=event[nuis$id], 
                              A=treat[nuis$id], 
                              fit.times=fit.times, 
                              eval.times=nuis$eval.times, 
                              S.hats=nuis$event.pred.1, 
                              G.hats=nuis$cens.pred.1, 
                              g.hats=nuis$prop.pred[nuis$id], 
                              obsWeights = nuis$w, 
                              id = nuis$id, 
                              isotonize = TRUE)
    
    surv.df.1 <- data.frame(time=unique(c(0,fit.times)), trt=1, surv=c(1, surv.1$surv.iso))
    result$IF.vals.1 <- surv.1$IF.vals
    
    q.01 <- quantile(time[event == 1 & treat == 1], .01)#min(fit.times[surv.1$surv <= .99])
    q.99 <- max(fit.times[surv.1$surv.iso >= .01])
    min.obs.time <- min(time[event == 1 & treat == 1])
    
    c.int <- .surv.confints(fit.times, 
                            surv.1$surv, 
                            surv.1$IF.vals, 
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

# Nuisance options for CFsurvival2

CFsurvival.nuisance.options <- function(cross.fit = TRUE, V = ifelse(cross.fit, 10, 1), folds = NULL, eval.times = NULL,
                                        event.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) c(alg, "survscreen.glmnet", "survscreen.marg", "All") ),  event.pred.0 = NULL, event.pred.1 = NULL,
                                        cens.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg", "survSL.weibreg", "survSL.loglogreg", "survSL.gam", "survSL.rfsrc"), function(alg) c(alg, "survscreen.glmnet", "survscreen.marg", "All") ),  cens.trunc=0, cens.pred.0 = NULL, cens.pred.1 = NULL,
                                        survSL.control = list(initWeightAlg = "survSL.rfsrc", verbose=FALSE), survSL.cvControl = list(V = 10), save.nuis.fits = FALSE,
                                        prop.SL.library = lapply(c("SL.mean", "SL.glm", "SL.gam", "SL.earth", "SL.ranger", "SL.xgboost"), function(alg) c(alg, "screen.glmnet", "screen.corRank", "All") ), prop.trunc=0, prop.pred = NULL,
                                        verbose=FALSE) {
  list(cross.fit = cross.fit, V = V, folds = folds, eval.times = eval.times,
       event.SL.library = event.SL.library, event.pred.0 = event.pred.0, event.pred.1 = event.pred.1,
       cens.SL.library = cens.SL.library, cens.trunc=cens.trunc, cens.pred.0 = cens.pred.0, cens.pred.1 = cens.pred.1,
       survSL.control = survSL.control,
       survSL.cvControl = survSL.cvControl,
       save.nuis.fits = save.nuis.fits,
       prop.SL.library = prop.SL.library,  prop.trunc=prop.trunc, prop.pred = prop.pred, verbose=verbose)
}

# Modified function to check input for CFsurvival2 function

.check.input <- function(time, event, treat, W_event, W_cens, fit.times, fit.treat, nuisance.options, conf.band, conf.level, contrasts, verbose, ...) {
  if(any(time < 0)) stop("Only non-negative event/censoring times allowed!")
  if(any(time == 0 & event == 1)) stop("Events at time zero not allowed.")
  if(any(!(event %in% c(0,1)))) stop("Event must be binary.")
  if(any(!(treat %in% c(0,1)))) stop("Treatment must be binary.")
  if(length(time) != length(event) | length(time) != length(treat)) stop("time, event, and treat must be n x 1 vectors")
  if(any(!(fit.treat %in% c(0,1)))) stop("fit.treat must be a subset of c(0,1).")
  if(!is.null(nuisance.options$event.pred.1)) {
    if(!all.equal(dim(nuisance.options$event.pred.1), c(length(time), length(nuisance.options$eval.times)))) {
      stop("event.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
    }
  }
  if(!is.null(nuisance.options$event.pred.0)) {
    if(!all.equal(dim(nuisance.options$event.pred.0), c(length(time), length(nuisance.options$eval.times)))) {
      stop("event.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
    }
  }
  if(!is.null(nuisance.options$cens.pred.1)) {
    if(!all.equal(dim(nuisance.options$cens.pred.1), c(length(time), length(nuisance.options$eval.times)))) {
      stop("cens.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
    }
  }
  if(!is.null(nuisance.options$cens.pred.0)) {
    if(!all.equal(dim(nuisance.options$cens.pred.0), c(length(time), length(nuisance.options$eval.times)))) {
      stop("cens.pred must be an n x k matrix (n is number of observations, k is length of eval.times)")
    }
  }
  if(is.null(nuisance.options$event.SL.library)) {
    if(0 %in% fit.treat & is.null(nuisance.options$event.pred.0)) {
      stop("event.pred.0 must be provided if event.SL.library is not specified and 0 is a treatment of interest.")
    }
    if(1 %in% fit.treat & is.null(nuisance.options$event.pred.1)) {
      stop("event.pred.1 must be provided if event.SL.library is not specified and 1 is a treatment of interest.")
    }
  }
  if(is.null(nuisance.options$cens.SL.library)) {
    if(0 %in% fit.treat & is.null(nuisance.options$cens.pred.0)) {
      stop("cens.pred.0 must be provided if cens.SL.library is not specified and 0 is a treatment of interest.")
    }
    if(1 %in% fit.treat & is.null(nuisance.options$cens.pred.1)) {
      stop("cens.pred.1 must be provided if cens.SL.library is not specified and 1 is a treatment of interest.")
    }
  }
  if(length(contrasts) > 0 & !identical(sort(fit.treat), c(0,1))) {
    warning("contrast specified but both treatment regimens not requested -- contrasts will not be provided. Re-run with fit.treat = c(0,1) for survival contrasts.")
  }
  if(any(is.na(time) | is.na(event))) {
    stop("Missing time or event detected; missing data not allowed.")
  }
  if(any(is.na(W_event)) | any(is.na(W_cens)) | any(is.na(treat))) {
    warning("Missing confounders or exposures detected; missing data not allowed.")
  }
}

#' Estimate counterfactual survival functions based on nuisance parameter fits to augmented dataset
#'
#' This function computes the one-step estimator using the EIF based on coarsening of the right-censored negative control predictor. It is applied to the "augmented" dataset where the weight of observations with right-censored predictor were redistributed to the right according to the conditional survival function. These weights are included in the argument 'obsWeights'.

.get.survival.2 <- function(Y, Delta, A, fit.times, eval.times, S.hats, G.hats, g.hats, obsWeights, id, isotonize=TRUE) {
  
  # Y=time[nuis$id];
  # Delta=event[nuis$id];
  # A=1-treat[nuis$id];
  # fit.times=c(0, fit.times); 
  # eval.times=nuis$eval.times[-1]; 
  # S.hats=nuis$event.pred.0;
  # G.hats=nuis$cens.pred.0;
  # g.hats=1-nuis$prop.pred[nuis$id]; 
  # obsWeights = nuis$w;
  # id = nuis$id;
  # isotonize = TRUE
  
  # Collect fitting times, order them, and order (conditional) survival and censoring probabilities
  fit.times <- fit.times[fit.times > 0]
  n <- length(Y)
  N <- length(unique(id))
  ord <- order(eval.times)
  eval.times <- eval.times[ord]
  S.hats <- S.hats[,ord]
  G.hats <- G.hats[,ord]
  
  # Discretized compensator -- used to approximate the integral term in the EIF
  # loops through rows of S.hats and G.hats
  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S.hats[i,]) * 1/ G.hats[i,-ncol(G.hats)] #corrects for censoring
    if(any(eval.times[-1] > Y[i])) vals[eval.times[-1] > Y[i]] <- 0
    c(0,cumsum(vals))
  }))
  
  # Estimate survival/censoring at observed times Y
  S.hats.Y <- sapply(1:n, function(i) stepfun(eval.times, c(1,S.hats[i,]), right = FALSE)(Y[i]))
  G.hats.Y <- sapply(1:n, function(i) stepfun(eval.times, c(1,G.hats[i,]), right = TRUE)(Y[i]))
  
  IF.vals <- matrix(NA, nrow=N, ncol=length(fit.times))
  surv <- rep(NA, length(fit.times))
  
  # loop through target times t0
  for(t0 in fit.times) {
    k <- min(which(eval.times >= t0))
    #S.hats.t0 <- S.hats[,k]
    S.hats.t0 <- sapply(1:n, function(i){S.hats[i,k]})
    # This corresponds to the I(Y ≤ t0, \Delta=1)/{S(Y|X)G(Y|X)}, the first interior term in EIF
    inner.func.1 <- ifelse(Y <= t0 & Delta == 1, 1/(S.hats.Y * G.hats.Y), 0 )
    # This corresponds to the compensator, the second interior term of the EIF
    inner.func.2 <- int.vals[,k]
    # Assemble EIF -- with weights that correspond to redistribution weights
    if.func <- as.numeric(A == 1) * S.hats.t0 * ( -inner.func.1 + inner.func.2) / g.hats + S.hats.t0
    if.func <- obsWeights * if.func
    if.func <- tapply(if.func, id, sum)
    k1 <- which(fit.times == t0)
    # survival is estimated as mean of EIF, as this is the uncentered version
    surv[k1] <- mean(if.func)
    # the "centered" if values are the if minus the survival probabilities (NOW mean zero)
    IF.vals[,k1] <- if.func - surv[k1]
  }
  res <- list(times=fit.times, surv=pmin(1,pmax(0,surv)), IF.vals=IF.vals)
  if(isotonize) {
    res$surv.iso <- NA
    res$surv.iso[!is.na(res$surv)] <- 1 - isoreg(res$times[!is.na(res$surv)], 1-res$surv[!is.na(res$surv)])$yf
  }
  
  return(res)
}

# Same as from CFSurvival -- function for computing confidence intervals based on IF values

.surv.confints <- function(times, est, IF.vals, isotonize=TRUE, conf.band=TRUE, band.end.pts=c(0,Inf), conf.level=.95) {
  logit <- function(x) log(x / (1-x))
  logit.prime <- function(x) 1/(x * (1-x))
  expit <- function(x) 1/(1 + exp(-x))
  
  n <- nrow(IF.vals)
  in.bounds <- est > 0 & est < 1
  IF.vals.logit <- IF.vals
  for(j in 1:length(est)) IF.vals.logit[,j] <- IF.vals[,j] * logit.prime(est[j])
  
  res <- NULL
  res$se <- sqrt(colMeans(IF.vals^2)) / sqrt(n)
  res$se.logit <- sqrt(colMeans(IF.vals.logit^2)) / sqrt(n)
  
  res$se[res$se == 0] <- NA
  res$se.logit[is.infinite(res$se.logit)] <- NA
  quant <- qt(1-(1-conf.level)/2, n - 1) #qnorm(1-(1-conf.level)/2)
  
  #Intervals and bands based on un-isotonized survival ests
  res$ptwise.lower <- pmax(est - quant * res$se, 0)
  res$ptwise.upper <- pmin(est + quant * res$se, 1)
  res$ptwise.logit.lower <- expit(logit(est) - quant * res$se.logit) #pmax(est - quant * res$se, 0)
  res$ptwise.logit.upper <- expit(logit(est) + quant * res$se.logit) #pmin(est + quant * res$se, 1)
  #
  # F.vals <-  isoreg(times, 1-est)$yf
  # F.inv.hat <- stepfun(F.vals[-1], times)
  #
  # ln.q.ptwise <- as.numeric(F.inv.hat(expit(logit(F.vals) - quant * res$se.logit )))
  # un.q.ptwise <- as.numeric(F.inv.hat(expit(logit(F.vals) + quant * res$se.logit )))
  
  # res$ptwise.inv.upper <- 1-sapply(times, function(t) {
  #     if(any(un.q.ptwise <= t,na.rm=TRUE)) F.vals[max(which(un.q.ptwise <= t), na.rm=TRUE)]
  #     else 0
  # })
  # res$ptwise.inv.lower <- 1-sapply(times, function(t) {
  #     if(any(ln.q.ptwise >= t,na.rm=TRUE)) F.vals[min(which(ln.q.ptwise >= t),na.rm=TRUE)]
  #     else min(ln.q.ptwise, na.rm=TRUE)
  # })
  
  # Isotonized intervals and bands
  out <- NULL
  if(conf.band) {
    if(any(!is.na(res$se))) {
      unif.vals <- .estimate.uniform.quantile(IF.vals[,!is.na(res$se), drop=FALSE], conf.level, scale = FALSE)
      unif.quant <- unif.vals$quantile
      out$ew.sim.maxes <- unif.vals$maxes
      out$unif.ew.quant <- unif.quant
      res$unif.ew.lower <- pmax(est - unif.quant / sqrt(n), 0)
      res$unif.ew.upper <- pmin(est + unif.quant / sqrt(n), 1)
      
      if(isotonize) {
        res$unif.ew.lower[!is.na(res$unif.ew.lower)] <- 1 - isoreg(times[!is.na(res$unif.ew.lower)], 1-res$unif.ew.lower[!is.na(res$unif.ew.lower)])$yf
        res$unif.ew.upper[!is.na(res$unif.ew.upper)] <- 1 - isoreg(times[!is.na(res$unif.ew.upper)], 1-res$unif.ew.upper[!is.na(res$unif.ew.upper)])$yf
      }
    } else {
      res$unif.ew.lower <- rep(NA, length(res$se))
      res$unif.ew.upper <- rep(NA, length(res$se))
    }
    
    if(any(!is.na(res$se.logit) & times >= band.end.pts[1] & times <= band.end.pts[2]))  {
      unif.logit.vals <- .estimate.uniform.quantile(IF.vals.logit[,!is.na(res$se.logit) & times >= band.end.pts[1] & times <= band.end.pts[2], drop=FALSE], conf.level)
      unif.logit.quant <- unif.logit.vals$quantile
      out$logit.sim.maxes <- unif.logit.vals$maxes
      out$unif.logit.quant <- unif.logit.quant
      res$unif.logit.lower <- expit(logit(est) - unif.logit.quant * res$se.logit) #pmax(est - unif.quant * res$se, 0)
      res$unif.logit.upper <- expit(logit(est) + unif.logit.quant * res$se.logit) # pmin(est + unif.quant * res$se, 1)
      res$unif.logit.lower[times < band.end.pts[1] | times > band.end.pts[2]] <- NA
      res$unif.logit.upper[times < band.end.pts[1] | times > band.end.pts[2]] <- NA
      # res$unif.lower[!in.bounds] <- pmin(approx(times[in.bounds], res$unif.lower[in.bounds], xout=times[!in.bounds], rule=2)$y, est[!in.bounds])
      # res$unif.upper[!in.bounds] <- pmax(approx(times[in.bounds], res$unif.lower[in.bounds], xout=times[!in.bounds], rule=2)$y, est[!in.bounds])
      if(isotonize) {
        res$unif.logit.lower[!is.na(res$unif.logit.lower)] <- 1 - isoreg(times[!is.na(res$unif.logit.lower)], 1-res$unif.logit.lower[!is.na(res$unif.logit.lower)])$yf
        res$unif.logit.upper[!is.na(res$unif.logit.upper)] <- 1 - isoreg(times[!is.na(res$unif.logit.upper)], 1-res$unif.logit.upper[!is.na(res$unif.logit.upper)])$yf
      }
    } else {
      res$unif.logit.lower <-  rep(NA, length(res$se))
      res$unif.logit.upper <-  rep(NA, length(res$se))
    }
  }
  out$res <- res
  return(out)
}

# Same as from CFSurvival -- function for computing confidence band

.estimate.uniform.quantile <- function(IF.vals, conf.level=.95, scale=TRUE) {
  IF.vals <- cbind(IF.vals)
  n <- nrow(IF.vals)
  if(scale) IF.vals <- scale(IF.vals)
  maxes <- replicate(1e4, max(abs(rbind(rt(n, df = n - 1)/sqrt(n)) %*% IF.vals)))
  return(list(quantile=quantile(maxes, conf.level), maxes=maxes))
}

# Adapations of nuisance estimators for our problem

.estimate.conditional.survival <- function(Y, Delta, A, W, newW, fit.times, fit.treat, event.SL.library, cens.SL.library, survSL.control, survSL.cvControl, cens.trunc, verbose, save.fit) {
  ret <- list(fit.times=fit.times)
  AW <- cbind(A, W)
  if(0 %in% fit.treat & 1 %in% fit.treat) {
    newAW <- rbind(cbind(A=0, newW), cbind(A=1, newW))
  } else {
    newAW <- cbind(A=fit.treat, newW)
  }
  res <- require(survSuperLearner)
  if(!res) stop("Please install the package survSuperLearner via:\n devtools::install_github('tedwestling/survSuperLearner')")
  if(is.null(survSL.control)) survSL.control <- list(saveFitLibrary = save.fit)
  
  fit <- survSuperLearner(time = Y, event = Delta,  X = AW, newX = newAW, new.times = fit.times, event.SL.library = event.SL.library, cens.SL.library = cens.SL.library, verbose=verbose, control = survSL.control, cvControl = survSL.cvControl)
  
  if(save.fit) ret$surv.fit <- fit
  if(0 %in% fit.treat) {
    ret$event.pred.0 <- fit$event.SL.predict[1:nrow(newW),]
    if(any(ret$event.pred.0 == 0)) ret$event.pred.0[ret$event.pred.0 == 0] <- min(ret$event.pred.0[ret$event.pred.0 > 0])
    ret$cens.pred.0 <- fit$cens.SL.predict[1:nrow(newW),]
    ret$cens.pred.0 <- pmax(ret$cens.pred.0, cens.trunc)
    if(any(ret$cens.pred.0 == 0)) ret$cens.pred.0[ret$cens.pred.0 == 0] <- min(ret$cens.pred.0[ret$cens.pred.0 > 0])
    if(1 %in% fit.treat) {
      ret$event.pred.1 <- fit$event.SL.predict[-(1:nrow(newW)),]
      if(any(ret$event.pred.1 == 0)) ret$event.pred.1[ret$event.pred.1 == 0] <- min(ret$event.pred.1[ret$event.pred.1 > 0])
      
      ret$cens.pred.1 <- fit$cens.SL.predict[-(1:nrow(newW)),]
      ret$cens.pred.1 <- pmax(ret$cens.pred.1, cens.trunc)
      if(any(ret$cens.pred.1 == 0)) ret$cens.pred.1[ret$cens.pred.1 == 0] <- min(ret$cens.pred.1[ret$cens.pred.1 > 0])
    }
  } else {
    ret$event.pred.1 <- fit$event.SL.predict
    if(any(ret$event.pred.1 == 0)) ret$event.pred.1[ret$event.pred.1 == 0] <- min(ret$event.pred.1[ret$event.pred.1 > 0])
    ret$cens.pred.1 <- fit$cens.SL.predict
    ret$cens.pred.1 <- pmax(ret$cens.pred.1, cens.trunc)
    if(any(ret$cens.pred.1 == 0)) ret$cens.pred.1[ret$cens.pred.1 == 0] <- min(ret$cens.pred.1[ret$cens.pred.1 > 0])
  }
  ret$event.coef <- fit$event.coef
  ret$cens.coef <- fit$cens.coef
  return(ret)
}

.estimate.propensity <- function(A, W, newW, SL.library, fit.treat, prop.trunc, save.fit, verbose) {
  ret <- list()
  library(SuperLearner)
  if (length(SL.library) == 1) {
    if(length(unlist(SL.library)) == 2 & ncol(W) > 1) {
      screen <- get(SL.library[[1]][2])
      whichScreen <- screen(Y = A, X = W, family = 'binomial')
    } else {
      whichScreen <- rep(TRUE, ncol(W))
    }
    learner <- get(SL.library[[1]][1])
    prop.fit <- learner(Y = A, X = W[,whichScreen, drop=FALSE], newX = newW[,whichScreen, drop=FALSE], family='binomial')
    ret$prop.pred <- prop.fit$pred
    if(save.fit) {
      ret$prop.fit <- list(whichScreen = whichScreen, pred.alg = prop.fit)
    }
  } else {
    prop.fit <- SuperLearner(Y=A, X=W, newX=newW, family='binomial',
                             SL.library=SL.library, method = "method.NNloglik", verbose = verbose)
    ret$prop.pred <- c(prop.fit$SL.predict)
    if(save.fit) ret$prop.fit <- prop.fit
  }
  if(1 %in% fit.treat) ret$prop.pred <- pmax(ret$prop.pred, prop.trunc)
  if(0 %in% fit.treat) ret$prop.pred <- pmin(ret$prop.pred, 1-prop.trunc)
  return(ret)
}

.estimate.conditional.survival.2 <- function(Y, Delta, A, W_event, W_cens, newW_event, newW_cens, fit.times, fit.treat, event.SL.library, cens.SL.library, survSL.control, survSL.cvControl, cens.trunc, verbose, save.fit, obsWeights, id){
  
  # Y=aug_train[["Y"]];
  # Delta=aug_train[["DeltaY"]];
  # A=aug_train[["A"]];
  # W_event=aug_train %>% dplyr::select(c("N", "DeltaN", colnames(W)));
  # W_cens=aug_train %>% dplyr::select(colnames(W));
  # newW_event = aug_train %>% dplyr::select(c("N", "DeltaN", colnames(W))) #%>% slice(1:10);
  # newW_cens = aug_train %>% dplyr::select(colnames(W)) #%>% slice(1:10);
  # event.SL.library=nuis$event.SL.library;
  # fit.times=nuis$eval.times;
  # fit.treat=fit.treat;
  # cens.SL.library=nuis$cens.SL.library;
  # survSL.control=nuis$survSL.control;
  # survSL.cvControl = nuis$survSL.cvControl;
  # cens.trunc=nuis$cens.trunc;
  # save.fit = nuis$save.nuis.fits;
  # verbose = nuis$verbose;
  # obsWeights = aug_train$w;
  # id = aug_train$id
  
  ret <- list(fit.times=fit.times)
  AW_event <- cbind(A, W_event)
  AW_cens <- cbind(A, W_cens)
  if(0 %in% fit.treat & 1 %in% fit.treat) {
    newAW_event <- rbind(cbind(A=0, newW_event), cbind(A=1, newW_event))
    newAW_cens <- rbind(cbind(A=0, newW_cens), cbind(A=1, newW_cens))
  } else {
    newAW_event <- cbind(A=fit.treat, newW_event)
    newAW_cens <- cbind(A=fit.treat, newW_cens)
  }
  #res <- require(survSuperLearner)
  #if(!res) stop("Please install the package survSuperLearner via:\n devtools::install_github('tedwestling/survSuperLearner')")
  if(is.null(survSL.control)) survSL.control <- list(saveFitLibrary = save.fit)
  
  fit <- survSuperLearner2(time = Y, 
                           event = Delta, 
                           X_event = AW_event, 
                           X_cens = AW_cens, 
                           newX_event = newAW_event,
                           newX_cens = newAW_cens, 
                           new.times = fit.times, 
                           event.SL.library = event.SL.library, 
                           cens.SL.library = cens.SL.library, 
                           verbose=verbose, 
                           control = survSL.control, 
                           cvControl = survSL.cvControl,
                           obsWeights = obsWeights,
                           id = id)
  
  if(save.fit) ret$surv.fit <- fit
  if(0 %in% fit.treat) {
    ret$event.pred.0 <- fit$event.SL.predict[1:nrow(newW_event),]
    if(any(ret$event.pred.0 == 0)) ret$event.pred.0[ret$event.pred.0 == 0] <- min(ret$event.pred.0[ret$event.pred.0 > 0])
    ret$cens.pred.0 <- fit$cens.SL.predict[1:nrow(newW_cens),]
    ret$cens.pred.0 <- pmax(ret$cens.pred.0, cens.trunc)
    if(any(ret$cens.pred.0 == 0)) ret$cens.pred.0[ret$cens.pred.0 == 0] <- min(ret$cens.pred.0[ret$cens.pred.0 > 0])
    if(1 %in% fit.treat) {
      ret$event.pred.1 <- fit$event.SL.predict[-(1:nrow(newW_event)),]
      if(any(ret$event.pred.1 == 0)) ret$event.pred.1[ret$event.pred.1 == 0] <- min(ret$event.pred.1[ret$event.pred.1 > 0])
      
      ret$cens.pred.1 <- fit$cens.SL.predict[-(1:nrow(newW_cens)),]
      ret$cens.pred.1 <- pmax(ret$cens.pred.1, cens.trunc)
      if(any(ret$cens.pred.1 == 0)) ret$cens.pred.1[ret$cens.pred.1 == 0] <- min(ret$cens.pred.1[ret$cens.pred.1 > 0])
    }
  } else {
    ret$event.pred.1 <- fit$event.SL.predict
    if(any(ret$event.pred.1 == 0)) ret$event.pred.1[ret$event.pred.1 == 0] <- min(ret$event.pred.1[ret$event.pred.1 > 0])
    ret$cens.pred.1 <- fit$cens.SL.predict
    ret$cens.pred.1 <- pmax(ret$cens.pred.1, cens.trunc)
    if(any(ret$cens.pred.1 == 0)) ret$cens.pred.1[ret$cens.pred.1 == 0] <- min(ret$cens.pred.1[ret$cens.pred.1 > 0])
  }
  ret$event.coef <- fit$event.coef
  ret$cens.coef <- fit$cens.coef
  return(ret)
}

.estimate.propensity.2 <- function(A, W, newW, SL.library, fit.treat, prop.trunc, save.fit, id, obsWeights, verbose) {
  ret <- list()
  library(SuperLearner)
  if (length(SL.library) == 1) {
    if(length(unlist(SL.library)) == 2 & ncol(W) > 1) {
      screen <- get(SL.library[[1]][2])
      whichScreen <- screen(Y = A, X = W, family = 'binomial')
    } else {
      whichScreen <- rep(TRUE, ncol(W))
    }
    learner <- get(SL.library[[1]][1])
    prop.fit <- learner(Y = A, X = W[,whichScreen, drop=FALSE], newX = newW[,whichScreen, drop=FALSE], family='binomial', obsWeights = obsWeights)
    ret$prop.pred <- prop.fit$pred
    if(save.fit) {
      ret$prop.fit <- list(whichScreen = whichScreen, pred.alg = prop.fit)
    }
  } else {
    prop.fit <- SuperLearner(Y=A, X=W, newX=newW, family='binomial',
                             SL.library=SL.library, method = "method.NNloglik", verbose = verbose, id=id, obsWeights = obsWeights)
    ret$prop.pred <- c(prop.fit$SL.predict)
    if(save.fit) ret$prop.fit <- prop.fit
  }
  if(1 %in% fit.treat) ret$prop.pred <- pmax(ret$prop.pred, prop.trunc)
  if(0 %in% fit.treat) ret$prop.pred <- pmin(ret$prop.pred, 1-prop.trunc)
  return(ret)
}

# Changes the predict method to align better with inputs -- SUPERIOR while other does not work for some odd reason
predict.survSuperLearner.2 <- function (object, newdata, new.times, X = NULL, time = NULL, event = NULL, onlySL = FALSE, threshold = 1e-4, ...) {
  if (missing(newdata)) {
    out <- list(survpred = object$event.SL.predict, cens.pred = object$cens.SL.predict,
                survlibrary.predict = object$event.library.predict, cens.library.predict = object$cens.library.predict)
    return(out)
  }
  if (!object$control$saveFitLibrary) {
    stop("This SuperLearner fit was created using control$saveFitLibrary = FALSE, so new predictions cannot be made.")
  }
  event.k <- nrow(object$event.libraryNames)
  event.pred <- array(0, dim=c(nrow(newdata), length(new.times), event.k))
  dimnames(event.pred)[[3]] <- apply(object$event.libraryNames, 1, paste, collapse = '_')
  dimnames(event.pred)[[2]] <- new.times
  
  cens.k <- nrow(object$cens.libraryNames)
  cens.pred <- array(0, dim=c(nrow(newdata), length(new.times), cens.k))
  dimnames(cens.pred)[[3]] <- apply(object$cens.libraryNames, 1, paste, collapse = '_')
  dimnames(cens.pred)[[2]] <- new.times
  if (onlySL) {
    event.whichLibrary <- which(object$event.coef > threshold)
    event.coef <- object$event.coef
    event.coef[-event.whichLibrary] <- 0
    event.coef <- event.coef / sum(event.coef)
    
    cens.whichLibrary <- which(object$cens.coef > threshold)
    cens.coef <- object$cens.coef
    cens.coef[-cens.whichLibrary] <- 0
    cens.coef <- cens.coef / sum(cens.coef)
  } else {
    event.whichLibrary <- seq(event.k)
    event.coef <- object$event.coef
    cens.whichLibrary <- seq(cens.k)
    cens.coef <- object$cens.coef
  }
  for (mm in event.whichLibrary) {
    newdataMM <- subset(newdata,
                        select = object$event.whichScreen[object$event.SL.library$library[mm, 2], ])
    XMM <- if (is.null(X)) {
      NULL
    } else {
      subset(X, select = object$event.whichScreen[object$event.SL.library$library[mm, 2], ])
    }
    event.pred[, ,mm] <- do.call(paste0("predict.", attr(object$event.fitLibrary[[mm]], "class")), list(object = object$event.fitLibrary[[mm]],
                                                                                                        newX = newdataMM, new.times = new.times,
                                                                                                        X = XMM, time = time, event = event, ...))
  }
  event.SL.predict <- matrix(NA, nrow = nrow(newdata), ncol=length(new.times))
  for (j in seq(length(new.times))) {
    if(length(event.coef)==1){
      event.SL.predict[,j] <- event.pred[,j,] * event.coef
    }else{
      event.SL.predict[,j] <- event.pred[,j,] %*% event.coef
    }
  }
  
  for (mm in cens.whichLibrary) {
    newdataMM <- subset(newdata,
                        select = object$cens.whichScreen[object$cens.SL.library$library[mm, 2], ])
    XMM <- if (is.null(X)) {
      NULL
    } else {
      subset(X, select = object$cens.whichScreen[object$cens.SL.library$library[mm, 2], ])
    }
    cens.pred[, ,mm] <- do.call(paste0("predict.", attr(object$cens.fitLibrary[[mm]], "class")), list(object = object$cens.fitLibrary[[mm]],
                                                                                                      newX = newdataMM, new.times = new.times,
                                                                                                      X = XMM, time = time, event = event, ...))
  }
  cens.SL.predict <- matrix(NA, nrow = nrow(newdata), ncol=length(new.times))
  for (j in seq(length(new.times))) {
    if(length(cens.coef)==1){
      cens.SL.predict[,j] <- cens.pred[,j,] * cens.coef
    }else{
      cens.SL.predict[,j] <- cens.pred[,j,] %*% cens.coef
    }
  }
  
  out <- list(event.SL.predict = event.SL.predict, event.library.predict = event.pred,
              cens.SL.predict = cens.SL.predict, cens.library.predict = cens.pred)
  return(out)
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
      
      out<-predict.survSuperLearner.2(object = FN_1, 
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
      
      out<-predict.survSuperLearner.2(object = FN_1, 
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

# Function to compute contrast
.risk.ratio <- function(fit.times, surv.0, surv.1, IF.vals.0, IF.vals.1, conf.band=TRUE, band.end.pts=c(0,Inf), conf.level=.95) {
  n <- nrow(IF.vals.0)
  quant <- qt(1-(1-conf.level)/2, df = n - 1)
  risk.0 <- 1-surv.0
  risk.1 <- 1-surv.1
  IF.vals.0 <- -IF.vals.0
  IF.vals.1 <- -IF.vals.1
  df <- data.frame(time=c(0,fit.times), log.risk.ratio=c(NA,log(risk.1) - log(risk.0)), risk.ratio=c(NA,risk.1 / risk.0))
  
  IF.log.ratio <- IF.vals.1 / matrix(risk.1, nrow=nrow(IF.vals.1), ncol=ncol(IF.vals.1), byrow=TRUE) -
    IF.vals.0 / matrix(risk.0, nrow=nrow(IF.vals.0), ncol=ncol(IF.vals.0), byrow=TRUE)
  se.log.ratio <- sqrt(colMeans(IF.log.ratio^2))
  se.log.ratio[se.log.ratio == 0] <- NA
  se.log.ratio[is.infinite(se.log.ratio)] <- NA
  df$se.log.ratio <- c(NA,se.log.ratio)/sqrt(n)
  df$ptwise.lower.log <- df$log.risk.ratio - quant * df$se.log.ratio
  df$ptwise.upper.log <- df$log.risk.ratio + quant * df$se.log.ratio
  df$ptwise.lower <- exp(df$ptwise.lower.log)
  df$ptwise.upper <- exp(df$ptwise.upper.log)
  df$ptwise.pval <- c(1,pchisq((df$log.risk.ratio[-1]/(df$se.log.ratio[-1]))^2, df=1, lower.tail = FALSE))
  
  
  if(conf.band) {
    if(any(!is.na(se.log.ratio) & fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2])) {
      log.unif.info <- .estimate.uniform.quantile(IF.log.ratio[,!is.na(se.log.ratio) & fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2],drop=FALSE], conf.level)
      log.unif.quant <- log.unif.info$quantile
      df$unif.lower.log <- df$log.risk.ratio - log.unif.quant * df$se.log.ratio
      df$unif.upper.log <- df$log.risk.ratio + log.unif.quant * df$se.log.ratio
      df$unif.lower <- exp(df$unif.lower.log)
      df$unif.upper <- exp(df$unif.upper)
      
      df$unif.lower[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA
      df$unif.upper[fit.times < band.end.pts[1] | fit.times > band.end.pts[2]] <- NA
      # df$unif.pval <- c(1, pchisq((df$log.risk.ratio[-1]/( log.unif.quant * df$se.log.ratio[-1]))^2, df=1, lower.tail = FALSE))
    } else {
      df$unif.lower.log <- df$unif.upper.log <- df$unif.lower <- df$unif.upper <- NA
    }
    
  }
  
  res <- list(risk.ratio.df=df)
  if(conf.band & any(!is.na(se.log.ratio) & fit.times >= band.end.pts[1] & fit.times <= band.end.pts[2])) {
    res$log.risk.ratio.unif.quant <- log.unif.quant
    res$log.risk.ratio.sim.maxes <- log.unif.info$sim.maxes
  }
  
  return(res)
}
