# AIPWCC version of CFSurvival

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
  # treat = data$A;
  # W = data %>% dplyr::select("X");
  # time_N = data$N1;
  # event_N = data$DeltaN1;
  # admin_t = 730;
  # prop.fit<-predict(glm(A~X, data=data, family='binomial'), type='response')
  # id = 1:nrow(data);
  # fit.times = sort(unique(c(data$Y[data$Y>0], data$N[data$N>0])));
  # fit.treat = c(0,1);
  # nuisance.options = CFsurvival.nuisance.options(
  #   cross.fit = TRUE,
  #   V = 20,
  #   folds = NULL,
  #   eval.times = NULL,
  #   event.SL.library = c("survSL.km", "survSL.coxph", "survSL.expreg"),
  #   event.pred.0 = NULL,
  #   event.pred.1 = NULL,
  #   cens.SL.library = c("survSL.km", "survSL.coxph", "survSL.expreg"),
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
    nuis$fN <- matrix(NA, nrow=0, ncol=k)
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
        
        if(verbose) message("Step 1: Estimate nuisances via 'One-Way' Factorization")
        
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
        
        # Estimate N | A,X
        if(verbose) message("\tStep 1(a): Estimate distn of N | A, X")
        SN<-estimate_conditional_survival_N(dat=dat_train, 
                                            N_name="N", 
                                            DeltaN_name = "DeltaN", 
                                            A_name="A", 
                                            X_names = colnames(W), 
                                            admin_t = dat_train$admin_t, 
                                            fit.treat=c(0,1),
                                            event.SL.library = nuis$event.SL.library,
                                            cens.SL.library = nuis$cens.SL.library,
                                            survSL.control = nuis$survSL.control,
                                            survSL.cvControl = nuis$survSL.cvControl,
                                            verbose=FALSE)
        # Output -- list containing two functions, survest_1 and survest_0 which 
        # take as input covariates x and time t to output the predicted survivor function for N
        # estimated using SuperLearner ensemble
        
        # Dataset augmentation
        if(verbose) message("\tStep 1(b): Augmenting train fold w/ redistribute-to-the-right weights")
        
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
        
        if(verbose){
          pb <- progress_bar$new(total = nrow(dat_train),
                                 format = "Train fold [:bar] :percent ETA: :eta")
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
                                               SN = SN)
          if(nuis$verbose){pb$tick()}
        }
        aug_train <- do.call(rbind, aug_train_list)
        rm(aug_train_list)
        
        if(verbose) message("\tStep 1(c): Augmenting test along N grid")
        if(verbose){
          pb <- progress_bar$new(total = nrow(dat_test),
                                 format = "Test fold [:bar] :percent ETA: :eta")
        }
        
        aug_test_list <- vector(mode="list", length=nrow(dat_test))
        for(i in 1:nrow(dat_test)){
          aug_test_list[[i]] <- bind_cols(dat_test[i,] %>% dplyr::select(-c("N", "DeltaN")), 
                    data.frame("N"=fit.times[fit.times<=dat_test[i,]$admin_t], 
                               "DeltaN"=c(rep(1, length(fit.times[fit.times<=dat_test[i,]$admin_t])-1), 0)))

          if(nuis$verbose){pb$tick()}
        }
        aug_test <- do.call(rbind, aug_test_list)
        rm(aug_test_list)
        
        if(verbose) message("\tStep 1(c): Weighted estimation of survivor function Y | A, X, N on augmented dataset")
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
        
        # Desired output
        # 4 x array
        # row i, column j, slice k refers to the predicted survival of Y and C for ith participant (a, X_i, N_i)
        # at time j, for learner k in the augmented test dataset.
        # key summary is what the superlearner ensemble says -- 4 x matrix
        surv.fit <- .estimate.conditional.survival.2(Y=aug_train[["Y"]], 
                                                     Delta=aug_train[["DeltaY"]], 
                                                     A=aug_train[["A"]], 
                                                     W_event=aug_train %>% 
                                                       dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                     W_cens=aug_train %>% 
                                                       dplyr::select(colnames(W)),
                                                     newW_event = aug_test %>% 
                                                       dplyr::select(c("N", "DeltaN", colnames(W))), 
                                                     newW_cens = aug_test %>% 
                                                       dplyr::select(colnames(W)) %>%
                                                       distinct(),
                                                     fit.times=fit.times, 
                                                     fit.treat=fit.treat, 
                                                     event.SL.library=nuis$event.SL.library, 
                                                     cens.SL.library=nuis$cens.SL.library, 
                                                     survSL.control=nuis$survSL.control, 
                                                     survSL.cvControl = nuis$survSL.cvControl, 
                                                     cens.trunc=nuis$cens.trunc, 
                                                     save.fit = nuis$save.nuis.fits, 
                                                     verbose = nuis$verbose,
                                                     obsWeights = aug_train$w, 
                                                     id = aug_train$id)
        
        gc()
        # Now I have test-fold per-participant predicted survival of 
        # 1) Y conditional on A, X, N (over a huge grid) 
        # 2) C on A, X
        
        if(verbose) message("Step 2: Derive Backwards Conditional Density via Bayes rule on Test Fold")
        
        # Next I need to derive the reverse (posterior) density of N at each grid point
        
        fN <- matrix(data=NA, nrow=nrow(dat_test), ncol=length(fit.times))
        
        ### Loop through each row of augmented dataset (repeated participants, unique gridpoint N)
        
        if(verbose){
          pb <- progress_bar$new(total = nrow(dat_test),
                                 format = "Test fold [:bar] :percent ETA: :eta")
        }
        for(i in 1:nrow(dat_test)){

          Ai = dat_test[i,][["A"]]
          
          # return vector of survival probabilities for N
          if(Ai==1){
          SNi <- SN$survest_1(t=fit.times, x=dat_test[i,] %>% select(colnames(W)))
          }else{
          SNi <- SN$survest_0(t=fit.times, x=dat_test[i,] %>% select(colnames(W)))
          }
          
          # conditional probability of irrelevant infection failure at each gridpoint Nj
          LNi <- diff(1-c(1, SNi))
          
          Yi=dat_test[i,][["Y"]]
          DeltaYi=dat_test[i,][["DeltaY"]]
          idi = dat_test[i,][["id"]]
          
          #Ni=_test[i,][["N"]]
          #idi <- aug_test[i,][["id"]]
          indx_aug <- which(aug_test[["id"]]==idi)
          
          if(Ai==1 & DeltaYi==0){
            #Li = surv.fit$event.pred.1[i,which(fit.times==Yi)] #P(Y > \tilde{y} \mid A, X, N=n)
            Lvec = surv.fit$event.pred.1[indx_aug, which(fit.times==Yi)] #probs of failure at different N for fixed ppt and Y
            fN[i,] <- sapply(1:length(LNi), FUN=function(j){(Lvec[j] *LNi[j]) / sum(Lvec * LNi)})
          }else if(Ai==1 & DeltaYi==1){
            #P(Y=k \mid X) = P(Y > k-1 \mid X) - P(Y > k \mid X)
            Lvec = cbind(1, surv.fit$event.pred.1)[indx_aug,which(fit.times==Yi)] - cbind(1, surv.fit$event.pred.1)[indx_aug,which(fit.times==Yi)]
            fN[i,] <- sapply(1:length(LNi), FUN=function(j){(Lvec[j] *LNi[j]) / sum(Lvec * LNi)})
           }else if(Ai==0 & DeltaYi==0){
             Lvec = surv.fit$event.pred.0[indx_aug, which(fit.times==Yi)] #probs of failure at different N for fixed ppt and Y
             fN[i,] <- sapply(1:length(LNi), FUN=function(j){(Lvec[j] *LNi[j]) / sum(Lvec * LNi)})
          }else{
            Lvec = cbind(1, surv.fit$event.pred.0)[indx_aug,which(fit.times==Yi)] - cbind(1, surv.fit$event.pred.0)[indx_aug,which(fit.times==Yi)]
            fN[i,] <- sapply(1:length(LNi), FUN=function(j){(Lvec[j] *LNi[j]) / sum(Lvec * LNi)})
          }
          if(nuis$verbose){pb$tick()}
        }
        
        if(do.event.pred.0) nuis$event.pred.0 <- rbind(nuis$event.pred.0, surv.fit$event.pred.0)
        if(do.event.pred.1) nuis$event.pred.1 <- rbind(nuis$event.pred.1, surv.fit$event.pred.1)
        if(do.cens.pred.0) nuis$cens.pred.0 <- rbind(nuis$cens.pred.0, surv.fit$cens.pred.0)
        if(do.cens.pred.1) nuis$cens.pred.1 <- rbind(nuis$cens.pred.1, surv.fit$cens.pred.1)
        nuis$fN <- rbind(nuis$fN, fN)
        nuis$id <- c(nuis$id, aug_test$id)
        if(nuis$save.nuis.fits) result$surv.fits[[v]] <- surv.fit$surv.fit
        nuis$event.coef[[v]] <- surv.fit$event.coef
        nuis$cens.coef[[v]] <- surv.fit$cens.coef
        
        rm(surv.fit)
        gc()
        if(verbose){message(paste("Fold ", v, "Complete"))}
      }}else{
        stop("Please specify folds > 1 to ensure cross-fitting is used.")
      }
  }
  
  #### ESTIMATE CF SURVIVALS ####
  if(verbose) message("Computing counterfactual survivals...")
  
  # Key function required to compute survival
  # .get.survival.AIPWCC <- function(Y, DeltaY, N, DeltaN, A, fit.times, eval.times, S.hats, G.hats, g.hats, obsWeights=NULL, id=NULL, isotonize=TRUE, verbose=TRUE) {
  #   
  #   # Y=time
  #   # Delta = event
  #   # N = time_N
  #   # DeltaN = event_N
  #   # A = treat
  #   # eval.times = fit.times
  #   # S.hats = nuis$event.pred.1
  #   # G.hats = nuis$cens.pred.1
  #   # nuis$fN_1 <- matrix(nuis$fN, nrow = nrow(data), ncol=length(fit.times), byrow = FALSE)
  #   # nuis$fN_1 <- nuis$fN_1 / rowSums(nuis$fN_1)
  #   # g.hats <- prop.fit
  #   
  #   # Collect fitting times, order them, and order (conditional) survival and censoring probabilities
  #   fit.times <- fit.times[fit.times > 0]
  #   n <- length(Y)
  #   #N <- length(unique(id))
  #   ord <- order(fit.times)
  #   fit.times <- fit.times[ord]
  #   S.hats <- S.hats[,ord]
  #   G.hats <- G.hats[,ord]
  #   
  #   # This holds the per-id IF contribution
  #   IF.vals <- matrix(NA, nrow=n, ncol=length(fit.times))
  #   # This holds the final output of influence function contributions.
  #   surv <- rep(NA, length(fit.times))
  #   
  #   Y_ind = match(Y, fit.times)
  #   N_ind = match(N, fit.times)
  #   
  #   # loop through target times t0
  # 
  #   if(verbose){
  #     pb <- progress_bar$new(total = length(Y) * length(fit.times),
  #                            format = "Computing IFs [:bar] :percent ETA: :eta")
  #   }
  #   
  #   for(i in 1:length(Y)){
  #       #LOOP OVER INDIVIDUALS
  #       fN.i <- nuis$fN_1[i,] #(k x 1)conditional probabilities of N|A_i, X_i
  #       indx.hats <- which(nuis$id==unique(nuis$id)[i]) #(k x 1) index of S.hats
  #       
  #       S.hats.i <- S.hats[indx.hats,] #(kxk) rows index potential values of N, columns index Y
  #       G.hats.i <- G.hats[i,] #(k x 1) single row for single ppt, columns index times t
  #       
  #       # list of length k
  #       # each element contains (k x 1) vector of integral term
  #       # each list element indexes a value of N
  #       # each vector indexes values of Y
  #       vals.i <- t(sapply(1:nrow(S.hats.i), function(j){
  #         vals <- diff(1/S.hats.i[j,]) * (1/G.hats.i[-length(G.hats.i)]) #(k-1 x 1) vector
  #         if(any(fit.times[-1] > Y[i])) vals[fit.times[-1] > Y[i]] <- 0
  #         c(0,cumsum(vals)) #(kx1) vector
  #       }))
  #       # else return first inner term over all values over grid of N
  #       inner.func.2 <- vals.i
  #       
  #       # convert to kxk matrix
  #       # each row corresponds to value of N
  #       # each column corresponds to value of Y
  #       
  #       S.hats.y<-S.hats.i[,Y_ind[i]] #(k x 1) S.hats.Y over grid of N at Y
  #       G.hats.y<-G.hats.i[Y_ind[i]] #(scalar) G.hats.Y at Y
  #       
  #       for(t0 in eval.times) {
  #       # return index of smallest time greater than t0 (in fit times)
  #       k <- min(which(fit.times >= t0))
  # 
  #       S.hats.t0<-S.hats.i[,k] #(k x 1) S.hats.Y over grid of N at t0
  #       G.hats.t0<-G.hats.i[k] #(scalar) G.hats.Y at t0
  #       
  #       inner.func.1 <- ifelse(Y[i] <= t0 & Delta[i] == 1, 1/(S.hats.y * G.hats.y), 0) #return 0 if no event before t0
  #       if(length(inner.func.1)==1){
  #       inner.func.1 <- rep(inner.func.1, length(fit.times))
  #       }
  #       
  #       #list of length k 
  #       #each list element indexes values for N
  #       #each vector element indexes 
  #       IF.vals.i <- S.hats.t0 + as.numeric(A[i]==1) * S.hats.t0 * (-inner.func.1 + inner.func.2[,k])/g.hats[i]
  #       # S.hats.t0 denotes FIXED Y and variable N
  #       # inner.func.1 denotes fixed Y and variable N
  #       # inner.func.2 denotes fixed Y (column k) and variable n
  #       
  #       ########
  #       # NOW MOVE ON TO OBS IF calculation
  #       
  #       term1 <- (DeltaN[i]/G.hats.i[N_ind[i]]) * IF.vals.i[N_ind[i]]
  #       Q0N <- (sum(I(fit.times >= N[i]) * IF.vals.i * fN.i)/sum(I(fit.times >= N[i]) * fN.i))
  #       term2 <- (1-DeltaN[i])/G.hats.i[N_ind[i]] * Q0N
  #       
  #       term3 <- sum(sapply(
  #         1:length(fit.times[fit.times<=N[i]]), FUN=function(j){
  #         (sum(I(fit.times >= fit.times[j]) * IF.vals.i * fN.i)/sum(I(fit.times >= fit.times[j]) * fN.i)) * diff(1/c(1,G.hats.i))[j]
  #       }))
  #       
  #     IF.vals[i,k] <- term1 + term2 - term3
  #     if(verbose){pb$tick()}
  #     }
  #   }
  #   
  #   surv <- colMeans(IF.vals)
  #   
  #   IF.vals.center <- IF.vals - surv
  # 
  #   res <- list(times=eval.times, surv=pmin(1,pmax(0,surv)), IF.vals=IF.vals, IF.vals.center = IF.vals.center)
  #   if(isotonize) {
  #     res$surv.iso <- NA
  #     res$surv.iso[!is.na(res$surv)] <- 1 - isoreg(res$times[!is.na(res$surv)], 1-res$surv[!is.na(res$surv)])$yf
  #   }
  #   
  #   return(res)
  # }
  
  #### THIS WAS FAST BUT RETURNED SOMETHING WRONG
  # NEED TO REVISIT
  .get.survival.AIPWCC_fast <- function(
    Y, DeltaY, N, DeltaN, A,
    fit.times, eval.times,
    S.hats, G.hats,
    g.hats, obsWeights,
    nuis,
    isotonize = TRUE,
    verbose = TRUE
  ) {
    
    # Y=time;
    # DeltaY = event;
    # N = time_N;
    # A = 1-treat;
    # DeltaN = event_N;
    # fit.times = fit.times;
    # eval.times = fit.times;
    # S.hats = nuis$event.pred.0;
    # G.hats = nuis$cens.pred.0;
    # g.hats = 1-prop.fit;
    # nuis=nuis;
    # verbose=TRUE
    
    fit.times <- fit.times[fit.times > 0]
    n  <- length(Y)
    K  <- length(fit.times)
    
    if(verbose){
      pb <- progress_bar$new(total = n,
                             format = "Computing survival [:bar] :percent ETA: :eta")
    }
    
    ord <- order(fit.times)
    eval.times <- fit.times[ord]
    S.hats <- S.hats[, ord, drop = FALSE]
    G.hats <- G.hats[, ord, drop = FALSE]
    
    Y_ind <- sapply(Y, function(y)
      which(fit.times >= y)[1]
    )
    N_ind <- sapply(N, function(n)
      which(fit.times >= n)[1]
    )
    
    IF.vals <- matrix(NA_real_, n, K)
    
    unique_ids <- unique(nuis$id)
    
    for(i in seq_len(n)) {
      
      fN.i <- nuis$fN[i, ] #select conditional distribution of N for participant i
      
      indx.hats <- which(nuis$id == unique_ids[i]) #select rows that correspond to unique id
      S.hats.i  <- S.hats[indx.hats, , drop = FALSE] #select survival of Y (kxk matrix)
      # each row refers to a grid point of N, each column refers to a new value for Y
      G.hats.i  <- G.hats[i, ] #select censoring survival
      
      ## =============================
      ## Discretized compensator
      ## =============================
      
      S_inv  <- 1 / S.hats.i
      dS_inv <- S_inv[, -1] - S_inv[, -ncol(S_inv)] #selects all columns but first minus all columns minus last
      # check for first part of int.vals
      
      vals.i <- sweep(dS_inv, 2, (1 / G.hats.i[-length(G.hats.i)]), FUN="*") #I fixed and validated this code chunk
      #takes each row of dS_inv and multiplies it by inverse of G.hats.i **elementwise**
      
      # zero beyond Y[i]
      vals.i[, eval.times[-1] > Y[i]] <- 0
      
      vals.i <- cbind(0, t(apply(vals.i, 1, cumsum)))
      # inner command `apply(vals.i, 1, cumsum)` calculates cumulative sum over each row, 
      # then transposes to ensure that the rowwise cumsums are represented per row
      # then cbinds 0 to ensure margin
      
      #for particular value of Y=Yi, select 
      S_y <- S.hats.i[, Y_ind[i]] #(Kx1) vector associated with each N
      G_y <- G.hats.i[Y_ind[i]] #scalar
      
      dG_inv <- diff(1 / c(1, G.hats.i)) #discrete derivative of 1/G.hats (measure to integrate against)
      
      #loop through length of fit.times
      for(k in seq_len(K)) {
        
        S_t0 <- S.hats.i[, k] #return predicted survival at that time Y
        
        ## CORRECT time-dependent inner.func.1
        if (Y[i] <= fit.times[k] && DeltaY[i] == 1) {
          inner1 <- 1 / (S_y * G_y) #return (Kx1) vector of inner terms if satisfied
        } else {
          inner1 <- rep(0, length(S_y)) #return (Kx1) vector of 0s if not satisfied
        }
        
        IF_k <- S_t0 - ((as.numeric(A[i] == 1)/ g.hats[i]) * S_t0 * (inner1 - vals.i[, k]) )
        
        ### Full-data IF complete (K x 1) vector
        
        Ni <- N_ind[i] #index for N value
        
        term1 <- (DeltaN[i] / G.hats.i[Ni]) * IF_k[Ni] #if observed, return IPW term
        
        ## fast cumulative Q sequence
        w <- IF_k * fN.i #product of fN.i and influence function
        denom_cum <- rev(cumsum(rev(fN.i)))
        num_cum   <- rev(cumsum(rev(w)))
        # these chunks above compute the tail sums
        # e.g., we want sum(I(t >= fit.times[j]) * fN)
        # rev(fN.i) reports f_4, f_3, f_2, f_1
        # cumsum(rev(fN.i)) reports f_4, f_4+f_3, f_4+f_3+f_2, f_4+f_3+f_2+f_1 (the tail sums from left)
        # rev(cumsum(rev(fN.i))) reports (f_4+f_3+f_2+f_1), (f_4+f_3+f_2), (f_4+f_3), f_4
        # this exactly reflects our desired target
        Q_seq <- num_cum / denom_cum ###this gives Q_r for all r
        
        term2 <- (1 - DeltaN[i]) / G.hats.i[Ni] * Q_seq[Ni]
        
        term3 <- sum(Q_seq[1:Ni] * dG_inv[1:Ni])
        
        IF.vals[i, k] <- term1 + term2 - term3
      }
      
      if(verbose) pb$tick()
    }
    
    surv <- colMeans(IF.vals)
    
    IF.vals.center <- sweep(IF.vals, 2, surv, "-")
    
    res <- list(
      times = fit.times,
      surv  = pmin(1, pmax(0, surv)),
      IF.vals = IF.vals,
      IF.vals.center = IF.vals.center
    )
    
    if(isotonize) {
      idx <- !is.na(res$surv)
      res$surv.iso <- res$surv
      res$surv.iso[idx] <-
        1 - isoreg(res$times[idx],
                   1 - res$surv[idx])$yf
    }
    
    return(res)
  }
  
  if(0 %in% fit.treat) {
    
    if(verbose){message("Computating Survival Probability for A=0")}
    
    surv.0 <- .get.survival.AIPWCC_fast(
      Y=time,
      DeltaY = event,
      N = time_N,
      A = 1-treat,
      DeltaN = event_N,
      fit.times = fit.times,
      eval.times = fit.times,
      S.hats = nuis$event.pred.0,
      G.hats = nuis$cens.pred.0,
      g.hats = 1-prop.fit,
      nuis=nuis, 
      verbose=TRUE
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
      Y=time,
      DeltaY = event,
      N = time_N,
      A = treat,
      DeltaN = event_N,
      fit.times = fit.times,
      eval.times = fit.times,
      S.hats = nuis$event.pred.1,
      G.hats = nuis$cens.pred.1,
      g.hats = prop.fit,
      nuis=nuis, 
      verbose=TRUE
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
  