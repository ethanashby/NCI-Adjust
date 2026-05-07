#' Super Learner for conditional survival functions (EDIT)
#'
#' This function estimates conditional survival functions for the event and censoring times from right-censored data. A key difference between this version and the version developed by Ted Westling is that the censoring and survival functions can depend on different covariates.
#'
#' The conditional survival function of the event at time \code{t} given covariates \code{X_1} is defined as the probability that the event occurs after time \code{t} given covariate values \code{x_1}. The conditional survival function of censoring is the probability that the censoring time occurs after \code{t} given covariates \code{x_2}. This function finds the optimal weighted combination, i.e. the Super Learner, of candidate learners for both of these functions simultaneously.
#'
#' @param time \code{n x 1} numeric vector of observed right-censored follow-up times; i.e. the minimum of the event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of whether an event was observed.
#' @param X_event \code{n x p1} data.frame of observed covariate values on which to train the SuperLearner for the survival function for the event process.
#' @param X_cens \code{n x p2} data.frame of observed covariate values on which to train the SuperLearner for the censoring process.
#' @param newX_event \code{m1 x p} data.frame of new observed covariate values at which to obtain predictions for the estimated algorithm. Must have the same names and structure as \code{X_event}.
#' @param newX_cens \code{m2 x p} data.frame of new observed covariate values at which to obtain predictions for the estimated algorithm. Must have the same names and structure as \code{X_cens}.
#' @param new.times \code{k x 1} numeric vector of times at which to obtain predicted conditional survivals.
#' @param event.SL.library Library of candidate learners to use to estimate the conditional survival of the event. Should have the same structure as the \code{SL.library} argument to the \code{SuperLearner} function in the \code{SuperLearner} package; see details below. Run \code{\link{survlistWrappers()}} to see a list of currently available prediction and screening algorithms.
#' #' @param cens.SL.library Library of candidate learners to use to estimate the conditional survival of censoring.
#' @param id Optional \code{n x 1} vector of observation clusters. If provided, cross-validation folds will respect clustering and \code{id} will be passed to every learner, though some learners may not make use of it. Default is every observation in its own cluster; i.e. iid observations.
#' @param verbose \code{TRUE/FALSE} indicating whether to print progress in the fitting process to the console.
#' @param control Named list of parameters controlling the fitting process. See \code{\link{survSuperLearner.control}} for details.
#' @param cvControl Named list of parameters controlling the cross-validation process. See \code{\link{survSuperLearner.cvControl}} for details.
#' @param obsWeights Optional \code{n x 1} vector of observation weights. If provided, these weights will be passed to each learner, which may or may not make use of them (or make use of them correctly), and will be used in the ensemble step to weight the empirical risk function.
#' @return \code{survSuperLearner} returns a named list with the following elements:
#' \item{call}{The matched call.}
#' \item{event.libraryNames, cens.libraryNames}{Parsed learner names.}
#' \item{event.SL.library, cens.SL.library}{Libraries used for fitting.}
#' \item{event.SL.predict, cens.SL.predict}{\code{m x k} matrices of SuperLearner predicted survival values. Rows index observations in \code{newX}; columns index times in \code{new.times.}}
#' \item{event.coef, cens.coef}{Fitted SuperLearner coefficients for the model for the conditional survival functions for the event and censoring times, respectively.}
#' \item{event.library.predict, cens.library.predict}{\code{m x k x p} predicted event and censoring survivals on \code{newX} and \code{new.times} from the candidate learners, where \code{p} is the number of candidated learners.}
#' \item{event.Z, cens.Z}{\code{n x l x p} cross-validated event and censoring survivals on the training data, where \code{l} is the number of elements in \code{control$event.t.grid} and \code{control$cens.t.grid}, respectively, and \code{p} is the number of candidate learners.}
#' \item{event.cvRisk, cens.cvRisk}{Cross-validated risks for the candidate conditional event and censoring survival functions.}
#' \item{event.fitLibrary, cens.fitLibrary}{Fitted conditional survival functions for all learners in the library on the full data.}
#' \item{varNames}{Variable names of the training data.}
#' \item{validRows}{Length \code{V} list containing the indices contained in each fold used for cross-validation.}
#' \item{event.whichScreen, cens.whichScreen}{Matrix indicating which variables were included in each screening algorithm in the full training data.}
#' \item{control, cvControl}{Parameters used for controlling the fitting and cross-validation processes, respectively.}
#' \item{event.errorsInCVLibrary, cens.errorsInCVLibrary}{Logical matrices indicating whether each learning algorithm encountered any errors in each cross-validation fold.}
#' \item{event.errorsInLibrary, cens.errorsInLibrary}{Logical vectors indicating whether each learning algorithm encountered any errors on the full data.}
#' \item{times}{Timing data.}
#' @references van der Laan, M. J., Polley, E. C., & Hubbard, A. E. (2007). Super learner. \emph{Statistical Applications in Genetics and Molecular Biology}, 6(1).
#' @references van der Laan, M. J., and Rose, S. (2011). \emph{Targeted Learning: Causal inference for observational and experimental data}. Springer-Verlag New York.
#' @examples

survSuperLearner2 <- function(time, event, X_event, X_cens, newX_event = NULL, newX_cens = NULL, new.times, event.SL.library, cens.SL.library, id = NULL, verbose = FALSE, control = list(), cvControl = list(), obsWeights = NULL)  {
  
  # Check to see if required packages are installed
  # packages <- c("survival", "mgcv", "randomForestSRC", "glmnet", "SuperLearner",
  #               "nnls", "Rsolnp")
  # for (pkg in packages) {
  #   if (!requireNamespace(pkg, quietly=TRUE)) {
  #     stop(paste0("Package '", pkg, "' is required by survSuperLearner."))
  #   }
  # }
  
  time_start <- proc.time()
  call <- match.call(expand.dots = TRUE)
  if (!is.null(dim(time)) && ncol(time) > 0) stop("time must be an (n x 1)  numeric vector.")
  if (!is.null(dim(event)) && ncol(event) > 0) stop("event must be an (n x 1) numeric vector.")
  time <- as.numeric(time)
  event <- as.numeric(event)
  if (is.null(obsWeights)) obsWeights <- rep(1, length(time))
  if (is.null(newX_event)) newX_event <- X_event
  if (is.null(newX_cens)) newX_cens <- X_cens
  
  varNames_event <- colnames(X_event)
  varNames_cens <- colnames(X_cens)
  N <- dim(X_event)[1L]
  p_event <- dim(X_event)[2L]
  p_cens <- dim(X_cens)[2L]
  
  .checkInputs(time=time, event=event, X_event = X_event, X_cens = X_cens, newX_event=newX_event, newX_cens=newX_cens, id=id, obsWeights=obsWeights, verbose=verbose)
  
  if (is.null(control$event.t.grid)) {
    control$event.t.grid <- seq(0, max(time[event == 1]), length.out = 250)
  } else {
    if (!is.null(dim(control$event.t.grid)) && ncol(control$event.t.grid) > 0) stop("t.grid must be an (k x 1)  numeric vector.")
    control$event.t.grid <- sort(unique(as.numeric(control$event.t.grid)))
    if(any(is.na(control$event.t.grid))) stop("No missing values allowed in t.grid")
    if (any(control$event.t.grid < 0)) stop("Values in t.grid must be non-negative")
    if (any(control$event.t.grid > max(time))) stop("Values in t.grid must be no large than max(time)")
  }
  
  if (is.null(control$cens.t.grid)) {
    control$cens.t.grid <- seq(0, max(time[event == 0]), length.out = 250)
    control$cens.t.grid <- control$cens.t.grid - max(min(diff(sort(unique(time)))) / 2, 1e-5)
    control$cens.t.grid <- c(0, control$cens.t.grid[control$cens.t.grid > 0])
  } else {
    if (!is.null(dim(control$cens.t.grid)) && ncol(control$cens.t.grid) > 0) stop("t.grid must be an (k x 1)  numeric vector.")
    control$cens.t.grid <- sort(unique(as.numeric(control$cens.t.grid)))
    if(any(is.na(control$cens.t.grid))) stop("No missing values allowed in t.grid")
    if (any(control$cens.t.grid < 0)) stop("Values in t.grid must be non-negative")
    if (any(control$cens.t.grid > max(time))) stop("Values in t.grid must be no large than max(time)")
  }
  
  control <- do.call("survSuperLearner.control", control)
  
  cvControl <- do.call("survSuperLearner.CV.control", cvControl)
  
  event.library <- .createLibrary(event.SL.library)
  cens.library <- .createLibrary(cens.SL.library)
  
  event.k <- nrow(event.library$library)
  cens.k <- nrow(cens.library$library)
  event.kScreen <- length(event.library$screenAlgorithm)
  cens.kScreen <- length(cens.library$screenAlgorithm)
  ### These are the steps that are leading to computational exhaustion
  event.Z <- array(NA, dim = c(N, length(control$event.t.grid), event.k))
  cens.Z <- array(NA, dim = c(N, length(control$cens.t.grid), cens.k))
  event.libraryNames <- cbind(predAlgorithm = event.library$library$predAlgorithm,
                              screenAlgorithm = event.library$screenAlgorithm[event.library$library$rowScreen])
  cens.libraryNames <- cbind(predAlgorithm = cens.library$library$predAlgorithm,
                             screenAlgorithm = cens.library$screenAlgorithm[cens.library$library$rowScreen])
  if (p_event < 2 & !identical(event.library$screenAlgorithm, "All")) {
    warning("Screening algorithms specified in combination with single-column X.")
  }
  
  if (p_cens < 2 & !identical(cens.library$screenAlgorithm, "All")) {
    warning("Screening algorithms specified in combination with single-column X.")
  }
  
  event.errorsInCVLibrary <- event.errorsInLibrary <- rep(0, event.k)
  cens.errorsInCVLibrary <- cens.errorsInLibrary <- rep(0, cens.k)
  
  validRows <- .survCVFolds2(N = N, id = id, event = event, cvControl = cvControl)
  
  if (is.null(id)) id <- seq(N)
  
  ### Get cross-validated survivals for S and T
  
  .crossValFUN <- function(valid, time, event, dataX, id, obsWeights, t.grid,
                           library, kScreen, k, p) {
    foldNum <- as.numeric(which(unlist(lapply(validRows, function(v) all.equal(v, valid) == TRUE))))
    tempLearn <- dataX[-valid, , drop = FALSE]
    tempTime <- time[-valid]
    tempEvent <- event[-valid]
    tempValid <- dataX[valid, , drop = FALSE]
    tempWhichScreen <- matrix(NA, nrow = kScreen, ncol = p)
    tempId <- id[-valid]
    tempObsWeights <- obsWeights[-valid]
    for (s in seq(kScreen)) {
      if(verbose) message(paste("CV ", library$screenAlgorithm[s],
                                ", fold ", foldNum, sep = ""))
      screen_fn <- get(library$screenAlgorithm[s])
      testScreen <- try(do.call(screen_fn, list(time = tempTime,
                                                event = tempEvent,
                                                X = tempLearn, id = tempId,
                                                obsWeights = tempObsWeights)))
      if (inherits(testScreen, "try-error")) {
        warning(paste("replacing failed screening algorithm,",
                      library$screenAlgorithm[s], ", with All()",
                      "\n "))
        tempWhichScreen[s, ] <- TRUE
      }
      else {
        tempWhichScreen[s, ] <- testScreen
      }
      if (verbose) {
        message(paste("Number of covariates in ", library$screenAlgorithm[s],
                      " is: ", sum(tempWhichScreen[s, ]), sep = ""))
      }
    }
    
    uniqueScreen <- unique(tempWhichScreen)
    screenMap <- apply(uniqueScreen, 1, function(row) which(apply(tempWhichScreen, 1, function(row2) all.equal(row, row2) == TRUE)))
    
    out <- array(NA, dim = c(nrow(tempValid), length(t.grid), k))
    
    for (predAlg in unique(library$library$predAlgorithm)) {
      if (verbose) message(paste("CV ", predAlg, ", fold ", foldNum, sep = ""))
      pred_fn <- get(predAlg)
      for(j in seq(nrow(uniqueScreen))) {
        testAlg <- try(do.call(pred_fn, list(time = tempTime, event = tempEvent,
                                             X = subset(tempLearn, select = uniqueScreen[j,], drop = FALSE),
                                             newX = subset(tempValid, select = uniqueScreen[j,], drop = FALSE),
                                             new.times = t.grid,
                                             id = tempId,
                                             obsWeights = tempObsWeights)))
        if (inherits(testAlg, "try-error")) {
          warning(paste("Error in algorithm", predAlg,
                        "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
        } else {
          libraryRows <- which(library$library$predAlgorithm == predAlg & library$library$rowScreen %in% unlist(screenMap[j]))
          for (row in libraryRows) {
            out[,,row] <- testAlg$pred
          }
        }
      }
    }
    
    
    invisible(list(out = out))
  }
  
  time_train_start <- proc.time()
  event.crossValFUN_out <- lapply(validRows, FUN = .crossValFUN,
                                  time = time, event = event, dataX = X_event, id = id, obsWeights = obsWeights,
                                  t.grid = control$event.t.grid, library = event.library,
                                  kScreen = event.kScreen, k = event.k, p = p_event)
  
  for(v in 1:cvControl$V) {
    event.Z[validRows[[v]], ,] <- event.crossValFUN_out[[v]]$out
  }
  
  event.errorsInCVLibrary <- apply(event.Z, 3, anyNA)
  if (sum(event.errorsInCVLibrary) > 0) event.Z[, , as.logical(event.errorsInCVLibrary)] <- 0
  if (all(event.Z == 0))  stop("All algorithms dropped from event library")
  
  cens.crossValFUN_out <- lapply(validRows, FUN = .crossValFUN,
                                 time = time, event = 1-event, dataX = X_cens, id = id, obsWeights = obsWeights,
                                 t.grid = control$cens.t.grid, library = cens.library,
                                 kScreen = cens.kScreen, k = cens.k, p = p_cens)
  
  for(v in 1:cvControl$V) {
    cens.Z[validRows[[v]], ,] <- cens.crossValFUN_out[[v]]$out
  }
  
  cens.errorsInCVLibrary <- apply(cens.Z, 3, anyNA)
  if (sum(cens.errorsInCVLibrary) > 0) cens.Z[, , as.logical(cens.errorsInCVLibrary)] <- 0
  if (all(cens.Z == 0))  stop("All algorithms dropped from censoring library")
  
  ###  Do iterative SuperLearning to get coefficients
  # This is what originally exhausted vector memory
  # Had to edit to stream in data rather than include everything in memory
  getCoef <- .surviterativeSL2(event.Z = event.Z, cens.Z = cens.Z, time = time, event = event, X_event = X_event,
                               X_cens = X_cens, obsWeights = obsWeights, id = id, control = control, verbose = verbose,
                               event.errorsInLibrary = event.errorsInCVLibrary, cens.errorsInLibrary = cens.errorsInCVLibrary)
  event.coef <- getCoef$event.coef
  cens.coef <- getCoef$cens.coef
  
  event.cvRisks <- getCoef$event.cvRisks
  cens.cvRisks <- getCoef$cens.cvRisks
  
  names(event.coef) <- names(event.cvRisks) <- apply(event.libraryNames, 1, paste, collapse = "_")
  names(cens.coef) <- names(cens.cvRisks) <- apply(cens.libraryNames, 1, paste, collapse = "_")
  
  time_train <- proc.time() - time_train_start
  
  ### Get predictions
  
  time_predict_start <- proc.time()
  .screenFun <- function(fun, list) {
    screen_fn <- get(fun)
    testScreen <- try(do.call(screen_fn, list))
    if (inherits(testScreen, "try-error")) {
      warning(paste("replacing failed screening algorithm,",
                    fun, ", with All() in full data", "\n "))
      out <- rep(TRUE, ncol(list$X))
    }
    else {
      out <- testScreen
    }
    return(out)
  }
  
  event.whichScreen <- sapply(event.library$screenAlgorithm, FUN = .screenFun,
                              list = list(time = time, event = event, id = id, X = X_event, obsWeights = obsWeights),
                              simplify = FALSE)
  event.whichScreen <- do.call(rbind, event.whichScreen)
  
  cens.whichScreen <- sapply(cens.library$screenAlgorithm, FUN = .screenFun,
                             list = list(time = time, event = 1-event, id = id, X = X_cens, obsWeights = obsWeights),
                             simplify = FALSE)
  cens.whichScreen <- do.call(rbind, cens.whichScreen)
  
  assign("event.fitLibrary", vector("list", length = event.k))
  names(event.fitLibrary) <- apply(event.libraryNames, 1, paste, collapse = "_")
  
  assign("cens.fitLibrary", vector("list", length = cens.k))
  names(cens.fitLibrary) <- apply(cens.libraryNames, 1, paste, collapse = "_")
  
  .predFun <- function(index, lib, time, event, dataX, newX, whichScreen, t.grid,
                       family, id, obsWeights, verbose, control, libraryNames) {
    if (verbose) {
      message(paste("full", libraryNames[index]))
    }
    pred_fn <- get(lib$predAlgorithm[index])
    testAlg <- try(do.call(pred_fn, list(time = time, event = event,
                                         X = subset(dataX, select = whichScreen[lib$rowScreen[index], ],
                                                    drop = FALSE),
                                         newX = subset(newX, select = whichScreen[lib$rowScreen[index],
                                         ], drop = FALSE), id = id,
                                         obsWeights = obsWeights, new.times = t.grid)))
    if (inherits(testAlg, "try-error")) {
      warning(paste("Error in algorithm", lib$predAlgorithm[index],
                    " on full data", "\n  The Algorithm will be removed from the Super Learner (i.e. given weight 0) \n"))
      out <- rep.int(NA, times = nrow(newX))
      model_out <- NULL
    }
    else {
      out <- testAlg$pred
      if(control$saveFitLibrary) model_out <- testAlg$fit
      else model_out <- NULL
    }
    
    invisible(list(out = out, model_out = model_out))
  }
  
  event.pred <- lapply(seq(event.k), FUN = .predFun,
                       lib = event.library$library, time = time, event = event,
                       dataX = X_event, newX = newX_event, t.grid = new.times,
                       whichScreen = event.whichScreen, id = id,
                       obsWeights = obsWeights, verbose = verbose, control = control,
                       libraryNames = event.libraryNames)
  
  event.libraryPred <-  array(NA, dim = c(nrow(newX_event), length(new.times), event.k))
  
  for(j in 1:event.k) {
    event.libraryPred[, , j] <- event.pred[[j]]$out
  }
  
  event.errorsInLibrary <- apply(event.libraryPred, 3, function(algorithm) anyNA(algorithm))
  
  if(control$saveFitLibrary) {
    event.fitLibrary <- lapply(event.pred, "[[", "model_out")
    names(event.fitLibrary) <- apply(event.libraryNames, 1, paste, collapse = "_")
  } else {
    event.fitLibrary <- NULL
  }
  
  cens.pred <- lapply(seq(cens.k), FUN = .predFun,
                      lib = cens.library$library, time = time, event = 1 - event,
                      dataX = X_cens, newX = newX_cens, t.grid = new.times,
                      whichScreen = cens.whichScreen, id = id,
                      obsWeights = obsWeights, verbose = verbose, control = control,
                      libraryNames = cens.libraryNames)
  
  cens.libraryPred <-  array(NA, dim = c(nrow(newX_cens), length(new.times), cens.k))
  
  for(j in 1:cens.k) {
    cens.libraryPred[, , j] <- cens.pred[[j]]$out
  }
  
  cens.errorsInLibrary <- apply(cens.libraryPred, 3, function(algorithm) anyNA(algorithm))
  
  if(control$saveFitLibrary) {
    cens.fitLibrary <- lapply(cens.pred, "[[", "model_out")
    names(cens.fitLibrary) <- apply(cens.libraryNames, 1, paste, collapse = "_")
  } else {
    cens.fitLibrary <- NULL
  }
  
  
  if(event.k == 1) {
    event.SL.predict <- event.libraryPred[,,1]
  } else {
    event.SL.predict <- matrix(NA, nrow = nrow(newX_event), ncol = length(new.times))
    for(j in seq(length(new.times))) event.SL.predict[,j] <- event.libraryPred[,j,!event.errorsInLibrary] %*% event.coef[!event.errorsInLibrary]
  }
  
  if(cens.k == 1) {
    cens.SL.predict <- cens.libraryPred[,,1]
  } else {
    cens.SL.predict <- matrix(NA, nrow = nrow(newX_cens), ncol = length(new.times))
    for(j in seq(length(new.times))) cens.SL.predict[,j] <- cens.libraryPred[,j,!cens.errorsInLibrary] %*% cens.coef[!cens.errorsInLibrary]
  }
  
  time_predict <- proc.time() - time_predict_start
  
  time_end <- proc.time()
  times <- list(everything = time_end - time_start, train = time_train,
                predict = time_predict)
  
  out <- list(call = call,
              event.libraryNames = event.libraryNames, cens.libraryNames = cens.libraryNames,
              event.SL.library = event.library, cens.SL.library = cens.library,
              event.SL.predict = event.SL.predict, cens.SL.predict = cens.SL.predict,
              event.coef = event.coef, cens.coef = cens.coef,
              event.library.predict = event.libraryPred, cens.library.predict = cens.libraryPred,
              event.Z = event.Z, cens.Z = cens.Z,
              event.cvRisk = event.cvRisks, cens.cvRisk = cens.cvRisks,
              event.fitLibrary = event.fitLibrary, cens.fitLibrary = cens.fitLibrary,
              varNames_event = varNames_event, varNames_cens = varNames_cens, validRows = validRows,
              event.whichScreen = event.whichScreen, cens.whichScreen = cens.whichScreen,
              control = control, cvControl = cvControl,
              event.errorsInCVLibrary = event.errorsInCVLibrary, cens.errorsInCVLibrary = cens.errorsInCVLibrary,
              event.errorsInLibrary = event.errorsInLibrary, cens.errorsInLibrary = cens.errorsInLibrary,
              times = times)
  class(out) <- c("survSuperLearner")
  return(out)
}

#' Prediction function for survival Super Learner
#'
#' This function predicts the fitted survival Super Learner on new data.
#'
#' This is the prediction function for \code{\link{survSuperLearner2}}.
#'
#' @param object Fitted \code{survSuperLearner} object.
#' @param newdata \code{r x p} data.frame of new covariates for which to obtain predicted survivals. Must have same names and structure as original covariates \code{X} used to train the SuperLearner.
#' @param new.times \code{q x 1} vector of times for which to obtain predicted survivals.
#' @param X_event Original training covariates, which may be needed for obtaining predictions for some learners. Defaults to \code{NULL}.
#' @param X_cens Original training covariates, which may be needed for obtaining predictions for some learners. Defaults to \code{NULL}.
#' @param time Original training follow-up times, which may be needed for obtaining predictions for some learners. Defaults to \code{NULL}.
#' @param event Original training event indicators, which may be needed for obtaining predictions for some learners. Defaults to \code{NULL}.
#' @param onlySL Logical indicating whether to only fit the learners with non-zero SuperLearner coefficient. Defaults to \code{FALSE}.
#' @param threshold Minimum coefficient weight for a learner to be inclued in the super learner prediction. Helps reduce computation time if there are many learners with small weights. Defaults to \code{10e-4}.
#' @return \code{survSuperLearner} returns a named list with the following elements:
#' \item{call}{The matched call.}
#' \item{event.libraryNames, cens.libraryNames}{Parsed learner names.}
#' \item{event.SL.library, cens.SL.library}{Libraries used for fitting.}
#' \item{event.SL.predict, cens.SL.predict}{\code{m x k} matrices of SuperLearner predicted survival values. Rows index observations in \code{newX}; columns index times in \code{new.times.}}
#' \item{event.coef, cens.coef}{Fitted SuperLearner coefficients for the model for the conditional survival functions for the event and censoring times, respectively.}

# predict.survSuperLearner <- function (object, newdata_event, newdata_cens, new.times, X_event = NULL, X_cens = NULL, time = NULL, event = NULL, onlySL = FALSE, threshold = 1e-4, ...) {
#   if (missing(newdata)) {
#     out <- list(survpred = object$survSL.predict, cens.pred = object$cens.SL.predict,
#                 survlibrary.predict = object$survlibrary.predict, cens.library.predict = object$cens.library.predict)
#     return(out)
#   }
#   if (!object$control$saveFitLibrary) {
#     stop("This SuperLearner fit was created using control$saveFitLibrary = FALSE, so new predictions cannot be made.")
#   }
#   event.k <- nrow(object$event.libraryNames)
#   event.pred <- array(0, dim=c(nrow(newdata), length(new.times), event.k))
#   dimnames(event.pred)[[3]] <- apply(object$event.libraryNames, 1, paste, collapse = '_')
#   dimnames(event.pred)[[2]] <- new.times
#   
#   cens.k <- nrow(object$cens.libraryNames)
#   cens.pred <- array(0, dim=c(nrow(newdata), length(new.times), cens.k))
#   dimnames(cens.pred)[[3]] <- apply(object$cens.libraryNames, 1, paste, collapse = '_')
#   dimnames(cens.pred)[[2]] <- new.times
#   if (onlySL) {
#     event.whichLibrary <- which(object$event.coef > threshold)
#     event.coef <- object$event.coef
#     event.coef[-event.whichLibrary] <- 0
#     event.coef <- event.coef / sum(event.coef)
#     
#     cens.whichLibrary <- which(object$cens.coef > threshold)
#     cens.coef <- object$cens.coef
#     cens.coef[-cens.whichLibrary] <- 0
#     cens.coef <- cens.coef / sum(cens.coef)
#   } else {
#     event.whichLibrary <- seq(event.k)
#     event.coef <- object$event.coef
#     cens.whichLibrary <- seq(cens.k)
#     cens.coef <- object$cens.coef
#   }
#   for (mm in event.whichLibrary) {
#     newdataMM <- subset(newdata_event,
#                         select = object$event.whichScreen[object$event.SL.library$library[mm, 2], ])
#     XMM <- if (is.null(X_event)) {
#       NULL
#     } else {
#       subset(X_event, select = object$event.whichScreen[object$event.SL.library$library[mm, 2], ])
#     }
#     event.pred[, ,mm] <- do.call("predict", list(object = object$event.fitLibrary[[mm]],
#                                                  newX = newdataMM, new.times = new.times,
#                                                  X = XMM, time = time, event = event, ...))
#   }
#   event.SL.predict <- matrix(NA, nrow = nrow(newdata_event), ncol=length(new.times))
#   for (j in seq(length(new.times))) {
#     event.SL.predict[,j] <- event.pred[,j,] %*% event.coef
#   }
#   
#   for (mm in cens.whichLibrary) {
#     newdataMM <- subset(newdata_cens,
#                         select = object$cens.whichScreen[object$cens.SL.library$library[mm, 2], ])
#     XMM <- if (is.null(X_cens)) {
#       NULL
#     } else {
#       subset(X_cens, select = object$cens.whichScreen[object$cens.SL.library$library[mm, 2], ])
#     }
#     cens.pred[, ,mm] <- do.call("predict", list(object = object$cens.fitLibrary[[mm]],
#                                                 newX = newdataMM, new.times = new.times,
#                                                 X = XMM, time = time, event = event, ...))
#   }
#   cens.SL.predict <- matrix(NA, nrow = nrow(newdata), ncol=length(new.times))
#   for (j in seq(length(new.times))) {
#     cens.SL.predict[,j] <- cens.pred[,j,] %*% cens.coef
#   }
#   
#   out <- list(event.SL.predict = event.SL.predict, event.library.predict = event.pred,
#               cens.SL.predict = cens.SL.predict, cens.library.predict = cens.pred)
#   return(out)
# }

# Functionn which checks inputs to survSuperLearner2
.checkInputs <- function(time, event, X_event, newX_event, X_cens, newX_cens, id, obsWeights, verbose) {
  if(any(time < 0)) stop("Only non-negative event/censoring times allowed!")
  if(any(!(event %in% c(0,1)))) stop("Event must be binary.")
  if(any(is.na(time)) | any(is.na(event))) stop("No missing values allowed in time or event.")
  if(any(is.na(X_event)) | any(is.na(newX_event)) | any(is.na(X_cens)) | any(is.na(newX_cens))) stop("No missing values allowed in X_event, X_cens, new X_event, or new X_cens.")
  if(length(time) != length(event) | length(time) != nrow(X_event) | length(time) != nrow(X_cens)) stop("time and event must be n x 1 vectors and X_event & X_cens must have n rows.")
  if(!is.data.frame(X_event) | !is.data.frame(newX_event)) stop("X_event and newX_event must be data frames.")
  if(!is.data.frame(X_cens) | !is.data.frame(newX_cens)) stop("X_cens and newX_cens must be data frames.")
  if(!identical(names(X_event), names(newX_event))) stop("X_event and newX_event must have the same features.")
  if(!identical(names(X_cens), names(newX_cens))) stop("X_event and newX_event must have the same features.")
  if(any(obsWeights < 0)) stop("obsWeights < 0 not allowed.")
  if(!(verbose %in% c(TRUE, FALSE))) stop("verbose must be TRUE/FALSE.")
  if(!is.null(id) && !identical(length(id), length(time))) stop("id vector must have the same dimension as time")
}

#' Control parameters for the survival Super Learner
#'
#' This function initiates control parameters for the \code{\link{survSuperLearner}} function.
#'
#' @param initWeightAlg Algorithm to use for the first iteration of the iterative SuperLearner algorithm. Defaults to \code{survSL.rfsrc}
#' @param initWeight Whether to start the iterative SuperLearner by fitting censoring weights (\code{"censoring"}, default) or by fitting event weights (\code{"event"}).
#' @param max.SL.iter Maximum iterations of the iterative SuperLearner algorithm. Defaults to 20.
#' @param event.t.grid Grid of times to use to approximate the integral in the risk function for the conditional survival function of the event. Defaults to 250 points equally spaced between 0 and the last uncensored follow-up time.
#' @param cens.t.grid Grid of times to use to approximate the integral in the risk function for the conditional censoring survival function. Defaults to 250 points equally spaced between 0 and the last censored follow-up time, minus a small constant in order to approximate left-continuous survivals.
#' @param saveFitLibrary Logical indicating whether to save the fit library on the full data. Defaults to \code{TRUE}. If \code{FALSE}, cannot obtain predicted values on new data later.
#' @return Returns a named list with control parameters.

survSuperLearner.control <- function (initWeightAlg = "survSL.rfsrc", initWeight = "censoring", max.SL.iter = 20, event.t.grid, cens.t.grid, saveFitLibrary = TRUE, ...) {
  if(!(initWeight %in% c("event", "censoring")))
    stop("initWeight must be one of event or censoring.")
  
  list(initWeightAlg = initWeightAlg, initWeight = initWeight, max.SL.iter = max.SL.iter, event.t.grid = event.t.grid, cens.t.grid = cens.t.grid, saveFitLibrary = saveFitLibrary)
}

#' Control parameters for the cross validation steps in survival Super Learner
#'
#' This function initiates control parameters for the cross-validation in \code{\link{survSuperLearner}} function.
#'
#' @param V Number of cross-validation folds. Defaults to 10.
#' @param stratifyCV Logical indicating whether to balance number of observed events across folds. Defaults to \code{TRUE}.
#' @param shuffle Logical indicating whether to shuffle the indices, or to simply assign sequentially. Defaults to \code{TRUE}. Should almost always be set to \code{TRUE} unless it is explicitly desired to assign sequentially.
#' @param validRows Optional custom list of indices for validation folds.
#' @return Returns a list of length \code{V} with validation indices for each of the folds.

survSuperLearner.CV.control <- function (V = 10L, stratifyCV = TRUE, shuffle = TRUE, validRows = NULL) {
  V <- as.integer(V)
  if (!is.null(validRows)) {
    if (!is.list(validRows)) {
      stop("validRows must be a list of length V containing the row numbers for the corresponding validation set")
    }
    if (!identical(V, length(validRows))) {
      stop("V and length(validRows) must be identical")
    }
  }
  list(V = V, stratifyCV = stratifyCV, shuffle = shuffle, validRows = validRows)
}

.createLibrary <- function (survSL.library)  {
  if (is.character(survSL.library)) {
    k <- length(survSL.library)
    whichScreen <- matrix(1, nrow = 1, ncol = k)
    screenAlgorithm <- "All"
    library <- data.frame(predAlgorithm = survSL.library, rowScreen = 1,
                          stringsAsFactors = FALSE)
  }
  else if (is.list(survSL.library)) {
    predNames <- sapply(survSL.library, FUN = "[", 1)
    NumberScreen <- (sapply(survSL.library, FUN = length) - 1)
    if (sum(NumberScreen == 0) > 0) {
      for (ii in which(NumberScreen == 0)) {
        SL.library[[ii]] <- c(SL.library[[ii]], "All")
        NumberScreen[ii] <- 1
      }
    }
    screenAlgorithmFull <- unlist(lapply(survSL.library, FUN = "[", -1))
    screenAlgorithm <- unique(screenAlgorithmFull)
    library <- data.frame(predAlgorithm = rep(predNames,
                                              times = NumberScreen), rowScreen = match(screenAlgorithmFull,
                                                                                       screenAlgorithm), stringsAsFactors = FALSE)
  }
  else {
    stop("format for survSL.library is not recognized")
  }
  out <- list(library = library, screenAlgorithm = screenAlgorithm)
  return(out)
}

.survCVFolds2 <- function (N, id, event, cvControl) {
  if (!is.null(cvControl$validRows)) return(cvControl$validRows)
  stratifyCV <- cvControl$stratifyCV
  shuffle <- cvControl$shuffle
  V <- cvControl$V
  if (!stratifyCV) {
    if (shuffle) {
      if (is.null(id)) {
        validRows <- split(sample(1:N), rep(1:V, length = N))
      }
      else {
        n.id <- length(unique(id))
        id.split <- split(sample(1:n.id), rep(1:V, length = n.id))
        validRows <- vector("list", V)
        for (v in seq(V)) {
          validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
        }
      }
    }
    else {
      if (is.null(id)) {
        validRows <- split(1:N, rep(1:V, length = N))
      }
      else {
        n.id <- length(unique(id))
        id.split <- split(1:n.id, rep(1:V, length = n.id))
        validRows <- vector("list", V)
        for (v in seq(V)) {
          validRows[[v]] <- which(id %in% unique(id)[id.split[[v]]])
        }
      }
    }
  } else {
    # if (sum(event) < V | sum(1-event) < V) {
    #   stop("number of (event = 1) or (event = 0) is less than the number of folds")
    # }
    if (shuffle) {
      if (is.null(id)) {
        event.0 <- which(event == 0)
        event.1 <- which(event == 1)
        rows.0 <- split(sample(event.0), rep(1:V, length = length(event.0)))
        rows.1 <- split(sample(event.1), rep(1:V, length = length(event.1)))
        validRows <- vector("list", length = V)
        names(validRows) <- paste(seq(V))
        for (vv in seq(V)) {
          if (length(rows.0) >= vv) {
            if (length(rows.1) >= vv) validRows[[vv]] <- c(rows.0[[vv]], rows.1[[vv]])
            else validRows[[vv]] <- rows.0[[vv]]
          } else {
            validRows[[vv]] <- rows.1[[vv]]
          }
        }
      } else {
        ## Stratified CV WITH id (NEW IMPLEMENTATION)
        
        # unique ids
        uid <- unique(id)
        n.id <- length(uid)
        
        # subject-level event indicator
        event_id <- tapply(event, id, max)
        
        # ids by event
        id0 <- names(event_id)[event_id == 0]
        id1 <- names(event_id)[event_id == 1]
        
        if (length(id0) < V || length(id1) < V) {
          stop("Not enough ids in one of the strata to create V folds.")
        }
        
        # shuffle ids within strata
        id0.split <- split(sample(id0), rep(1:V, length = length(id0)))
        id1.split <- split(sample(id1), rep(1:V, length = length(id1)))
        
        validRows <- vector("list", V)
        names(validRows) <- paste(seq(V))
        
        for (v in seq(V)) {
          ids_v <- c(id0.split[[v]], id1.split[[v]])
          validRows[[v]] <- which(id %in% ids_v)
        }
      }
    }else {
      if (is.null(id)) {
        within.split <- suppressWarnings(tapply(1:N,
                                                INDEX = event, FUN = split, rep(1:V)))
        validRows <- vector("list", length = V)
        names(validRows) <- paste(seq(V))
        for (vv in seq(V)) {
          validRows[[vv]] <- c(within.split[[1]][[vv]],
                               within.split[[2]][[vv]])
        }
      } else {
        ## Stratified CV WITH id (no shuffle)
        
        uid <- unique(id)
        event_id <- tapply(event, id, max)
        
        id0 <- names(event_id)[event_id == 0]
        id1 <- names(event_id)[event_id == 1]
        
        id0.split <- split(id0, rep(1:V, length = length(id0)))
        id1.split <- split(id1, rep(1:V, length = length(id1)))
        
        validRows <- vector("list", V)
        names(validRows) <- paste(seq(V))
        
        for (v in seq(V)) {
          ids_v <- c(id0.split[[v]], id1.split[[v]])
          validRows[[v]] <- which(id %in% ids_v)
        }
      }
    }
  }
  return(validRows)
}

.surviterativeSL2 <- function(event.Z, cens.Z, time, event, X_event, X_cens, obsWeights, id, control, verbose, event.errorsInLibrary, cens.errorsInLibrary) {
  if(verbose) message("Performing iterative SuperLearner...")
  event.k <- dim(event.Z)[3]
  cens.k <- dim(cens.Z)[3]
  N <- length(time)
  event.n.time <- length(control$event.t.grid)
  cens.n.time <- length(control$cens.t.grid)
  
  epsilon <- max(min(diff(sort(unique(time)))), 1e-5)
  
  event.idx <- findInterval(time, control$event.t.grid, rightmost.closed = TRUE)
  event.idx[event.idx == 0] <- 1
  
  event.Z.obs <- matrix(NA_real_, N, event.k)
  for (j in seq_len(event.k)) {
    event.Z.obs[, j] <- event.Z[cbind(seq_len(N), event.idx, j)]
  }
  
  cens.idx <- findInterval(time, control$cens.t.grid, rightmost.closed = TRUE)
  cens.idx[cens.idx == 0] <- 1
  
  cens.Z.obs <- matrix(NA_real_, N, cens.k)
  for (j in seq_len(cens.k)) {
    cens.Z.obs[, j] <- cens.Z[cbind(seq_len(N), cens.idx, j)]
  }
  
  initWeightAlg <- get(control$initWeightAlg)
  if (control$initWeight == "censoring") {
    # initFit <- try({initWeightAlg(time = time, event = 1 - event, X = X_cens, newX = X_cens,
    #                          new.times = time - epsilon,
    #                          obsWeights = obsWeights, id = id)})
    # 
    # if(inherits(initFit, "try-error")){
    initFit <- survSL.km(
      time = time,
      event = 1 - event,
      X = X_cens,
      newX = X_cens,
      new.times = control$cens.t.grid,
      obsWeights = obsWeights,
      id = id
    )
    #}else{
    #  
    #}
    
    # Vector to store subject-specific values
    obs.cens.vals <- numeric(N)
    
    # Loop over subjects (memory-efficient)
    for (i in seq_len(N)) {
      obs.cens.vals[i] <- approx(
        x = control$cens.t.grid,    # or the grid used in initFit
        y = initFit$pred[i, ],
        xout = time[i] - epsilon,   # for censoring, subtract epsilon
        method = "constant",
        rule = 2
      )$y
    }
    
    #### obs.event.vals -- where does it come from?
    S.coef <- rep(0, event.k)
    S.coef[!event.errorsInLibrary] <- .survcomputeCoef_stream(
      time = time, event = 1 - event, t.grid = control$cens.t.grid, cens.vals = obs.cens.vals,
      Z= cens.Z[, , !cens.errorsInLibrary, drop = FALSE], obsWeights = obsWeights)
    
    
    # Vector to store subject-specific values
    obs.event.vals <- numeric(N)
    
    obs.event.vals <- as.numeric(event.Z.obs %*% S.coef)
    
    #obs.event.vals <- rep(c(event.Z.obs %*% S.coef), length(control$cens.t.grid))
  } else {
    initFit <- initWeightAlg(time = time, event = event, X = X_event, newX = X_event,
                             new.times = time,
                             obsWeights = obsWeights, id = id)
    obs.event.vals <- rep(diag(initFit$pred), length(control$cens.t.grid))
  }
  obs.event.vals[obs.event.vals == 0] <- min(obs.event.vals[obs.event.vals > 0])
  
  iter <- 1
  while(TRUE) {
    #print(iter)
    #gc()
    if(iter > control$max.SL.iter) {
      warning("Did not converge in ", control$max.SL.iter, " iterations")
      break
    }
    if(!is.null(obs.cens.vals)) obs.cens.vals.old <- obs.cens.vals
    if(!is.null(obs.event.vals)) obs.event.vals.old <- obs.event.vals
    
    G.coef <- rep(0, cens.k)
    G.coef[!cens.errorsInLibrary] <- .survcomputeCoef_stream(
      time = time, 
      event = 1 - event, 
      t.grid = control$cens.t.grid, 
      cens.vals = obs.event.vals,
      Z= cens.Z[, , !cens.errorsInLibrary, drop = FALSE], 
      obsWeights = obsWeights)
    
    obs.cens.vals <- as.numeric(cens.Z.obs %*% G.coef)
    #obs.cens.vals <- rep(c(cens.Z.obs %*% G.coef), length(control$event.t.grid))
    obs.cens.vals[obs.cens.vals == 0] <- min(obs.cens.vals[obs.cens.vals > 0])
    
    S.coef[!event.errorsInLibrary] <- .survcomputeCoef_stream(
      time = time, event = event, t.grid = control$event.t.grid, cens.vals = obs.cens.vals,
      Z= event.Z[, , !event.errorsInLibrary, drop = FALSE], obsWeights = obsWeights)
    
    obs.event.vals <- as.numeric(event.Z.obs %*% S.coef)
    obs.event.vals[obs.event.vals == 0] <- min(obs.event.vals[obs.event.vals > 0])
    
    if(!is.null(obs.cens.vals.old) & !is.null(obs.event.vals.old)) {
      cens.delta <- max(abs(obs.cens.vals - obs.cens.vals.old))
      event.delta <- max(abs(obs.event.vals - obs.event.vals.old))
      if(cens.delta + event.delta < 1e-5) {
        if(verbose) message("Converged in ", iter, " iterations.")
        break
      }
    }
    iter <- iter + 1
  }
  # event.cvRisks <- apply(event.Z.long, 2, function(col) {
  #   mean((obsWeights.event.long * event.event.long / obs.cens.vals) * (as.numeric(time.event.long > event.t.grid.long) - col)^2)
  # })
  
  event.cvRisks <- numeric(event.k)
  for (j in seq_len(event.k)) {
    risk <- 0
    for (tt in seq_along(control$event.t.grid)) {
      y <- 1 - as.numeric(time <= control$event.t.grid[tt]) *
        event / obs.cens.vals
      res <- y - event.Z[, tt, j]
      risk <- risk + sum(obsWeights * res * res)
    }
    event.cvRisks[j] <- risk / (length(control$event.t.grid) * sum(obsWeights))
  }
  
  # event.cvRisks <- apply(event.Z.long, 2, function(col) {
  #   -mean(obsWeights.event.long * ifelse(time.event.long <= event.t.grid.long & event.event.long == 1,
  #                                            (1 - 1/obs.cens.vals) * log(col) + (1 / obs.cens.vals) * log(1 - col),
  #                                            log(col)))
  # })
  # cens.cvRisks <- apply(cens.Z.long, 2, function(col) {
  #   mean((obsWeights.cens.long * (1-event.cens.long) / obs.event.vals) * (as.numeric(time.cens.long > cens.t.grid.long) - col)^2)
  # })
  cens.cvRisks <- numeric(cens.k)
  for (j in seq_len(cens.k)) {
    risk <- 0
    for (tt in seq_along(control$event.t.grid)) {
      y <- 1 - as.numeric(time <= control$cens.t.grid[tt]) *
        event / obs.cens.vals
      res <- y - cens.Z[, tt, j]
      risk <- risk + sum(obsWeights * res * res)
    }
    cens.cvRisks[j] <- risk / (length(control$cens.t.grid) * sum(obsWeights))
  }
  
  # cens.cvRisks <- apply(cens.Z.long, 2, function(col) {
  #   -mean(obsWeights.cens.long * ifelse(time.cens.long <= cens.t.grid.long & event.cens.long == 1,
  #                                        (1 - 1/obs.event.vals) * log(col) + (1 / obs.event.vals) * log(1 - col),
  #                                        log(col)))
  # })
  return(list(event.coef = S.coef, cens.coef = G.coef, event.cvRisks = event.cvRisks, cens.cvRisks = cens.cvRisks))
}

.survcomputeCoef_stream <- function(time, event, t.grid, cens.vals, Z, obsWeights) {
  
  k <- dim(Z)[3]
  if (k == 1) return(1)
  
  cens.vals[cens.vals < 1e-4] <- 1e-4
  sw <- sqrt(obsWeights)
  
  XtWX <- matrix(0, k, k)
  XtWy <- numeric(k)
  
  for (tt in seq_along(t.grid)) {
    t0 <- t.grid[tt]
    y  <- 1 - as.numeric(time <= t0) * event / cens.vals
    
    Zt <- Z[, tt, ]        # N × k
    WZ <- sw * Zt
    Wy <- sw * y
    
    XtWX <- XtWX + crossprod(WZ)
    XtWy <- XtWy + crossprod(WZ, Wy)
  }
  
  fit <- nnls::nnls(XtWX, XtWy)
  coef <- coef(fit)
  if (sum(coef) == 0) coef <- rep(1, k)
  coef / sum(coef)
}

.alogb <- function(a,b) {
  ifelse(a == 0, )
}