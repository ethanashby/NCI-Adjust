# Counterfactual survival utilities used by the simulation analysis.

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

.estimate.uniform.quantile <- function(IF.vals, conf.level=.95, scale=TRUE) {
  IF.vals <- cbind(IF.vals)
  n <- nrow(IF.vals)
  if(scale) IF.vals <- scale(IF.vals)
  maxes <- replicate(1e4, max(abs(rbind(rt(n, df = n - 1)/sqrt(n)) %*% IF.vals)))
  return(list(quantile=quantile(maxes, conf.level), maxes=maxes))
}

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
