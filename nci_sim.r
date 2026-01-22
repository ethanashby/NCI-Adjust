library(CFsurvival,
        lib.loc=c("/home/eashby/R_lib"))
library(SimEngine,
        lib.loc=c("/home/eashby/R_lib"))
library(survSuperLearner,
        lib.loc=c("/home/eashby/R_lib"))
library(copula,
        lib.loc=c("/home/eashby/R_lib"))
library(tidyverse)
library(survival)
library(msm)
library(parallel)
library(progress)

if (Sys.getenv("SLURM_JOB_ID") == "") { n_cores <- parallel::detectCores() } else { n_cores <- 600 }

run_on_cluster(
  
  first={
    
    sim<-new_sim()
    
    sim %<>% set_levels(
      n=c(1200, 2400, 4800),
      trt_prob=0.5,
      RR = c(0.75), #25% efficacy
      RR_C = c(1),
      per_year_incidence_Y = c(0.05),
      per_year_incidence_N = c(0.12),
      per_year_incidence_C = c(0.08),
      J=1 ,
      rho=c(0, 0.5, 0.75, 0.9),
      X = TRUE,
      rho_X = c(0.1),
      admin_date = 2
    )
    
    sim %<>% set_script(function(){
      
      source("create_multivar_survival.r")
      source("CFSurvival_funcs.r")
      source("survSuperLearner_funcs.r")
      
      data <- create_multivar_survival(n=L$n, 
                                       trt_prob=L$trt_prob, 
                                       RR = L$RR, 
                                       RR_C = L$RR_C, 
                                       per_year_incidence_Y = L$per_year_incidence_Y, 
                                       per_year_incidence_N = L$per_year_incidence_N, 
                                       per_year_incidence_C = L$per_year_incidence_C,
                                       J = L$J, 
                                       rho = L$rho, 
                                       admin_date = L$admin_date,
                                       X=L$X,
                                       rho_X = L$rho_X, 
                                       convert_days = TRUE
      )
      
      km<-summary(survfit(Surv(Y, DeltaY)~A, data=data), times=max(data$Y[data$DeltaY==1]))
      res_KM <- list(
        "logRR"=log((1-km$surv[2])/(1-km$surv[1])),
        "logRR_var" = as.numeric(msm::deltamethod(~log((1-x2)/(1-x1)), mean=km$surv, diag(km$std.err^2), ses = FALSE))
      )
      
      fit.treat <- glm(A ~ X, data=data, family="binomial")
      prop.fit <- predict(fit.treat, type = "response")
      
      est0 <- CFsurvival::CFsurvival(
        time = data$Y, 
        event = data$DeltaY, 
        treat = data$A, 
        confounders = data.frame(X=data$X), 
        fit.times = sort(unique(data$Y)), 
        fit.treat = c(0,1), 
        nuisance.options = CFsurvival.nuisance.options(
          cross.fit = TRUE,
          V = 5,
          folds = NULL,
          eval.times = NULL,
          event.SL.library = c("survSL.km", 
                               "survSL.coxph", 
                               "survSL.expreg"),
          event.pred.0 = NULL,
          event.pred.1 = NULL,
          cens.SL.library = c("survSL.km", 
                              "survSL.coxph", 
                              "survSL.expreg"),
          cens.trunc = 0,
          cens.pred.0 = NULL,
          cens.pred.1 = NULL,
          survSL.control = list(initWeightAlg = "survSL.km", verbose = FALSE),
          survSL.cvControl = list(V = 5),
          save.nuis.fits = FALSE,
          prop.SL.library = NULL,
          prop.trunc = 0,
          prop.pred = prop.fit,
          verbose = FALSE
        ), conf.band = FALSE, conf.level = 0.95, contrasts = c("risk.ratio"), verbose = TRUE
      )
      
      est <- CFsurvival2(
        time = data$Y, 
        event = data$DeltaY, 
        time_N = data$N1, 
        event_N = data$DeltaN1, 
        admin_t = round(365*L$admin_date), 
        treat = data$A, 
        W = data.frame(X=data$X), 
        id = 1:nrow(data), 
        fit.times = sort(unique(data$Y)), 
        fit.treat = c(0,1), 
        nuisance.options = CFsurvival.nuisance.options(
          cross.fit = TRUE,
          V = 5,
          folds = NULL,
          eval.times = NULL,
          event.SL.library = c("survSL.km", 
                               "survSL.coxph", 
                               "survSL.expreg"),
          event.pred.0 = NULL,
          event.pred.1 = NULL,
          cens.SL.library = c("survSL.km", 
                              "survSL.coxph", 
                              "survSL.expreg"),
          cens.trunc = 0,
          cens.pred.0 = NULL,
          cens.pred.1 = NULL,
          survSL.control = list(initWeightAlg = "survSL.km", verbose = FALSE),
          survSL.cvControl = list(V = 5),
          save.nuis.fits = FALSE,
          prop.SL.library = NULL,
          prop.trunc = 0,
          prop.pred = prop.fit,
          verbose = TRUE
        ), 
        conf.band = FALSE, 
        conf.level = 0.95, 
        contrasts = c("risk.ratio"), 
        verbose = TRUE
      )
      
      results <- bind_rows(
        data.frame("logRR"=res_KM$logRR, "logRR_var"=res_KM$logRR_var, "type"="KM"),
        data.frame("logRR"=tail(est0$risk.ratio.df$log.risk.ratio, 1), 
                   "logRR_var"=tail(est0$risk.ratio.df$se.log.ratio^2, 1), 
                   "type"="X-Adj"),
        data.frame("logRR"=tail(est$risk.ratio.df$log.risk.ratio, 1), 
                   "logRR_var"=tail(est$risk.ratio.df$se.log.ratio^2, 1), 
                   "type"="NCI-Adj")
      )                            
      
      r <- results %>% 
        dplyr::mutate(
          bias = logRR - log(L$RR),
          releff = logRR_var/logRR_var[type=="KM"],
          cover = as.numeric(between(rep(log(L$RR), nrow(results)), 
                                     logRR - qnorm(0.975) * sqrt(logRR_var), 
                                     logRR + qnorm(0.975) * sqrt(logRR_var)))
        )
      
      return(
        list(
          ".complex"= list(
            "output" = r
          )))
      
    })
    
    
    sim %<>% set_config(num_sim=500, 
                        seed=47,
                        stop_at_error=TRUE,
                        n_cores = n_cores)
    
  },
  
  main={
    sim %<>% run()
  },
  
  last = {
    
  },
  
  cluster_config = list(js="slurm")
)
