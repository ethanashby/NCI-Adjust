#####
# This is code to submit to the FH rhino cluster
# install.packages('name', lib="/home/users/eashby/R_lib")

library(SimEngine,
        lib.loc=c("/home/eashby/R_lib"))
library(survival,
        lib.loc=c("/home/eashby/R_lib"))
library(Matrix,
        lib.loc=c("/home/eashby/R_lib"))
library(glmnet,
        lib.loc=c("/home/eashby/R_lib"))
library(mgcv,
        lib.loc=c("/home/eashby/R_lib"))
library(randomForestSRC,
        lib.loc=c("/home/eashby/R_lib"))
library(gam)
library(nnls)
library(Rsolnp,
        lib.loc=c("/home/eashby/R_lib"))
library(SuperLearner, lib.loc=c("/home/eashby/R_lib"))
library(survSuperLearner,
        lib.loc=c("/home/eashby/R_lib"))
library(CFsurvival,
        lib.loc=c("/home/eashby/R_lib"))
library(copula,
        lib.loc=c("/home/eashby/R_lib"))
library(msm,
        lib.loc=c("/home/eashby/R_lib"))
library(tidyverse)
library(parallel)
library(progress)

if (Sys.getenv("SLURM_JOB_ID") == "") { n_cores <- parallel::detectCores() } else { n_cores <- 800 }

run_on_cluster(
  
  first={
    
    sim<-new_sim()
    
    sim %<>% set_levels(
      #n=c(1200, 2400, 4800),
      n=c(1500),
      trt_prob=0.5,
      RR = c(0.75), #25% efficacy
      RR_C = c(1),
      per_year_incidence_Y = c(0.06),
      per_year_incidence_N = c(0.12),
      #per_year_incidence_C = c(0.05, 0.10, 0.20, 0.30),
      per_year_incidence_C = c(0.10),
      J=1,
      rho=c(0, 0.5, 0.8),
      X = TRUE,
      #rho_X = c(0.1, 0.5),
      rho_X = c(0),
      admin_date = 2
    )
    
    # sim %<>% set_levels(.keep=
    #                       unique(
    #                         c(
    #                           which(sim$levels_grid$per_year_incidence_C==0.10), # Help me understand influence of prognostic N versus X
    #                           which(sim$levels_grid$rho_X==0.1 & sim$levels_grid$rho==0.75 & sim$levels_grid$n==2000) # Help me understand impact of coarsening versus prognostic N
    #                           )
    #                         ))
    
    sim %<>% set_script(function(){
      
      source("./create_multivar_survival.r")
      source("./CFSurvival_funcs.r")
      rm(CFsurvival2)
      rm(CFsurvival3)
      source("./survSuperLearner_funcs.r")
      
      message("[1] Generating survival data")
      
      L<-sim$levels_grid[3,]
      
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
      
      message("[2] Running unadjusted (KM) analysis")
      
      km_fit<-survfit(Surv(Y, DeltaY)~A, data=data) 
      
      #km<-summary(km_fit, times=((L$admin_date*365)))
      
      km<-summary(km_fit, times=max(data$Y), extend=TRUE)
      
      res_KM <- list(
        "logRR"=log((1-km$surv[2])/(1-km$surv[1])),
        "logRR_var" = as.numeric(msm::deltamethod(~log((1-x2)/(1-x1)), mean=km$surv, diag(km$std.err^2), ses = FALSE))
      )
      
      res_KM[["logRR_lci"]] = res_KM[["logRR"]] - qt(1-(1-0.95)/2, df=L$n) * sqrt(res_KM[["logRR_var"]])
      res_KM[["logRR_uci"]] = res_KM[["logRR"]] + qt(1-(1-0.95)/2, df=L$n) * sqrt(res_KM[["logRR_var"]])
      
      results <- bind_rows(
        data.frame("logRR"=res_KM$logRR, 
                   "logRR_lci" = res_KM$logRR_lci,
                   "logRR_uci" = res_KM$logRR_uci,
                   "logRR_var"=res_KM$logRR_var, 
                   "type"="KM")
      )
      rm(km)
      
      message("[3] Running X-Adj Analysis")
      
      fit.treat <- glm(A ~ X, data=data, family="binomial")
      prop.fit <- predict(fit.treat, type = "response")
      rm(fit.treat)
      
      est0 <- CFsurvival::CFsurvival(
        time = data$Y, 
        event = data$DeltaY, 
        treat = data$A, 
        confounders = data.frame(X=data$X), 
        fit.times = sort(unique(c(data$Y[data$DeltaY==1], ((L$admin_date*365))))), 
        fit.treat = c(0,1), 
        nuisance.options = CFsurvival.nuisance.options(
          cross.fit = TRUE,
          V = 5,
          folds = NULL,
          eval.times = NULL,
          event.SL.library = c("survSL.coxph", 
                               "survSL.km", 
                               "survSL.expreg"),
          event.pred.0 = NULL,
          event.pred.1 = NULL,
          cens.SL.library = c("survSL.coxph", 
                              "survSL.km", 
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
        ), conf.band = FALSE, 
        conf.level = 0.95, 
        contrasts = c("risk.ratio"), 
        verbose = TRUE
      )
      
      results <- bind_rows(results,
                           data.frame("logRR"=tail(est0$risk.ratio.df$log.risk.ratio, 1), 
                                      "logRR_lci" = tail(est0$risk.ratio.df$ptwise.lower.log, 1),
                                      "logRR_uci" = tail(est0$risk.ratio.df$ptwise.upper.log, 1),
                                      "logRR_var"=tail(est0$risk.ratio.df$se.log.ratio^2, 1), 
                                      "type"="X-Adj")
      )
      
      rm(est0)
      gc()
      
      source("./oracle_attempt_2.R")

      oracle <- AIPWCC_newEIF_oracle(
       time=data$Y,
       event=data$DeltaY,
       time_N=data$N1,
       event_N=data$DeltaN1,
       admin_t=730,
       treat=data$A,
       W=data.frame(X=data$X),
       id=1:nrow(data),
       fit.times        = sort(unique(data$Y)),
       eval.times       = seq(0, 730, by=1),
       fit.treat        = c(0, 1),
       conf.band        = FALSE,
       conf.level       = 0.95,
       contrasts        = c("risk.ratio"),
       verbose          = TRUE,
       dgm_params       = list(
         n = 5000000,
         trt_prop = L$trt_prob,
         RR = L$RR,
         RR_C = L$RR_C,
         per_year_incidence_Y = L$per_year_incidence_Y,
         per_year_incidence_N = L$per_year_incidence_N,
         per_year_incidence_C = L$per_year_incidence_C,
         rho = L$rho,
         admin_date = L$admin_date,
         rho_X = L$rho_X)
      )

      results <- bind_rows(results,
                          data.frame("logRR"=tail(oracle$risk.ratio$log.risk.ratio, 1),
                                     "logRR_lci" = tail(oracle$risk.ratio$ptwise.lower.log, 1),
                                     "logRR_uci" = tail(oracle$risk.ratio$ptwise.upper.log, 1),
                                     "logRR_var"=tail(oracle$risk.ratio$se.log.ratio^2, 1),
                                     "type"="NCI-adj (Oracle)")
      )
      
      rm(oracle)
      
      message("[4] Running NCI-Adj Analysis (IF-based variance)")
      source("./AIPWCC_est_fullaug2.R")
      
      eval.times <- sort(unique(c(0, data$N1)))
      
      est1 <- AIPWCC_newEIF_em(
        time = data$Y,
        event = data$DeltaY,
        time_N = data$N1,
        event_N = data$DeltaN1,
        admin_t = round(365*L$admin_date),
        treat = data$A,
        W = data.frame(X=data$X),
        id=1:nrow(data),
        fit.times=sort(unique(c(data$Y, (L$admin_date*365)))),
        fit.treat=c(0,1),
        conf.band=FALSE,
        conf.level=.95,
        contrasts = c("risk.ratio"),
        verbose=TRUE,
        nuisance.options = .CFsurvival.nuis.opt(
          cross.fit = TRUE,
          V = 5,
          folds = NULL,
          eval.times = eval.times,
          event.SL.library = c("survSL.coxph", "survSL.expreg", "survSL.km"),
          event.SL.library_SY = c("survSL.coxph", "survSL.coxph.ps", "survSL.rfsrc"),
          event.pred.0 = NULL,
          event.pred.1 = NULL,
          cens.SL.library = c("survSL.coxph", "survSL.expreg", "survSL.km", "survSL.rfsrc"),
          cens.trunc = 0,
          cens.pred.0 = NULL,
          cens.pred.1 = NULL,
          survSL.control = list(initWeightAlg = "survSL.coxph", verbose = FALSE),
          survSL.cvControl = list(V = 5),
          save.nuis.fits = FALSE,
          prop.SL.library = NULL,
          prop.trunc = 0,
          prop.pred = prop.fit,
          verbose = TRUE
        ), 
        single_cens_mechanism = TRUE,
        em = TRUE, 
        em_fast = TRUE,
        max_iter = 5,
        tol  = 5E-4
      )
      
      tmp<-tail(est1$risk.ratio,1)
      
      results <- bind_rows(results,
                           data.frame(
                             "logRR" = tmp$log.risk.ratio,
                             "logRR_lci" = tmp$log.risk.ratio - qt(1-(1-0.95)/2, L$n) * tmp$se.log.ratio,
                             "logRR_uci" = tmp$log.risk.ratio + qt(1-(1-0.95)/2, L$n) * tmp$se.log.ratio,
                             "logRR_var" = (tmp$se.log.ratio)^2,
                             "type" = "NCI-adj (CF)"
                           )
      ) 
      
      rm(est1)
      gc()
      
      
      r <- results %>% 
        dplyr::mutate(
          bias = logRR - log(L$RR),
          releff = logRR_var/logRR_var[type=="KM"],
          cover = as.numeric(between(rep(log(L$RR), nrow(results)), 
                                     logRR_lci, 
                                     logRR_uci))
        )
      
      return(
        list(
          ".complex"= list(
            "output" = r
          )))
      
    })
    
    
    sim %<>% set_config(num_sim=1000, 
                        seed=474,
                        stop_at_error=TRUE,
                        n_cores = n_cores)
    
  },
  
  main={
    sim %<>% run()
  },
  
  last = {
    sim
  },
  
  cluster_config = list(js="slurm")
)