library(tidyverse)
library(survival)
library(parallel)
library(SimEngine,
        lib.loc=c("/home/eashby/R_lib"))
library(copula,
        lib.loc=c("/home/eashby/R_lib"))
library(RobinCar, 
        lib.loc=c("/home/eashby/R_lib"))

if (Sys.getenv("SLURM_JOB_ID") == "") { n_cores <- parallel::detectCores() } else { n_cores <- 600 }

run_on_cluster(
  
  first={
    
    sim<-new_sim()
    
    sim %<>% set_levels(
      n=c(1200, 2400, 4800),
      trt_prob=0.5,
      RR = c(0.8, 0.2), #20% and 80% efficacy
      RR_C = c(1),
      per_year_incidence_Y = c(0.05),
      per_year_incidence_N = c(0.04, 0.08, 0.16),
      per_year_incidence_C = c(0.10),
      J=3,
      rho=c(0, 0.3, 0.5, 0.7, 0.9),
      admin_date = 2,
      burn_in = 20,
      MC_samples = 20
    )
    
    source("/home/eashby/NCI_adj/create_multivar_survival.r")
    source("/home/eashby/NCI_adj/MI_infections.r")
    
    sim %<>% set_script(function(){
      
      #L<-sim$levels_grid[60,]
      
      L<-list(n=5000, 
              trt_prob=0.50, 
              RR=0.5, 
              RR_C=2.5, 
              per_year_incidence_Y = 0.05, 
              per_year_incidence_N = 0.15,
              per_year_incidence_C = 0.25,
              J = 3,
              rho = 0.85,
              admin_date=2, burn_in =15, MC_samples = 15)
      
      data <- create_multivar_survival(n=L$n, 
                                       trt_prob=L$trt_prob, 
                                       RR = L$RR, 
                                       RR_C = L$RR_C, 
                                       per_year_incidence_Y = L$per_year_incidence_Y, 
                                       per_year_incidence_N = L$per_year_incidence_N, 
                                       per_year_incidence_C = L$per_year_incidence_C, 
                                       J = L$J, 
                                       rho = L$rho, 
                                       admin_date=L$admin_date)
      
      res_KM <- KM_est(dat=data, T.name = "Y", Delta.name = "DeltaY", A.name = "A", eval_t = L$admin_date)
      
      X <- data %>% transmute(
        N = N1,
        DeltaN = N1
      )
        
        transmute(N=pmin(N1, N2, N3), 
                              DeltaN=case_when(N==N1 ~ DeltaN1, 
                                               N==N2 ~ DeltaN2, 
                                               N==N3 ~ DeltaN3)) %>% as.matrix()
      #%>% dplyr::select(starts_with("N"), starts_with("DeltaN"))
      
      fit<-CFsurvival(
        time = data$Y,
        event = data$DeltaY,
        treat = data$A,
        confounders = X,
        fit.times = sort(unique(data$Y)),
        contrasts = "risk.ratio",
        conf.band=FALSE,
        nuisance.options = list(
          cross.fit=TRUE,
          V=3,
          event.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg",
                             "survSL.weibreg", "survSL.loglogreg"), function(alg)
                             c(alg, "survscreen.glmnet", "survscreen.marg", "All")),
          #event.SL.library = "survSL.coxph",
          # cens.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg",
          #                            "survSL.weibreg", "survSL.loglogreg"), function(alg)
          #                              c(alg, "survscreen.glmnet", "survscreen.marg", "All")),
          prop.pred  = rep(0.5, length(data$Y))
        )
      )
      
      gs_single <- multivariate_gibbs_sampling(dat=data %>% mutate(N=pmin(N1, N2, N3), 
                                        DeltaN=case_when(N==N1 ~ DeltaN1, N==N2 ~ DeltaN2, N==N3 ~ DeltaN3)), 
                                        T.names=c("Y", "N"), 
                                        Delta.names = c("DeltaY", "DeltaN"), 
                                        primary_outcome_indx = 1, 
                                        A.name = "A", 
                                        burn_in = L$burn_in, 
                                        MC_samples = L$MC_samples)
      
      res_NCI_single <- NCI_adj(imputed_datasets = gs_single, 
                         indices_to_analyze = tail(seq_along(1:length(gs_single)), L$MC_samples), 
                         Y.name = "Y", 
                         A.name = "A", 
                         adj_vars = c("N", "DeltaN"),
                         eval_t = L$admin_date, calibrate=TRUE)
      
      gs <- multivariate_gibbs_sampling(dat=data, 
                                  T.names=c("Y", "N1", "N2", "N3"), 
                                  Delta.names = c("DeltaY", "DeltaN1", "DeltaN2", "DeltaN3"), 
                                  primary_outcome_indx = 1, 
                                  A.name = "A", 
                                  burn_in = L$burn_in, 
                                  MC_samples = L$MC_samples)
      
      res_NCI <- NCI_adj(imputed_datasets = gs, 
              indices_to_analyze = tail(seq_along(1:length(gs)), L$MC_samples), 
              Y.name = "Y", 
              A.name = "A", 
              adj_vars = c("N1", "N2", "N3", "DeltaN1", "DeltaN2", "DeltaN3"),
              eval_t = L$admin_date,
              calibrate = TRUE)
              
      results <- bind_rows(
        data.frame("logRR"=res_KM$logRR, "logRR_var"=res_KM$logRR_var, "type"="KM"),
        data.frame("logRR"=res_NCI_single$logRR, "logRR_var"=res_NCI_single$var_logRR, "type"="NCI-Adj (J=1)"),
        data.frame("logRR"=res_NCI$logRR, "logRR_var"=res_NCI$var_logRR, "type"="NCI-Adj (J=3)")
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
    
    
    sim %<>% set_config(num_sim=1000, 
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
