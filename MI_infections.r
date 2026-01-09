
KM_est <- function(dat, T.name, Delta.name, A.name, eval_t){
  
  #dat=test
  #T.name = "Y"
  #Delta.name = "DeltaY"
  #A.name="A"
  
  fits<-survfit(as.formula(paste0("Surv(", T.name, ", ", Delta.name, ") ~", A.name)), data=dat)
  cumincs<-1-summary(fits, t=eval_t)$surv
  ses <- summary(fits, t=eval_t)$std.err
  
  logRR <- log(cumincs[2]/cumincs[1])
  logRR_var <- msm::deltamethod(g= ~ log(x2/x1), mean=cumincs, cov=diag(ses^2), ses = FALSE)
  #RR <- cumincs[2]/cumincs[1]
  #RR_var <- msm::deltamethod(g= ~ x2/x1, mean=cumincs, cov=diag(ses^2), ses = FALSE)
  
  return(list(
    logRR = logRR,
    logRR_var = logRR_var
    #RR = RR,
    #RR_var = RR
  ))
}

multivariate_gibbs_sampling <- function(dat, 
                         T.names, 
                         Delta.names, 
                         primary_outcome_indx=1, 
                         A.name, 
                         burn_in = 25, 
                         MC_samples = 15,
                         verbose = TRUE, admin_t=2){

  # dat=dat_704_assess %>% dplyr::select(A, fu_year_wk104, statuswk104, ends_with("_TTE"),starts_with("delta")); 
  # T.names = T.names;
  # Delta.names = D.names;
  # primary_outcome_indx = 1;
  # A.name = "A";
  # burn_in = 20;
  # MC_samples = 20; 
  # verbose = TRUE; 
  # admin_t=2
  
  #dat = test
  #A.name = "A"
  #T.names = c("Y", "N1", "N2", "N3")
  #Delta.names = c("DeltaY", "DeltaN1", "DeltaN2", "DeltaN3")
  #primary_outcome_indx=1
  Y.name = T.names[primary_outcome_indx]
  DeltaY.name= Delta.names[primary_outcome_indx]
  
  num_outcomes=length(T.names)
  
  #Step 1: imputation Y and N (marginal)
  imputed_datasets = replicate(burn_in+MC_samples+1, dat %>% dplyr::select(c(A.name, T.names, Delta.names)), simplify=FALSE)
  
  if(verbose){message("Step 0: Initialization (marginal) KM imputation")}
  # Iterate through participants
  for(i in 1:nrow(dat)){
    Ai=dat[[A.name]][i]
    
    for(j in 1:num_outcomes){
      D.name = Delta.names[j]
      T.name = T.names[j]
      
      Ti=dat[[T.name]][i]
      Deltai = dat[[D.name]][i]
      if(is.null(admin_t) | is.na(admin_t)){admin_t = max(dat[[T.name]], na.rm=TRUE)}
    
      if(Ti==admin_t & Deltai == 0 | Deltai ==1){
        imputed_datasets[[1]][[T.name]][i]<-Ti
        imputed_datasets[[1]][[D.name]][i]<-Deltai
      }else{
      
        to_fit <- dat %>% filter(get(A.name)==Ai & get(T.name) >= Ti)
      
        if(nrow(to_fit)==0){
          imputed_datasets[[1]][[T.name]][i]<-admin_t
          imputed_datasets[[1]][[D.name]][i]<-0
        }else{
          km<-survfit(Surv(get(T.name), get(D.name)) ~ 1, data=to_fit)
          u = runif(1)
          Tistar = quantile(km, probs=u)$quantile
          imputed_datasets[[1]][[T.name]][i]<-ifelse(is.na(Tistar), admin_t, Tistar)
          imputed_datasets[[1]][[D.name]][i]<-ifelse(is.na(Tistar), 0, 1)
        }
      }
    }
  
  }
  if(verbose){
    message("Step 1: Gibbs sampler to impute from conditional survival functions")
  }
  if(verbose){
    pb1 <- progress::progress_bar$new(format = "Burn-in [:bar] :percent eta: :eta", 
                                      total = burn_in, 
                                      clear = FALSE, 
                                      width = 100,
                                      show_after = 0,
                                      force=TRUE)
    pb2 <- progress::progress_bar$new(format = "MC samples [:bar] :percent eta: :eta",
                                      total = MC_samples, 
                                      clear = FALSE, 
                                      width = 100,
                                      show_after = 0,
                                      force=TRUE)
  }
  
  for(iter in 1:(burn_in + MC_samples)){

  if(iter==1){
    pb1$tick(0)
  }
    
  if(iter==burn_in+1){
    pb2$tick(0)
  }
    
    # Initialize values from previous iteration
    imputed_datasets[[(iter+1)]] <- imputed_datasets[[(iter)]]
  
    n = nrow(imputed_datasets[[(iter+1)]])
    p = ncol(imputed_datasets[[(iter+1)]])
    
    for(j in 1:num_outcomes){
      
    D.name = Delta.names[j]
    T.name = T.names[j]
  
    #message(paste0("\tImputing censored ", T.name, " values from distn ", T.name, "|",  paste0(setdiff(T.names, T.name), collapse=", ")))
    
    admin_t = max(dat[[T.name]])
    
    for(i in 1:n){ 
  
      Ai=dat[[A.name]][i]
      Ti = dat[[T.name]][i]
      Deltai = dat[[D.name]][i]
      
      # If observed, do not impute
      if(Ti==admin_t & Deltai == 0 | Deltai==1){
        
        imputed_datasets[[(iter+1)]][[T.name]][i]<-Ti
        imputed_datasets[[(iter+1)]][[D.name]][i]<-Deltai
      
      # If censored, impute
      }else{
        
        rs_indx <- which(dat[[T.name]]>Ti & dat[[A.name]] == Ai)
        
        risk_set <- dat %>% filter(get(T.name) > Ti, get(A.name)==Ai)
        
        if(nrow(risk_set)==0){
        
          imputed_datasets[[(iter+1)]][[T.name]][i]<-admin_t
          imputed_datasets[[(iter+1)]][[D.name]][i]<-0
          
        }else{
        
          w = matrix(0, nrow=nrow(risk_set), ncol=num_outcomes)
        
          for(k in 1:num_outcomes){
          
            if(k==j){
              w[,k]<-1
            }else{
              # Return 
              bw = 2.34 * stats::sd(imputed_datasets[[(iter+1)]][[T.names[k]]][rs_indx]) * (nrow(risk_set))^(-1/(num_outcomes - 1 + 4))
              z <- (imputed_datasets[[(iter+1)]][[T.names[k]]][rs_indx] - imputed_datasets[[(iter+1)]][[T.names[k]]][i]) / bw
              w_tmp <- ifelse(abs(z) <= 1, 0.75 * (1 - z^2), 0)
              w_tmp <- w_tmp / sum(w_tmp)
              w[,k]<-w_tmp
          }
        }
        
        w=apply(w, MARGIN=1, prod)
        w = w/sum(w)
        
        if(all(is.nan(w)) | all(is.na(w))){
          w<-rep(1, length=nrow(risk_set))
        }
        
        km<- survfit(Surv(get(T.name), get(D.name))~1, 
                     data=risk_set, 
                     weights = w)
        
        #plot(km, fun="event", mark.time=FALSE)
        u = runif(1)
        Tistar = quantile(km, probs=u)$quantile
        imputed_datasets[[(iter+1)]][[T.name]][i]<-ifelse(is.na(Tistar), admin_t, Tistar)
        imputed_datasets[[(iter+1)]][[D.name]][i]<-ifelse(is.na(Tistar), 0, 1)
        }
      }
      
    }
    }
    
    if(verbose){
      if(iter <= burn_in){
        pb1$tick()
      }else{
        pb2$tick()
      }
    }
  }
  
  return(
    imputed_datasets
  )
}


library(doParallel)
library(foreach)

multivariate_gibbs_sampling_speedup <- function(dat, 
                                        T.names, 
                                        Delta.names, 
                                        primary_outcome_indx=1, 
                                        A.name, 
                                        burn_in = 25, 
                                        MC_samples = 15,
                                        verbose = TRUE, 
                                        admin_t=2,
                                        predictor_transform = c("identity", "log", "prob"),
                                        kern_shape = c("epanechnikov", "gaussian"),
                                        bw = 1){
  
  # dat=dat_704_assess %>% dplyr::select(A, fu_year_wk104, statuswk104, ends_with("_TTE"),starts_with("delta"));
  # T.names = T.names;
  # Delta.names = D.names;
  # primary_outcome_indx = 1;
  # A.name = "A";
  # burn_in = 10;
  # MC_samples = 10;
  # verbose = TRUE;
  # admin_t=2 + 1/365.25
  # kern_shape = "gaussian"
  # predictor_transform = "identity"
  
  #dat = test
  #A.name = "A"
  #T.names = c("Y", "N1", "N2", "N3")
  #Delta.names = c("DeltaY", "DeltaN1", "DeltaN2", "DeltaN3")
  #primary_outcome_indx=1
  Y.name = T.names[primary_outcome_indx]
  DeltaY.name= Delta.names[primary_outcome_indx]
  
  num_outcomes=length(T.names)
  
  #Step 1: imputation Y and N (marginal)
  
  imputed_datasets = replicate(burn_in+MC_samples+1, dat %>% dplyr::select(c(A.name, T.names, Delta.names)), simplify=FALSE)
  
  if(verbose){message("Step 0: Initialization (marginal) KM imputation")}
  
  # Iterate through participants
  # ~ 30 fold faster initialization step by computing KM estimator *ONCE* at start of loop rather than at each participant
  #tic()
  for(i in 1:nrow(dat)){
    #print(i)
    Ai=dat[[A.name]][i]
    
    for(j in 1:num_outcomes){
      
      #Collect time and event indicator names
      D.name = Delta.names[j]
      T.name = T.names[j]
      
      # fit (unweighted) KM estimators to outcome j
      if(i==1){
        km1<-survfit(Surv(get(T.name), get(D.name)) ~ 1, data=dat %>% filter(get(A.name)==1))
        Fhat1 <- stepfun(km1$time, c(0, 1 - km1$surv))
        km0<-survfit(Surv(get(T.name), get(D.name)) ~ 1, data=dat %>% filter(get(A.name)==0))
        Fhat0 <- stepfun(km0$time, c(0, 1 - km0$surv))
      }
      
      # Collect time and event indicator for participant i
      Ti=dat[[T.name]][i]
      Deltai = dat[[D.name]][i]
      
      if(is.null(admin_t) | is.na(admin_t)){admin_t = max(dat[[T.name]], na.rm=TRUE)}
      
      if(Ti==admin_t & Deltai == 0 | Deltai ==1){
        imputed_datasets[[1]][[T.name]][i]<-Ti
        imputed_datasets[[1]][[D.name]][i]<-Deltai
      }else{
        
        if(Ai==1){
          q <- Fhat1(Ti)
          u <- runif(1, min = q, max = 1)
          if(!any(Fhat1(km1$time) >= u)){
            Tistar <- admin_t
            imputed_datasets[[1]][[T.name]][i]<- admin_t
            imputed_datasets[[1]][[D.name]][i]<- 0
          }else{
            Tistar <- min(km1$time[Fhat1(km1$time) >= u])
            imputed_datasets[[1]][[T.name]][i]<-ifelse(is.infinite(Tistar) | is.na(Tistar), admin_t, Tistar)
            imputed_datasets[[1]][[D.name]][i]<-ifelse(is.infinite(Tistar) | is.na(Tistar), 0, 1)
          }
        }
        
        if(Ai==0){
          q <- Fhat0(Ti)
          u <- runif(1, min = q, max = 1)
          if(!any(Fhat0(km0$time) >= u)){
            Tistar <- admin_t
            imputed_datasets[[1]][[T.name]][i]<-ifelse(is.infinite(Tistar) | is.na(Tistar), admin_t, Tistar)
            imputed_datasets[[1]][[D.name]][i]<-ifelse(is.infinite(Tistar) | is.na(Tistar), 0, 1)
          }else{
            Tistar <- min(km0$time[Fhat0(km0$time) >= u])
            imputed_datasets[[1]][[T.name]][i]<-ifelse(is.infinite(Tistar) | is.na(Tistar), admin_t, Tistar)
            imputed_datasets[[1]][[D.name]][i]<-ifelse(is.infinite(Tistar) | is.na(Tistar), 0, 1)
          }
        }
      }
    }
    
  }
  #toc()
  
  if(verbose){
    message("Step 1: Gibbs sampler to impute from conditional survival functions")
  }
  if(verbose){
    pb1 <- progress::progress_bar$new(format = "Burn-in [:bar] :percent eta: :eta", 
                                      total = burn_in, 
                                      clear = FALSE, 
                                      width = 100,
                                      show_after = 0,
                                      force=TRUE)
    pb2 <- progress::progress_bar$new(format = "MC samples [:bar] :percent eta: :eta",
                                      total = MC_samples, 
                                      clear = FALSE, 
                                      width = 100,
                                      show_after = 0,
                                      force=TRUE)
  }
  
  # Helper function for weight KM estimator
  wkm <- function(time, status, w){
    
    ord <- order(time)
    time   <- time[ord]
    status <- status[ord]
    
    # Unique event times
    event_times <- sort(unique(time[status == 1]))
    
    # Indices of events at each event time
    event_idx <- split(seq_along(time)[status == 1], time[status == 1])
    
    w <- w[ord]
    
    # Risk set weights (reverse cumulative sum)
    Y <- rev(cumsum(rev(w)))
    
    # Event weights at each event time
    d <- vapply(event_idx,
                function(idx) sum(w[idx]),
                numeric(1))
    
    # Risk set weights at event times
    Y_event <- Y[match(names(d), time)]
    
    # KM estimate
    S <- cumprod(1 - d / Y_event)
    
    list(time = as.numeric(names(d)), surv = S)
  }
  
  for(iter in 1:(burn_in + MC_samples)){
    
    #print(iter)
    if(iter==1){
      pb1$tick(0)
    }
    
    if(iter==burn_in+1){
      pb2$tick(0)
    }
    
    # Initialize values from previous iteration
    imputed_datasets[[(iter+1)]] <- imputed_datasets[[(iter)]]
    
    n = nrow(imputed_datasets[[(iter+1)]])
    p = ncol(imputed_datasets[[(iter+1)]])
    
    # Precompute kernel distance matrix
    # w = vector(mode='list', length=num_outcomes)
    # for(k in 1:num_outcomes){
    #   
    #   if(predictor_transform=="identity"){
    #     if(kern_shape=="epanechnikov"){
    #     bw = 2.34 * stats::sd(imputed_datasets[[(iter+1)]][[T.names[k]]]) * (n)^(-1/(num_outcomes - 1 + 4))}
    #     if(kern_shape=="gaussian"){
    #     bw = 1.06 * stats::sd(imputed_datasets[[(iter+1)]][[T.names[k]]]) * (n)^(-1/(num_outcomes - 1 + 4))}
    #     
    #     z <- outer(imputed_datasets[[(iter+1)]][[T.names[k]]], imputed_datasets[[(iter+1)]][[T.names[k]]], `-`) / bw
    #     z<-as(z, "sparseMatrix")
    #   }else if(predictor_transform=="log"){
    #     
    #     if(kern_shape=="epanechnikov"){
    #       bw = 2.34 * stats::sd(log(imputed_datasets[[(iter+1)]][[T.names[k]]]), na.rm=TRUE) * (n)^(-1/(num_outcomes - 1 + 4))}
    #     if(kern_shape=="gaussian"){
    #       bw = 1.06 * stats::sd(log(imputed_datasets[[(iter+1)]][[T.names[k]]]), na.rm=TRUE) * (n)^(-1/(num_outcomes - 1 + 4))}
    #     z <- outer(log(imputed_datasets[[(iter+1)]][[T.names[k]]]), log(imputed_datasets[[(iter+1)]][[T.names[k]]]), `-`) / bw
    #     z<-as(z, "sparseMatrix")
    #   }else if(predictor_transform=="prob"){
    #       
    #     km1<-survfit(Surv(get(T.names[k]), get(D.names[k])) ~ 1, data=dat %>% filter(get(A.name)==1))
    #     Fhat1 <- stepfun(km1$time, c(0, 1 - km1$surv))
    #     km0<-survfit(Surv(get(T.names[k]), get(D.names[k])) ~ 1, data=dat %>% filter(get(A.name)==0))
    #     Fhat0 <- stepfun(km0$time, c(0, 1 - km0$surv))
    #       
    #     q <- ifelse(imputed_datasets[[(iter+1)]][[A.name]] == 0, 
    #                 Fhat0(imputed_datasets[[(iter+1)]][[T.names[k]]]), 
    #                 Fhat1(imputed_datasets[[(iter+1)]][[T.names[k]]]))
    #     if(kern_shape=="epanechnikov"){
    #       bw = 2.34 * stats::sd(q, na.rm=TRUE) * (n)^(-1/(num_outcomes - 1 + 4))}
    #     if(kern_shape=="gaussian"){
    #       bw = 1.06 * stats::sd(q, na.rm=TRUE) * (n)^(-1/(num_outcomes - 1 + 4))}
    #     z <- outer(q, q, `-`) / bw 
    #     z<-as(z, "sparseMatrix")
    #   }
    # 
    #   if(kern_shape=="epanechnikov"){
    #   w_tmp <- MatrixExtra::mapSparse(z, function(x){
    #     ifelse(abs(x) <= 1, 0.75 * (1 - x^2), 0)}) 
    #   }
    #   if(kern_shape=="gaussian"){
    #   w_tmp <- MatrixExtra::mapSparse(z, function(x){1/bw * dnorm(x)})
    #   }
    #   #w_tmp <- w_tmp / sum(w_tmp)
    #   w[[k]]<-w_tmp
    #   rm(z)
    # }
    
    # Precompute kernel distance matrix
    # w = matrix(1, nrow=n, ncol=n)
    z_list = vector(mode='list', length=num_outcomes)
    for(k in 1:num_outcomes){
      
      if(predictor_transform=="identity"){
        #bw = 1.06 * stats::sd(imputed_datasets[[(iter+1)]][[T.names[k]]]) * (n)^(-1/(2*2 + num_outcomes))
        #print(bw)
        z <- outer(imputed_datasets[[(iter+1)]][[T.names[k]]], imputed_datasets[[(iter+1)]][[T.names[k]]], `-`) / bw
        z<-as(z, "sparseMatrix")
      }else if(predictor_transform=="log"){
        bw = 1.06 * stats::sd(log(imputed_datasets[[(iter+1)]][[T.names[k]]]), na.rm=TRUE) * (n)^(-1/(num_outcomes - 1 + 4))
        z <- outer(log(imputed_datasets[[(iter+1)]][[T.names[k]]]), log(imputed_datasets[[(iter+1)]][[T.names[k]]]), `-`) / bw
        z<-as(z, "sparseMatrix")
      }else if(predictor_transform=="prob"){
        
        km1<-survfit(Surv(get(T.names[k]), get(D.names[k])) ~ 1, data=dat %>% filter(get(A.name)==1))
        Fhat1 <- stepfun(km1$time, c(0, 1 - km1$surv))
        km0<-survfit(Surv(get(T.names[k]), get(D.names[k])) ~ 1, data=dat %>% filter(get(A.name)==0))
        Fhat0 <- stepfun(km0$time, c(0, 1 - km0$surv))
        
        q <- ifelse(imputed_datasets[[(iter+1)]][[A.name]] == 0, 
                    Fhat0(imputed_datasets[[(iter+1)]][[T.names[k]]]), 
                    Fhat1(imputed_datasets[[(iter+1)]][[T.names[k]]]))
          bw = 1.06 * stats::sd(q, na.rm=TRUE) * (n)^(-1/(num_outcomes - 1 + 4))
        z <- outer(q, q, `-`) / bw 
        z<-as(z, "sparseMatrix")
      }
      z_list[[k]] <- z
    }
    
    for(j in 1:num_outcomes){
      
      D.name <- Delta.names[j]
      T.name <- T.names[j]
      
      Tvec <- dat[[T.name]]
      Dvec <- dat[[D.name]]
      
      if (is.null(admin_t) || is.na(admin_t)) {
        admin_t_local <- max(Tvec, na.rm = TRUE)
      } else {
        admin_t_local <- admin_t
      }
      
      indices_to_impute <- which(Tvec < admin_t_local & Dvec == 0)
      
      w <- Reduce(`+`, lapply(z_list[-j], function(z) z^2))
      w  <- exp(-0.5 * w)
      
      for (i in indices_to_impute) {
        
        Ai <- dat[[A.name]][i]
        Ti <- Tvec[i]
        
        weights <- w[, i]
        
        km <- wkm(
          dat[[T.name]][dat[[A.name]] == Ai],
          dat[[D.name]][dat[[A.name]] == Ai],
          w = weights[dat[[A.name]] == Ai]
        )
        
        Fhat <- stepfun(km$time, c(0, 1 - km$surv))
        
        q <- Fhat(Ti)
        u <- runif(1, min = q, max = 1)
        
        if (!any(Fhat(km$time) >= u)) {
          Tvec[i] <- admin_t_local
          Dvec[i] <- 0
        } else {
          Tistar <- min(km$time[Fhat(km$time) >= u])
          if (is.infinite(Tistar) || is.na(Tistar)) {
            Tvec[i] <- admin_t_local
            Dvec[i] <- 0
          } else {
            Tvec[i] <- Tistar
            Dvec[i] <- 1
          }
        }
      }
      
      imputed_datasets[[iter + 1]][[T.names[j]]] <- Tvec
      imputed_datasets[[iter + 1]][[Delta.names[j]]] <- Dvec
      
      # update kernel weights
      # recompute kernel matrix
      w = matrix(1, nrow=n, ncol=n)
      #z_list = vector(mode='list', length=num_outcomes)
      #for(k in 1:num_outcomes){
      if(predictor_transform=="identity"){
        #bw = 1.06 * stats::sd(imputed_datasets[[(iter+1)]][[T.names[k]]]) * (n)^(-1/(2*2 + num_outcomes))
        #print(bw)
        z <- outer(imputed_datasets[[(iter+1)]][[T.names[j]]], 
                   imputed_datasets[[(iter+1)]][[T.names[j]]], `-`) / bw
        z<-as(z, "sparseMatrix")
      }else if(predictor_transform=="log"){
        bw = 1.06 * stats::sd(log(imputed_datasets[[(iter+1)]][[T.names[k]]]), na.rm=TRUE) * (n)^(-1/(num_outcomes - 1 + 4))
        z <- outer(log(imputed_datasets[[(iter+1)]][[T.names[j]]]), log(imputed_datasets[[(iter+1)]][[T.names[j]]]), `-`) / bw
        z<-as(z, "sparseMatrix")
      }else if(predictor_transform=="prob"){
          
        km1<-survfit(Surv(get(T.names[j]), get(D.names[j])) ~ 1, data=dat %>% filter(get(A.name)==1))
        Fhat1 <- stepfun(km1$time, c(0, 1 - km1$surv))
        km0<-survfit(Surv(get(T.names[j]), get(D.names[j])) ~ 1, data=dat %>% filter(get(A.name)==0))
        Fhat0 <- stepfun(km0$time, c(0, 1 - km0$surv))
          
          q <- ifelse(imputed_datasets[[(iter+1)]][[A.name]] == 0, 
                      Fhat0(imputed_datasets[[(iter+1)]][[T.names[j]]]), 
                      Fhat1(imputed_datasets[[(iter+1)]][[T.names[j]]]))
          bw = 1.06 * stats::sd(q, na.rm=TRUE) * (n)^(-1/(num_outcomes - 1 + 4))
          z <- outer(q, q, `-`) / bw 
          z<-as(z, "sparseMatrix")
        }
      
      z_list[[j]] <- z
    }
    
    if(iter<=burn_in){
      pb1$tick()
    }else{
      pb2$tick()
    }
    
  }
  
  return(
    imputed_datasets
  )
}
  
#library(RobinCar)

NCI_adj <- function(imputed_datasets, indices_to_analyze, Y.name, A.name, adj_vars, eval_t, calibrate = TRUE, SL=FALSE){
  
  #imputed_datasets = imputed_datasets
  #indices_to_analyze = tail(1:length(imputed_datasets), 5)
  #eval_t=10
  #A.name="A"
  #T.names = c("Y", "N1", "N2", "N3")
  #Delta.names = c("DeltaY", "DeltaN1", "DeltaN2", "DeltaN3")
  #Y.name = T.names[1]
  #DeltaY.name = T.names[1]
  
  formula=as.formula(paste0("V ~ ", A.name, "* (", paste0(adj_vars, collapse="+"), ")"))
  
  #res_RR<-rep(0, length(indices_to_analyze))
  #res_RR_var<-rep(0, length(indices_to_analyze))
  res_logRR<-rep(0, length(indices_to_analyze))
  res_logRR_var<-rep(0, length(indices_to_analyze))
  
  # log_ratio <- function(est){
  #   base <- est[1] / (1-est[1])
  #   cont <- log(est[-1] / (1-est[-1]) / base)
  #   return(cont)
  # }
  # 
  # jacobian.LOGRATIO <- function(settings, est){
  #   col1 <- -1 / est[1]
  #   mat1 <- diag(1/est[-1], length(est)-1)
  #   jacob <- cbind(col1, mat1)
  #   return(jacob)
  # }
  
  message("Pooling MI analysis results using Rubin's Rules")
  if(calibrate){message("Using calibration to guarantee efficiency gain")}
  for(i in 1:length(indices_to_analyze)){
    
    indx = indices_to_analyze[i]
    
    if(SL==FALSE){
    rc_glm <-imputed_datasets[[indx]]  %>%
      mutate(A = factor(get(A.name), levels=c(0,1)),
             V = ifelse(get(Y.name) < eval_t, 1, 0)) %>%
      RobinCar::robincar_glm(df = .,
                             response_col="V",
                             treat_col=A.name,
                             formula=formula,
                             car_scheme="simple",
                             g_family=binomial(link="logit"),
                             g_accuracy=7)
    }else{
      
      rc_glm <- imputed_datasets[[indx]]  %>%
        mutate(A = factor(get(A.name), levels=c(0,1)),
               V = ifelse(get(Y.name) < eval_t, 1, 0)) %>%
        robincar_SL(df=.,
                    response_col="V",
                    treat_col=A.name,
                    covariate_cols = adj_vars,
                    SL_libraries=c("SL.ranger", "SL.glm"),
                    car_scheme="simple",
                    covariate_to_include_strata=FALSE)
      
    }
    
    if(calibrate){
    rc_glm <- robincar_calibrate(
      result=rc_glm,
      joint=FALSE
    )
    }
    #res_tmp<-RobinCar::robincar_contrast(result=rc_glm, contrast_h = "ratio")
    #res_RR[i]<-res_tmp$result$estimate
    #res_RR_var[i]<-res_tmp$varcov
    
    res_tmp<-RobinCar::robincar_contrast(result=rc_glm, contrast_h = "log_ratio")
    res_logRR[i]<-res_tmp$result$estimate
    res_logRR_var[i]<-res_tmp$varcov
  }
  
  logRR_MI = mean(res_logRR)
  V_logRR_MI = mean(res_logRR_var) + (length(indices_to_analyze)+1)/length(indices_to_analyze) * var(res_logRR)
  
  return(
  list(
    logRR = logRR_MI,
    var_logRR = V_logRR_MI,
    M = length(indices_to_analyze),
    within_imputation_variance = mean(res_logRR_var),
    between_imputation_variance = (length(indices_to_analyze)+1)/length(indices_to_analyze) * var(res_logRR)
  )
  )
  
}


NCI_adj_2 <- function(imputed_datasets, indices_to_analyze, original_dat, Y.name, D.name, A.name, adj_vars, eval_t, calibrate = TRUE){
  
  # imputed_datasets = gs_single
  # eval_t=2
  # A.name="A"
  # adj.vars = c("N1", "N2", "N3", "DeltaN1", "DeltaN2", "DeltaN3")
  # indices_to_analyze <- tail(seq(1, L$burn_in + 1 + L$MC_samples), L$MC_samples)
  
  #formula=as.formula(paste0("Surv(", Y.name, ",", D.name, ")~ ", A.name, "* (", paste0(adj_vars, collapse="+"), ")"))
  
  #res_RR<-rep(0, length(indices_to_analyze))
  #res_RR_var<-rep(0, length(indices_to_analyze))
  res_logRR<-rep(0, length(indices_to_analyze))
  res_logRR_var<-rep(0, length(indices_to_analyze))
  
  # log_ratio <- function(est){
  #   base <- est[1] / (1-est[1])
  #   cont <- log(est[-1] / (1-est[-1]) / base)
  #   return(cont)
  # }
  # 
  # jacobian.LOGRATIO <- function(settings, est){
  #   col1 <- -1 / est[1]
  #   mat1 <- diag(1/est[-1], length(est)-1)
  #   jacob <- cbind(col1, mat1)
  #   return(jacob)
  # }
  
  message("Pooling MI analysis results using Rubin's Rules")
  for(i in 1:length(indices_to_analyze)){
    print(i)
    indx = indices_to_analyze[i]
    
    rc <-imputed_datasets[[indx]]  #%>%
      #mutate(A = factor(get(A.name), levels=c(0,1)))
    
    rc<-bind_cols(data %>% dplyr::select(Y, DeltaY), rc %>% dplyr::select(-c("Y", "DeltaY")))
    
    fit <- CFsurvival(
          time = rc$Y,
          event = data$DeltaY,
          treat = data$A,
          confounders = data %>% dplyr::select(starts_with("N"), starts_with("DeltaN")) %>% as.matrix(),
          fit.times = c(sort(unique(data$Y[data$DeltaY==1])), 2),
          contrasts = "risk.ratio",
          conf.band=FALSE,
        nuisance.options = list(
          cross.fit=TRUE,
          V=3,
          event.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg",
                             "survSL.weibreg", "survSL.loglogreg"), function(alg)
                             c(alg, "All")),
          #event.SL.library = "survSL.coxph",
          cens.SL.library = lapply(c("survSL.km", "survSL.coxph", "survSL.expreg",
                                      "survSL.weibreg", "survSL.loglogreg"), function(alg)
                                        c(alg, "All")),
          prop.SL.library=NULL,
          prop.pred  = rep(0.5, length(data$Y))
        )
      )
    
    res_tmp<-fit$risk.ratio.df[nrow(fit$risk.ratio.df),]
    res_logRR[i]<-res_tmp$risk.ratio
    res_logRR_var[i]<-res_tmp$se.log.ratio^2
  }
  
  logRR_MI = mean(res_logRR)
  V_logRR_MI = mean(res_logRR_var) + (length(indices_to_analyze)+1)/length(indices_to_analyze) * var(res_logRR)
  
  return(
    list(
      logRR = logRR_MI,
      var_logRR = V_logRR_MI,
      M = length(indices_to_analyze))
  )
  
}
# res <- NCI_adj(imputed_datasets, 
#         indices_to_analyze = tail(1:length(imputed_datasets), 10),
#         Y.name = "Y",
#         A.name = "A",
#         adj_vars = c("N1", "N2", "N3", "DeltaN1", "DeltaN2", "DeltaN3"),
#         eval_t = 2
#         )
# 
# KM_est(test, T.name = "Y", Delta.name = "DeltaY", A.name = "A", eval_t = 2)

# df_all <- bind_rows(
#   lapply(1:(burn_in+1), function(i) {
#     bind_cols(impute_T[[i]] |> mutate(frame = i), impute_Delta[[i]] %>% dplyr::select(-"A"))
#   })
# )
# 
# require(gganimate)
# p <- ggplot(df_all, aes(x=N1, y=Y, color=interaction(DeltaN1, DeltaY))) +
#   geom_point() +
#   transition_states(frame, transition_length = 1, state_length = 3) +
#   labs(title = "Iteration {closest_state}")
# 
# animate(p, fps = 5, duration = 10)
# 
# df_all %>% 
#   pivot_longer(cols=c("Y", "N1", 'N2', "N3")) %>%
#   group_by(frame, name) %>%
#   summarise(
#     prob = mean(value < 2)
#   ) %>%
#   ggplot(aes(x=frame, y=prob, color=name))+
#   geom_line()+
#   scale_x_continuous(breaks=seq(0, burn_in+1,5))