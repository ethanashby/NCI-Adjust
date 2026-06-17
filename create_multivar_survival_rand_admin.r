generate_lambdas<-function(eval_time=1, per_year_incidence=0.05, RR){
  
  #Derived from CDF of exponential distn
  lambda_0 = - log(1-per_year_incidence)
  lambda_1 = -log(1-(per_year_incidence * RR))
  
  return(c("lambda_0"=lambda_0, "lambda_1"=lambda_1))
}

create_multivar_survival_var_admin <- function(n,
                                     trt_prob,
                                     RR,
                                     RR_C,
                                     per_year_incidence_Y,
                                     per_year_incidence_N,
                                     per_year_incidence_C,
                                     J=1,
                                     rho,
                                     admin_date=2,
                                     convert_days = TRUE, 
                                     X = TRUE,
                                     rho_X = 0.3
){
  
  require(tidyverse)
  require(copula)
  
  lambdas_Y=generate_lambdas(eval_time=1, per_year_incidence=per_year_incidence_Y, RR=RR)
  lambdas_N=generate_lambdas(eval_time=1, per_year_incidence = per_year_incidence_N, RR=1)
  lambdas_C=generate_lambdas(eval_time=1, per_year_incidence = per_year_incidence_C, RR = RR_C)
  lambdas<-c(lambdas_Y, lambdas_C, lambdas_N)
  
  if(X==TRUE){
    Sigma <- diag(1, 1+1+1+1+1+J) #Y0, Y1, C0, C1, X, J
    Sigma[1,2]<-Sigma[2,1]<-rho
    Sigma[1,4+seq(J)]<-Sigma[4+seq(J),1]<-Sigma[2,4+seq(J)]<-Sigma[4+seq(J),2]<-rho
    Sigma[ncol(Sigma),1:(ncol(Sigma)-1)]<-Sigma[1:(ncol(Sigma)-1),ncol(Sigma)]<-rho_X
    
    Sigma = Matrix::nearPD(Sigma, corr=FALSE)$mat
    
    cp = normalCopula(param=P2p(Sigma), dim = 5+J, dispstr = "un")
    
    v<-vector(mode="list", length=5+J)
    for(i in 1:length(v)){
      if(i %in% c(1,2)){
        v[[i]]<-c("rate"=as.numeric(lambdas[i]))
      }else if(i %in% c(3,4)){
        v[[i]]<-c("rate"=as.numeric(lambdas[i]))
      }else if(i > 4 & i < length(v)){
        v[[i]]<-c("rate"=as.numeric(lambdas[i]))
      }else{
        v[[i]]<-c('mean'=0, 'sd'=1)
      }
    }
    
    myMvd = mvdc(copula=cp, margins=c(rep(c("exp"), J+4), "norm"),
                 paramMargins=v, marginsIdentical=FALSE)
    
    oracle_vars = rMvdc(n=n, myMvd)
  }else{
    Sigma <- diag(1, 1+1+1+1+J) #Y0, Y1, C0, C1, J
    Sigma[1,2]<-Sigma[2,1]<-rho
    Sigma[1,4+seq(J)]<-Sigma[4+seq(J),1]<-Sigma[2,4+seq(J)]<-Sigma[4+seq(J),2]<-rho
    #Sigma[ncol(Sigma),1:(ncol(Sigma)-1)]<-Sigma[1:(ncol(Sigma)-1),ncol(Sigma)]<-rho_X
    
    Sigma = Matrix::nearPD(Sigma, corr=FALSE)$mat
    
    cp = normalCopula(param=P2p(Sigma), dim = 4+J, dispstr = "un")
    
    v<-vector(mode="list", length=4+J)
    for(i in 1:length(v)){
      if(i %in% c(1,2)){
        v[[i]]<-c("rate"=as.numeric(lambdas[i]))
      }else if(i %in% c(3,4)){
        v[[i]]<-c("rate"=as.numeric(lambdas[i]))
      }else{
        v[[i]]<-c("rate"=as.numeric(lambdas[i]))
      }
    }
    
    myMvd = mvdc(copula=cp, margins=rep(c("exp"), J+4),
                 paramMargins=v, marginsIdentical=FALSE)
    
    oracle_vars = rMvdc(n=n, myMvd)
  }
  
  A<-c(rep(1, n*trt_prob), rep(0, n-c(n*trt_prob)))
  
  Y = ifelse(
    A==0, oracle_vars[,1], oracle_vars[,2]
  )
  
  C = ifelse(
    A==0, oracle_vars[,3], oracle_vars[,4]
  )
  
  admin_dates = runif(length(Y), admin_date-0.05*admin_date, admin_date+0.05*admin_date)
  observed_Y = pmin(pmin(Y, C), admin_dates)
  DeltaY = as.numeric(observed_Y == Y)
  
  if(X){
    if(convert_days){
      observed_vars<- data.frame(A=A, 
                                 Y=pmax(1, round(observed_Y * 365)),
                                 DeltaY = DeltaY,
                                 X=oracle_vars[,ncol(oracle_vars)])
    }else{
      observed_vars<- data.frame(A=A, 
                                 Y=observed_Y, 
                                 DeltaY = DeltaY,
                                 X=oracle_vars[,ncol(oracle_vars)])
    }
  }else{
    if(convert_days){
      observed_vars<- data.frame(A=A, Y=pmax(1, round(observed_Y * 365)), DeltaY=DeltaY)
    }else{
      observed_vars<- data.frame(A=A, Y=observed_Y, DeltaY=DeltaY)
    }
  }
  
  for(i in 1:J){
    if(convert_days){
      observed_vars[[paste0("N", i)]]<-pmax(1, round(pmin(C, oracle_vars[,i+4], admin_dates) * 365))
      observed_vars[[paste0("DeltaN", i)]]<-as.numeric(pmin(C, oracle_vars[,i+4], admin_dates)==oracle_vars[,i+4])
      observed_vars[['admin_t']]=round(admin_dates*365)
    }else{
      observed_vars[[paste0("N", i)]]<-pmin(C, oracle_vars[,i+4], admin_dates)
      observed_vars[[paste0("DeltaN", i)]]<-as.numeric(pmin(C, oracle_vars[,i+4], admin_dates)==oracle_vars[,i+4])
      observed_vars[['admin_t']]=admin_dates
    }
  }
  
  return(observed_vars)
}

# set.seed(2)
# test <- create_multivar_survival_var_admin(n=3000,
#                     trt_prob=0.5,
#                     RR=0.7,
#                     RR_C = 1,
#                     per_year_incidence_Y = 0.06,
#                     per_year_incidence_N = 0.15,
#                     per_year_incidence_C = 0.10,
#                     J=1,
#                     rho=0,
#                     admin_date=2)
