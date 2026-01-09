#library(copula)

generate_lambdas<-function(eval_time=1, per_year_incidence=0.05, RR){
  
  #Derived from CDF of exponential distn
  lambda_0 = - log(1-per_year_incidence)
  lambda_1 = -log(1-(per_year_incidence * RR))
  
  return(c("lambda_0"=lambda_0, "lambda_1"=lambda_1))
}

create_multivar_survival <- function(n,
                        trt_prob,
                        RR,
                        RR_C,
                        per_year_incidence_Y,
                        per_year_incidence_N,
                        per_year_incidence_C,
                        J=1,
                        rho,
                        admin_date=2,
                        convert_days = TRUE
){
  
  require(tidyverse)
  require(copula)
  # n=1000
  # per_year_incidence_Y=0.05
  # RR=0.5
  # per_year_incidence_N=0.10
  # per_year_incidence_C=0.10
  # RR_C= 4
  # admin_date=2
  # J=1
  # rho=0.7
  # trt_prob=0.5
  
  lambdas_Y=generate_lambdas(eval_time=1, per_year_incidence=per_year_incidence_Y, RR=RR)
  lambdas_N=generate_lambdas(eval_time=1, per_year_incidence = per_year_incidence_N, RR=1)
  lambdas_C=generate_lambdas(eval_time=1, per_year_incidence = per_year_incidence_C, RR = RR_C)
  
  Sigma <- diag(1, J+2)
  #Sigma[1,2]<-Sigma[2,1]<-1
  Sigma[row(Sigma) != col(Sigma)] <- rho
  #Sigma[c(1,3:ncol(Sigma)),2]<-Sigma[2:ncol(Sigma),1]<-Sigma[2,c(1,3:ncol(Sigma))]<-Sigma[1,2:ncol(Sigma)]<-rho
  
  Sigma = Matrix::nearPD(Sigma, corr=FALSE)$mat
  
  cp = normalCopula(param=P2p(Sigma), dim = J+2, dispstr = "un")
  
  v<-c(lambdas_Y[1], lambdas_Y[2], rep(lambdas_N[1], J))
  names(v)<-rep("rate", J+2)
  
  params=split(v, factor(c("Y0", "Y1", paste0("N", 1:J)), levels=c("Y0", "Y1", paste0("N", 1:J))))
  
  myMvd = mvdc(copula=cp, margins=rep(c("exp"), J+2),
               paramMargins=params, marginsIdentical=FALSE)
  
  oracle_vars = rMvdc(n=n, myMvd)
  
  A<-c(rep(1, n*trt_prob), rep(0, n-c(n*trt_prob)))
  
  observed_vars = ifelse(
    A==0, oracle_vars[,1], oracle_vars[,2]
  )
  
  if(convert_days){
    observed_vars<- data.frame(A=A, Y=pmax(1, round(observed_vars * 365)))
  }else{
    observed_vars<- data.frame(A=A, Y=observed_vars)
  }
  
  for(i in 1:J){
    if(convert_days){
    observed_vars[[paste0("N", i)]]<-pmax(1, round(oracle_vars[,i+2] * 365))
    }else{
    observed_vars[[paste0("N", i)]]<-oracle_vars[,i+2]
    }
  }
  
  #censoring
  Sigma<-rbind(c(1,1), c(1,1))
  cp = normalCopula(param=P2p(Sigma), dim = 2, dispstr = "un")
  myMvd = mvdc(copula=cp, margins=rep(c("exp"), 2),
               paramMargins=list(lambdas_C[1], lambdas_C[2]), marginsIdentical=FALSE)
  cens_times = rMvdc(n=n, myMvd)
  observed_cens = ifelse(
    A==0, cens_times[,1], cens_times[,2]
  )
  
  if(convert_days){
    observed_cens <- pmax(1, round(pmin(observed_cens, admin_date) * 365))
  }else{
    observed_cens <- pmin(observed_cens, admin_date)
  }
  
  #observed_vars$DeltaY<-ifelse(observed_vars$Y<=observed_cens, observed_vars$DeltaY, 0)
  #observed_vars$Y <- ifelse(observed_vars$Y<=observed_cens, observed_vars$Y, observed_cens)
  for(i in 1:J){
    
    for(j in 1:nrow(observed_vars)){
      
      if(i==1){
        observed_vars$DeltaY[j]<-ifelse(
          observed_vars$Y[j] > observed_cens[j],
          0,
          1)
        
        observed_vars$Y[j] <- ifelse(
          observed_vars$Y[j] > observed_cens[j],
          observed_cens[j],
          observed_vars$Y[j])
      }
      
      
      observed_vars[[paste0("DeltaN", i)]][j]<- ifelse(
        observed_vars[[paste0("N", i)]][j] > observed_cens[j],
        0,
        1)
      
      observed_vars[[paste0("N", i)]][j]<- ifelse(
        observed_vars[[paste0("N", i)]][j] > observed_cens[j],
        observed_cens[j],
        observed_vars[[paste0("N", i)]][j])
      
    }
  }
  
  return(observed_vars)
}

# set.seed(2)
# test <- create_multivar_survival(n=3000,
#                     trt_prob=0.5,
#                     RR=0.7,
#                     RR_C = 1.7,
#                     per_year_incidence_Y = 0.04,
#                     per_year_incidence_N = 0.15,
#                     per_year_incidence_C = 0.10,
#                     J=J,
#                     rho=0.75,
#                     admin_date=2)