library(kableExtra)
library(knitr)
library(survival)

source("./create_multivar_survival.r")

lRR <- rep(NA, 300)

# POI cannot be exactly controlled through simulations
# Simulate several large datasets and obtain approximate truth

set.seed(47)

for(i in 1:length(lRR)){
print(i)

rho_tmp <- case_when(
  i<=100 ~ 0,
  i>=101 & i<=200 ~ 0.5,
  i>=201 ~ 0.8
)
  
data <- create_multivar_survival(n=1000000,
                                 trt_prob=0.5,
                                 RR = 0.75,
                                 RR_C = 1,
                                 per_year_incidence_Y = 0.06,
                                 per_year_incidence_N = 0.12,
                                 per_year_incidence_C = 0.10,
                                 J = 1,
                                 rho = rho_tmp,
                                 admin_date = 2,
                                 X=TRUE,
                                 rho_X = 0,
                                 convert_days = TRUE
)

l<-summary(survfit(Surv(Y, DeltaY) ~ A, data), times=365*2)
lRR[i] = log((1-l$surv[2])/(1-l$surv[1]))
#se<-msm::deltamethod(~(1-x2)/(1-x1), mean=l$surv, cov = diag(l$std.err^2), ses=TRUE)
}
lRRs <- c(mean(lRR[1:100]), mean(lRR[101:200]), mean(lRR[201:300]))


lRR <- rep(NA, 3000)

# POI cannot be exactly controlled through simulations
# Simulate several large datasets and obtain approximate truth

set.seed(47)

for(i in 1:length(lRR)){
  print(i)
  
  rho_tmp <- case_when(
    i<=1000 ~ 0,
    i>=1001 & i<=2000 ~ 0.5,
    i>=2001 ~ 0.8
  )
  
  data <- create_multivar_survival(n=50000,
                                   trt_prob=0.5,
                                   RR = 0.75,
                                   RR_C = 1,
                                   per_year_incidence_Y = 0.06,
                                   per_year_incidence_N = 0.12,
                                   per_year_incidence_C = 0.10,
                                   J = 1,
                                   rho = rho_tmp,
                                   admin_date = 2,
                                   X=TRUE,
                                   rho_X = 0,
                                   convert_days = TRUE
  )
  
  l<-summary(survfit(Surv(Y, DeltaY) ~ A, data), times=365*2)
  lRR[i] = log((1-l$surv[2])/(1-l$surv[1]))
  #se<-msm::deltamethod(~(1-x2)/(1-x1), mean=l$surv, cov = diag(l$std.err^2), ses=TRUE)
}

lRRs <- c(mean(lRR[1:1000]), mean(lRR[1001:2000]), mean(lRR[2001:3000]))
lRRs <- mean(lRR)

# files <- list.files("~/Desktop/sims/sim_results/")
# files<-files[grepl(".rds", files) & grepl("r_", files)]
# ids<-as.numeric(gsub("sim_uid_", "", names(sim$results_complex)))

# for(i in 1:length(files)){
#   print(i)
#   tmp<-readRDS(paste0("~/Desktop/sims/sim_results/", files[i]))
#   
#   for(j in 1:length(tmp$results_complex)){
#     
#     if(i==1 & j==1){
#       res <- bind_cols(tmp$results[j,], tmp$results_complex[[j]]$output) |>
#              mutate(
#                bias = logRR - lRR,
#                cover = as.numeric(
#                         logRR_lci <= lRR & lRR <= logRR_uci)
#              )
#     }else{
#       res <- bind_rows(res, bind_cols(tmp$results[j,], tmp$results_complex[[j]]$output)) |>
#       mutate(
#         bias = logRR - lRR,
#         cover = as.numeric(
#           logRR_lci <= lRR & lRR <= logRR_uci)
#       )
#     }
#   }
#   
# }

#res<-left_join(as.data.frame(res), sim$results, by='sim_uid')
  

### Start here

files <- list.files("~/Desktop/sims/sim_results/")
files<-files[grepl(".rds", files) & grepl("r_", files)]

#lRR<-mean(lRR)

for(i in 1:length(files)){
  print(i)
  tmp<-readRDS(paste0("~/Desktop/sims/sim_results/", files[i]))

  for(j in 1:length(tmp$results_complex)){

    if(i==1 & j==1){
      
    #lRR_tmp <- ifelse(tmp$results[j,]$rho==0, lRRs[1], ifelse(tmp$results[j,]$rho==0.5, lRRs[2], lRRs[3]))
    lRR_tmp <- lRRs
       
    res <- bind_cols(tmp$results[j,], tmp$results_complex[[j]]$output) |>
           mutate(
             bias = (logRR - lRR_tmp)/lRR_tmp,
             cover = as.numeric(
                      logRR_lci <= lRR_tmp & lRR_tmp <= logRR_uci)
           )
    }else{
    
    #lRR_tmp <- ifelse(tmp$results[j,]$rho==0, lRRs[1], ifelse(tmp$results[j,]$rho==0.5, lRRs[2], lRRs[3]))  
    lRR_tmp <- lRRs
      
    res <- bind_rows(res, bind_cols(tmp$results[j,], tmp$results_complex[[j]]$output)) |>
    mutate(
      bias = (logRR - lRR_tmp)/lRR_tmp,
      cover = as.numeric(
        logRR_lci <= lRR_tmp & lRR_tmp <= logRR_uci)
    )
    }
  }

}

library(knitr)

res %>%
  filter(per_year_incidence_C==0.10) %>%
  mutate(type=factor(type, levels=c("KM", 
                                    "X-Adj", 
                                    "NCI-adj (Oracle)",
                                    "NCI-adj (CF)"))) %>%
  group_by(type, n, rho) %>%
  summarise(
    #n_sim = n(),
    bias = mean(bias, na.rm=TRUE),
    releff = median(releff, na.rm=TRUE),
    coverage = mean(cover, na.rm=TRUE)
  ) %>%
  ungroup() %>%
  pivot_longer(cols=c("bias", "releff", "coverage"), names_to="param", values_to="val") %>%
  pivot_wider(names_from=c("rho", "param"), values_from="val") %>%
  arrange(type) %>%
  #dplyr::select(-"n") %>%
  mutate(across(c(3,5,6,8,9,11), function(x){signif(x*100,3)})) %>%
  mutate(across(c(4,7,10), function(x){signif(x,3)})) %>%
  arrange(n) %>%
  kable(format='latex', digits=3)

# res_1 %>%
#   filter(per_year_incidence_C==0.10) %>%
#   mutate(type=factor(type, levels=c("KM", 
#                                     "X-Adj", 
#                                     "NCI-adj (IF)", 
#                                     "NCI-adj (DXF IF)"))) %>%
#   group_by(type, n, rho) %>%
#   summarise(
#     bias = mean(bias),
#     releff = mean(releff, na.rm=TRUE),
#     coverage = mean(cover, na.rm=TRUE)
#   ) %>%
#   ungroup() %>%
#   pivot_longer(cols=c("bias", "releff", "coverage"), names_to="param", values_to="val") %>%
#   pivot_wider(names_from=c("rho", "param"), values_from="val") %>%
#   arrange(type) %>%
#   #dplyr::select(-"n") %>%
#   mutate(across(c(5,8,11), function(x){signif(x*100,3)})) %>%
#   mutate(across(c(3,4,6,7,9,10), function(x){signif(x,3)})) %>%
#   kable(format='latex')

sim<-readRDS("~/Desktop/sims/sim_earlier.rds")

ids<-as.numeric(gsub("sim_uid_", "", names(sim$results_complex)))

for(i in 1:length(sim$results_complex)){
  print(i)
  if(i==1){
    res <- cbind(
      sim$results_complex[[i]]$output #|>
        # mutate(
        #   bias = logRR - mean(log(RR)),
        #   cover = as.numeric(
        #     logRR_lci <= mean(log(RR)) & mean(log(RR)) <= logRR_uci
        #   )
        #)
      , 'sim_uid'=ids[i]
    )
  }else{
    
    res <- rbind(res,
                 cbind(
                   sim$results_complex[[i]]$output #|>
                     # mutate(
                     #   bias = logRR - log(mean(RR)),
                     #   cover = as.numeric(
                     #     logRR_lci <= log(mean(RR)) & log(mean(RR)) <= logRR_uci
                     #   )
                     #)
                   , 'sim_uid'=ids[i]
                 )
    )
    
  }
}

res<-left_join(as.data.frame(res), sim$results, by='sim_uid')

# sim <- readRDS("~/Desktop/sims/sim.rds")
# ids<-as.numeric(gsub("sim_uid_", "", names(sim$results_complex)))
# 
# for(i in 1:length(sim$results_complex)){
#   print(i)
#   if(i==1){
#     res_2 <- cbind(
#       sim$results_complex[[i]]$output #|>
#       # mutate(
#       #   bias = logRR - mean(log(RR)),
#       #   cover = as.numeric(
#       #     logRR_lci <= mean(log(RR)) & mean(log(RR)) <= logRR_uci
#       #   )
#       #)
#       , 'sim_uid'=ids[i]
#     )
#   }else{
#     
#     res_2 <- rbind(res_2,
#                  cbind(
#                    sim$results_complex[[i]]$output #|>
#                    # mutate(
#                    #   bias = logRR - log(mean(RR)),
#                    #   cover = as.numeric(
#                    #     logRR_lci <= log(mean(RR)) & log(mean(RR)) <= logRR_uci
#                    #   )
#                    #)
#                    , 'sim_uid'=ids[i]
#                  )
#     )
#     
#   }
# }
# 
# res_2<-left_join(as.data.frame(res_2), sim$results, by='sim_uid')
# 
# res <- bind_rows(res, res_2)

res %>%
  filter(per_year_incidence_C==0.10) %>%
  mutate(type=factor(type, levels=c("KM", 
                                    "X-Adj",
                                    "NCI-adj (IF)",
                                    "NCI-adj (DXF IF)"))) %>%
  group_by(type, n, rho) %>%
  summarise(
    #nsim = n(),
    bias = mean(bias),
    releff = median(releff, na.rm=TRUE),
    coverage = mean(cover, na.rm=TRUE)
  ) %>%
  ungroup() %>%
  pivot_longer(cols=c("bias", "releff", "coverage"), names_to="param", values_to="val") %>%
  pivot_wider(names_from=c("rho", "param"), values_from="val") %>%
  arrange(type) %>%
  #dplyr::select(-"n") %>%
  mutate(across(c(5,8,11), function(x){signif(x*100,3)})) %>%
  mutate(across(c(3,4,6,7,9,10), function(x){signif(x,3)})) %>%
  arrange(n) %>%
  kable(format="latex")



sim<-readRDS("~/Desktop/sims/sim_MOSTRECENT_dxf.rds")

ids<-as.numeric(gsub("sim_uid_", "", names(sim$results_complex)))

for(i in 1:length(sim$results_complex)){
  print(i)
  if(i==1){
    res <- cbind(
      sim$results_complex[[i]]$output |>
      mutate(
        bias = logRR - mean(log(RR)),
        cover = as.numeric(
          logRR_lci <= mean(log(RR)) & mean(log(RR)) <= logRR_uci
        )
      )
      , 'sim_uid'=ids[i]
    )
  }else{
    
    res <- rbind(res,
                 cbind(
                   sim$results_complex[[i]]$output |>
                   mutate(
                     bias = logRR - log(mean(RR)),
                     cover = as.numeric(
                       logRR_lci <= log(mean(RR)) & log(mean(RR)) <= logRR_uci
                     )
                   )
                   , 'sim_uid'=ids[i]
                 )
    )
    
  }
}

res<-left_join(as.data.frame(res), sim$results, by='sim_uid')

res %>%
  filter(per_year_incidence_C==0.10) %>%
  mutate(type=factor(type, levels=c("KM", 
                                    "X-Adj", 
                                    "NCI-adj (IF)", 
                                    "NCI-adj (DXF IF)"))) %>%
  group_by(type, n, rho) %>%
  summarise(
    nsim = n(),
    bias = mean(bias),
    releff = mean(releff, na.rm=TRUE),
    coverage = mean(cover, na.rm=TRUE)
  ) %>%
  ungroup() %>%
  pivot_longer(cols=c("bias", "releff", "coverage"), names_to="param", values_to="val") %>%
  pivot_wider(names_from=c("rho", "param"), values_from="val") %>%
  arrange(type) %>%
  #dplyr::select(-"n") %>%
  mutate(across(c(5,8,11), function(x){signif(x*100,3)})) %>%
  mutate(across(c(3,4,6,7,9,10), function(x){signif(x,3)}))
