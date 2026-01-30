sim<-readRDS("~/Desktop/sim.rds")

files<-list.files("~/Desktop/res")
for(i in 1:length(files)){
  
  tmp<-readRDS(paste0("~/Desktop/res/", files[i]))
  uids<-as.numeric(gsub("sim_uid_", "", names(tmp$results_complex)))
  #match(uids, tmp$results$sim_uid)
  for(j in 1:length(tmp$results_complex)){
    if(i==1 & j==1){
    out<-bind_cols(tmp$results[match(uids[j], tmp$results$sim_uid),], tmp$results_complex[[j]]$output)
    }else{
    out <- bind_rows(out, bind_cols(tmp$results[match(uids[j], tmp$results$sim_uid),], tmp$results_complex[[j]]$output))
    }
  }
}

# ids<-as.numeric(gsub("sim_uid_", "", names(sim$results_complex)))
# 
# for(i in 1:length(sim$results_complex)){
#   
#   if(i==1){
#   res<-cbind(sim$results_complex[[i]]$output, 'sim_uid'=ids[i])
#   }else{
#   res<-rbind(res, cbind(sim$results_complex[[i]]$output, 'sim_uid'=ids[i]))
#   }
# }
# 
# res<-left_join(as.data.frame(res), sim$results, by='sim_uid')

# res %>%
#   group_by(type, n, rho, RR_C) %>%
#   summarise(
#     MSE = mean(bias^2)
#   ) %>%
#   pivot_wider(names_from=c("type", "RR_C"), values_from="MSE") %>%
#   dplyr::select(c(3,7,5,4,8,6)) %>%
#   kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
#   kable_styling(
#     latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
#     full_width = FALSE # Prevent the table from stretching to full page width
#   )

out %>%
  filter(rho_X==0.1) %>%
  group_by(type, n, per_year_incidence_C) %>%
  summarise(
    #n = n(),
    #prop_is_na = mean(is.na(cover)),
    cover = mean(cover, na.rm=TRUE)
  ) %>%
  pivot_wider(names_from="type", values_from="cover") %>%
  dplyr::select(c(1,2,3,5,4)) %>%
  kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
  kable_styling(
    latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
    full_width = FALSE # Prevent the table from stretching to full page width
  )

out %>%
  filter(rho_X==0.1) %>%
  group_by(type, n, per_year_incidence_C, rho) %>%
  summarise(
    MSE = mean(bias^2, na.rm=TRUE)
  ) %>%
  pivot_wider(names_from="type", values_from="MSE") %>%
  dplyr::select(c(1,2,3,5,4)) %>%
  kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
  kable_styling(
    latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
    full_width = FALSE # Prevent the table from stretching to full page width
  )

out %>%
  filter(rho_X==0.1) %>%
  group_by(type, n, per_year_incidence_C, rho) %>%
  summarise(
    releff = mean(releff, na.rm=TRUE)
  ) %>%
  pivot_wider(names_from="type", values_from="releff") %>%
  dplyr::select(c(1,2,3,5,4)) %>%
  kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
  kable_styling(
    latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
    full_width = FALSE # Prevent the table from stretching to full page width
  )

out %>%
  filter(per_year_incidence_C==0.08) %>%
  group_by(type, n, rho, rho_X) %>%
  summarise(
    releff = mean(releff, na.rm=TRUE)
  ) %>%
  pivot_wider(names_from="type", values_from="releff") %>%
  

# res %>%
#   group_by(type, n, rho, RR_C) %>%
#   summarise(
#     MSE = mean(bias^2)
#   ) %>%
#   pivot_wider(names_from=c("type"), values_from="MSE") %>%
#   dplyr::select(c(1,2,4,6,5)) %>%
#   kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
#   kable_styling(
#     latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
#     full_width = FALSE # Prevent the table from stretching to full page width
#   )
# 
# res %>%
#   group_by(type, n, rho) %>%
#   summarise(
#     cover = mean(cover)
#   ) %>%
#   pivot_wider(names_from="type", values_from="cover") %>%
#   dplyr::select(c(1,2,3,5,4)) %>%
#   kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
#   kable_styling(
#     latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
#     full_width = FALSE # Prevent the table from stretching to full page width
#   )
# 
# res %>%
#   group_by(type, rho) %>%
#   summarise(
#     med_releff = mean(releff)
#   ) %>%
#   pivot_wider(names_from="type", values_from="med_releff") %>%
#   dplyr::select(1, 2, 4, 3) %>%
#   kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
#   kable_styling(
#     latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
#     full_width = FALSE # Prevent the table from stretching to full page width
#   )
