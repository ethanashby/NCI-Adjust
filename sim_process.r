sim<-readRDS("~/Desktop/NCI_adj/sim.rds")

library(kableExtra)
library(knitr)

ids<-as.numeric(gsub("sim_uid_", "", names(sim$results_complex)))

for(i in 1:length(sim$results_complex)){

  if(i==1){
  res<-cbind(sim$results_complex[[i]]$output, 'sim_uid'=ids[i])
  }else{
  res<-rbind(res, cbind(sim$results_complex[[i]]$output, 'sim_uid'=ids[i]))
  }
}

res<-left_join(as.data.frame(res), sim$results, by='sim_uid')

res %>%
  filter(per_year_incidence_C==0.10) %>%
  group_by(type, n, rho, rho_X) %>%
  summarise(
    MSE = mean(bias^2, na.rm=TRUE)
  ) %>%
  pivot_wider(names_from=c("type"), values_from="MSE") %>%
  arrange(n, rho_X, rho) %>%
  dplyr::select(c(1, 3, 2, 4, 9, 5, 6, 7, 8)) %>%
  mutate(
    across(where(is.numeric), ~ signif(.x, 3))
  ) %>%
  kbl(., format='latex', booktabs=TRUE) %>%
  kable_styling(
    latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
    full_width = FALSE # Prevent the table from stretching to full page width
  )

res %>%
  filter(per_year_incidence_C==0.10) %>%
  group_by(type, n, rho, rho_X) %>%
  summarise(
    releff = mean(releff, na.rm=TRUE)
  ) %>%
  pivot_wider(names_from="type", values_from="releff") %>%
  arrange(n, rho_X, rho) %>%
  dplyr::select(c(1, 3, 2, 4, 9, 5, 6, 7, 8)) %>%
  kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
  kable_styling(
    latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
    full_width = FALSE # Prevent the table from stretching to full page width
  )

res %>%
  filter(n==2000 & rho==0.75 & rho_X==0.1) %>%
  group_by(type, n, per_year_incidence_C) %>%
  summarise(
    releff = mean(releff, na.rm=TRUE)
  ) %>%
  pivot_wider(names_from="type", values_from="releff") %>%
  arrange(n, rho_X, rho) %>%
  dplyr::select(c(1, 3, 2, 4, 9, 5, 6, 7, 8)) %>%
  kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
  kable_styling(
    latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
    full_width = FALSE # Prevent the table from stretching to full page width
  )

res %>%
  filter(per_year_incidence_C==0.10) %>%
  group_by(type, n, rho) %>%
  summarise(
    coverage = mean(cover, na.rm=TRUE)
  ) %>%
  pivot_wider(names_from="type", values_from="coverage") %>%
  arrange(n, rho) %>%
  dplyr::select(c(1, 2, 3, 8, 4, 5, 6, 7)) %>%
  kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
  kable_styling(
    latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
    full_width = FALSE # Prevent the table from stretching to full page width
  )

# res %>%
#   filter(per_year_incidence_C==0.10) %>%
#   group_by(type, n, rho, rho_X) %>%
#   summarise(
#     coverage = mean(cover, na.rm=TRUE)
#   ) %>%
#   pivot_wider(names_from="type", values_from="coverage") %>%
#   arrange(n, rho_X, rho) %>%
#   dplyr::select(c(1, 3, 2, 4, 9, 5, 6, 7, 8)) %>%
#   kbl(., format='latex', digits = 3, booktabs=TRUE) %>%
#   kable_styling(
#     latex_options = c("basic", "striped", "hold_position"), # Apply classic styling options
#     full_width = FALSE # Prevent the table from stretching to full page width
#   )
  

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
