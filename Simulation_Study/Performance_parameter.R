## ########################################################################## ##
#
# Joint Framework of estimands and missing values in diagnostic studies    #####
# Authors: Katharina Stahlmann and Alexander Fierenz                           
# Date: October 2025                                                             
# Purpose of program: Calculation of performance parameters
#
## ########################################################################## ##

##                                Info                                     #####

# perfromance parameters needed to be estimated
## relative and absolute bias
## coverage probability
## root mean squared error
## all Monte Carlo standard errors
## all parameters estimated using estix_sens and estix_spec 
## all parameters estimated using the input se/sp


##                  Data preparation: Long format                          #####

# convert all data to long format for plots
se_sp_vars <- list(
  c("esti1_sens", "esti2_sens", "esti3_sens"),
  c("sens.e1", "sens.e2", "sens.e3"),
  c("sens.cil.e1", "sens.cil.e2", "sens.cil.e3"),
  c("sens.ciu.e1", "sens.ciu.e2", "sens.ciu.e3"),
  c("sens.ci.length.e1", "sens.ci.length.e2", "sens.ci.length.e3"),
  c("esti1_spec", "esti2_spec", "esti3_spec"),
  c("spec.e1", "spec.e2", "spec.e3"),
  c("spec.cil.e1", "spec.cil.e2", "spec.cil.e3"),
  c("spec.ciu.e1", "spec.ciu.e2", "spec.ciu.e3"),
  c("spec.ci.length.e1", "spec.ci.length.e2", "spec.ci.length.e3")
)
names <- c("esti_sens", "sens", "sens.cil", "sens.ciu", "sens.ci.length", "esti_spec", "spec", "spec.cil", "spec.ciu", "spec.ci.length")
res_long <- reshape(res_all[,!(names(res_all) %in% c("mice_meth", "mice_out"))], varying=se_sp_vars, v.names = names, times = c("Estimand 1", "Estimand 2", "Estimand 3"), timevar = "Estimand", idvar = "id_var", direction = "long")
res_long <- res_long %>%
  dplyr::select(scenario, sim, id_var, everything()) %>%
  arrange(scenario, id_var, sim)



##                  Performance Parameter dynamically                     ##### 

# using the rsimsum package for calculating mean of estimates across iterations,
# absolute bias, relative bias, mse and the MCSE

# performance parameter based on esti_sen/_spec 
dta_sens <- perform_param(res_long, varname="sens", var_label="sensitivity", true_var="esti_sens")
dta_spec <- perform_param(res_long, varname="spec", var_label="specificity", true_var="esti_spec")
# performance parameter based on input se/sp values
dta_sens2 <- perform_param(res_long, varname="sens", var_label="sensitivity (based on input se)", true_var="se")
dta_spec2 <- perform_param(res_long, varname="spec", var_label="specificity (based on input sp)", true_var="sp")




##                  Additional performance Parameter                       ##### 

# additionally calculate coverage probability, mean of ci lengths
nsim <- max(res_long$sim)
dta_other <- res_long %>%
  mutate(
    sens_cov = 1*(sens.cil<=esti_sens & esti_sens<=sens.ciu),
    spec_cov = 1*(spec.cil<=esti_spec & esti_spec<=spec.ciu),
    sensinput_cov = 1*(sens.cil<=se & se<=sens.ciu),
    specinput_cov = 1*(spec.cil<=sp & sp<=spec.ciu)
  ) %>%
  group_by(scenario, Estimand) %>%
  summarise(
    sens_ci_length_mean = mean(sens.ci.length, na.rm=TRUE),
    spec_ci_length_mean = mean(spec.ci.length, na.rm=TRUE),
    data_pm_mean = mean(data_pm, na.rm=TRUE),
    dis_pm_mean = mean(dis_pm, na.rm=TRUE),
    # the bias variables only for checking the rsimsum results
    bias_sens_check = mean(sens, na.rm=TRUE)-mean(esti_sens, na.rm=TRUE),
    relbias_sens_check = (mean(sens, na.rm=TRUE)-mean(esti_sens, na.rm=TRUE))/mean(esti_sens, na.rm=TRUE),
    # coverage
    sens_cov = mean(sens_cov, na.rm=TRUE),
    spec_cov = mean(spec_cov, na.rm=TRUE),
    sensinput_cov = mean(sensinput_cov, na.rm=TRUE),
    specinput_cov = mean(specinput_cov, na.rm=TRUE)
  ) %>%
  mutate(
    sens_mcse_cov = sqrt(sens_cov*(1-sens_cov)/nsim),
    spec_mcse_cov = sqrt(spec_cov*(1-spec_cov)/nsim),
    sensinput_mcse_cov = sqrt(sensinput_cov*(1-sensinput_cov)/nsim),
    specinput_mcse_cov = sqrt(specinput_cov*(1-specinput_cov)/nsim)
  ) %>%
  sjlabelled::var_labels(
    sens_ci_length_mean = "Mean of CI length for sensitivity ",
    spec_ci_length_mean = "Mean of CI length for specificity",
    data_pm_mean = "Mean of the proportion of missing/non-existent values in the data",
    dis_pm_mean = "Mean of the proportion of missing/non-existent values in the sample with target condition",
    sens_cov = "Coverage probability sensitivity",
    spec_cov = "Coverage probability specificity",
    sensinput_cov = "Coverage probability sensitivity (based on input se)",
    specinput_cov = "Coverage probability specificity (based on input sp)",
    sens_mcse_cov = "MCSE for coverage probability sensitivity",
    spec_mcse_cov = "MCSE for coverage probability specificity",
    sensinput_mcse_cov = "MCSE for coverage probability sensitivity (based on input se)",
    specinput_mcse_cov = "MCSE for coverage probability specificity (based on input sp)"
  )
  


##                  Merge performance parameter results                    ##### 

res_perform <- res_long %>%
  filter(sim==1) %>%
  dplyr::select(all_of(c("scenario", "Estimand", simulation_parameter))) %>%
  full_join(dta_sens, by = c("scenario", "Estimand")) %>%
  full_join(dta_sens2, by = c("scenario", "Estimand")) %>%
  full_join(dta_spec, by = c("scenario", "Estimand")) %>%
  full_join(dta_spec2, by = c("scenario", "Estimand")) %>%
  full_join(dta_other, by = c("scenario", "Estimand")) %>%
  mutate(
    sens_est_rmse = sqrt(sens_est_mse),
    spec_est_rmse = sqrt(spec_est_mse),
    sensinput_est_rmse = sqrt(sensinput_est_mse),
    specinput_est_rmse = sqrt(specinput_est_mse)
  ) %>%
  sjlabelled::var_labels(
    sens_est_rmse = "Root mean squared error sensitivity",
    spec_est_rmse = "Root mean squared error specificity",
    sensinput_est_rmse = "Root mean squared error sensitivity (based on input se)",
    specinput_est_rmse = "Root mean squared error specificity (based on input sp)"
  )

# check bias variables 
check1 <- round(res_perform$bias_sens_check, 2)==round(res_perform$sens_est_bias, 2)
table(check1)
check2 <- round(res_perform$relbias_sens_check, 2)==round(res_perform$sens_est_rbias, 2)
table(check2)
# --> bias and relative bias are correct at two decimal points


rm(dta_sens, dta_sens2, dta_spec, dta_spec2, dta_other, check1, check2, se_sp_vars, names)
