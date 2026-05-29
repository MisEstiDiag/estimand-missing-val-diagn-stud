## ########################################################################## ##
#
# Joint Framework of estimands and missing values in diagnostic studies    #####
# Authors: Katharina Stahlmann and Alexander Fierenz                           
# Date: April 2025                                                             
# Purpose of program: Functions
#
## ########################################################################## ##

##                         Data generation                                 #####

# Description: function to generate simulated data
# Simulation parameters:
# N = sample size
# p = prevalence of target condition
# se = sensitivity
# sp = specificity
# r_cov = correlation between index test and covariates & among covariates 
# r_ie_cov = correlation between Interfering events and covariates
# r_ie_ind = correlation between IEs and index test
# r_ie_ie = correlation between interfering events
# mech = missingness mechanism (MCAR, MAR, MNAR)
# pm = proportion of missing values
# mypattern = pattern of missing values (only in the index test)
# weights_mar = pattern of weights for MAR
# p_ie1 = proportion of IE 1 
# eff_ie1 = difference in proportion of occurring IE1 between diseased and healhy subjects (IEdiseased - IEhealthy)
# p_ie2 = proportion of IE 2
# eff_ie2 = difference in proportion of occurring IE2 between diseased and healhy subjects
# mu_ie1 = mean for IE 1 for diseased group
# mu_ie2 = mean for IE 2 for diseased group
# pr_ind = probability that the index test is not affected if IE1 occured
# pos_ie1 = probability that the index test is positive if IE1 occured and affected the index test
# ie10 = proportion of healthy subjects with IE1 
# ie11 = proportion of disesed subjects with IE1
# ie20 = proportion of healthy subjects with IE2
# ie21 = proportion of diseased subjects with IE2
# seed = seed for data generation

data_fun <- function(N, p, se, sp, r_cov, r_ie_ind, r_ie_cov, r_ie_ie, p_ie1, eff_ie1, pr_ind, pos_ie1, 
                     p_ie2, eff_ie2, ie10, ie11, ie20, ie21, mech, pm, mypattern, weights_mar, seed){
 
  
  sigma <- matrix(c(1,r_ie_ind,r_ie_ind,r_cov,r_cov,r_cov, 
                    r_ie_ind,1,r_ie_ie,r_ie_cov,r_ie_cov,r_ie_cov,
                    r_ie_ind,r_ie_ie,1,r_ie_cov,r_ie_cov,r_ie_cov,
                    r_cov, r_ie_cov,r_ie_cov,1,r_cov,r_cov,
                    r_cov,r_ie_cov,r_ie_cov,r_cov,1,r_cov,
                    r_cov,r_ie_cov,r_ie_cov,r_cov,r_cov,1), ncol=6,nrow=6, byrow = TRUE)
  
  n <- rbinom(N, 1, prob=p)
  n1 <- sum(n) # with target condition / diseased
  n0 <- length(n)-n1 # without target condition
  mu_1=qnorm(sp)-qnorm(1-se) # calculate mean of index test for diseased sample
  
  mu_ie1=qnorm(1-ie10)-qnorm(1-ie11) # qnorm(proportion of no IEs in healthy) - qnorm(proportion of no IEs in diseased) --> analog to sp/se
  #mu_ie1=qnorm(ie11)-qnorm(ie10) # alternative?
  mu_ie2=qnorm(1-ie20)-qnorm(1-ie21) 

  # generate data
  set.seed(seed)
  N0 <-as.data.frame(rmvnorm(n0, mean = c(0, 0, 0, 0, 5, 70), sigma = sigma,
                             method=c("eigen"), pre0.9_9994 = FALSE, checkSymmetry = TRUE)) # non-diseased population
  N1 <- as.data.frame(rmvnorm(n1, mean = c(mu_1, mu_ie1, mu_ie2, 0, 5, 70), sigma = sigma, 
                              method=c("eigen"), pre0.9_9994 = FALSE, checkSymmetry = TRUE)) # diseased population; V1 is index test
  N0$D <- 0 # reference test non-diseased
  N1$D <- 1 # reference test diseased
  
  data <- N0 %>%
    bind_rows(N1) %>%
    mutate(
      V4 = if_else(V4<=median(V4), 1, 2), # e.g. sex
      ind = if_else(V1<qnorm(sp), 0, 1), # 0=neg, 1=pos
      IE1 = if_else(V2<qnorm(1-ie10), 0, 1), # 1 = present, 0 = not present
      IE2 = if_else(V3<qnorm(1-ie20), 0, 1) # 1 = present, 0 = not present 
    ) %>%
   select(ind, D, IE1, IE2, V4, V5, V6) # ind bin, reftest, IE1 bin, IE2 bin, cov1, cov2, cov3
  
  
  # insert missing values with ampute function (mice package)
  # data must be numeric
  if (mech=="MAR") {
    amp <- ampute(data, prop = pm, patterns = mypattern, mech = mech, weights = weights_mar)
    data_mis <- amp$amp
  } else {
    amp <- ampute(data, prop = pm, patterns = mypattern, mech = mech)
    data_mis <- amp$amp
  }
  
  # change ind value if IE1 occured, set ind value for non-existent IE to NA & label
  set.seed(seed+34)
  data_mis <- data_mis %>%
    rowwise() %>%
    mutate(
      ind = if_else(IE1 == 1 & !is.na(ind), sample(c(ind,rbinom(1,1,pos_ie1)),1, prob=c(pr_ind,1-pr_ind)), ind) #if IE1 occurred chose with a 1-pr_ind% chance if the value is affected and with a pos_IE% chance that the affected value is positive (=1)
    ) %>%
    ungroup %>%
    mutate(
      ind = if_else(IE2 == 1, NA, ind)
    ) %>%
    mutate(
      V4 = as.factor(V4), # e.g. sex
      D = factor(D, levels = c(1,0), labels = c("diseased", "healthy")),
      ind = factor(ind, levels = c(1,0),
                   labels = c("positive", "negative")),
      IE1 = factor(IE1, levels = c(1,0),
                   labels=c("IE present", "IE not present")),
      IE2 = factor(IE2, levels = c(1,0),
                   labels=c("IE present", "IE not present"))
    ) %>%
    sjlabelled::var_labels(
      D = "Reference test",
      ind = "Index test",
      V4 = "Covariate z1 (binary)",
      IE1 = "Interfering event (existent)",
      IE2 = "Interfering event (non-existent)",
      V5 = "Covariate z2 (continuous)",
      V6 = "Covariate z3 (continuous)"
    )
  
  return(data_mis)
  
}





##                       Transform imputed data                            #####

# Description: transform imputed data depending for calculating sensitivity and specificity, respectively
# Function parameters:
# data = individual imputed dataset
# a = string indicating "sens" or "spec" depending on what should be calculating in the following

transform_dat_cp <- function(data, a) {
  
  if (a=="sens"){
    data <- data %>%
      mutate(
        ind_cp = case_when(ind=="positive" ~ 1,
                           ind=="negative" ~ 0)
      ) %>%
      filter(
        D == "diseased"
      )
    
  } else if (a=="spec") {
    
    data <- data %>%
      mutate(
        ind_cp = case_when(ind=="positive" ~ 0,
                           ind=="negative" ~ 1)
      ) %>%
      filter(
        D == "healthy"
      )
    
  }
  
  return(data)
  
}



##                   Estimate performance parameter                        #####

# Description: Function for estimating performance parameter using the simsum function (rsimsum package)
# dynamically for sensitivity and specificity and for input se/sp and esti_sens/esti_spec as true values
# Function parameters:
# data: data with simulation results in long format (methods compared in one variable)
# varname: varname of the variable providing the estimated value (theta)
# var_label: label of the estimated variable 
# true_var: varname of the column providing the "true" value needed for estimating the performance parameter

perform_param <- function(data, varname, var_label, true_var) {
  
  results <- simsum(data, estvarname = varname, true = true_var, methodvar = "Estimand",
                    ref = "Estimand 1", by = "scenario", 
                    ci.limits = c(paste0(varname, ".cil"), paste0(varname, ".ciu")),
                    df = NULL, dropbig = FALSE, x = FALSE, control = list())
  
  # !!! coverage should also be estimated using the ci.limits, why does it not work? !!!

  ### write results into dataframe  
  dta <- results[["summ"]] %>%
    filter(stat %in% c("thetamean", "bias", "rbias", "mse", "nsim")) %>%
    pivot_wider(names_from = stat, values_from = c(est, mcse)) %>%
    dplyr::select(-mcse_nsim, -mcse_thetamean) %>%
    sjlabelled::var_labels(
      est_nsim = !!sym(paste0("Number of estimated ", var_label, " estimates")),
      est_thetamean = !!sym(paste0("Mean of ", var_label, " estimates")),
      est_bias = !!sym(paste0("Absolute bias ", var_label)),
      est_rbias = !!sym(paste0("Relative bias ", var_label)),
      est_mse = !!sym(paste0("Mean squared error ", var_label)),
      mcse_bias = !!sym(paste0("MCSE for absolute bias ", var_label)),
      mcse_rbias = !!sym(paste0("MCSE for relative bias ", var_label)),
      mcse_mse = !!sym(paste0("MCSE for mean squared error ", var_label))
    ) %>%
    mutate(
      scenario = as.numeric(scenario),
      Estimand = as.character(Estimand)
    )
  
  # if performance parameters are estimated based on the input se/sp values, rename variable names
  if (true_var=="se" | true_var=="sp") {
    dta <- dta %>%
      dplyr::select(-est_nsim, -est_thetamean) %>%
      rename_with(~ paste0(varname, "input_", .x), starts_with("est_")) %>%
      rename_with(~ paste0(varname, "input_", .x), starts_with("mcse")) 
  } else {
    dta <- dta %>%
      rename_with(~ paste0(varname, "_", .x), starts_with("est_")) %>%
      rename_with(~ paste0(varname, "_", .x), starts_with("mcse")) 
  }
  
  return(dta)
  
}



##                       Metamodels 1st approach                           #####

# Description: Function for metamodels (stratified by Estimand and missingness mechanism)

mm1 <- function(data, performance, mechanism, estimand) {
  
  data <- data %>%
    filter(mech == mechanism & Estimand==estimand) %>%
    mutate(
      sens_est_rbias = abs(sens_est_rbias),
      spec_est_rbias = abs(spec_est_rbias),
      across(all_of(sim_param), ~as.factor(.x))
    )
  # mixed linear model
  form <- paste0(performance, "*100", " ~ (", paste(sim_param, collapse=" + "), ")^2")
  mm <- lm(as.formula(form), data=data)
  invisible(capture.output(backward <- stats::step(mm, direction='backward')))
  r2 <- summary(backward)$adj.r.squared
  names(r2) <- "Adjusted R squared"
  
  dta <- broom::tidy(backward, conf.int=T) %>%
    dplyr::select(term, estimate, conf.low, conf.high, p.value) %>%
    mutate(
      across(contains("estimate") | contains("conf"), ~round(.x, 2)),
      p.value = round(p.value, 3),
      Mechanism = mechanism,
      Estimand = estimand,
      r2 = paste0("Adjusted R-squared: ", round(r2, 3))
    )
  
  return(dta)
}



##                       Metamodels 2nd approach                           #####

# Description: Function for metamodels (stratified by missingness mechanism)

mm2 <- function(data, performance, mechanism) {
  
  data <- data %>%
    filter(mech == mechanism) %>%
    mutate(
      Estimand = as.factor(Estimand),
      sens_est_rbias = abs(sens_est_rbias),
      spec_est_rbias = abs(spec_est_rbias),
      across(all_of(sim_param), ~as.factor(.x))
    )
  # mixed linear model
  predictors <- paste("Estimand", sim_param, sep = "*")
  form <- paste0(performance, "*100", " ~ ", paste(predictors, collapse=" + "), " + (1|scenario)")
  mm <- lme4::lmer(as.formula(form), data=data)
  
  dta <- broom.mixed::tidy(mm, effects="fixed", conf.int=T) %>%
    dplyr::select(term, estimate, conf.low, conf.high) %>%
    mutate(
      across(where(is.numeric), ~round(.x, 2)),
      Mechanism = mechanism,
      CI_sign = case_when((conf.low>0 & conf.high>0) | 
                            (conf.low<0 & conf.high<0) ~ 1,
                          TRUE ~0)
    )

  return(dta)
}
