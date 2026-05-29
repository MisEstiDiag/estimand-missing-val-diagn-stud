## ########################################################################## ##
#
# Joint Framework of estimands and missing values in diagnostic studies    #####
# Authors: Katharina Stahlmann and Alexander Fierenz                           
# Date: August 2025                                                             
# Purpose of program: parallel execution of the simulation study
#
## ########################################################################## ##

# print title to SLURM output
print("Joint Framework: Simulation")


## Settings #####
# Paths for Hummel2
path_log <- "/beegfs/u/bbb2836/jf/log"   # the respective folders must be generated in Hummel2
path_res <- "/beegfs/u/bbb2836/jf/results" 
path_home <- "/home/bbb2836/jf"
# Paths for locale
#path_log <- path_res <- path_home <- getwd()

# packages path
.libPaths("/usw/bbb2836/R-library") # only if executed on Hummel2

# list of necessary packages
pckn <- c("mvtnorm", # data generation
          "foreach", # necessary for simulation in parallel
          "doFuture", # necessary for simulation in parallel
          "future",   # necessary for simulation in parallel
          "parallelly",# necessary for simulation in parallel with doFuture
          "dplyr",   # necessary for data manipulation
          "mice", # multiple imputation
          "tidyr",
          "mitools",
          "miceafter",
          "sjlabelled",
          "mitml",
          "DescTools",
          "parallel"
)

# load all packages
invisible(lapply(pckn,library, character.only = TRUE))

# info
sessionInfo()

# load functions
source(paste0(path_home,"/Functions_jf.R"))

## Set up simulation #####
nsim <- 1000
nsim

# 768 different scenarios
# extract pm and mech scalars from batch file
pm_0   <- as.numeric(Sys.getenv("PM_0"))
mech_0 <- Sys.getenv("MECH_0")
pm_0
mech_0

grid = expand.grid(
  sim = 1:nsim 
  , N = c(100, 500) 
  , p = c(0.1, 0.3)
  , se = c( 0.8, 0.9)
  , sp = c(0.9, 0.95)
  , r_cov = c(0.4, 0.7)
  , r_ie_cov = c(0.2)
  , r_ie_ind = c(0) # r_ie_ie defined in data function
  , p_ie1 = c(0.1)
  , eff_ie1 = c(0, 0.1)
  , pr_ind = c(0) # pos_ie1 also defined in function
  , p_ie2 = c(0.1)
  , eff_ie2 = c(0, 0.1)
  , pm = pm_0
  , mech = mech_0 
)

# calculation of the true estimand

# Estimand 1
grid <- grid %>% mutate(pos_ie1 = ifelse(eff_ie1==0,0.5,0.8),
                        ie10 = p_ie1 - (eff_ie1 * p),
                        ie11 = ie10 + eff_ie1,
                        ie20 = p_ie2 - (eff_ie2 * p),
                        ie21 = ie20 + eff_ie2,
                        esti1_sens = se * (1 - ie11) + pos_ie1 * ie11,
                        esti1_spec = sp * (1 - ie10) + (1 - pos_ie1) * ie10)

# Estimand 2
grid <- grid %>% mutate(r_ie_ie = ifelse(eff_ie1!=0 & eff_ie2!=0,0.2,0),
                        p_ie1a2_d = ie11 * ie21 + r_ie_ie*sqrt(ie11*(1-ie11)*ie21*(1-ie21)),  #joint probability of IE1 and IE2
                        p_ie_d = p_ie1a2_d + (ie11 - p_ie1a2_d) + (ie21 - p_ie1a2_d),  # Probability of at least one IE = joint probability plus probability of one IE without the other
                        p_ie1a2_h = ie10 * ie20 + r_ie_ie*sqrt(ie10*(1-ie10)*ie21*(1-ie20)),
                        p_ie_h = p_ie1a2_h + (ie10 - p_ie1a2_h) + (ie20 - p_ie1a2_h),
                        esti2_sens = se * (1 - p_ie_d),
                        esti2_spec = sp * (1 - p_ie_h)) 

# Estimand 3
# grid <- grid %>% mutate(esti3_sens = se * (1 - p_ie_d),
#                         esti3_sens = ifelse(eff_ie1==0, esti3_sens + (pos_ie1 * (ie11 - p_ie1a2_d)), esti3_sens + (se * (ie11 - p_ie1a2_d))),
#                         esti3_sens = ifelse(eff_ie2!=0, esti3_sens + (ie21 - p_ie1a2_d), esti3_sens + (se * (ie21 - p_ie1a2_d))),
#                         esti3_sens = ifelse(eff_ie1==0 & eff_ie2!=0, esti3_sens + p_ie1a2_d, esti3_sens + (se * p_ie1a2_d)),
#                         )

grid <- grid %>% mutate(esti3_sens = ifelse(eff_ie1!=0 & eff_ie2==0, se, se * (1 - p_ie_d)),
                        esti3_sens = ifelse(eff_ie1!=0 & eff_ie2!=0, se * ((1 - p_ie_d)/(1 - ie11)) + (ie21 - p_ie1a2_d)/(1 - ie11), esti3_sens),
                        esti3_sens = ifelse(eff_ie1==0 & eff_ie2==0, esti3_sens + (pos_ie1*(ie11 - p_ie1a2_d)) + (se * ie21), esti3_sens),
                        esti3_sens = ifelse(eff_ie1==0 & eff_ie2!=0, esti3_sens + (pos_ie1*(ie11 - p_ie1a2_d)) + ie21, esti3_sens),
                        esti3_spec = ifelse(eff_ie1!=0 & eff_ie2==0, sp, sp * (1 - p_ie_h)),
                        esti3_spec = ifelse(eff_ie1!=0 & eff_ie2!=0, sp * ((1 - p_ie_h)/(1 - ie10)), esti3_spec),
                        esti3_spec = ifelse(eff_ie1==0 & eff_ie2==0, esti3_spec + ((1-pos_ie1)*(ie10 - p_ie1a2_h)) + (sp * ie20), esti3_spec),
                        esti3_spec = ifelse(eff_ie1==0 & eff_ie2!=0, esti3_spec + ((1-pos_ie1)*(ie10 - p_ie1a2_h)), esti3_spec),)



grid$scen_id <- rep(1:(nrow(grid)/nsim),each=nsim)

mypattern <- matrix(c(0,1,1,1,1,1,1), nrow=1, ncol=7, byrow=TRUE)
weights_mar <- matrix(c(0,0,0,0,1,1,1), nrow=1, ncol=7, byrow=TRUE)


# document progress: log file
writeLines(c("Log of iterations"), paste0(path_log, "/log.txt"))

# assign nodes - only for Hummel2
hosts <- system('srun hostname', intern = TRUE)
cat(paste("available nodes:",hosts,"\n"))

hosts_core <- rep(hosts, each=8)
cl <- parallelly::makeClusterPSOCK(hosts_core, rscript_libs = .libPaths())
plan(cluster, workers=cl)

# if on local machine
#plan(multisession) 

cat(paste("available cores:",availableCores(),"\n"))
cat(paste("number of workers:",nbrOfWorkers(),"\n"))


## Start simulation #####
suppressMessages({

  set.seed(5273)
  
  result <- foreach(i = 1:nrow(grid), .combine = rbind, .options.future=list(seed=TRUE, packages=pckn)) %dofuture% {
    
          
      sink(paste0(path_log, "/log.txt"), append=TRUE)
      on.exit(sink())
      cat(paste("Starting iteration",i,Sys.time(),Sys.info()[['nodename']],"\n"))
      
      scen_nr <- grid$scen_id[i]
      sim <- grid$sim[i]
      
      ### generate data #####
      # use only data where there is at least 1 observed value for the index test of diseased and non-diseased subjects, respectively
      data <- data.frame(matrix(ncol = 7, nrow = 1))
      colnames(data) <- c("ind", "D", "IE1", "IE2", "V4", "V5", "V6")
      data$D <- 0
      seed <- 0
      while (sum(!is.na(data[which(data$D=="diseased"), "ind"]))==0 | sum(!is.na(data[which(data$D=="healthy"), "ind"]))==0) {
        data <- data_fun(N=grid$N[i], p=grid$p[i], se=grid$se[i], sp=grid$sp[i], r_cov=grid$r_cov[i], 
                         r_ie_cov=grid$r_ie_cov[i], r_ie_ind=grid$r_ie_ind[i], r_ie_ie=grid$r_ie_ie[i],
                         p_ie1=grid$p_ie1[i], eff_ie1 = grid$eff_ie1[i], pr_ind = grid$pr_ind[i], 
                         pos_ie1=grid$pos_ie1[i], p_ie2=grid$p_ie2[i], eff_ie2 = grid$eff_ie2[i], 
                         ie10=grid$ie10[i], ie11=grid$ie11[i], ie20=grid$ie20[i], ie21=grid$ie21[i],
                         mech=grid$mech[i], pm=grid$pm[i], mypattern=mypattern, weights_mar = weights_mar, 
                         seed=(473+seed)*i)
        seed <- seed+1
      }
      
      
      # save some general info on dataset to simulation output
      data_pm <- sum(is.na(data$ind))/nrow(data) # overall proportion of missing/non-existent results
      dis_pm <- sum(is.na(data[which(data$D=="diseased"),"ind"]))/nrow(data[which(data$D=="diseased"),]) # proportion of NAs in the index test of diseased sample
      ndis_pm <- sum(is.na(data[which(data$D=="healthy"),"ind"]))/nrow(data[which(data$D=="healthy"),]) # proportion of NAs in the index test of healthy sample
      
      result <- tryCatch({
        
       # if (i %in% c(100674, 100838, 104121, 104202, 104360, 104433, 108753, 116461, 120484, 
        #             120644,  32072,  32811,  36347,  48886,  96136,  96701)) {
         # save(data, file = paste0(path_res, "/data0_", i, ".RData"))
          
      ### Estimand 1: naive approach #####
    
      ## Se/Sp and CI 
      testdta <- data %>%
        group_by(ind, D) %>%
        summarise(n = n()) %>%
        drop_na()
      # Wilson score CI
      est <- DescTools::BinomCI(sum(testdta[which(testdta$ind=="positive" & testdta$D=="diseased"),"n"]), sum(testdta[which(testdta$D=="diseased"), "n"]),
              conf.level=0.95, sides="two.sided", method="wilson") 
      sens.e1 <- est[1, "est"]
      sens.cil.e1 <- est[1,"lwr.ci"]
      sens.ciu.e1 <- est[1,"upr.ci"]
      est <- DescTools::BinomCI(sum(testdta[which(testdta$ind=="negative" & testdta$D=="healthy"),"n"]), sum(testdta[which(testdta$D=="healthy"), "n"]),
                                conf.level=0.95, sides="two.sided", method="wilson")
      spec.e1 <- est[1,"est"]
      spec.cil.e1 <- est[1,"lwr.ci"]
      spec.ciu.e1 <- est[1,"upr.ci"]
      
      
      
      
      
      ### Estimand 2: worst case scenario #####
      ## set index test values if IE1 or IE2 occured and missing values to false test decision (false negative or false positive)
      data <- data %>%
        mutate(
          ind_e2 = factor(case_when((IE1=="IE present" | IE2=="IE present" | is.na(ind)) & D=="healthy" ~ "positive",
                             (IE1=="IE present" | IE2=="IE present" | is.na(ind)) & D=="diseased" ~ "negative",
                             TRUE ~ ind), 
                          levels = c("positive", "negative"))
        )
  
      ## Se/Sp and CI with epi.tests
      testdta <- data %>%
        group_by(ind_e2, D) %>%
        summarise(n = n()) %>%
        drop_na()
      
      # wilson score ci
      est <- DescTools::BinomCI(sum(testdta[which(testdta$ind_e2=="positive" & testdta$D=="diseased"),"n"]), sum(testdta[which(testdta$D=="diseased"), "n"]),
                                conf.level=0.95, sides="two.sided", method="wilson")
      sens.e2 <- est[1,"est"]
      sens.cil.e2 <- est[1,"lwr.ci"]
      sens.ciu.e2 <- est[1,"upr.ci"]
      est <- DescTools::BinomCI(sum(testdta[which(testdta$ind_e2=="negative" & testdta$D=="healthy"),"n"]), sum(testdta[which(testdta$D=="healthy"), "n"]),
                                conf.level=0.95, sides="two.sided", method="wilson")
      spec.e2 <- est[1,"est"]
      spec.cil.e2 <- est[1,"lwr.ci"]
      spec.ciu.e2 <- est[1,"upr.ci"]
      
      
      
      
      ### Estimand 3: targeted approach #####
      
      
      #### check IE1 #####
      if (grid$eff_ie1[i] != 0) { # IE1 informative > PS 
        
        # PS 1 - IE present
        # we do not perform any analyses for this group as it is too small
        
        
        # PS 2- IE not present
        data <- data %>%
          filter(IE1=="IE not present") %>%
          select(-ind_e2) # keep only original index test
        
        predMatrix <- make.predictorMatrix(data) # in some cases predMatrix should be manually altered
        predMatrix[c("IE1", "IE2"),] <- 0
        predMatrix[,c("IE1", "IE2")] <- 0
        
        
      } else {
        
        data <- data %>%
          select(-ind_e2) # keep only original index test
        
        predMatrix <- make.predictorMatrix(data) # in some cases predMatrix should be manually altered
        predMatrix[c("IE2"),] <- 0
        predMatrix[,c( "IE2")] <- 0
        
      }
      
      
      #### mice for missing values and non-existent results for HYP strategy #####
      set.seed(56*i)
      
      # MI
      impMethod <- make.method(data) 
      impMethod[] <- ""
      impMethod["ind"] <- "logreg" 
      invisible(capture.output(imp <- mice(data, m=30, maxit=10, method = impMethod, predictorMatrix = predMatrix)))
      implist <- mitml::mids2mitml.list(imp) 
      
      # document logged events
      if (length(unique(as.character(imp$loggedEvents[, "out"])))==0) {
        mice_out <- "none"
        mice_meth <- "none"
      } else {
        mice_out <- unique(as.character(imp$loggedEvents[, "out"])) # which variables were removed from the imputation model 
        mice_meth <- unique(as.character(imp$loggedEvents[, "meth"])) # which problem arose 
      }
     
      
      #### check IE2 (non-existent) #####
      if (grid$eff_ie2[i] != 0) { # IE2 informative > IE+
        
        implist <- lapply(implist, function(data)
          data <- data %>%
            mutate(
              ind = case_when(IE2=="IE present" ~ "positive",
                                 TRUE ~ ind)
            ) 
        )

      } 
      
      
      # transform imputed data for later calculation of sens and spec
      datlist_sens <- lapply(implist, transform_dat_cp, a="sens")
      datlist_spec <- lapply(implist, transform_dat_cp, a="spec")
      
      ## with miceafter package wilson CI; Sensitivity
      # check whether all diseased subject have NAs --> sensitivity cannot be calculated 
      #save(data, file = paste0(path_res, "/data1_", i, ".RData"))
      
      err_test <- data %>% filter(D=="diseased") %>% select(ind)
      if(sum(!is.na(err_test))==0){
        sens.e3 <- NA
        sens.cil.e3 <- NA
        sens.ciu.e3 <- NA
      } else if (sum(is.na(err_test))==0 | 
          (sum(is.na(err_test))==sum(is.na(data[which(data$IE2=="IE present" & data$D=="diseased"), "ind"])) & grid$eff_ie2[i]!=0)) { # index test in the diseased has no missing values, 
        # estimator and CI will be estimated based on the non-imputed original data
        testdta <- data %>%
          mutate( # indicative event strategy for IE2 present
            ind = factor(case_when(IE2=="IE present" ~ "positive",
                                   TRUE ~ ind), 
                         levels = c("positive", "negative"))
          ) %>%
          group_by(ind, D) %>%
          summarise(n = n()) %>%
          drop_na()
        est <- DescTools::BinomCI(sum(testdta[which(testdta$ind=="positive" & testdta$D=="diseased"),"n"]), sum(testdta[which(testdta$D=="diseased"), "n"]),
                                  conf.level=0.95, sides="two.sided", method="wilson") 
        sens.e3 <- est[1, "est"]
        sens.cil.e3 <- est[1,"lwr.ci"]
        sens.ciu.e3 <- est[1,"upr.ci"]
        
      } else if (sum(!is.na(err_test))!=0) {
        imp_dat <- list2milist(datlist_sens)
        ra <- with(imp_dat, expr=prop_wald(ind_cp ~ 1))
        res <- pool_prop_wilson(ra)  #--> reference to Lott & Reiter 2018 (simulation study)!
        sens.e3 <- res[1]
        sens.cil.e3 <- res[2]
        sens.ciu.e3 <- res[3]
      }
      
      ## Specificity
      if (sum(!is.na(data[which(data$D=="healthy"),"ind"]))==0) {
        spec.e3 <- NA
        spec.cil.e3 <- NA
        spec.ciu.e3 <- NA
      } else if (sum(is.na(data[which(data$D=="healthy"),"ind"]))==0 | 
          (sum(is.na(data[which(data$D=="healthy"),"ind"]))==sum(is.na(data[which(data$IE2=="IE present" & data$D=="healthy"), "ind"])) & grid$eff_ie2[i]!=0)) {
        
        testdta <- data %>%
          mutate( # indicative event strategy for IE2 present
            ind = factor(case_when(IE2=="IE present" ~ "positive",
                                   TRUE ~ ind), 
                         levels = c("positive", "negative"))
          ) %>%
          group_by(ind, D) %>%
          summarise(n = n()) %>%
          drop_na()
        est <- DescTools::BinomCI(sum(testdta[which(testdta$ind=="negative" & testdta$D=="healthy"),"n"]), sum(testdta[which(testdta$D=="healthy"), "n"]),
                                  conf.level=0.95, sides="two.sided", method="wilson")
        spec.e3 <- est[1,"est"]
        spec.cil.e3 <- est[1,"lwr.ci"]
        spec.ciu.e3 <- est[1,"upr.ci"]
      } else if (sum(!is.na(data[which(data$D=="healthy"),"ind"]))!=0) {
        imp_dat <- list2milist(datlist_spec)
        ra <- with(imp_dat, expr=prop_wald(ind_cp ~ 1))
        res2 <- pool_prop_wilson(ra)  #--> reference to Lott & Reiter 2018 (simulation study)!
        spec.e3 <- res2[1]
        spec.cil.e3 <- res2[2]
        spec.ciu.e3 <- res2[3]
      }
      
        #} else {
         # sens.e1 <- sens.cil.e1 <- sens.ciu.e1 <- spec.e1 <- spec.cil.e1 <- spec.ciu.e1 <- NA # Estimand 1
        #  sens.e2 <- sens.cil.e2 <- sens.ciu.e2 <- spec.e2 <- spec.cil.e2 <- spec.ciu.e2 <- NA # Estimand 2
         # sens.e3 <- sens.cil.e3 <- sens.ciu.e3 <- spec.e3 <- spec.cil.e3 <- spec.ciu.e3 <- NA # Estimand 3
        #  mice_meth <- mice_out <- "no info"
        #}
      
      # write results to output
      res <- c(i, scen_nr, sim, data_pm, dis_pm, ndis_pm, 
               sens.e1, sens.cil.e1, sens.ciu.e1, spec.e1, spec.cil.e1, spec.ciu.e1, # Estimand 1
               sens.e2, sens.cil.e2, sens.ciu.e2, spec.e2, spec.cil.e2, spec.ciu.e2, # Estimand 2
               sens.e3, sens.cil.e3, sens.ciu.e3, spec.e3, spec.cil.e3, spec.ciu.e3, # Estimand 3
               mice_meth, mice_out) # info on problems with mice
      return(res)
      
    }, error = function(err){
      
      print(err)
      save(data, file = paste0(path_res, "/data_err_", i, "_part1.RData"))
      sens.e1 <- sens.cil.e1 <- sens.ciu.e1 <- spec.e1 <- spec.cil.e1 <- spec.ciu.e1 <- NA # Estimand 1
      sens.e2 <- sens.cil.e2 <- sens.ciu.e2 <- spec.e2 <- spec.cil.e2 <- spec.ciu.e2 <- NA # Estimand 2
      sens.e3 <- sens.cil.e3 <- sens.ciu.e3 <- spec.e3 <- spec.cil.e3 <- spec.ciu.e3 <- NA # Estimand 3
      mice_meth <- mice_out <- "no info"
      
      res <- c(i, scen_nr, sim, data_pm, dis_pm, ndis_pm, 
               sens.e1, sens.cil.e1, sens.ciu.e1, spec.e1, spec.cil.e1, spec.ciu.e1, # Estimand 1
               sens.e2, sens.cil.e2, sens.ciu.e2, spec.e2, spec.cil.e2, spec.ciu.e2, # Estimand 2
               sens.e3, sens.cil.e3, sens.ciu.e3, spec.e3, spec.cil.e3, spec.ciu.e3, # Estimand 3
               mice_meth, mice_out) # info on problems with mice
      return(res)
      
      
    })
    
    return(result)
  }
  
})

result <- as.data.frame(result)

colnames(result) <- c("iteration", "scen_id", "sim", "data_pm", "dis_pm", "ndis_pm", 
                      "sens.e1", "sens.cil.e1", "sens.ciu.e1", "spec.e1", "spec.cil.e1", "spec.ciu.e1", 
                      "sens.e2", "sens.cil.e2", "sens.ciu.e2", "spec.e2", "spec.cil.e2", "spec.ciu.e2",
                      "sens.e3", "sens.cil.e3", "sens.ciu.e3", "spec.e3", "spec.cil.e3", "spec.ciu.e3",
                      "mice_meth", "mice_out")

result <- result %>% 
   mutate(
     across(!starts_with("mice"), as.numeric)
   )
  
# combine results with grid 
res <- merge(grid, result, by = c("scen_id", "sim"))


## save results #####
save(res, file = paste0(path_res, "/res_", nsim, "_", pm_0, "_", mech_0, ".RData"))


# stop cluster
parallel::stopCluster(cl)
