######################################################################################################
#           Code to run a Bayesian adaptive design with interim decision rules described in 
#           Section 2.4 of Maleyeff et al. (2024). This outer shell is used with all three
#                       models (proposed FK-BMA, PLTY, and LKR) in the manuscript.
#                                   Contact: laramaleyeff@gmail.com                                  
#                                       Last updated: April 2024                                       
######################################################################################################

# First, load the code for each model fitting procedure
source("mcmc_liu.R")
source("mcmc_park.R")
source("mcmc_maleyeff.R")
#
# function interimAnalysis
#
# Author: Lara Maleyeff
#
# Function:       interimAnalysis
# Description:    Perform one interim analysis: fit the respective model based on currently
#                 available data and assess whether the trial should be stopped for efficacy or futility
# Parameters:     trialfunction     function to fit the appropriate model
#                 curdata           A data frame with the observations arranged by row, and including the columns:
#                                   - trt: the group indicator (=1 for experimental group; =0 for control group)
#                                   - Y: outcome
#                                   - One column for each of candsplinevars, named appropriately (i.e. if
#                                     candsplinevars = c("X_1","X_2) then there are two respective columns 
#                                     named "X_1" and "X_2)
#                                   - One column for each of candbinaryvars
#                 candpslinevars    Vector with names of continuous variables
#                 candbinaryvars    Vector with names of binary variables
#                 last              Boolean indicating whether this is the last interim analysis (if so, an 
#                                   assessment of futility is unnecessary)
#                 family            Family and link function of outcome, defaults to continuous
#                 B                 The number of posterior samples
#                 burnin            The number of burn-in samples
#                 thin              The thinning parameter
#                 alpha             Cutoff for effective subspace
#                 B_1               Efficacy cutoff (significance)
#                 B_2               Futility cutoff (significance)
#                 b_1               Efficacy cutoff (magnitude)
#                 b_2               Futility cutoff (magnitude)
#                 e_1               Cutoff for effective subspace magnitude
#                 pi                Cutoff for prevalence of the effective subspace
#
# Returns:        If the fitting procedure is successful, the function returns a list with:
#                 - success: Indicates whether the procedure was successful based on geweke convergence
#                 - included_vars: Selected tailoring variables based on cutoff pi 
#                 - stop_efficacy: Boolean indicating whether the trial is to be stopped for efficacy
#                 - stop_futility: Boolean indicating whether the trial is to be stopped for futility
#                 - mse: Mean squared error of each individual's treatment effect estimate
#                 - prop_eff: The prevalence of the effective subspace
#                 - trial_results: Results from the given interim analysis, contains a list with method-specific
#                   entries. All include the posterior distribution of treatment effect for each individual
#                 - mean_subgroup_ate: The average treatment effect in the effective subspace (used for debugging)
#                 - mean_subgroup_ate_gr_cutoff: Proportion of effective subspace who met the efficacy criteria (used for debugging)
#                 - mean_subgroup_ate_l_cutoff: Proportion of effective subspace who met the futility criteria  (used for debugging)
#                 If the procedure is not successful, the function returns list(success = FALSE)
interimAnalysis <- function(trialfunction,
                            curdata,
                            candsplinevars,
                            candbinaryvars,
                            last,
                            family,
                            B,
                            burnin,
                            thin,
                            alpha,
                            B_1,
                            B_2,
                            b_1,
                            b_2,
                            e_1,
                            pi,
                            ...) {
  
  for (i in 1:5) {
    trial_results = trialfunction(curdata,
                                  candsplinevars,
                                  candbinaryvars,
                                  family,
                                  B,
                                  burnin,
                                  thin,
                                  pi,
                                  ...)
    if (trial_results$success) {
      break
    }
  }
  
  if (trial_results$success) {
    trt_eff_posterior = trial_results$trt_eff_posterior
    quantiles = rowQuantiles(trt_eff_posterior,probs = c(alpha))
    curdata_subgroup = curdata[which(quantiles>e_1),]

    trt_eff_posterior_subgroup = trt_eff_posterior[which(quantiles>0),]
    
    trt_eff_per_person = rowMeans(trt_eff_posterior)
    mse = mean((trt_eff_per_person  - curdata$truth)^2)
    print(paste("mse: ", mse))
    print(paste("trial success: ", trial_results$geweke.conv))
   
    prop_eff = nrow(curdata_subgroup)/nrow(curdata)
    print(paste("prop_eff:", prop_eff))
    print(paste("included vars: ", paste(trial_results$included_vars,collapse=",")))

    mean_subgroup_ate = 0
    mean_subgroup_ate_gr_cutoff = 0
    mean_subgroup_ate_l_cutoff = 0
    stop_efficacy = FALSE
    stop_futility = FALSE

    if (prop_eff > 0.1) {
      # Compute the average treatment effect for the sensitive subgroup
      combined_subgroup_ate = colMeans(trt_eff_posterior_subgroup)
      mean_subgroup_ate = mean(combined_subgroup_ate)
      mean_subgroup_ate_gr_cutoff = mean(combined_subgroup_ate > b_1)
      mean_subgroup_ate_l_cutoff = mean(combined_subgroup_ate < b_2)

      if (mean(combined_subgroup_ate > b_1) > B_1) {
        stop_efficacy = TRUE
      } 
      
      if (((mean(combined_subgroup_ate < b_2) > B_2)) & !last & !stop_efficacy) {
        stop_futility = TRUE
      }
    } else {
      stop_futility = TRUE
    }
    
    print(paste("stop_futility",stop_futility))
    print(paste("stop_efficacy",stop_efficacy))
    
    return(list(success = TRUE,
                included_vars = trial_results$included_vars,
                stop_efficacy = stop_efficacy, 
                stop_futility = stop_futility,
                mse = mse,
                prop_eff = prop_eff,
                trial_results = trial_results,
                mean_subgroup_ate = mean_subgroup_ate,
                mean_subgroup_ate_gr_cutoff = mean_subgroup_ate_gr_cutoff,
                mean_subgroup_ate_l_cutoff = mean_subgroup_ate_l_cutoff
               )
    )
  } else {
    return(list(success=FALSE))
  }
 

}

#
# function runTrial
#
# Author: Lara Maleyeff
#
# Function:       runTrial
# Description:    Run complete Bayesian adaptive trial
# Parameters:     data_pool         A population data frame with the observations arranged by row, and including the columns:
#                                   - trt: the group indicator (=1 for experimental group; =0 for control group)
#                                   - Y: outcome
#                                   - One column for each of candsplinevars, named appropriately (i.e. if
#                                     candsplinevars = c("X_1","X_2) then there are two respective columns 
#                                     named "X_1" and "X_2)
#                                   - One column for each of candbinaryvars
#                 data_test         A large, external testing data frame with the observations arranged by row, and including the columns:
#                                   - One column for each of candsplinevars, named appropriately (i.e. if
#                                     candsplinevars = c("X_1","X_2) then there are two respective columns 
#                                     named "X_1" and "X_2)
#                                   - One column for each of candbinaryvars
#                                   - true_trt: = 1 if treatment is truly effective and = 0 if not
#                 candpslinevars    Vector with names of continuous variables
#                 candbinaryvars    Vector with names of binary variables
#                 trialfunction     function to fit the appropriate model
#                 quantilefunction  function to compute the alpha-row quantile for each individual
#                 family            Family and link function of outcome, defaults to continuous
#                 true_tailoring vars   Names of the true tailoring variables
#                 B                 The number of posterior samples
#                 burnin            The number of burn-in samples
#                 thin              The thinning parameter
#                 alpha             Cutoff for effective subspace
#                 B_1               Efficacy cutoff (significance)
#                 B_2               Futility cutoff (significance)
#                 b_1               Efficacy cutoff (magnitude); defaults to 0
#                 b_2               Futility cutoff (magnitude); defaults to 0
#                 e_1               Cutoff for effective subspace magnitude; defaults to 0
#                 pi                Cutoff for prevalence of the effective subspace (defaults to 0.1)
#                 enrich            Boolean indicating whether we should perform adaptive enrichment (defaults to TRUE)
#
# Returns:        If the fitting procedure is successful, the function returns a data.frame with:
#                 - success: Indicates whether the procedure was successful based on geweke convergence
#                 - B: The number of posterior samples
#                 - burnin: The number of burn-in samples
#                 - thin: The thinning parameter
#                 - final_stop_efficacy: Indicates whether the efficacy criteria was met in final interim analysis
#                 - final_stop_futility: Indicates whether the futility criteria was met in final interim analysis
#                 - prop_eff_final: Prevalence of effective subspace in the final analysis
#                 - final_trial_size: Sample size at end of trial
#                 - effect_and_subgroup: Boolean indicating whether efficacy criteria was met for only a subset of the entire sample
#                 - effect_overall: Boolean indicating whether efficacy criteria was met for the entire sample
#                 - subgroup_spec: Boolean indicating whether a subgroup-specific analysis was performed
#                 - subgroup: String of selected tailoring variables, separated by a comma
#                 - correct_subgroup: Boolean indicating whether ONLY the correct tailoring variables were selected
#                 - include_correct_subgroup:  Boolean indicating whether the correct tailoring variables were selected
#                 - accuracy: % of correct treatment decisions based on data_test
#                 - mse: Mean squared error of each individual's treatment effect estimate
#                 - mean_subgroup_ate: The average treatment effect in the effective subspace (used for debugging)
#                 - mean_subgroup_ate_gr_cutoff: Proportion of effective subspace who met the efficacy criteria (used for debugging)
#                 - mean_subgroup_ate_l_cutoff: Proportion of effective subspace who met the futility criteria  (used for debugging)
#                 - for each interim analysis i: 
#                   - stop_efficacy_i: Boolean indicating whether the trial stopped for efficacy at interim analysis i
#                   - stop_futility_i: Boolean indicating whether the trial stopped for futility at interim analysis i
#                   - prop_eff_i: Prevalence of the effective subspace in interim analysis i
#                   - true_prop_eff_i: True proportion of individuals expected to benefit from treatment at interim i
#                 If the procedure is not successful, the function returns data.frame(success = FALSE)
runTrial <- function(data_pool,
                     data_test,
                     candsplinevars,
                     candbinaryvars,
                     trialfunction,
                     quantilefunction,
                     family = gaussian(),
                     true_tailoring_vars,
                     B,
                     burnin,
                     thin,
                     interim_n,
                     alpha,
                     B_1,
                     B_2,
                     b_1 = 0,
                     b_2 = 0,
                     e_1 = 0,
                     pi = 0.1,
                     enrich = TRUE,
                     ...
) {
  library(matrixStats)
  library(coda)
  library(splines)
  library(dplyr)
  
  subgroup_spec = FALSE
  for (interim in 1:length(interim_n)) {
    max_interim = interim
    print(paste("interim: ", interim))
    # If we are adaptively enriching trial population, after the first interim analysis
    # we restrict enrollment to individuals expected to benefit from treatment, i.e.
    # the effective subspace criteron from previous interim analysis
    if (interim > 1 & enrich) {
      quantile_interim <- quantilefunction(data_pool,
                                        candsplinevars,
                                        candbinaryvars,
                                        interim_analysis_results$trial_results,
                                        alpha)
      data_indx = sample(which(quantile_interim>e_1),interim_n[interim])
      curdata = rbind(curdata, data_pool[data_indx,])
    } 
    
    if (interim == 1 | !enrich) {
      data_indx = sample(nrow(data_pool), interim_n[interim])
      curdata = data_pool[data_indx,]
    } 
    
    assign(paste0("true_prop_eff_",interim),mean(curdata$truth>0))
    
    data_pool = data_pool[-data_indx,]
    
    interim_analysis_results = interimAnalysis(trialfunction,
                                               curdata,
                                               candsplinevars,
                                               candbinaryvars,
                                               last = (interim == length(interim_n)),
                                               family,
                                               B,
                                               burnin,
                                               thin,
                                               alpha,
                                               B_1,
                                               B_2,
                                               b_1,
                                               b_2,
                                               e_1,
                                               pi,
                                               ...)
    
    if (!interim_analysis_results$success) {
      return(data.frame(success = FALSE))
    }
    
    assign(paste0("prop_eff_",interim),interim_analysis_results$prop_eff)
    assign(paste0("stop_efficacy_",interim),interim_analysis_results$stop_efficacy)
    assign(paste0("stop_futility_",interim),interim_analysis_results$stop_futility)

    # If the proportion effective is ever less than one,
    # we are performing a subgroup-specific analysis.
    
    # For example, if prop_eff = 0.5 in interim analysis 1, 
    # then we restrict the entry criteria and prop_eff may be 
    # equal to 1 in the second cohort. An overall analysis is
    # performed if prop_eff=1 for all interim analyses
    if (interim_analysis_results$prop_eff < 1) {
      subgroup_spec = TRUE
    }
    
    # If we have reached the last interim analysis, or the trial is being 
    # stopped early for efficacy or futility: compute external accuracy and break
    # out of the outer for loop
    if (interim == length(interim_n) | 
        interim_analysis_results$stop_efficacy | 
        interim_analysis_results$stop_futility) {
      final_trial_size = sum(interim_n[1:interim])
      # Assess accuracy based on large, external dataset
      quantile_test = quantilefunction(data_test,
                                    candsplinevars,
                                    candbinaryvars,
                                    interim_analysis_results$trial_results,
                                    alpha
                                    )
      # Here, "trt" refers to the treatment decisions that would be made based on the observed
      # trial data
      data_test$trt = as.numeric(quantile_test > e_1)
      accuracy = sum(as.numeric(data_test$trt == data_test$true_trt))/nrow(data_test)
      break
    } 
  }
  
  # correct_subgroup tells us if the selecting tailoring variables are identical to the true
  # tailoring variables
  correct_subgroup = identical(true_tailoring_vars,interim_analysis_results$included_vars)
  # identical() doesn't handle length 0 vectors, so we adjust for this
  if (length(true_tailoring_vars) == 0) {
    correct_subgroup = (length(interim_analysis_results$included_vars)==0)
  }
  
  # check if the selected tailoring variables contain the true tailoring variables; i.e.
  # = TRUE even if we have selected extra 
  include_correct_subgroup = all(true_tailoring_vars %in% interim_analysis_results$included_vars)
  
  returned = data.frame(success = TRUE,
                        B = B,
                        burnin = burnin,
                        thin = thin,
                        final_stop_efficacy = interim_analysis_results$stop_efficacy, 
                        final_stop_futility = interim_analysis_results$stop_futility,
                        prop_eff_final = interim_analysis_results$prop_eff,
                        final_trial_size = final_trial_size, 
                        effect_and_subgroup = (interim_analysis_results$stop_efficacy & subgroup_spec),
                        effect_overall = (interim_analysis_results$stop_efficacy & !subgroup_spec),
                        subgroup_spec = subgroup_spec,
                        subgroup = paste(interim_analysis_results$included_vars,collapse=","),
                        correct_subgroup = correct_subgroup,
                        include_correct_subgroup = include_correct_subgroup,
                        accuracy = accuracy,
                        mse = interim_analysis_results$mse,
                        mean_subgroup_ate = interim_analysis_results$mean_subgroup_ate,
                        mean_subgroup_ate_gr_cutoff = interim_analysis_results$mean_subgroup_ate_gr_cutoff,
                        mean_subgroup_ate_l_cutoff = interim_analysis_results$mean_subgroup_ate_l_cutoff
  )
  
  for (i in 1:max_interim) {
    returned[[paste0("stop_efficacy_",i)]] = eval(parse(text=paste0("stop_efficacy_",i)))
    returned[[paste0("stop_futility_",i)]] = eval(parse(text=paste0("stop_futility_",i)))
    returned[[paste0("prop_eff_",i)]] = eval(parse(text=paste0("prop_eff_",i)))
    returned[[paste0("true_prop_eff_",i)]] = eval(parse(text=paste0("true_prop_eff_",i)))
  }
  
  return(returned)
}
