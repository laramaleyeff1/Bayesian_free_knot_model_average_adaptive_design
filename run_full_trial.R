source("mcmc_liu.R")
source("mcmc_park.R")
source("mcmc_maleyeff.R")
interimAnalysis <- function(trialfunction,
                            curdata,
                            candsplinevars,
                            candbinaryvars,
                            last,
                            family,
                            B = 10000,
                            burnin = 5000,
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
      
      print(paste("combined_subgroup_ate", mean(combined_subgroup_ate)))
      print(paste("mean(combined_subgroup_ate > b_1)", mean(combined_subgroup_ate > b_1) ))
      
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
                geweke.conv = trial_results$geweke.conv,
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

runTrial <- function(data_pool,
                     data_test,
                     candsplinevars,
                     candbinaryvars,
                     trialfunction,
                     quantilefunction,
                     family = gaussian(),
                     true_tailoring_vars,
                     alpha = 0.5,
                     B_1 = 0.975,
                     B_2 = 0.8,
                     B,
                     burnin,
                     thin,
                     pi = 0.1,
                     interim_n,
                     b_1 = 0,
                     b_2 = 0,
                     e_1 = 0,
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
    # If the proportion effective is ever less than one,
    # we are performing a subgroup-specific analysis.
    
    # For example, if prop_eff = 0.5 in interim analysis 1, 
    # then we restrict the entry criteria and prop_eff may be 
    # equal to 1 in the second cohort. An overall analysis is
    # performed if prop_eff=1 for all interim analyses
    if (!interim_analysis_results$success) {
      return(data.frame(success = FALSE))
    }
    
    assign(paste0("prop_eff_",interim),interim_analysis_results$prop_eff)
    assign(paste0("stop_efficacy_",interim),interim_analysis_results$stop_efficacy)
    assign(paste0("stop_futility_",interim),interim_analysis_results$stop_futility)

    if (interim_analysis_results$prop_eff < 1) {
      subgroup_spec = TRUE
    }
    if (interim == length(interim_n) | 
        interim_analysis_results$stop_efficacy | 
        interim_analysis_results$stop_futility) {
      final_trial_size = sum(interim_n[1:interim])
      # Assess accuracy based on external dataset
      quantile_test = quantilefunction(data_test,
                                    candsplinevars,
                                    candbinaryvars,
                                    interim_analysis_results$trial_results,
                                    alpha
                                    )
      data_test$trt = as.numeric(quantile_test > e_1)
      accuracy = sum(as.numeric(data_test$trt == data_test$true_trt))/nrow(data_test)
      break
    } 
  }
  
  correct_subgroup = identical(true_tailoring_vars,interim_analysis_results$included_vars)
  if (length(true_tailoring_vars) == 0) {
    correct_subgroup = (length(interim_analysis_results$included_vars)==0)
  }
  
  include_correct_subgroup = all(true_tailoring_vars %in% interim_analysis_results$included_vars)
  
  returned = data.frame(success = TRUE,
                        B = B,
                        burnin = burnin,
                        thin = thin,
                        B_1 = B_1,
                        B_2 = B_2,
                        b_1 = b_1,
                        b_2 = b_2,
                        e_1 = e_1,
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
    returned[[paste0("stop_futility_",i)]] = eval(parse(text=paste0("stop_futility_",i)))
    returned[[paste0("stop_efficacy_",i)]] = eval(parse(text=paste0("stop_efficacy_",i)))
    returned[[paste0("prop_eff_",i)]] = eval(parse(text=paste0("prop_eff_",i)))
    returned[[paste0("true_prop_eff_",i)]] = eval(parse(text=paste0("true_prop_eff_",i)))
  }
  
  return(returned)
}
