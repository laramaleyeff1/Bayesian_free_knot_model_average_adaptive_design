######################################################################################################
#                       Implementation of proposed MCMC procedure of Maleyeff et al. (2024)  
#                                   Contact: laramaleyeff@gmail.com                                  
#                                       Last updated: April 2024                                       
######################################################################################################

source("mcmc_maleyeff_helper.R")

#
# function runMHMaley
#
# Author: Lara Maleyeff
#
# Function:       runMHMaley
# Description:    Perform MCMC proceure to generate posterior distribution, using Bayesian model averaging and
#                 free-knot B-splines
# Parameters:     data              A data frame with the observations arranged by row, and including the column:
#                                   - trt: the group indicator (=1 for experimental group; =0 for control group)
#                                   - Y: continuous-valued outcome
#                                   - All the variables described in candsplinevars and candbinaryvars
#                 candsplinevars    Vector with names of continuous predictive candidate variables
#                 candbinaryvars    Vector with names of binary predictive candidate variables
#                 candinter         Which of candpslinevars and candbinaryvars are tailoring
#                 B                 The number of posterior samples
#                 burnin            The number of burn-in samples
#                 thin              The thinning parameter
#                 family            A description of the error distribution and link function to be used in the model
#                 pi                Cutoff for variable inclusion (defaults to 0.1)
#                 degree            Degree of B-splines (defaults to 3)
#                 k_max             Maximum number of knots for each spline term (defaults to 9)
#                 sigma_v           Heterogeneity parameter for "jump" terms
#                 lambda_1          Prior parameter for number of terms in the model (defaults to 0.1)
#                 lambda_2          Prior parameter for number of knots in each spline (defaults to 1)
#                 a_0               For continuous outcomes, prior parameter on individual-level heterogeneity
#                 b_0               For continuous outcomes, prior parameter on individual-level heterogeneity
#                 sigma_B           Prior variance for model coefficients
#                 w                 Window to propose knot location changes
#                 sigma_epsilon     Proposal variance for coefficient updates
#
# Returns:        If the procedure is successful, a list with the following elements is returned:
#                 - success: Indicates whether the procedure was successful based on geweke convergence
#                 - accept_var: A matrix with two columns and one row per iteration indicating whether the (1)
#                   proposed variable addition, or (2) the proposed variable removal was accepted
#                 - accept_add_knot: A matrix with a column for each spline term and a row for each iteration 
#                   indicating (1/0) whether a proposed knot addition was accepted
#                 - accept_remove_knot: A matrix with a column for each spline term and a row for each iteration 
#                   indicating (1/0) whether a proposed knot removal was accepted
#                 - accept_move_knot: A matrix with a column for each spline term and a row for each iteration 
#                   indicating (1/0) whether a proposed knot position change was accepted
#                 - accept_inter_trt: A matrix with a column for the intercept (1) and main effect of treatment (2),
#                   and a row for each iteration indicating (1/0) whether a proposed coefficient change was accepted
#                 - accept_binary: A matrix with a column for each binary term and a row for each iteration indicating 
#                   (1/0) whether a proposed coefficient change was accepted
#                 - accept_spline: A matrix with a column for each spline term and a row for each iteration containing
#                   the average acceptance rate for each coefficient within a term in that iteration
#                 - trt_eff_posterior: A matrix with one individual per row, and B columns per individual, describing the 
#                   complete posterior distribution of the treatment effect for each individual, based on their tailoring v
#                   variables
#                 - splines_fitted: A list with one element for each interaction spline term. For each element, there is a matrix 
#                   with B rows and a column for each individual containing the fitted value for the given spline term. This helps 
#                   us estimate fitted values for new data in the getCutoffMaley function
#                 - binary_param: A matrix with B rows containing the posterior distribution for each binary term
#                 - inter_trt_param: A matrix with B rows containing the posterior distribution for the intercept and main effect
#                   of treatment
#                 - sd: Posterior distribution of model standard deviation
#                 - k: Posterior distribution of number of knots for each spline term
#                 - included_vars: Selected tailoring variables based on cutoff pi 
#                 - candsplineinter: Candidate tailoring spline variables
#                 - candbinaryinter: Candidate tailoring binary variables
#                 - prop_incl: Posterior inclusion probability for each term
#                 If the procedure is not successful, the function returns list(success = FALSE)

runMHMaley <- function(data, 
                       candsplinevars, 
                       candbinaryvars, 
                       candinter, 
                       B,
                       burnin,
                       thin,
                       family = gaussian(link = "identity"),
                       pi = 0.1,
                       degree = 3,
                       k_max = 9,  
                       sigma_v = 0.1,
                       lambda_1 = 0.1,
                       lambda_2 = 1,
                       a_0 = 0.01,
                       b_0 = 0.01,
                       sigma_B = 10,
                       w = 1,
                       sigma_epsilon = 0.1
) {
  
  if(is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if(is.function(family)) family <- family()
  if(is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  iterations = B*thin + burnin
  pb <- txtProgressBar(min=1,max=iterations+1,style=3,width=50,char="=")
  
  
  inter_trt_param = array(dim = c(iterations+1,2))
  colnames(inter_trt_param) = c("intercept", "trt")
  spline_param = list()
  sd = array(dim = c(iterations+1,1))
  accept_var = matrix(nrow = iterations,
                   ncol = 2)
  accept_inter_trt = matrix(nrow = iterations,
                  ncol = 2)
  colnames(accept_inter_trt) = colnames(inter_trt_param)
  splines_fitted = list()

  candsplinevars_ext = c(paste0(candsplinevars,"_main", recycle0 = T),
                 paste0(intersect(candsplinevars, candinter),"_inter", recycle0 = T))
  candbinaryvars_ext = c(paste0(candbinaryvars,"_main", recycle0 = T),
                         paste0(intersect(candbinaryvars, candinter),"_inter", recycle0 = T))
  
  curspline_ext = candsplinevars_ext
  curbinary_ext = candbinaryvars_ext
  n_cand_vars = length(candsplinevars_ext) + length(candbinaryvars_ext)
  candmain = c(candsplinevars,candbinaryvars)
  curinter = candinter
  curmain = candmain
  
  accept_binary = matrix(nrow = iterations,
                  ncol = length(candbinaryvars_ext))
  colnames(accept_binary) = candbinaryvars_ext
  accept_spline = matrix(nrow = iterations,
                         ncol = length(candsplinevars_ext))
  colnames(accept_spline) = candsplinevars_ext


  if (length(candbinaryvars_ext) > 0) {
    binary_param = array(dim = c(iterations+1,length(candbinaryvars_ext)))
    colnames(binary_param) = candbinaryvars_ext
    if (length(intersect(candbinaryvars, candinter))== 0) {
      binary_mod_mat = data[candbinaryvars]
    } else {
      binary_mod_mat = cbind(data[candbinaryvars],
                             data[intersect(candbinaryvars, candinter)]*data$trt)
    }
   
    colnames(binary_mod_mat) = candbinaryvars_ext
  } else {
    binary_param = as.matrix(rep(0,iterations+1))
    binary_mod_mat = 0
  }
 
  
  # List of all variables in the model in their extended form
  candvars_ext = c(paste0(candmain,"_main"),
                  paste0(candinter,"_inter"))
  
  curvars_ext = candvars_ext
  
  nvars_max = length(candvars_ext)
  
  k = matrix(nrow = iterations + 1,
             ncol = length(candsplinevars_ext))
  accept_move_knot = matrix(nrow = iterations + 1,
                            ncol = length(candsplinevars_ext))
  accept_add_knot = matrix(nrow = iterations + 1,
                            ncol = length(candsplinevars_ext))
  accept_remove_knot = matrix(nrow = iterations + 1,
                            ncol = length(candsplinevars_ext))
  colnames(accept_move_knot) = candsplinevars_ext
  colnames(accept_add_knot) = candsplinevars_ext
  colnames(accept_remove_knot) = candsplinevars_ext
  
  vars_prop = matrix(nrow = iterations+1,
                         ncol = nvars_max)
  
  vars_prop[1,] <- rep(1,nvars_max)
  
  if (length(candsplinevars_ext) > 0) {
    knotscand = list()
    knotscur_idx = vector("list", length = length(candsplinevars_ext)) 
    names(knotscur_idx) = candsplinevars_ext
    
    for (i in 1:length(candsplinevars_ext)) {
      knotscand_i = quantile(data[[sub("_[^_]+$", "", candsplinevars_ext[i])]], 
                             seq(0,1,length.out=11))[-c(1,11)]
      
      knotscand[[candsplinevars_ext[i]]] = knotscand_i
    }
    
    
    spline_mod_mat = list()
    spline_mod_mat_raw = list()
    for (i in 1:length(candsplinevars_ext)) {
      mod_mat_i = bs(data[[sub("_[^_]+$", "", candsplinevars_ext[i])]],
                     degree = degree,
                     knots = knotscand[[candsplinevars_ext[i]]][knotscur_idx[[candsplinevars_ext[i]]]],
                     intercept = F)
      if (length(grep("inter",candsplinevars_ext[i]))>0) {
        spline_mod_mat_raw[[candsplinevars_ext[i]]] = mod_mat_i
        mod_mat_i = mod_mat_i*data$trt
      }
      
      spline_mod_mat[[candsplinevars_ext[i]]] = mod_mat_i
    }
    
    
    k[1,] <- unlist(lapply(knotscur_idx,length))
  } else {
    spline_mod_mat = 0
  }
  if (length(candbinaryvars_ext) == 0 & length(candsplinevars_ext) == 0) {
    mod_start <- glm(Y ~ trt, family=family)
  } else if (length(candbinaryvars_ext) == 0) {
    mod_start <- glm(Y ~ trt + 
                       do.call(cbind, spline_mod_mat),data=data,family=family)
  } else if (length(candsplinevars_ext) == 0) {
    mod_start <- glm(Y ~ trt + do.call(cbind,binary_mod_mat), data=data,family=family)
    binary_param[1,] = coef(mod_start)[-c(1,2)]
  } else {
    mod_start <- glm(Y ~ trt + do.call(cbind,binary_mod_mat) + 
                       do.call(cbind, spline_mod_mat),data=data,family=family)
    binary_param[1,] = coef(mod_start)[3:(length(candbinaryvars_ext)+2)]
  }

  
  inter_trt_param[1,] = c(coef(mod_start)[1], coef(mod_start)[2])
  if (family$family == "gaussian") {
    sd[1] = sigma(mod_start)
  }
  candsplineinter = intersect(candsplinevars, candinter)
  if (length(candsplinevars_ext) > 0) {
    ncoef_perx = k[1,] + degree
    coefs = split(coef(mod_start)[-c(1:(length(candbinaryvars_ext)+2))],rep(1:(length(ncoef_perx)),ncoef_perx))
    
    spline_ols_param = list()
    for (i in 1:length(candsplinevars_ext)) {
      spline_ols_param[[candsplinevars_ext[i]]] = unlist(coefs[i])
    }
    
    spline_param[[1]] = spline_ols_param
    
    for (l in 1:length(candsplineinter)) {
      splines_fitted[[candsplineinter[l]]] = matrix(nrow=iterations+1,
                                                    ncol=nrow(data))
      splines_fitted[[candsplineinter[l]]][1,] = splinesFitted(0,
                                                 spline_param[[1]][[paste0(candsplineinter[l],"_inter")]], 
                                                 spline_mod_mat_raw[[paste0(candsplineinter[l],"_inter")]])

    }
  }
  

  for (i in 1:iterations){
    setTxtProgressBar(pb, i)
    if (length(candsplinevars_ext) > 0) {
      spline_param[[i+1]] = spline_param[[i]]
      k[i+1,] = k[i,]
      for (j in 1:length(candsplinevars_ext)) {
        j_name = candsplinevars_ext[j]
        if (j_name %in% curspline_ext) {
          knotscand_x = knotscand[[j_name]]
          knotscur_idx_x = knotscur_idx[[j_name]]
          knotscur_x = knotscand_x[knotscur_idx_x]


          # move knot
          move_knot = moveKnot(data$Y, data, w, j_name, j, knotscand_x, knotscur_x,
                               knotscur_idx_x, binary_param[i,],
                               binary_mod_mat, spline_mod_mat, spline_mod_mat_raw,
                               spline_ols_param[[j_name]],
                               inter_trt_param[i,], spline_param[[i+1]], sd[i], family)

          knotscur_idx_x = move_knot$knotscur_idx_x
          spline_ols_param[[j_name]] = move_knot$spline_ols_param_x
          spline_mod_mat[[j_name]] = move_knot$spline_mod_mat_x
          if (length(grep("inter",j_name))>0) {
            spline_mod_mat_raw[[j_name]] = move_knot$spline_mod_mat_raw_x
          }
          knotscur_x = move_knot$knotscur_x
          accept_move_knot[i,] = move_knot$accept

          v <- rnorm(1,0,sigma_v)
          u_1 <- runif(1)

          # Add a knot
          if (u_1 < 0.5) {
            add_remove_knot = addKnot(data$Y, data, v, sigma_v, sigma_B, lambda_2, j_name, j,
                                      k[i,j], k_max, knotscand_x, knotscur_x,
                                      knotscur_idx_x,binary_param[i,],
                                      binary_mod_mat, spline_mod_mat, spline_mod_mat_raw,
                                      spline_ols_param[[j_name]],
                                      inter_trt_param[i,], spline_param[[i+1]], sd[i], family)
            accept_add_knot[i,] = add_remove_knot$accept
          }

          # remove knot
          if (u_1 >= 0.5) {
            add_remove_knot = removeKnot(data$Y, data, v, sigma_v, sigma_B, lambda_2, j_name, j,
                                         k[i,j], k_max, knotscand_x, knotscur_x,
                                         knotscur_idx_x, binary_param[i,],
                                         binary_mod_mat, spline_mod_mat, spline_mod_mat_raw,
                                         spline_ols_param[[j_name]],
                                         inter_trt_param[i,], spline_param[[i+1]], sd[i], family)
            accept_remove_knot[i,] = add_remove_knot$accept
          }

          # Update spline parameters
          spline_param[[i+1]][[j_name]] = add_remove_knot$spline_param_x

          # Update k parameters
          k[i+1,j] = add_remove_knot$k

          # Update knot indices
          knotscur_idx_x = add_remove_knot$knotscur_idx_x

          # Update knot values
          knotscur_x = add_remove_knot$knotscur_x

          # Update current ols spline parameters (used to scale coefficients)
          spline_ols_param[[j_name]] = add_remove_knot$spline_ols_param_x

          # Update spline model matrix
          spline_mod_mat[[j_name]] = add_remove_knot$spline_mod_mat_x
          if (length(grep("inter",j_name))>0) {
            spline_mod_mat_raw[[j_name]] = add_remove_knot$spline_mod_mat_raw_x
          }
          
          
          accept_spline_coef = vector(length=length(spline_param[[i+1]][[j_name]]))
          # Update all splines
          for (coef in 1:(length(spline_param[[i+1]][[j_name]]))) {
            spline_change_propose = spline_param[[i+1]][[j_name]][coef] +
              rnorm(1,0,sigma_epsilon)
            
            splines_propose = spline_param[[i+1]]
            splines_propose[[j_name]][coef] = spline_change_propose
            
            log_prob = logLikelihoodCustom(data$Y,
                                           data$trt,
                                           inter_trt_param[i,],
                                           binary_param[i,],
                                           binary_mod_mat,
                                           splines_propose,
                                           spline_mod_mat,
                                           sd[i], family) -
              logLikelihoodCustom(data$Y,
                                  data$trt,
                                  inter_trt_param[i,],
                                  binary_param[i,],
                                  binary_mod_mat,
                                  spline_param[[i+1]],
                                  spline_mod_mat,
                                  sd[i], family) +
              dnorm(spline_change_propose, mean = 0, sd = sigma_B, log = T) -
              dnorm(spline_param[[i+1]][[j_name]][coef], mean = 0,
                        sd = sigma_B, log = T)
            gamma <- runif(1)
            if (gamma < min(1,exp(log_prob))) {
              spline_param[[i+1]][[j_name]][coef] = spline_change_propose
              accept_spline_coef[coef] = 1
            } else {
              accept_spline_coef[coef] = 0
            }
          }
          accept_spline[i,j] = mean(accept_spline_coef)
          knotscur_idx[[j_name]] = knotscur_idx_x
        }
      } 
    } else {
      spline_param[[i+1]] = 0
    }
   
    # Next, we add or remove one term (with equal probabibility). First, we 
    # compute the set of terms that are eligible to be added (eligibletoadd) and
    # the set eligible to be removed (eligibletoremove)
    internotinmodel = c(setdiff(candinter,curinter))
    eligibletoremove = c(paste_(setdiff(curmain, curinter), "_main"),
                         paste_(curinter, "_inter"))

    eligibletoadd = c(paste_(c(setdiff(candmain,curmain)),"_main"),
                      paste_(internotinmodel[internotinmodel %in% curmain],"_inter"))

    n_cur_vars = length(curmain) + length(curinter)

    u_2 = runif(1)
    if (u_2 < 0.5) {
      add_remove_var = addVar(data$Y,
                              data$trt,
                              n_cand_vars,
                              n_cur_vars,
                               sigma_v,
                               sigma_B,
                               lambda_1,
                               eligibletoremove,
                               eligibletoadd,
                               curvars_ext,
                               curmain,
                               curinter,
                              curbinary_ext,
                              curspline_ext,
                              candbinaryvars_ext,
                               inter_trt_param[i,],
                               binary_param[i,],
                               binary_mod_mat,
                               spline_param[[i+1]],
                               spline_mod_mat,
                               sd[i],
                               family)
      accept_var[i,1] = add_remove_var$accept
    }

    if (u_2 >= 0.5) {
      add_remove_var = removeVar(data$Y,
                                 data$trt,
                                 n_cand_vars,
                                 n_cur_vars,
                                sigma_v,
                                sigma_B,
                                lambda_1,
                                eligibletoremove,
                                eligibletoadd,
                                curvars_ext,
                                curmain,
                                curinter,
                                curbinary_ext,
                                curspline_ext,
                                candbinaryvars_ext,
                                inter_trt_param[i,],
                                binary_param[i,],
                                binary_mod_mat,
                                spline_param[[i+1]],
                                spline_mod_mat,
                                sd[i],
                                family)
      accept_var[i,2] = add_remove_var$accept

    }

    spline_param[[i+1]] = add_remove_var$spline_param_i
    inter_trt_param[i+1,] = add_remove_var$inter_trt_param_i
    binary_param[i+1,] = add_remove_var$binary_param_i
    curmain = add_remove_var$curmain
    curinter = add_remove_var$curinter
    curvars_ext = add_remove_var$curvars_ext
    curbinary_ext = add_remove_var$curbinary_ext
    curspline_ext = add_remove_var$curspline_ext
    vars_prop[i+1,] = as.numeric(candvars_ext %in% curvars_ext)
    
    # Update 
    for (j in 1:2) {
      inter_trt_propose = inter_trt_param[i+1,]
      single_inter_trt_propose = inter_trt_param[i+1,j] +
        rnorm(1,0,sigma_epsilon)
      inter_trt_propose[j] = single_inter_trt_propose
      log_prob = logLikelihoodCustom(data$Y,
                                     data$trt,
                                     inter_trt_propose,
                                     binary_param[i+1,],
                                     binary_mod_mat,
                                     spline_param[[i+1]],
                                     spline_mod_mat,
                                     sd[i],
                                     family) -
        logLikelihoodCustom(data$Y,
                            data$trt,
                            inter_trt_param[i+1,],
                            binary_param[i+1,],
                            binary_mod_mat,
                            spline_param[[i+1]],
                            spline_mod_mat,
                            sd[i],
                            family) +
        sum(dnorm(single_inter_trt_propose, mean = 0, sd = sigma_B, log = T)) -
        sum(dnorm(inter_trt_param[i+1,j], mean = 0,
                  sd = sigma_B, log = T))
      
      gamma <- runif(1)
      if (gamma < min(1,exp(log_prob))) {
        inter_trt_param[i+1,] = inter_trt_propose
        accept_inter_trt[i,j] = 1
      } else {
        accept_inter_trt[i,j] = 0
      }
    }
    
    if (length(curbinary_ext) > 0) {
      for (j in 1:ncol(binary_param)) {
        if (candbinaryvars_ext[j] %in% curbinary_ext) {
          binary_param_propose = binary_param[i+1,]
          single_binary_param_propose = binary_param[i+1,j] +
            rnorm(1,0,sigma_epsilon)
          binary_param_propose[j] = single_binary_param_propose
          log_prob = logLikelihoodCustom(data$Y,
                                         data$trt,
                                         inter_trt_param[i+1,],
                                         binary_param_propose,
                                         binary_mod_mat,
                                         spline_param[[i+1]],
                                         spline_mod_mat,
                                         sd[i],
                                         family) -
            logLikelihoodCustom(data$Y,
                                data$trt,
                                inter_trt_param[i+1,],
                                binary_param[i+1,],
                                binary_mod_mat,
                                spline_param[[i+1]],
                                spline_mod_mat,
                                sd[i],
                                family) +
            sum(dnorm(single_binary_param_propose, mean = 0, sd = sigma_B, log = T)) -
            sum(dnorm(binary_param[i+1,j], mean = 0,
                      sd = sigma_B, log = T))

          gamma <- runif(1)
          if (gamma < min(1,exp(log_prob))) {
            binary_param[i+1,] = binary_param_propose
            accept_binary[i,j] = 1
          } else {
            accept_binary[i,j] = 0
          }
        }
      }
    }
  
    # Update sigma^2
    if (length(grep("gaussian",family))>0) {
      resids = data$Y - computeFittedValues(data$trt,
                                           inter_trt_param[i+1,],
                                           binary_param[i+1,],
                                           binary_mod_mat,
                                           spline_param[[i+1]],
                                           spline_mod_mat
        
      ) 
      sd[i+1] <- 1/rgamma(1, shape=nrow(data)/2+a_0, rate=sum(resids^2)/2+b_0)
    }
    
    # 
    # Update treatment effect, use the "raw" model.matrix which contains the 
    # spline matrix for the interaction effect, before being multiplied by the
    # treatment indicator
    if (length(candsplinevars_ext) > 0) {
      for (l in 1:length(candsplineinter)) {
        splines_fitted[[candsplineinter[l]]][i+1,] = splinesFitted(0,
                                                     spline_param[[i+1]][[paste0(candsplineinter[l],"_inter")]],
                                                     spline_mod_mat_raw[[paste0(candsplineinter[l],"_inter")]])
      
      }
    }
    
  }
  
  close(pb)
  colnames(accept_var) = c("add var", "remove var")
  colnames(accept_move_knot)
  colnames(accept_add_knot)
  colnames(accept_remove_knot)
  
  colnames(vars_prop) = candvars_ext
  colnames(k) = candsplinevars_ext
  
  ############################################ MCMC diagnostics ###############################################
  
  vars_prop_res = vars_prop[-(1:burnin),]
  names(vars_prop_res) = names(vars_prop)
  vars_prop_res = vars_prop_res[seq(1,nrow(vars_prop_res),thin),]
  splines_fitted_res = 0
  if (length(candsplinevars_ext)>0) {
    splines_fitted_res = list()
    for(l in 1:length(candsplineinter)){
      splines_fitted_res[[candsplineinter[l]]] = splines_fitted[[candsplineinter[l]]][-(1:burnin),]
      splines_fitted_res[[candsplineinter[l]]] = splines_fitted_res[[candsplineinter[l]]][seq(1,nrow(splines_fitted_res[[candsplineinter[l]]]),thin),]
    }
  }

  binary_param_res = 0
  if (length(candbinaryvars_ext) > 0) {
    binary_param_res = binary_param[-(1:burnin),]
    names(binary_param_res) = names(binary_param)
    binary_param_res = binary_param_res[seq(1,nrow(binary_param_res),thin),]
  }
   
  inter_trt_param_res = inter_trt_param[-(1:burnin),]
  names(binary_param_res) = names(inter_trt_param)
  inter_trt_param_res = inter_trt_param_res[seq(1,nrow(inter_trt_param_res),thin),]
  
  k_res = k[-(1:burnin),]
  names(k_res) = names(k)
  k_res = k_res[seq(1,nrow(k_res),thin),]
  
  sd_res = sd[-(1:burnin)]
  sd_res = sd_res[seq(1,length(sd_res),thin)]
  
  
  prop_incl = colMeans(vars_prop_res)
  included_vars = names(which(prop_incl > pi))
  included_vars = included_vars[grep("inter",included_vars)]
  included_vars = sub("_[^_]+$", "", included_vars)

  trt_eff_posterior = matrix(1,nrow=nrow(data)) %*% t(matrix(inter_trt_param_res[,2]))
  if (length(candsplineinter) > 0) {
    for (m in 1:length(candsplineinter)) {
      trt_eff_posterior = trt_eff_posterior + t(splines_fitted_res[[candsplineinter[m]]])
    }
  } 
  
  candbinaryinter = intersect(candbinaryvars,candinter)
  binary_param_inter_res = 0
  if (length(candbinaryinter) > 0) {
    binary_param_inter_res = binary_param_res[,paste0(candbinaryinter,"_inter")]
    trt_eff_posterior = trt_eff_posterior + 
      as.matrix(data[,candbinaryinter]) %*% t(binary_param_inter_res)
  }
  
  geweke.trt_eff_posterior <- rep(NA, nrow(data))
  names(geweke.trt_eff_posterior) = 1:nrow(data)
  for (t in 1:nrow(data)) {
    geweke.trt_eff_posterior[t] <- geweke.diag(trt_eff_posterior[t,], frac1=0.25, frac2=0.25)[[1]]
  }
  
  geweke.trt_param_res <- geweke.diag(inter_trt_param_res[,2], frac1=0.25, frac2=0.25)[[1]]
  
  
  if (family$family == "gaussian") {
    geweke.sd_res <- geweke.diag(sd_res, frac1=0.25, frac2=0.25)[[1]]
  } else {
    geweke.sd_res = 0
  }

  # Assess convergence
  geweke.conv <- !(max(abs(geweke.trt_eff_posterior))>4 | max(abs(geweke.sd_res))>4)

  if (geweke.conv) {
    return(list(
      success = TRUE,
      accept_var = accept_var,
      accept_add_knot = accept_add_knot,
      accept_remove_knot = accept_remove_knot,
      accept_move_knot = accept_move_knot,
      accept_inter_trt = accept_inter_trt,
      accept_binary = accept_binary,
      accept_spline = accept_spline,
      trt_eff_posterior = trt_eff_posterior,
      splines_fitted = splines_fitted_res,
      binary_param = binary_param_res,
      inter_trt_param = inter_trt_param_res,
      sd = sd_res,
      vars_prop = vars_prop_res,
      k = k_res,
      included_vars = included_vars,
      candsplineinter = candsplineinter,
      candbinaryinter = candbinaryinter,
      prop_incl = prop_incl))
  } else {
    return(list(success = FALSE))
  }

}

#
# function parseMHMaley
#
# Author: Lara Maleyeff
#
# Function:       parseMHMaley
# Description:    Function to parse MCMC results of proposed method
# Parameters:     curdata           A data frame with the observations arranged by row, and including the column:
#                                   - trt: the group indicator (=1 for experimental group; =0 for control group)
#                                   - Y: outcome
#                                   - variables (candsplinevars, candbinaryvars) in the columns
#                 candpslinevars    Vector with names of continuous variables
#                 candbinaryvars    Vector with names of binary variables
#                 family            Family and link function of outcome, defaults to gaussian()
#                 B                 The number of posterior samples
#                 burnin            The number of burn-in samples
#                 thin              The thinning parameter
#                 pi                Not used
#
# Returns:        If the fitting procedure was successful, the function returns a list with:
#                 - unparsed_results: The results of MCMC procedure
#                 - trt_eff_posterior: The posterior distribution of treatment effect for each individual (no. ind x B)
#                 - prop_incl: Posterior probability of inclusion for each variable
#                 - data_fit: Data used to fit the model, used later for linear interpolation on an external dataset
#                 - included_vars: Selected tailoring variables
#                 - candsplineinter: candidate continuous tailoring variables (because runMHMaley allows for set of predictive
#                   variables to differ from tailoring variables)
#                 - candbinaryinter: Candidate binary tailoring variables
#                 - trt_param: Posterior distribution of treamtent effect
#                 - binary_param: Posterior distribution of the binary tailoring effects
#                 - splines_fitted: List of matrices containing the posterior distribution of the fitted values for 
#                   inidividual for each spline term
#                 - success: Indicates whether the procedure was successful based on geweke convergence
#                 If the fitting procedure failed, returns list(success = FALSE))

parseMHMaley <- function(curdata, 
                         candsplinevars, 
                         candbinaryvars,
                         family,
                         B, 
                         burnin,
                         thin,
                         pi,
                         ...) {
  
  candinter = c(candsplinevars, candbinaryvars)
  unparsed_results = runMHMaley(curdata, 
                                candsplinevars, 
                                candbinaryvars, 
                                candinter, 
                                B,
                                burnin,
                                thin,
                                family,
                                pi,
                                ...)
  
  if (!unparsed_results$success) {
    return(list(success = FALSE))
  } else {
    return(list(unparsed_results = unparsed_results,
                trt_eff_posterior = unparsed_results$trt_eff_posterior,
                prop_incl = unparsed_results$prop_incl,
                data_fit = curdata,
                included_vars = unparsed_results$included_vars,
                candsplineinter = unparsed_results$candsplineinter,
                candbinaryinter = unparsed_results$candbinaryinter,
                trt_param = t(matrix(unparsed_results$inter_trt_param[,2])),
                binary_param = unparsed_results$binary_param,
                splines_fitted = unparsed_results$splines_fitted,
                success = TRUE
  ))
  }
}
#
# function getCutoffMaley
#
# Author: Lara Maleyeff
#
# Function:       getCutoffMaley
# Description:    Internal function to find the alpha-quantile of the posterior treatment effect distribution for
#                 each individual (row). It uses the model that was fit in the most recent interim analysis (original data) to 
#                 find the alpha-quantiles of external data (used for either the adaptive enrichment steps or 
#                 accuracy calculations)
# Parameters:     data              A data frame with one individual per row and information on candsplineinter and candbinaryinter
#                 candsplinevars    Vector with names of continuous variables (not used)
#                 candbinaryvars    Vector with names of binary variables (not used)
#                 trial_results     Results from most recent interim analysis, a list containing:
#                                   - included_vars: selected tailoring variables
#                                   - candsplineinter: candidate spline tailoring variables
#                                   - candbinaryinter: candidate binary tailoring variables
#                                   - trt_param: posterior distribution of treament effect
#                                   - data_fit: original data used to fit the model, used for interpolation
#                                   - splines_fitted: posterior distribution of fitted values for each spline
#                                   - binary_param: posterior distribution of binary coefficients
#                 alpha             Cutoff for the effective subspace
#
# Returns:        alpha-row quantile for each individuals to then be compared with e_1. If
#                 the alpha-quantile of the treatment effect is > e_1 then that individual's
#                 variable combination is in the effective subspace
getCutoffMaley <- function(data,
                           candsplinevars,
                           candbinaryvars,
                           trial_results,
                           alpha
                           ) {

  included_vars = trial_results$included_vars
  candsplineinter = trial_results$candsplineinter
  candbinaryinter = trial_results$candbinaryinter
  
  combined_posterior = matrix(1,nrow=nrow(data)) %*% trial_results$trt_param
  
  if (length(candsplineinter) > 0) {
    for (m in 1:length(candsplineinter)) {
      interpolated_splines_fitted <- sapply(1:nrow(trial_results$splines_fitted[[candsplineinter[m]]]), function(j) {
        approx(x = trial_results$data_fit[,candsplineinter[m]], 
               y = trial_results$splines_fitted[[candsplineinter[m]]][j,], xout = data[,candsplineinter[m]],rule=2)$y
      })
      combined_posterior = combined_posterior + interpolated_splines_fitted
    }
  }
  if (length(candbinaryinter) > 0) {
    combined_posterior = combined_posterior + 
      as.matrix(data[,candbinaryinter]) %*% t(trial_results$binary_param[,paste0(candbinaryinter,"_inter")])
  }
  
  return(rowQuantiles(combined_posterior,probs = c(alpha)))

}


