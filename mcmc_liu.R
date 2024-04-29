source("mcmc_liu_helper.R")
######################################################################################################
#   Implementation of MCMC procedure of Liu et al. (2022), described in Maleyeff et al. (2024)      
#                           Adapted from Liu et al. (2022), by Lara Maleyeff
#         Original code: https://github.com/YushaLiu/Continuous_marker_adaptive_trial_design
#                                   Contact: laramaleyeff@gmail.com                                  
#                                       Last updated: April 2024                                       
######################################################################################################

#
# function runMHLiu
#
# Author: Liu et al. (2022), with adaptations by Lara Maleyeff
#
# Function:       runMHLiu
# Description:    Perform MCMC to generate posterior samples of penalized spline coefficients and regularization 
#                 parameters
# Parameters:     curdata           A data frame with the observations arranged by row, and including the column:
#                                   - trt: the group indicator (=1 for experimental group; =0 for control group)
#                                   - Y: outcome
#                 candpslinevars    Vector with names of continuous variables
#                 family            Family and link function of outcome, defaults to continuous
#                 B                 The number of posterior samples
#                 burnin            The number of burn-in samples
#                 thin              The thinning parameter
#                 numIntKnots       The number of interior knots for each spline term
#                 multiplier1       Factor to adaptively adjust var_g and var_h proposals (default 0.5)
#                 multiplier2       Factor to adaptively adjust var_jh and var_jg proposals (default 0.5)
#                 multiplier3       Factor to adaptively adjust var_jh and var_jg proposals (default 0.25)
#                 prior             List with prior parameter information on alpha, beta, and sigma_B, 
#                                   defaults to list(alpha = 0.01, beta = 0.01, sigma_B = 10)
#
# Returns:        Results of MCMC procedure
runMHLiu <- function(curdata, 
                     candsplinevars,
                     family = gaussian(link = "identity"),
                     B = 10000,
                     burnin = 5000,
                     thin = 5,
                     numIntKnots = 5,
                     multiplier1 = 0.5,
                     multiplier2 = 0.5,
                     multiplier3 = 0.25,
                     prior = list(alpha = 0.01, beta = 0.01, sigma_B = 10)
) {
  
  library(MASS)
  library(Matrix)
  
  X_cont = curdata[,candsplinevars, drop=FALSE]
  xrange = t(apply(as.matrix(X_cont), 2, range))
  splines = get_splines(curdata, 
                        X_cont, 
                        xrange, 
                        numIntKnots
  )
  
  if(any(splines$stabCheck)){
    return(list(success = FALSE))
  }
  
  data <- splines$data
  formula <- splines$formula
  ncolz <- splines$ncolz
  data.X <- data[, -c(1)]
  
  # For other exponential family outcomes, adjust the log-likelihood
  logLikelihood <- function(Y, betas, var) {
    pred = betas %*% t(data.X)
    if (family$family == "binomial") {
      singlelikelihoods = dbinom(Y, length(Y), prob = family$linkinv(pred), log = T)
    } else if (family$family == "gaussian") {
      singlelikelihoods = dnorm(Y, mean = pred, sd = sqrt(var), log = T)
    }
    sumll = sum(singlelikelihoods)
    return(sumll)   
  }
  
  logPrior <- function(betas, sigma_B) {
    prioreachbeta <- dnorm(betas,0,sigma_B, log=T) 
    priorbetas <- sum(prioreachbeta)
    return(priorbetas)
  }
  
  # calculate the number of markers from the input data
  J <- ncol(X_cont)
  
  ####################################### MCMC initialization ##############################
  # Determine the initial values of the MCMC using a simple regression analysis
  initfit <- glm(formula=formula, data=data, family=family)
  theta <- initfit$coef[-c(1)]
  theta_0 <- initfit$coef[1]
  sigma2_r = sigma(initfit)^2
  
  # determine the initial proposal variance for Metropolis random walk of the 
  # regression coefficients for each marker
  # Determine the initial value of sigma2_g and sigma2_h, i.e., 
  # regularization parameters for each marker
  var_g <- list(NULL)
  var_h <- list(NULL)
  sigma2_g <- rep(NA, J)
  sigma2_h <- rep(NA, J)
  
  # loop over each continuous marker
  for(j in 1:J){
    if(j==1){
      # prognostic marker effects
      idx.jg.start <- 1
      idx.jg.end <- idx.jg.start + ncolz
      data_jg <- data.X[, idx.jg.start:idx.jg.end]
      theta_jg <- theta[idx.jg.start:idx.jg.end]
      u_jg <- theta_jg[-c(1)]
      data_jg$offset <- as.numeric(as.matrix(data.X[,-c(idx.jg.start:idx.jg.end)]) %*% theta[-c(idx.jg.start:idx.jg.end)])
      
      # predictive marker effects
      idx.jh.start <- idx.jg.end + 1
      idx.jh.end <- idx.jh.start + ncolz + 1
      data_jh <- data.X[, idx.jh.start:idx.jh.end]
      theta_jh <- theta[idx.jh.start:idx.jh.end]
      u_jh <- theta_jh[-c(1,2)]
      data_jh$offset <- as.numeric(as.matrix(data.X[,-c(idx.jh.start:idx.jh.end)]) %*% theta[-c(idx.jh.start:idx.jh.end)])
    } else {
      # prognostic marker effects
      idx.jg.start <- idx.jh.end + 1
      idx.jg.end <- idx.jg.start + ncolz
      data_jg <- data.X[, idx.jg.start:idx.jg.end]
      theta_jg <- theta[idx.jg.start:idx.jg.end]
      u_jg <- theta_jg[-c(1)]
      data_jg$offset <- as.numeric(as.matrix(data.X[,-c(idx.jg.start:idx.jg.end)]) %*% theta[-c(idx.jg.start:idx.jg.end)])   
      
      # predictive marker effects
      idx.jh.start <- idx.jg.end + 1
      idx.jh.end <- idx.jh.start + ncolz
      data_jh <- data.X[, idx.jh.start:idx.jh.end]
      theta_jh <- theta[idx.jh.start:idx.jh.end]
      u_jh <- theta_jh[-c(1)]
      data_jh$offset <- as.numeric(as.matrix(data.X[,-c(idx.jh.start:idx.jh.end)]) %*% theta[-c(idx.jh.start:idx.jh.end)])
    }
    
    # determine the initial proposal variance for Metropolis random walk for coefficients related to prognostic marker effects
    formula1 <- as.formula(paste("Y ~ ", paste(colnames(data_jg), collapse="+")))
    initfit_jg <- glm(formula=formula1, data=cbind(Y=data[,1],data_jg),family=family)
    if(j==1){
      var_g[[j]] <- (multiplier1^2)*vcov(initfit_jg)[-c(ncol(data_jg)+1),-c(ncol(data_jg)+1)]
    }
    else{
      var_g[[j]] <- (multiplier1^2)*vcov(initfit_jg)[-c(1,ncol(data_jg)+1),-c(1,ncol(data_jg)+1)]
    }
    
    
    # determine the initial value of sigma2_jg
    sigma2_g[j] <- 1/rgamma(1, shape=ncolz/2+prior$alpha, rate=sum(u_jg^2)/2+prior$beta)
    
    # determine the initial proposal variance for Metropolis random walk for coefficients related to predictive marker effects
    formula2 <- as.formula(paste("Y ~ ", paste(colnames(data_jh), collapse="+") ))
    initfit_jh <- glm(formula=formula2, data=cbind(Y = data[,1],data_jh),family=family)
    var_h[[j]] <- (multiplier1^2)*vcov(initfit_jh)[-c(1, ncol(data_jh)+1),-c(1,ncol(data_jh)+1)]
    
    # determine the initial value of sigma2_jh
    sigma2_h[j] <- 1/rgamma(1, shape=ncolz/2+prior$alpha, rate=sum(u_jh^2)/2+prior$beta)
  } 
  
  
  data.X <- cbind(1, as.matrix(data.X))
  colnames(data.X) <- c("intercept", names(theta))
  theta = c(initfit$coef[1],theta)
  
  
  ########################################### MCMC starts #####################################
  L <- burnin+B*thin
  
  theta.all <- matrix(NA, nrow=L, ncol=length(theta))
  colnames(theta.all) <- colnames(data.X)
  theta.res <- matrix(NA, nrow=B, ncol=length(theta))
  colnames(theta.res) <- colnames(data.X)
  sigma2.g.all <- matrix(NA, nrow=L, ncol=J)
  sigma2.g.res <- matrix(NA, nrow=B, ncol=J)
  sigma2.h.all <- matrix(NA, nrow=L, ncol=J)
  sigma2.h.res <- matrix(NA, nrow=B, ncol=J)
  sigma2.r.all <- matrix(NA, nrow=L, ncol=1)
  sigma2.r.res <- matrix(NA, nrow=B, ncol=1)
  
  ll <- 0
  acpt_g <- rep(0, J)
  acpt_h <- rep(0, J) 
  
  for(l in 1:L){
    ### update theta_jg, sigma2_jg, theta_jh, sigma2_jh for each marker j
    for(j in 1:J){
      # prognostic marker effects 
      if(j==1){
        # indices for theta_jg
        idx.jg.start <- 1
        idx.jg.end <- idx.jg.start + ncolz + 1
        
        # current value of theta_jg
        theta_jg <- theta[idx.jg.start:idx.jg.end]
        u_jg <- theta_jg[-c(1,2)]
      }
      
      else{
        # indices for theta_jg
        idx.jg.start <- idx.jh.end + 1
        idx.jg.end <- idx.jg.start + ncolz
        
        # current value of theta_jg
        theta_jg <- theta[idx.jg.start:idx.jg.end]
        u_jg <- theta_jg[-c(1)]
      }
      
      data_jg <- data.X[, idx.jg.start:idx.jg.end]
      coef_jg <- as.numeric(data_jg %*% theta_jg)
      coef_all <- as.numeric(data.X %*% theta)
      
      # propose new value for theta_jg
      newtheta_jg <- mvrnorm(1, theta_jg, var_g[[j]])
      if(j==1){
        newu_jg <- newtheta_jg[-c(1,2)]
      }
      else{
        newu_jg <- newtheta_jg[-c(1)]
      }
      
      newtheta <- theta
      newtheta[idx.jg.start:idx.jg.end] <- newtheta_jg
      newcoef_jg <- as.numeric(data_jg %*% newtheta_jg)
      newcoef_all <- as.numeric(data.X %*% newtheta)
      
      # calculate acceptance probability
      log_acpt <- logLikelihood(data$Y, newtheta, sigma2_r)  - 0.5*sum(newu_jg^2)/sigma2_g[j] - 
        logLikelihood(data$Y, theta, sigma2_r)  + 0.5*sum(u_jg^2)/sigma2_g[j]
      
      # decide whether to accept the proposed value for theta_jg
      if(log(runif(1)) < log_acpt){
        theta <- newtheta
        u_jg <- newu_jg
        acpt_g[j] <- acpt_g[j] + 1
      }
      
      sigma2_g[j] <- 1/rgamma(1, shape=ncolz/2+prior$alpha, rate=sum(u_jg^2)/2+prior$beta)
      

      
      
      ### predictive marker effects
      if(j==1){
        # indices for theta_jh
        idx.jh.start <- idx.jg.end + 1
        idx.jh.end <- idx.jh.start + ncolz + 1
        
        # current values of theta_jh
        theta_jh <- theta[idx.jh.start:idx.jh.end]
        u_jh <- theta_jh[-c(1,2)]
      }
      
      else{
        # indices for theta_jg
        idx.jh.start <- idx.jg.end + 1
        idx.jh.end <- idx.jh.start + ncolz
        
        # current values of theta_jh
        theta_jh <- theta[idx.jh.start:idx.jh.end]
        u_jh <- theta_jh[-c(1)]
      }
      
      data_jh <- data.X[, idx.jh.start:idx.jh.end]
      coef_jh <- as.numeric(data_jh %*% theta_jh)
      coef_all <- as.numeric(data.X %*% theta)
      
      # propose new value for theta_jh
      newtheta_jh <- mvrnorm(1, theta_jh, var_h[[j]])
      if(j==1){
        newu_jh <- newtheta_jh[-c(1,2)]
      }
      else{
        newu_jh <- newtheta_jh[-c(1)]
      }
      
      newtheta <- theta
      newtheta[idx.jh.start:idx.jh.end] <- newtheta_jh
      newcoef_jh <- as.numeric(data_jh %*% newtheta_jh)
      newcoef_all <- as.numeric(data.X %*% newtheta)
      
      # calculate acceptance probability
      log_acpt <- logLikelihood(data$Y, newtheta, sigma2_r) - 0.5*sum(newu_jh^2)/sigma2_h[j] - 
        logLikelihood(data$Y, theta, sigma2_r)  + 0.5*sum(u_jh^2)/sigma2_h[j]
      
      # decide whether to accept the proposed value for theta_jh
      if(log(runif(1)) < log_acpt){
        theta <- newtheta
        u_jh <- newu_jh
        acpt_h[j] <- acpt_h[j] + 1
      }
      
      # update sigma2_jh
      sigma2_h[j] <- 1/rgamma(1, shape=ncolz/2+prior$alpha, rate=sum(u_jh^2)/2+prior$beta)
      
    }
    
    
    # update sd
    # update sigma2_jg
    if (family$family == "gaussian") {
      resids = data$Y - theta %*% t(data.X)
      sigma2_r <- 1/rgamma(1, shape=nrow(data)/2+prior$alpha, rate=sum(resids^2)/2+prior$beta)
    }
      ############################# save the values of the new iteration ###############################
    theta.all[l,] <- theta
    sigma2.g.all[l,] <- sigma2_g
    sigma2.h.all[l,] <- sigma2_h
    sigma2.r.all[l,] <- sigma2_r
    
    if( (l>burnin) && ((l-burnin)%%thin==0) ){
      ll <- ll+1
      theta.res[ll,] <- theta
      sigma2.g.res[ll,] <- sigma2_g
      sigma2.h.res[ll,] <- sigma2_h
      sigma2.r.res[ll,] <- sigma2_r
    }
    
    
    ############# adaptively adjust the proposal variance based on the cumulative acceptance probability #############
    if( (l>2000) && (l%%2000==1) ){
      for(j in 1:J){
        ### prognostic marker effects
        if(j==1){
          idx.jg.start <- 1
          idx.jg.end <- idx.jg.start + ncolz + 1
        }
        else{
          idx.jg.start <- idx.jh.end + 1
          idx.jg.end <- idx.jg.start + ncolz
        }
        
        # calculate the covariance matrix of theta_jg
        var_jg <- 2.4^2*(cov(theta.all[1:(l-1), idx.jg.start:idx.jg.end]) + 
                           1e-10*diag(rep(1,length(idx.jg.start:idx.jg.end))))/length(idx.jg.start:idx.jg.end)
        
        # adjust adaptively
        if(acpt_g[j]/l < 0.1){
          var_jg <- (multiplier3^2)*var_jg
        }
        
        else if(acpt_g[j]/l < 0.2){
          var_jg <- (multiplier2^2)*var_jg
        }
        
        var_g[[j]] <- var_jg
        
        
        ### predictive marker effects
        if(j==1){
          idx.jh.start <- idx.jg.end + 1
          idx.jh.end <- idx.jh.start + ncolz + 1
        }
        else{
          idx.jh.start <- idx.jg.end + 1
          idx.jh.end <- idx.jh.start + ncolz
        }
        
        # calculate the covariance matrix of theta_jh
        var_jh <- 2.4^2*(cov(theta.all[1:(l-1), idx.jh.start:idx.jh.end]) + 
                           1e-10*diag(rep(1,length(idx.jh.start:idx.jh.end))))/length(idx.jh.start:idx.jh.end)
        
        # adjust adaptively
        if(acpt_h[j]/l < 0.1){
          var_jh <- (multiplier3^2)*var_jh
        }
        
        else if(acpt_h[j]/l < 0.2){
          var_jh <- (multiplier2^2)*var_jh
        }
        var_h[[j]] <- var_jh
      }
    }
  }
  ############################################ MCMC diagnostics ###############################################
  geweke.theta <- rep(NA, ncol(theta.res))
  names(geweke.theta) <- colnames(theta.res)
  
  for(t in 1:ncol(theta.res)){
    geweke.theta[t] <- geweke.diag(theta.res[,t], frac1=0.25, frac2=0.25)[[1]]
  }
  
  geweke.sigma2.g <- rep(NA, J)
  names(geweke.sigma2.g) <- paste0("x", 1:J)
  
  for(t in 1:J){
    geweke.sigma2.g[t] <- geweke.diag(sigma2.g.res[,t], frac1=0.25, frac2=0.25)[[1]]
  }
  
  geweke.sigma2.h <- rep(NA, J)
  names(geweke.sigma2.h) <- paste0("x", 1:J)
  
  for(t in 1:J){
    geweke.sigma2.h[t] <- geweke.diag(sigma2.h.res[,t], frac1=0.25, frac2=0.25)[[1]]
  }
  
  
  ########################################### return the result ################################################
  geweke.conv <- !(max(abs(geweke.theta))>4 | max(abs(geweke.sigma2.g))>4 | max(abs(geweke.sigma2.h))>4)
  if(!is.na(geweke.conv) & geweke.conv){
    return(list(success=TRUE, 
                spline_formula = splines$formula,
                spline_data = splines$data,
                spline_intKnots = splines$intKnots,
                numIntKnots = numIntKnots,
                xrange = xrange,
                geweke.conv=geweke.conv, theta=theta.res, sigma2_g=sigma2.g.res, sigma2_h=sigma2.h.res,				# MCMC samples 
                acpt_g=acpt_g/L, acpt_h=acpt_h/L,							# acceptance rates
                geweke_theta=geweke.theta, geweke_sigma2_g=geweke.sigma2.g, geweke_sigma2_h=geweke.sigma2.h))			# Geweke test statistics 
  } else {
    return(list(success=FALSE))			# Geweke test statistics 
    
  }
}
#
# function getCutoffLiu
#
# Author: Lara Maleyeff
#
# Function:       getCutoffLiu
# Description:    Internal function to find the alpha-quantile of the posterior treatment effect distribution for
#                 each individual (row). It uses the model that was fit in the most recent interim analysis (original data) to 
#                 find the alpha-quantiles of external data (used for either the adaptive enrichment steps or 
#                 accuracy calculations)
#                 Note: some parameters (candbinaryvars, pi) are not used here. They are included
#                 so the parameters of this function matches the other cutoff functions (Maleyeff and Park)
# Parameters:     data              A data frame with one individual per row and information on the continuous 
#                                   variables (candsplinevars) in the columns
#                 candpslinevars    Vector with names of continuous variables
#                 candbinaryvars    Not used
#                 trial_results     Results from most recent interim analysis, a data frame containing:
#                                   - xrange: the range of the original continuous variables
#                                   - intKnots: the internal knots used for the original continuous variables
#                                   - numIntKnots: the number of internal knots for each continuous variable
#                 alpha             alpha used for the effective subspace in paper
#                 pi  Not used
#
# Returns:        alpha-row quantile for each individuals to then be compared with 0. If
#                 the alpha-quantile of the treatment effect is > 0 then that individual's
#                 variable combination is in the effective subspace
parseMHLiu <- function(curdata, 
                       candsplinevars, 
                       candbinaryvars,
                       family,
                       B = 10000, 
                       burnin = 5000,
                       thin,
                       pi,
                       ...
                      ) {
 
  # First, fit the model using the MCMC procedure of Liu et al. (2022)
  unparsed_results = runMHLiu(curdata, 
                              candsplinevars,
                              family,
                              B,
                              burnin,
                              thin,
                              ...)
  
  # If the procedure was unsuccessful, stop
  if (!unparsed_results$success) {
    return(list(success = FALSE))
  } else {
    mcmc_coefs = unparsed_results$theta
    mod_mat = model.matrix(as.formula(unparsed_results$spline_formula),
                           data=unparsed_results$spline_data)
    indices_main = setdiff(1:ncol(mod_mat),grep("trt", colnames(mod_mat)))
    
    # trt_eff_posterior is a matrix with one individual per row, and B columns 
    # per individual, describing the complete posterior distribution of the treatment
    # effect for each individual, based on their tailoring variables
    trt_eff_posterior = mod_mat[,indices_main] %*% t(mcmc_coefs[,-indices_main])
    
    return(list(unparsed_results = unparsed_results,
                mcmc_coefs = mcmc_coefs,
                trt_eff_posterior = trt_eff_posterior,
                indices_main = indices_main,
                xrange = unparsed_results$xrange,
                intKnots = unparsed_results$spline_intKnots,
                numIntKnots = unparsed_results$numIntKnots,
                included_vars = candsplinevars,
                success = TRUE
    ))
  }
}

#
# function getCutoffLiu
#
# Author: Lara Maleyeff
#
# Function:       getCutoffLiu
# Description:    Internal function to find the alpha-quantile of the posterior treatment effect distribution for
#                 each individual (row). It uses the model that was fit in the most recent interim analysis (original data) to 
#                 find the alpha-quantiles of external data (used for either the adaptive enrichment steps or 
#                 accuracy calculations)
#                 Note: Parameter candbinaryvars is not used here. It is included
#                 so the parameters of this function match the other cutoff functions (Maleyeff and Park)
# Parameters:     data              A data frame with one individual per row and information on the continuous 
#                                   variables (candsplinevars) in the columns
#                 candpslinevars    Vector with names of continuous variables
#                 candbinaryvars    Not used
#                 trial_results     Results from most recent interim analysis, a data frame containing:
#                                   - xrange: the range of the original continuous variables
#                                   - intKnots: the internal knots used for the original continuous variables
#                                   - numIntKnots: the number of internal knots for each continuous variable
#                 alpha             alpha used for the effective subspace in paper
#
# Returns:        alpha-row quantile for each individuals to then be compared with 0. If
#                 the alpha-quantile of the treatment effect is > 0 then that individual's
#                 variable combination is in the effective subspace
getCutoffLiu <- function(data,
                         candsplinevars,
                         candbinaryvars,
                         trial_results,
                         alpha
                        ) {
  mod_mat = get_trt_model_matrix(data[,candsplinevars, drop=F], 
                                         trial_results$xrange, 
                                         trial_results$intKnots,
                                         trial_results$numIntKnots)
  predeff = as.matrix(mod_mat) %*% t(trial_results$mcmc_coefs[,-trial_results$indices_main])
  return(rowQuantiles(predeff,probs = c(alpha)))
}

