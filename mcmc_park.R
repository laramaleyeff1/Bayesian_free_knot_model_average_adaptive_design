runMHPark <- function(data, 
                      startvalue, 
                      candvars,
                      mod_mat, 
                      B,
                      burnin,
                      thin,
                      family = gaussian(link = "identity"),
                      mu = 200,
                      tau = 0.1
                      ) {
  library(coda)
  iterations = B*thin + burnin
  n_params = length(startvalue) - 1
  n_vars = (n_params-2)/2
  accept = array(dim = c(iterations,n_params))
  sd = array(dim = c(iterations,1))
  chain = array(dim = c(iterations+1,n_params))
  vars_prop = array(dim = c(iterations+1,n_params - 1))
  p = array(dim = c(iterations+1,n_params - 1))
  chain[1,] = startvalue[-1]
  colnames(chain) = names(startvalue[-1])
  colnames(accept) = names(startvalue[-1])
  
  sd[1] = startvalue[1]
  

  # Always include the intercept
  vars_prop[1,] <- rep(1,n_params-1)
  p[1,] <- rep(0.5,n_params-1)
  tau_0 = 10
  a_0 = 0.01
  b_0 = 0.01

  for (i in 1:iterations){
    
    prob_main = p[i,1:(n_vars+1)]*dnorm(chain[i,2:(n_vars+2)],0,mu*tau)
    prob_main = prob_main/(prob_main +
                             (1-p[i,1:(n_vars+1)])*dnorm(chain[i,2:(n_vars+2)],0,tau))

    prob_main[n_vars + 1] = 1
    vars_prop[i+1,1:(n_vars+1)] = rbinom(n_vars + 1,1,prob_main)

    # probs_incl = p[i,]*dnorm(chain[i,-c(1)],0,mu*tau)
    # probs_incl = probs_incl/(probs_incl + (1-p[i,])*dnorm(chain[i,-c(1)],0,tau))
    # vars_prop[i+1,] = rbinom(n_params - 1 ,1, probs_incl)
    
    prob_inter = sapply((n_vars+2):(n_params-1), function(j) {
      return(prob_main[j-n_vars-1]*prob_main[n_vars + 1]*min(prob_main[j-n_vars-1],prob_main[n_vars + 1]))
    })

    vars_prop[i+1,(n_vars+2):(n_params-1)] = rbinom(n_vars,1,prob_inter)

    p[i+1,] = rbeta(n_params - 1, vars_prop[i+1,] + 1, 2 - vars_prop[i+1,])

    indicators = c(1,vars_prop[i+1,]) 
    # Set new chain to current value 
    
    # Now the indicators are given! If indicator == 1, then the 
    # coefficient value is 0!
    chain[i+1,] = chain[i,] * indicators
    for (j in 1:length(chain[i+1,])) {
      # chain[i+1,] may have changed in previous iterations, 
      # so we want to reset the proposal to reflect its current value
      if (indicators[j] == 1) {
        sigma2_beta_j <- ifelse((j == 1 | j == n_vars + 1), tau_0, tau*mu)
        V_beta_j <- 1 / (t(mod_mat[, j]) %*% mod_mat[, j] / sd[i] + 1 / sigma2_beta_j)
        beta_j_hat <- V_beta_j * t(mod_mat[, j]) %*% (data$Y - mod_mat[, -j] %*% chain[i+1, -j]) / sd[i]
        chain[i+1, j] <- rnorm(1, beta_j_hat, sqrt(V_beta_j))
      } 
    }
   
    
    resids = data$Y - chain[i+1,] %*% t(mod_mat)
    sd[i+1] <- 1/rgamma(1, shape=nrow(data)/2+a_0, rate=sum(resids^2)/2+b_0)
  

  }
  colnames(vars_prop) = names(startvalue[-c(1,2)])
  
  vars_prop_res = vars_prop[-(1:burnin),]
  names(vars_prop_res) = names(vars_prop)
  vars_prop_res = vars_prop_res[seq(1,nrow(vars_prop_res),thin),]
  
  chain_res = chain[-(1:burnin),]
  names(chain_res) = names(chain)
  chain_res = chain_res[seq(1,nrow(chain_res),thin),]
  
  p_res = p[-(1:burnin),]
  names(p_res) = names(p)
  p_res = p_res[seq(1,nrow(p_res),thin),]
  
  sd_res = sd[-(1:burnin)]
  sd_res = sd_res[seq(1,length(sd_res),thin)]
  
  indices = (2+length(candvars)):(2*length(candvars) + 2)
  trt_eff_posterior = mod_mat[,-indices] %*% t(chain_res[,indices])
  
  
  geweke.trt_eff_posterior <- rep(NA, nrow(data))
  names(geweke.trt_eff_posterior) = 1:nrow(data)
  for (t in 1:nrow(data)) {
    geweke.trt_eff_posterior[t] <- geweke.diag(trt_eff_posterior[t,], frac1=0.25, frac2=0.25)[[1]]
  }
  
  geweke.trt <- geweke.diag(colMeans(trt_eff_posterior), frac1=0.25, frac2=0.25)[[1]]
  geweke.sd_res <- geweke.diag(sd_res, frac1=0.25, frac2=0.25)[[1]]
  geweke.conv <- !(max(abs(geweke.sd_res))>4 | max(abs(geweke.trt))>4)
  
  if (!is.na(geweke.conv) & geweke.conv) {
    return(list(success = TRUE,
                geweke.conv = geweke.conv,
                chain_res = chain_res, 
                trt_eff_posterior = trt_eff_posterior,
                indices = indices,
                p_res = p_res, 
                vars_prop_res = vars_prop_res, 
                sd_res = sd_res, 
                accept = accept,
                geweke.sd_res = geweke.sd_res,
                geweke.trt = geweke.trt,
                geweke.trt_eff_posterior = geweke.trt_eff_posterior
                ))
  } else {
    return(list(success = FALSE))
  }
  
}

parseMHPark <- function(curdata, 
                        candsplinevars, 
                        candbinaryvars,
                        family,
                        B = 10000, 
                        burnin = 5000,
                        thin,
                        pi,
                        ...) {
  candvars = c(candsplinevars, candbinaryvars)
  formula_all = paste0("Y ~ ", paste(c(candvars,paste0(candvars, "*trt")),collapse="+"))
  if (family$family == "gaussian") {
    mod_ols = lm(as.formula(formula_all), data=curdata)
    startvalue = c(sigma(mod_ols),
                   coef(mod_ols)
    )
  } else {
    mod_ols = glm(as.formula(formula_all), data=curdata, family = family)
    startvalue = coef(mod_ols)
  }
  
  mod_mat = model.matrix(mod_ols)
  unparsed_results = runMHPark(curdata, 
                               startvalue,
                               candvars,
                               mod_mat,
                               B, 
                               burnin,
                               thin,
                               family,
                               ...
                               )
  if (!unparsed_results$success) {
    return(list(success = FALSE))
  } else {
    prop_incl = colMeans(unparsed_results$vars_prop_res > 0)
    names(prop_incl) = colnames(unparsed_results$vars_prop_res)
    included_vars = names(which(prop_incl > pi))
    included_vars = included_vars[grep(":trt",included_vars)]
    included_vars = sub(":trt", "", included_vars)
    return(list(unparsed_results = unparsed_results,
                mcmc_coefs = unparsed_results$chain,
                trt_eff_posterior = unparsed_results$trt_eff_posterior,
                prop_incl = prop_incl,
                indices = unparsed_results$indices,
                included_vars = included_vars,
                success = TRUE))
  }
}

getCutoffPark <- function(data,
                           candsplinevars,
                           candbinaryvars,
                           trial_results,
                           alpha) {
  candvars = c(candsplinevars, candbinaryvars)
  formula_trt = paste0("~ ", paste(candvars,collapse="+"))
  mod_mat = model.matrix(as.formula(formula_trt),data=data)
  predeff = as.matrix(mod_mat) %*% t(trial_results$mcmc_coefs[,trial_results$indices])
  return(rowQuantiles(predeff,probs = c(alpha)))
}