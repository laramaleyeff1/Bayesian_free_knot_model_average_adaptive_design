
expit <- function(x) {
  exp(x)/(1+exp(x))
}

sample_ <- function(x) {
  if (length(x) == 1) {
    return(x)
  } else {
    return(sample(x,1))
  }
}
computeFittedValues <- function(trt,
                                inter_trt_param_,
                             binary_param_,
                             binary_mod_mat_,
                             spline_param_minus_j, 
                             spline_mat_minus_j) {
  
  binary_eff = 0
  if (is.list(binary_mod_mat_)) {
    binary_eff = do.call(cbind,binary_mod_mat_) %*% binary_param_
  }
  
  spline_eff = 0
  if (is.list(spline_mat_minus_j)) {
    spline_eff = do.call(cbind, spline_mat_minus_j) %*% unlist(spline_param_minus_j)
  }
  
  return(cbind(1,trt) %*% inter_trt_param_ + binary_eff + spline_eff)
}

splinesFitted <- function(trt, spline, mat) {
  return(trt + mat %*% spline)
}

paste_ <- function(set, suffix) {
  if (length(set) == 0) {
    return(NULL)
  } else {
    return(paste0(set,suffix))
  }
  
}

logLikelihoodCustom <- function(Y,
                                trt,
                                inter_trt_param_,
                                binary_param_,
                                binary_mod_mat_,
                                spline_param_, 
                                spline_mat_,
                                sd,
                                family) {
  pred = computeFittedValues(trt,
                             inter_trt_param_,
                             binary_param_,
                             binary_mod_mat_,
                             spline_param_, 
                             spline_mat_)
  if (family$family == "binomial") {
    singlelikelihoods = dbinom(Y, length(Y), prob = family$linkinv(pred), log = T)
  } else if (family$family == "gaussian") {
    singlelikelihoods = dnorm(Y, mean = pred, sd = sd, log = T)
  }
  sumll = sum(singlelikelihoods)
  return(sumll)   
}

moveKnot <- function(Y, 
                     data,
                     w,
                     varcur,
                     j,
                     knotscand_x, 
                     knotscur_x,
                     knotscur_idx_x,
                     binary_param,
                     binary_mod_mat,
                     spline_mod_mat, 
                     spline_mod_mat_raw,
                     spline_ols_param_x, 
                     inter_trt_param_i,
                     spline_param_i, 
                     sd_i,
                     family) {
  spline_mod_mat_x = spline_mod_mat[[varcur]]
  spline_mod_mat_raw_x = spline_mod_mat_raw[[varcur]]
  
  if (length(knotscur_idx_x) > 0) {
    knot_tomove = knotscand_x[sample_(knotscur_idx_x)]
    knots_avail = setdiff(knotscand_x, knotscur_x)
    window_propose = knots_avail[knots_avail<=(knot_tomove+w) &
                                   knots_avail>=(knot_tomove-w)]
  } else {
    window_propose = NULL
  }
  
  if (length(window_propose) > 0) {
    knot_newvalue = sample_(window_propose)
    
    knots_propose = sort(setdiff(union(knotscur_x, knot_newvalue),knot_tomove))
    knotscur_idx_propose = sort(setdiff(union(which(knotscand_x == knot_newvalue),
                                              knotscur_idx_x),
                                        which(knotscand_x == knot_tomove)))
    
    knots_avail_reverse = setdiff(knotscand_x, knots_propose)
    window_reverse = knots_avail_reverse[knots_avail_reverse<=(knot_newvalue+w) &
                                           knots_avail_reverse>=(knot_newvalue-w)]
    
    
    offset = computeFittedValues(data$trt,
                                 inter_trt_param_i,
                              binary_param,
                              binary_mod_mat,
                              spline_param_i[-j],
                              spline_mod_mat[-j])
    
    new_spline = bs(data[[sub("_[^_]+$", "", varcur)]],
                    degree = 3,
                    knots = knots_propose)
    
    if (length(grep("inter",varcur))>0) {
      new_spline_raw = new_spline
      new_spline = new_spline*data$trt
    }
    
    mod_propose <- glm(Y ~ -1 + new_spline, 
                       data=data, 
                       family=family,
                       offset = offset)
    
    spline_mod_mat_propose = spline_mod_mat
    spline_mod_mat_propose[[varcur]] = new_spline

    log_prob = logLikelihoodCustom(Y,
                                   data$trt,
                                   inter_trt_param_i,
                                   binary_param,
                                   binary_mod_mat,
                                   spline_param_i,
                                   spline_mod_mat_propose,
                                   sd_i,
                                   family) -
      logLikelihoodCustom(Y,
                          data$trt,
                          inter_trt_param_i,
                          binary_param,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sd_i,
                          family) +
      log(length(window_propose)) -
      log(length(window_reverse))
    
    
    gamma <- runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of which knots
      knotscur_idx_x = knotscur_idx_propose
      knotscur_x = knots_propose
      # keep track of OLS model coefs
      spline_ols_param_x = coef(mod_propose)
      # keep track of model matrices
      spline_mod_mat_x = new_spline
      if (length(grep("inter",varcur))>0) {
        spline_mod_mat_raw_x = new_spline_raw
      }
      accept = 1
    } else {
      accept = 0
    }
  } else {
    accept = 0
  }
  
  return(list(knotscur_idx_x = knotscur_idx_x,
              knotscur_x = knotscur_x,
              spline_ols_param_x = spline_ols_param_x,
              spline_mod_mat_x = spline_mod_mat_x,
              spline_mod_mat_raw_x = spline_mod_mat_raw_x,
              accept = accept))
}

addKnot <- function(Y, 
                    data,
                    v,
                    sigma_v,
                    sigma_B,
                    lambda_2,
                    varcur,
                    j,
                    k_ij,
                    k_max,
                    knotscand_x, 
                    knotscur_x,
                    knotscur_idx_x,
                    binary_param,
                    binary_mod_mat,
                    spline_mod_mat, 
                    spline_mod_mat_raw,
                    spline_ols_param_x, 
                    inter_trt_param_i,
                    spline_param_i, 
                    sd_i,
                    family) {
  spline_mod_mat_x = spline_mod_mat[[varcur]]
  spline_mod_mat_raw_x = spline_mod_mat_raw[[varcur]]
  if (k_ij < k_max) {
    index_propose = sample_(setdiff(1:k_max,knotscur_idx_x))
    knot_propose = knotscand_x[index_propose]
    knots_propose = sort(union(knotscur_x, knot_propose))
    k_interval_propose = which(knots_propose == knot_propose)
    knotscur_idx_propose = sort(union(index_propose, knotscur_idx_x))
    k_cur_propose = k_ij + 1
    
    
    offset = computeFittedValues(data$trt,
                                 inter_trt_param_i,
                              binary_param,
                              binary_mod_mat,
                              spline_param_i[-j],
                              spline_mod_mat[-j])
    
    new_spline = bs(data[[sub("_[^_]+$", "", varcur)]],
                    degree = 3,
                    knots = knots_propose)
    
    if (length(grep("inter",varcur))>0) {
      new_spline_raw = new_spline
      new_spline = new_spline*data$trt
    }
    
    mod_propose <- glm(Y ~ -1 + new_spline, 
                       data=data, 
                       family=family,
                       offset = offset)
    
    spline_propose = coef(mod_propose)
    spline_propose[-(k_interval_propose+1)] = spline_propose[-(k_interval_propose+1)] +
      spline_param_i[[j]] -
      spline_ols_param_x
    
    spline_propose[(k_interval_propose+1)] = spline_propose[(k_interval_propose+1)] + v
    
    splines_propose = spline_param_i
    splines_propose[[varcur]] = spline_propose
    
    spline_mod_mat_propose = spline_mod_mat
    spline_mod_mat_propose[[varcur]] = new_spline
    
    log_prob = logLikelihoodCustom(Y,
                                   data$trt,
                                   inter_trt_param_i,
                                   binary_param,
                                   binary_mod_mat,
                                   splines_propose,
                                   spline_mod_mat_propose,
                                   sd_i,
                                   family) -
      logLikelihoodCustom(Y,
                          data$trt,
                          inter_trt_param_i,
                          binary_param,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sd_i,
                          family) +
      log((lambda_2*sigma_v)/(k_cur_propose*sigma_B)) -
      (v^2/2*sigma_v^2) + (1/(2*sigma_B^2))*(t(spline_param_i[[varcur]]) %*% spline_param_i[[varcur]] -
                                               t(spline_propose) %*% (spline_propose))
    
    gamma <- runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of spline coeffients
      spline_param_x = spline_propose
      # keep track of number of knots
      k = k_cur_propose
      # keep track of which knots
      knotscur_idx_x = knotscur_idx_propose
      # keep track of OLS model coeffieints
      spline_ols_param_x = coef(mod_propose)
      # keep track of model matrices
      spline_mod_mat_x = new_spline
      if (length(grep("inter",varcur))>0) {
        spline_mod_mat_raw_x = new_spline_raw
      }
      accept = 1
      
    } else {
      spline_param_x = spline_param_i[[varcur]]
      k = k_ij
      accept = 0
    }
  } else {
    spline_param_x = spline_param_i[[varcur]]
    k = k_ij
    accept = 0
  }
  
  return(list(spline_param_x = spline_param_x,
              k = k,
              knotscur_idx_x = knotscur_idx_x,
              knotscur_x = knotscur_x,
              spline_ols_param_x = spline_ols_param_x,
              spline_mod_mat_x = spline_mod_mat_x,
              spline_mod_mat_raw_x = spline_mod_mat_raw_x,
              accept = accept))
  
}

removeKnot <- function(Y, 
                       data,
                    v,
                    sigma_v,
                    sigma_B,
                    lambda_2,
                    varcur,
                    j,
                    k_ij,
                    k_max,
                    knotscand_x, 
                    knotscur_x,
                    knotscur_idx_x,
                    binary_param,
                    binary_mod_mat,
                    spline_mod_mat, 
                    spline_mod_mat_raw,
                    spline_ols_param_x, 
                    inter_trt_param_i,
                    spline_param_i, 
                    sd_i,
                    family) {
  spline_mod_mat_x = spline_mod_mat[[varcur]]
  spline_mod_mat_raw_x = spline_mod_mat_raw[[varcur]]
  
  if (k_ij > 0) {
    index_propose = sample_(knotscur_idx_x)
    knot_propose = knotscand_x[index_propose]
    knots_propose = setdiff(knotscur_x, knot_propose)
    
    k_interval_propose = which(knotscur_x == knot_propose)
    knotscur_idx_propose = setdiff(knotscur_idx_x, index_propose)
    k_cur_propose = k_ij - 1
    
    offset = computeFittedValues(data$trt,
                                 inter_trt_param_i,
                              binary_param,
                              binary_mod_mat, 
                              spline_param_i[-j],
                              spline_mod_mat[-j])
    
    new_spline = bs(data[[sub("_[^_]+$", "", varcur)]],
                    degree = 3,
                    knots = knots_propose)
    
    if (length(grep("inter",varcur))>0) {
      new_spline_raw = new_spline
      new_spline = new_spline*data$trt
    }
    
    mod_propose <- glm(Y ~ -1 + new_spline, 
                       data=data, 
                       family=family,
                       offset = offset)    
    spline_propose = coef(mod_propose)
    spline_propose = spline_propose +
      spline_param_i[[j]][-(k_interval_propose+1)] -
      spline_ols_param_x[-(k_interval_propose+1)]
    
    splines_propose = spline_param_i
    splines_propose[[varcur]] = spline_propose
    
    spline_mod_mat_propose = spline_mod_mat
    spline_mod_mat_propose[[varcur]] = new_spline

    log_prob = logLikelihoodCustom(Y,
                                   data$trt,
                                   inter_trt_param_i,
                                   binary_param,
                                   binary_mod_mat,
                                   splines_propose,
                                   spline_mod_mat_propose,
                                   sd_i,
                                   family) -
      logLikelihoodCustom(Y,
                          data$trt,
                          inter_trt_param_i,
                          binary_param,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sd_i,
                          family) +
      log(sigma_B*k_ij/(lambda_2*sigma_v)) -
      (v^2/2*sigma_v^2) + (1/(2*sigma_B^2)) *
      (t(spline_param_i[[varcur]]) %*% spline_param_i[[varcur]] -
         t(spline_propose) %*% (spline_propose))

    gamma <- runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of spline coefs
      spline_param_x = spline_propose
      # keep track of number of knots
      k = k_cur_propose
      # keep track of which knots
      knotscur_idx_x = knotscur_idx_propose
      # keep track of OLS model coefs
      spline_ols_param_x = coef(mod_propose)
      # keep track of model matrices
      spline_mod_mat_x = new_spline
      if (length(grep("inter",varcur))>0) {
        spline_mod_mat_raw_x = new_spline_raw
      }
      accept = 1
      
    } else {
      spline_param_x = spline_param_i[[varcur]]
      k = k_ij
      accept = 0
    }
  } else {
    spline_param_x = spline_param_i[[varcur]]
    k = k_ij
    accept = 0
  }
  
  return(list(spline_param_x = spline_param_x,
              k = k,
              knotscur_idx_x = knotscur_idx_x,
              knotscur_x = knotscur_x,
              spline_ols_param_x = spline_ols_param_x,
              spline_mod_mat_x = spline_mod_mat_x,
              spline_mod_mat_raw_x = spline_mod_mat_raw_x,
              accept = accept))
  
}

addVar <- function(Y, 
                   trt,
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
                   inter_trt_param_i, 
                   binary_param_i,
                   binary_mod_mat,
                   spline_param_i,
                   spline_mod_mat,
                   sd_i,
                   family) {
  if (length(eligibletoadd) > 0) {
    to_add = sample_(eligibletoadd)
    spline_param_propose = spline_param_i
    binary_param_propose = binary_param_i
    inter_trt_param_propose = inter_trt_param_i
    
    curbinary_propose = curbinary_ext
    curspline_propose = curspline_ext
    
    if (to_add %in% candbinaryvars_ext) {
      curbinary_propose =  c(curbinary_ext, to_add)
      index_add = which(to_add == curbinary_propose)

      dim_x = 1
      v <- rnorm(dim_x,0,sigma_v)
      binary_param_propose[to_add] = v
      if (length(curspline_ext) == 0) {
        mod_propose = glm(Y ~ trt + do.call(cbind, binary_mod_mat[curbinary_propose]), family=family)

        if (length(curbinary_ext) == 0) {
          mod_prev = glm(Y ~ trt, family=family)
        } else {
          mod_prev = glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]), family=family)
          binary_param_propose[curbinary_ext] = binary_param_i[curbinary_ext] -
            coef(mod_prev)[-c(1,2)] +
            coef(mod_propose)[setdiff(3:(length(curbinary_propose)+2),2+index_add)]
        }
      } else {
        mod_propose = glm(Y ~ trt + do.call(cbind, binary_mod_mat[curbinary_propose]) +
                           do.call(cbind, spline_mod_mat[curspline_ext]), family=family)

        if (length(curbinary_ext) == 0) {
          mod_prev = glm(Y ~ trt + do.call(cbind, spline_mod_mat[curspline_ext]), family=family)
        } else {
          mod_prev = glm(Y ~ trt + do.call(cbind, binary_mod_mat[curbinary_ext]) +
                             do.call(cbind, spline_mod_mat[curspline_ext]), family=family)
          binary_param_propose[curbinary_ext] = binary_param_i[curbinary_ext] -
            coef(mod_prev)[3:(length(curbinary_ext)+2)] +
            coef(mod_propose)[setdiff(3:(length(curbinary_ext)+2),2+index_add)]
        }

        spline_ols_prev = convertSplineCoefs(coef(mod_prev)[-c(1:(length(curbinary_ext)+2))],
                                             spline_mod_mat,
                                             curspline_ext)
        spline_ols_propose = convertSplineCoefs(coef(mod_propose)[-c(1:(length(curbinary_propose)+2))],
                                                spline_mod_mat,
                                                curspline_ext)

        for (m in 1:length(curspline_ext)) {
          spline_param_propose[[curspline_ext[m]]] =
            spline_param_i[[curspline_ext[m]]] -
            spline_ols_prev[[curspline_ext[m]]] + spline_ols_propose[[curspline_ext[m]]]
        }
      }

      binary_param_propose[to_add] = v + coef(mod_propose)[2+index_add]
    } else {
      curspline_propose =  c(curspline_ext, to_add)
      dim_x = length(spline_param_i[[to_add]])
      index_add = which(to_add == curspline_propose)
      
      v <- rnorm(dim_x,0,sigma_v)
      if (length(curbinary_ext) == 0) {
        mod_propose = glm(Y ~ trt + do.call(cbind, spline_mod_mat[curspline_propose]), family=family)

        spline_ols_propose = convertSplineCoefs(coef(mod_propose)[-c(1,2)],
                                             spline_mod_mat,
                                             curspline_propose)

        if (length(curspline_ext) == 0) {
          mod_prev = glm(Y ~ trt, family=family)
        } else {
          mod_prev = glm(Y ~ trt + do.call(cbind, spline_mod_mat[curspline_ext]), family=family)
          spline_ols_prev = convertSplineCoefs(coef(mod_prev)[-c(1,2)],
                                                  spline_mod_mat,
                                                  curspline_ext)

          for (m in 1:length(curspline_ext)) {
            spline_param_propose[[curspline_ext[m]]] =
              spline_param_i[[curspline_ext[m]]] -
              spline_ols_prev[[curspline_ext[m]]] + spline_ols_propose[[curspline_ext[m]]]
          }
        }

      } else {
        mod_propose = glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]) +
                           do.call(cbind, spline_mod_mat[curspline_propose]), family=family)

        spline_ols_propose = convertSplineCoefs(coef(mod_propose)[-c(1:(length(curbinary_ext)+2))],
                                             spline_mod_mat,
                                             curspline_propose)

        if (length(curspline_ext) == 0) {
          mod_prev = glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]), family=family)
        } else {
          mod_prev = glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]) +
                          do.call(cbind, spline_mod_mat[curspline_ext]), family=family)
          spline_ols_prev = convertSplineCoefs(coef(mod_prev)[-c(1:(length(curbinary_ext)+2))],
                                                  spline_mod_mat,
                                                  curspline_ext)

          for (m in 1:length(curspline_ext)) {
            spline_param_propose[[curspline_ext[m]]] =
              spline_param_i[[curspline_ext[m]]] -
              spline_ols_prev[[curspline_ext[m]]] + spline_ols_propose[[curspline_ext[m]]]
          }
        }
        binary_param_propose[curbinary_ext] = binary_param_i[curbinary_ext] -
          coef(mod_prev)[3:(length(curbinary_ext)+2)] +
          coef(mod_propose)[3:(length(curbinary_ext)+2)]
      }
      spline_param_propose[[to_add]] = v + spline_ols_propose[[to_add]]
      #spline_param_propose[[to_add]] = v
    }
    
    inter_trt_param_propose = inter_trt_param_i - coef(mod_prev)[c(1,2)] + coef(mod_propose)[c(1,2)]
    #inter_trt_param_propose = inter_trt_param_i
    
    log_prob = logLikelihoodCustom(Y,
                                   trt,
                                   inter_trt_param_propose,
                                   binary_param_propose,
                                   binary_mod_mat,
                                   spline_param_propose,
                                   spline_mod_mat,
                                   sd_i,
                                   family) -
      logLikelihoodCustom(Y,
                          trt,
                          inter_trt_param_i,
                          binary_param_i,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sd_i,
                          family) +
      sum(dnorm(c(inter_trt_param_propose, binary_param_propose, unlist(spline_param_propose)), mean = 0, sd = sigma_B, log = T)) -
      sum(dnorm(c(inter_trt_param_i, binary_param_i, unlist(spline_param_i)), mean = 0,
                sd = sigma_B, log = T)) +
      log(lambda_1/(n_cand_vars-n_cur_vars)) +
      log(length(eligibletoadd)) - log(length(eligibletoremove)+1) -
      sum(dnorm(v,0,sigma_v,log=T))
    
    curmain_propose = curmain
    curinter_propose = curinter
    if (length(grep("main", to_add)>0)) {
      curmain_propose = c(sub("_[^_]+$", "", to_add), curmain)
    } else {
      curinter_propose = c(sub("_[^_]+$", "", to_add), curinter)
    }
    
    gamma <- runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of spline coeffients
      spline_param_i = spline_param_propose
      binary_param_i = binary_param_propose
      inter_trt_param_i = inter_trt_param_propose
      curmain = curmain_propose
      curinter = curinter_propose
      curvars_ext = c(curbinary_propose, curspline_propose)
      curbinary_ext = curbinary_propose
      curspline_ext = curspline_propose
      accept = 1
    } else {
      accept = 0
    }
    
  } else {
    accept = 0
  }
  return(list(spline_param_i = spline_param_i,
              inter_trt_param_i = inter_trt_param_i,
              binary_param_i = binary_param_i,
              curmain = curmain,
              curinter = curinter,
              curvars_ext = curvars_ext,
              curbinary_ext = curbinary_ext,
              curspline_ext = curspline_ext,
              accept = accept))
}

convertSplineCoefs <- function(spline_coefs, 
                               spline_mat,
                               vars) {
  ncoef_perx = unlist(lapply(spline_mat[vars],ncol))

  coef_prop = split(spline_coefs,
                    rep(1:(length(ncoef_perx)),ncoef_perx))
  spline_ols = list()
  for (m in 1:length(vars)) {
    spline_ols[[vars[m]]] = unlist(coef_prop[m])
  }
  return(spline_ols)
}

removeVar <- function(Y, 
                      trt,
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
                   inter_trt_param_i, 
                   binary_param_i,
                   binary_mod_mat,
                   spline_param_i,
                   spline_mod_mat,
                   sd_i,
                   family) {
  if (length(eligibletoremove) > 0) {
    to_remove = sample_(eligibletoremove)
    spline_param_propose = spline_param_i
    
    binary_param_propose = binary_param_i
    inter_trt_param_propose = inter_trt_param_i
    curbinary_propose = curbinary_ext
    curspline_propose = curspline_ext
    
    if (to_remove %in% candbinaryvars_ext) {
      curbinary_propose =  setdiff(curbinary_ext, to_remove)
      dim_x = 1
      v <- rnorm(dim_x,0,sigma_v)
      binary_param_propose[to_remove] = 0
      
      if (length(curspline_ext) == 0) {
        binary_mod_mat_prev = binary_mod_mat[curbinary_ext]
        index_remove = which(to_remove == curbinary_ext)
        mod_prev = glm(Y ~ trt + do.call(cbind, binary_mod_mat_prev),family=family)

        if (length(curbinary_propose) == 0) {
          mod_propose = glm(Y ~ trt,family=family)
        } else {
          mod_propose = glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_propose]),family=family)
          binary_param_propose[curbinary_propose] = binary_param_i[curbinary_propose] -
            coef(mod_prev)[setdiff(3:(length(curbinary_ext)+2),2+index_remove)] +
            coef(mod_propose)[-c(1,2)]
        }
      } else {
        mod_prev = glm(Y ~ trt + do.call(cbind, binary_mod_mat[curbinary_ext]) + do.call(cbind, spline_mod_mat[curspline_ext]),family=family)

        if (length(curbinary_propose) == 0) {
          mod_propose = glm(Y ~ trt + do.call(cbind, spline_mod_mat[curspline_ext]),family=family)
        } else {
          mod_propose = glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_propose]) +
                             do.call(cbind, spline_mod_mat[curspline_ext]),family=family)
          index_remove = which(to_remove == curbinary_ext)
          binary_param_propose[curbinary_propose] = binary_param_i[curbinary_propose] -
            coef(mod_prev)[setdiff(3:(length(curbinary_ext)+2),2+index_remove)] +
            coef(mod_propose)[3:(length(curbinary_propose)+2)]
        }

        spline_ols_prev = convertSplineCoefs(coef(mod_prev)[-c(1:(length(curbinary_ext)+2))],
                                         spline_mod_mat,
                                         curspline_ext)
        spline_ols_propose = convertSplineCoefs(coef(mod_propose)[-c(1:(length(curbinary_propose)+2))],
                                         spline_mod_mat,
                                         curspline_ext)

        for (m in 1:length(curspline_ext)) {
          spline_param_propose[[curspline_ext[m]]] =
            spline_param_i[[curspline_ext[m]]] -
            spline_ols_prev[[curspline_ext[m]]] + spline_ols_propose[[curspline_ext[m]]]
        }
      }

      binary_param_propose[to_remove] = 0
    } else {
      dim_x = length(spline_param_i[[to_remove]])
      curspline_propose =  setdiff(curspline_ext, to_remove)
      
      v <- rnorm(dim_x,0,sigma_v)
      
      if (length(curbinary_ext) == 0) {
        mod_mat_prev = do.call(cbind, spline_mod_mat[curspline_ext])
        mod_prev = glm(Y ~ trt + mod_mat_prev, family=family)

        spline_ols_prev = convertSplineCoefs(coef(mod_prev)[-c(1,2)],
                                             spline_mod_mat,
                                             curspline_ext)

        if (length(curspline_propose) == 0) {
          mod_propose = glm(Y ~ trt, family=family)
        } else {
          mod_mat_propose = do.call(cbind, spline_mod_mat[curspline_propose])
          mod_propose = glm(Y ~ trt + mod_mat_propose, family=family)
          spline_ols_propose = convertSplineCoefs(coef(mod_propose)[-c(1,2)],
                                                  spline_mod_mat,
                                                  curspline_propose)

          for (m in 1:length(curspline_propose)) {
            spline_param_propose[[curspline_propose[m]]] =
              spline_param_i[[curspline_propose[m]]] -
              spline_ols_prev[[curspline_propose[m]]] + spline_ols_propose[[curspline_propose[m]]]
          }
        }

      } else {
        mod_mat_prev = do.call(cbind, spline_mod_mat[curspline_ext])

        mod_prev = glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]) + mod_mat_prev, family=family)

        spline_ols_prev = convertSplineCoefs(coef(mod_prev)[-c(1:(length(curbinary_ext)+2))],
                                             spline_mod_mat,
                                             curspline_ext)
        
        if (length(curspline_propose) == 0) {
          mod_propose = glm(Y ~ trt + do.call(cbind,binary_mod_mat[curbinary_ext]), family=family)
        } else {
          mod_propose = glm(Y ~ trt +  do.call(cbind,binary_mod_mat[curbinary_ext]) +
                             do.call(cbind, spline_mod_mat[curspline_propose]), family=family)
          spline_ols_propose = convertSplineCoefs(coef(mod_propose)[-c(1:(length(curbinary_ext)+2))],
                                                  spline_mod_mat,
                                                  curspline_propose)
          for (m in 1:length(curspline_propose)) {
            spline_param_propose[[curspline_propose[m]]] =
              spline_param_i[[curspline_propose[m]]] -
              spline_ols_prev[[curspline_propose[m]]] + spline_ols_propose[[curspline_propose[m]]]
          }

          
        }
        binary_param_propose[curbinary_ext] = binary_param_i[curbinary_ext] -
          coef(mod_prev)[3:(length(curbinary_ext)+2)] +
          coef(mod_propose)[3:(length(curbinary_ext)+2)]
      }
      spline_param_propose[[to_remove]] = rep(0,dim_x)
    }
    
    inter_trt_param_propose = inter_trt_param_i - coef(mod_prev)[c(1,2)] + coef(mod_propose)[c(1,2)]

    log_prob = logLikelihoodCustom(Y,
                                   trt,
                                   inter_trt_param_propose,
                                   binary_param_propose,
                                   binary_mod_mat,
                                   spline_param_propose,
                                   spline_mod_mat,
                                   sd_i,
                                   family) -
      logLikelihoodCustom(Y,
                          trt,
                          inter_trt_param_i,
                          binary_param_i,
                          binary_mod_mat,
                          spline_param_i,
                          spline_mod_mat,
                          sd_i,
                          family) +
      sum(dnorm(c(inter_trt_param_propose, binary_param_propose, unlist(spline_param_propose)), mean = 0, sd = sigma_B, log = T)) -
      sum(dnorm(c(inter_trt_param_i, binary_param_i, unlist(spline_param_i)), mean = 0,
                sd = sigma_B, log = T)) +
      log((n_cand_vars-n_cur_vars + 1)/lambda_1) +
      log(length(eligibletoremove)) - log(length(eligibletoadd)+1) +
      sum(dnorm(v,0,sigma_v,log=T))

    #print(log_prob)
    curmain_propose = curmain
    curinter_propose = curinter
    if (length(grep("main", to_remove)>0)) {
      curmain_propose = setdiff(curmain, sub("_[^_]+$", "", to_remove))
    } else {
      curinter_propose =  setdiff(curinter, sub("_[^_]+$", "", to_remove))
    }
    
    gamma <- runif(1,0,1)
    if (gamma < min(1,exp(log_prob))) {
      # keep track of spline coeffients
      spline_param_i = spline_param_propose
      binary_param_i = binary_param_propose
      inter_trt_param_i = inter_trt_param_propose
      curmain = curmain_propose
      curinter = curinter_propose
      curvars_ext = c(curbinary_propose, curspline_propose)
      curbinary_ext = curbinary_propose
      curspline_ext = curspline_propose
      accept = 1
    } else {
      accept = 0
    }
    
  } else {
    accept = 0
  }
  return(list(spline_param_i = spline_param_i,
              inter_trt_param_i = inter_trt_param_i,
              binary_param_i = binary_param_i,
              curmain = curmain,
              curinter = curinter,
              curvars_ext = curvars_ext,
              curbinary_ext = curbinary_ext,
              curspline_ext = curspline_ext,
              accept = accept))
}

