######################################################################################################
# Set of helper functions for methods of Liu et al. (2022), described in Maleyeff et al. (2024)      
#                           Adapted from Liu et al. (2022), by Lara Maleyeff
#         Original code: https://github.com/YushaLiu/Continuous_marker_adaptive_trial_design
#                                   Contact: laramaleyeff@gmail.com                                  
#                                       Last updated: April 2024                                       
######################################################################################################

#
# function get_trt_model_matrix
#
# Author: Liu et al. (2022), with adaptations by Lara Maleyeff
#
# Function:             get_trt_model_matrix
# Description:          Used to compute the design matrix for external data 
#                       (in adaptive enrichment steps and accuracy calculation).
#                       We need this function to ensure that the coefficient interpretation
#                       is identical between the original fitted model and treatment 
#                       effect estimates from external data
#                   
# Parameters:           X_cont        Matrix with continuous variable values in each column
#                       xrange        Range of each continuous variable
#                       intKnots      Location of internal knots
#                       numIntKnots   Number of internal knots for each spline
# Returns:              Spline design matrix for given continuous variables
get_trt_model_matrix <- function(X_cont, xrange, intKnots, numIntKnots) {
  J <- ncol(X_cont)
  for(j in 1:J){
    a <- xrange[j,1]
    b <- xrange[j,2]
    x <- as.numeric(X_cont[,j])
    
    if (numIntKnots == 0) {
      curKnots = c()
    } else if (numIntKnots == 1) {
      curKnots = intKnots[,j]
    } else {
      curKnots = intKnots[j,]
    }
    B <- bs(x, knots=curKnots, degree=3, Boundary.knots=c(a,b), intercept=TRUE)
    
    # Create the Omega matrix
    formOmega <- function(a,b,myKnots)
    {
      allKnots <- c(rep(a,4), myKnots, rep(b,4))
      K <- length(myKnots) 
      L <- 3*(K+8)
      xtilde <- (rep(allKnots, each=3)[-c(1,(L-1),L)] + rep(allKnots, each=3)[-c(1,2,L)])/2
      wts <- rep(diff(allKnots), each=3)*rep(c(1,4,1)/6, K+7)
      Bdd <- spline.des(allKnots, xtilde, derivs=rep(2,length(xtilde)), outer.ok=TRUE)$design 
      Omega <- t(Bdd*wts)%*%Bdd  
      return(Omega)
    }
    
    Omega <- formOmega(a,b,curKnots)
    
    # Obtain the spectral decomposition of Omega
    eigOmega <- eigen(Omega)
    
    # Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$
    indsZ <- 1:(numIntKnots+2)
    UZ <- eigOmega$vectors[,indsZ]
    LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
    
    # Form Z matrices
    Z <- B %*% LZ
    
    if(j==1){
      design.j <- cbind(x, Z)
      colnames(design.j) <- c(paste0("x", j), paste0("x", j, "_z", 1:ncol(Z)))
      design.spline <- cbind(data.frame(trt = 1), design.j)
    }
    else{
      design.j <- cbind(x, Z)
      colnames(design.j) <- c(paste0("x", j), paste0("x", j, "_z", 1:ncol(Z)))
      design.spline <- cbind(design.spline, design.j)
    }
  }
  return(design.spline)
}

#
# function get_splines
#
# Author: Liu et al. (2022), with adaptations by Lara Maleyeff
#
# Function:             get_splines
# Description:          Function to create a spline design matrix based on original data,
#                       including continuous variables, and number of knots
#                       Automatically selects knots to be at equally spaced quantiles
# Parameters:           data          A data frame with the observations arranged by row, and including the column:
#                                     - trt: the group indicator (=1 for experimental group; =0 for control group)
#                       X             Matrix with continuous variable values in each column, same patient order as data
#                       xrange        Matrix with 2 columns which gives the range of marker values, 
#                                     one marker per row, defaults to c(0,1) for each marker
#                       numIntKnots   Number of internal knots for each spline
# Returns:              data.full     Spline design matrix for given continuous variables
#                       formula       Formula describing mean relationship to splines, treatment, and 
#                                     their interactions
#                       stabCheck     Vector of booleans indicating whether the design matrix is stable for each 
#                                     continuous variable
#                       ncolz         Number of columns for each continuous variable
#                       intKnots      Location of the internal knots, to be used for adaptive
#                                     enrichment steps and accuracy calculations
get_splines <- function(data, X, xrange=cbind(rep(0, ncol(X)), rep(1, ncol(X))), numIntKnots){
  # Set up the design matrix and related quantities for each biomarker
  J <- ncol(X)
  stabCheck_idx <- rep(FALSE, J)
  intKnots=t(apply(as.matrix(X), 2, function(x) {
    quantile(unique(x), seq(0,1,length=(numIntKnots+2)))[-c(1,(numIntKnots+2))]
  }))
  for(j in 1:J){
    a <- xrange[j,1]
    b <- xrange[j,2]
    x <- as.numeric(X[,j])
    
    if (numIntKnots == 0) {
      curKnots = c()
    } else if (numIntKnots == 1) {
      curKnots = intKnots[,j]
    } else {
      curKnots = intKnots[j,]
    }
    B <- bs(x, knots=curKnots, degree=3, Boundary.knots=c(a,b), intercept=TRUE)
    
    # Create the Omega matrix
    formOmega <- function(a,b,intKnots)
    {
      allKnots <- c(rep(a,4), intKnots, rep(b,4))
      K <- length(intKnots) 
      L <- 3*(K+8)
      xtilde <- (rep(allKnots, each=3)[-c(1,(L-1),L)] + rep(allKnots, each=3)[-c(1,2,L)])/2
      wts <- rep(diff(allKnots), each=3)*rep(c(1,4,1)/6, K+7)
      Bdd <- spline.des(allKnots, xtilde, derivs=rep(2,length(xtilde)), outer.ok=TRUE)$design 
      Omega <- t(Bdd*wts)%*%Bdd  
      return(Omega)
    }
    
    Omega <- formOmega(a,b,curKnots)
    
    # Obtain the spectral decomposition of Omega
    eigOmega <- eigen(Omega)
    
    # Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$
    indsZ <- 1:(numIntKnots+2)
    UZ <- eigOmega$vectors[,indsZ]
    LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
    
    # Perform a stability check
    indsX <- (numIntKnots+3):(numIntKnots+4)
    UX <- eigOmega$vectors[,indsX]   
    L <- cbind(UX,LZ)
    stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))  
    stabCheck_idx[j] <- sum(stabCheck^2) > 1.0001*(numIntKnots+2)
    
    # Form Z matrices
    Z <- B%*%LZ
    
    # Interaction between trt and Z matrices
    trtZ <- (data$trt)*Z     
    
    if(j==1){
      design.j <- cbind(x, Z, data$trt, (data$trt)*x, trtZ)
      colnames(design.j) <- c(paste0("x", j), paste0("x", j, "_z", 1:ncol(Z)), "trt", paste0("x", j, "trt"), paste0("x", j, "_z", 1:ncol(Z), "trt"))
      data.full <- cbind(data.frame(Y = data$Y), design.j)
    }
    else{
      design.j <- cbind(x, Z, (data$trt)*x, trtZ)
      colnames(design.j) <- c(paste0("x", j), paste0("x", j, "_z", 1:ncol(Z)), paste0("x", j, "trt"), paste0("x", j, "_z", 1:ncol(Z), "trt"))
      data.full <- cbind(data.full, design.j)
    }
  }
  
  # Return the formula for the regression model
  formula <- as.formula(paste("Y ~ ", paste(colnames(data.full)[-c(1)], collapse="+"))) 
  
  return(list(data=data.full, formula=formula, stabCheck=stabCheck_idx, ncolz=numIntKnots+2, intKnots = intKnots))
}
 