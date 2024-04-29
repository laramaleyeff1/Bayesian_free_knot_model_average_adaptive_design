start_time = Sys.time()

params <- commandArgs(trailingOnly=TRUE)

scenario <- as.numeric(params[1])

B_1 <- as.numeric(params[2])

lambda_1 <- as.numeric(params[3])

lambda_2 <- as.numeric(params[4])

alpha <- as.numeric(params[5])

index <- as.numeric(params[6])


source("run_full_trial.R")

library(doParallel)
library(plyr)

numIntKnots = 5
no_cont = 2
no_bin = 0
sim_iters = 10
if (index == 0) {
  sim_iters = 1
}

# Sparsity parameters for Park et al. (2022), not used here
mu = 10
tau = 0.1

B_2 = 0.8

n_pool = 5000
n_test = 10000

candsplinevars = paste0("X_", 1:(no_cont))
candbinaryvars = c()

interim_n = c(300,200)

# Functions for each scenario
scen_1 = 0
scen_2 = 0.35
scen_3 <- function(x) {
  2.3*(x - 0.5)
}
scen_4 <- function(x) {
  cos(x*2*pi)
}
scen_5 <- function(x) {
  return((exp(25*(x-0.5))/(1+exp(25*(x-0.5))))*1.4 - 0.6)
}
scen_6 <- function(x) {
  return(I(x<=0.5)*(exp(30*(x-0.3))/(1+exp(30*(x-0.3))))*2 +
                 I(x>0.5)*(1/(1+exp(30*(x-0.7))))*2 - 1)
}
scen_7 <- function(x) {
  return(I(x<=0.5)*(1/(1+exp(30*(x-0.3))))*1.5 +
           I(x>0.5)*(exp(30*(x-0.7))/(1+exp(30*(x-0.7))))*1.5 - 0.75)
}
scen_8 <- function(x1,x2) {
 return(2.3*(x1 - 0.5) + cos(x2*2*pi)) 
}

data_test =  data.frame(X_1 = runif(n_test,0,1),
                        X_2 = runif(n_test,0,1),
                        trt = 1)

if (scenario == 1) {
  data_test$truth = scen_1
}

if (scenario == 2) {
  data_test$truth = scen_2
}

if (scenario == 3) {
  data_test$truth = scen_3(data_test$X_1)
}

if (scenario == 4) {
  data_test$truth = scen_4(data_test$X_1)
}

if (scenario == 5) {
  data_test$truth = scen_5(data_test$X_1)
}

if (scenario == 6) {
  data_test$truth = scen_6(data_test$X_1)
}

if (scenario == 7) {
  data_test$truth = scen_7(data_test$X_1)
}

if (scenario == 8) { 
  data_test$truth = scen_8(data_test$X_1,data_test$X_2)
}

data_test$true_trt = as.numeric(data_test$truth > 0)

generate_data <- function(scenario) {
  data_pool = data.frame(X_1 = runif(n_pool,0,1),
                         X_2 = runif(n_pool,0,1),
                         trt = rbinom(n_pool,1,0.5))
  
  data_pool$nu = 0.5*data_pool$X_1
  
  if (scenario == 1) {
    data_pool$truth = scen_1
    true_tailoring_vars = c()
  }
  
  if (scenario == 2) {
    # 0.5 good
    data_pool$truth = scen_2
    true_tailoring_vars = c()
  }
  
  if (scenario == 3) {
    data_pool$truth = scen_3(data_pool$X_1)
    true_tailoring_vars = c("X_1")
  }
  
  if (scenario == 4) {
    data_pool$truth = scen_4(data_pool$X_1)
    true_tailoring_vars = c("X_1")
  }
  
  if (scenario == 5) {
    data_pool$truth = scen_5(data_pool$X_1)
    true_tailoring_vars = c("X_1")
  }
  
  if (scenario == 6) {
    data_pool$truth = scen_6(data_pool$X_1)
    true_tailoring_vars = c("X_1")
    
  }
  
  if (scenario == 7) {
    data_pool$truth = scen_7(data_pool$X_1)
    true_tailoring_vars = c("X_1")
  }
  
  if (scenario == 8) { 
    data_pool$truth = scen_8(data_pool$X_1,data_pool$X_2)
    true_tailoring_vars = c("X_1", "X_2")
  }
  
  data_pool$nu = data_pool$nu + data_pool$truth*data_pool$trt
  data_pool$Y = data_pool$nu + rnorm(n_pool,0,1)
  
  return(list(data_pool = data_pool,
              true_tailoring_vars = true_tailoring_vars))
}

if (index != 0) {
  # Use the environment variable SLURM_CPUS_PER_TASK to set the number of cores.
  # This is for SLURM. Replace SLURM_CPUS_PER_TASK by the proper variable for your system.
  # Avoid manually setting a number of cores.
  ncores = Sys.getenv("SLURM_CPUS_PER_TASK")
  registerDoParallel(cores=ncores)# Shows the number of Parallel Workers to be used
  
}

results = foreach(i=1:sim_iters, .combine='rbind.fill') %dopar% {
  print(paste("iteration", i))
  datas = generate_data(scenario)
  family = gaussian()

  start_time_liu = Sys.time()
  print("Liu")
  results_liu = runTrial(datas$data_pool,
                         data_test,
                         candsplinevars,
                         candbinaryvars,
                         parseMHLiu,
                         getCutoffLiu,
                         alpha,
                         B_1 = B_1,
                         B_2 = B_2,
                         B = 2000,
                         burnin = 30000,
                         thin = 10,
                         family = family,
                         true_tailoring_vars = datas$true_tailoring_vars,
                         interim_n = interim_n,
                         numIntKnots = numIntKnots
  )
  end_time_liu = Sys.time()
  results_liu$time = difftime(end_time_liu, start_time_liu, units = "mins")
  colnames(results_liu) = paste0("liu_",colnames(results_liu))

  if (FALSE) {
    print("Park")
    start_time_park = Sys.time()
    results_park = runTrial(datas$data_pool,
                            data_test,
                            candsplinevars,
                            candbinaryvars,
                            parseMHPark,
                            getCutoffPark,
                            alpha,
                            B_1,
                            B_2,
                            B = 5000,
                            burnin = 10000,
                            thin = 1,
                            family = family,
                            true_tailoring_vars = datas$true_tailoring_vars,
                            interim_n = interim_n,
                            mu = mu,
                            tau = tau
    )
    end_time_park = Sys.time()
    results_park$time = difftime(end_time_park, start_time_park, units = "mins")
    colnames(results_park) = paste0("park_",colnames(results_park))
  } else {
    results_park = data.frame(park_success = FALSE)
  }
 
  print("Maleyeff")
  start_time_maley = Sys.time()
  results_maley = runTrial(datas$data_pool,
                           data_test,
                           candsplinevars,
                           candbinaryvars,
                           parseMHMaley,
                           getCutoffMaley,
                           alpha,
                           B_1,
                           B_2,
                           B = 2000,
                           burnin = 30000,
                           thin = 10,
                           family = family,
                           true_tailoring_vars = datas$true_tailoring_vars,
                           interim_n = interim_n,
                           lambda_1 = lambda_1,
                           lambda_2 = lambda_2
  )
  end_time_maley = Sys.time()
  results_maley$time = difftime(end_time_maley, start_time_maley, units = "mins")
  colnames(results_maley) = paste0("maley_",colnames(results_maley))
  return(cbind(results_liu,
               results_park,
               results_maley))
}

results$B_1 = B_1
results$B_2 = B_2
results$alpha = alpha
results$scenario = scenario
results$n_test = n_test
results$n_pool = n_pool
results$interim_n = paste(interim_n,collapse=",")
results$sim_iters = sim_iters
results$no_cont = no_cont
results$no_bin = no_bin
results$lambda_1 = lambda_1
results$lambda_2 = lambda_2
results$numIntKnots = numIntKnots
results$tau = tau
results$mu = mu

results$candbinaryvars = paste(candbinaryvars,collapse=",")
results$candsplinevars = paste(candsplinevars,collapse=",")

end_time = Sys.time()
results$time = difftime(end_time, start_time, units = "mins")
print(results$time)
setwd("results_two_interim")
out.name = paste0("all_trials_", paste( "scen", scenario,
                                        "B_1", B_1,
                                        "B_2", B_2,
                                        "no_cont", no_cont,
                                        "no_bin", no_bin,
                                        "mu", mu,
                                        "tau", tau,
                                        "lambda_1", lambda_1,
                                        "lambda_2", lambda_2,
                                        "alpha", alpha,
                                        "numIntKnots", numIntKnots,
                                        "index", index, sep="_"),".csv")

write.csv(results, out.name, row.names=FALSE)
