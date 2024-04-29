######################################################################################################
#                             Simulation Study I of Maleyeff et al. (2024)      
#           Compares MCMC procedures of Park et al. (2022) with proposed method, using proposed 
#           design for both. Assumes 5 binary and 1 continuous candidate tailoring variables
#                               Results are in Tables 2 and 3 of manuscript.
#                                   Contact: laramaleyeff@gmail.com                                  
#                                       Last updated: April 2024                                       
######################################################################################################

start_time = Sys.time()

params <- commandArgs(trailingOnly=TRUE)

scenario <- as.numeric(params[1])

# Cutoff for effectiveness
B_1 <- as.numeric(params[2])

# Sparsity parameter for Park et al. (2022)
mu <- as.numeric(params[3])

# Sparsity parameter for Park et al. (2022)
tau <- as.numeric(params[4])

# Prior parameter for number of model terms (higher = more terms)
lambda_1 <- as.numeric(params[5])

# Prior parameter for number of knots in each spline term (higher = more knots)
lambda_2 <- as.numeric(params[6])

# Controls the size of the effective subspace
alpha <- as.numeric(params[7])

# Trend in the control group, (=0: no trend; =1: trend described in Table S1)
background <- as.numeric(params[8])

index <- as.numeric(params[9])


source("run_full_trial.R")

library(doParallel)
library(plyr)

# no_cont: Number of continuous candidate tailoring variables
# no_bin: Number of binary candidate tailoring variables
no_cont = 1
no_bin = 5

# sim_iters: Number of simulation iterations, when index = 0 we are testing so
# only one iteration is run
sim_iters = 10
if (index == 0) {
  sim_iters = 1
}

# Cutoff for futility
B_2 = 0.8

# n_pool: Size of the population pool from which we sample during the trial
# n_test: Size of the external dataset on which accuracy is calcualted
n_pool = 5000
n_test = 10000

# Names of the candidate continuous and binary variables
candsplinevars = paste0("X_", 1:(no_cont))
candbinaryvars = paste0("Z_", 1:(no_bin))
if (no_cont == 0) {
  candsplinevars = c()
} else if (no_bin == 0) {
  candbinaryvars = c()
} 

# interim_n: vector with length equal to the number of interim analyses to be 
# performed, denoting how many individuals are added to the sample at each one.
# Here, the total sample size is 500 with one interim analysis at 300 individuals.
interim_n = c(300,200)

# Functions describing the functional form of the interaction term for each 
# true tailoring variable. See Table 1 for a description.
scen_1 <- 0
scen_2 <- 0.28
scen_3 <- function(x) {
  return(x - 0.3)
}
scen_4 <- function(x) {
  return(0.7*x - 0.14)
}
scen_5 <- function(x) {
  return(0.8*x - 0.3)
}
scen_6 <- function(x1,x2) {
  return(0.9*x1 + 0.9*x2 - 0.2)
}
scen_7 <- function(x) {
  return(2.3*x - 1.15)
}
scen_8 <- function(x) {
  return(cos(x*2*pi))
}

# Generate external dataset to assess accuracy of treatment recommendations
# at the end of each trial
data_test =  data.frame(X_1 = runif(n_test,0,1),
                        Z_1 = rbinom(n_test,1,0.35),
                        Z_2 = rbinom(n_test,1,0.5),
                        Z_3 = rbinom(n_test,1,0.65),
                        Z_4 = rbinom(n_test,1,0.2),
                        Z_5 = rbinom(n_test,1,0.35),
                        trt = 1)

if (scenario == 1) {
  data_test$truth = scen_1
}

if (scenario == 2) {
  data_test$truth = scen_2
}

if (scenario == 3) {
  data_test$truth = scen_3(data_test$Z_1)
}

if (scenario == 4) {
  data_test$truth = scen_4(data_test$Z_2)
}

if (scenario == 5) {
  data_test$truth = scen_5(data_test$Z_3)
}

if (scenario == 6) {
  data_test$truth = scen_6(data_test$Z_4, data_test$Z_1)
}

if (scenario == 7) {
  data_test$truth = scen_7(data_test$X_1)
}

if (scenario == 8) {
  data_test$truth = scen_8(data_test$X_1)
}

# "true_trt" = 1 if treatment is truly effective and = 0 if not
data_test$true_trt = as.numeric(data_test$truth > 0)

# Function to generate data pool based on scenario considered. Returns data
# and a vector "true_tailoring_vars" containing the names of the true tailoring variables.
# This will be compared with the selected tailoring variables to determine correct marker
# detection rates.
generate_data <- function(scenario) {
  data_pool = data.frame(X_1 = runif(n_pool,0,1),
                         Z_1 = rbinom(n_pool,1,0.35),
                         Z_2 = rbinom(n_pool,1,0.5),
                         Z_3 = rbinom(n_pool,1,0.65),
                         Z_4 = rbinom(n_pool,1,0.2),
                         Z_5 = rbinom(n_pool,1,0.35),
                         trt = rbinom(n_pool,1,0.5))
  
  if (background == 0) {
    data_pool$nu = 0 
  } else {
    if (scenario == 1) {
      data_pool$nu = 0
    } else if (scenario == 2) {
      data_pool$nu = 0
    } else if (scenario == 3) {
      data_pool$nu = 0.5*data_pool$Z_1
    } else if (scenario == 4) {
      data_pool$nu = 0.5*data_pool$Z_2
    } else if (scenario == 5) {
      data_pool$nu = 0.5*data_pool$Z_3
    } else if (scenario == 6) {
      data_pool$nu = 0.5*data_pool$Z_4 + 0.3*data_pool$Z_1
    } else if (scenario == 7) {
      data_pool$nu = 0.3*data_pool$X_1
    } else if (scenario == 8) {
      data_pool$nu = 0.3*data_pool$X_1
    }
  }
  
  # Overall null effect
  if (scenario == 1) {
    data_pool$truth = scen_1
    true_tailoring_vars = c()
  }
  
  # Overall positive effect
  if (scenario == 2) {
    data_pool$truth = scen_2
    true_tailoring_vars = c()
  }
  
  # Low prevalence binary
  if (scenario == 3) {
    data_pool$truth = scen_3(data_pool$Z_1)
    true_tailoring_vars = c("Z_1")
  }
  
  # Moderate prevalence binary
  if (scenario == 4) {
    data_pool$truth = scen_4(data_pool$Z_2)
    true_tailoring_vars = c("Z_2")
  }

  # High prevalence binary
  if (scenario == 5) {
    data_pool$truth = scen_5(data_pool$Z_3)
    true_tailoring_vars = c("Z_3")
  }
  
  # Two binary
  if (scenario == 6) {
    data_pool$truth = scen_6(data_pool$Z_4,data_pool$Z_1)
    true_tailoring_vars = c("Z_4","Z_1")
  }

  # Linear continuous
  if (scenario == 7) {
    data_pool$truth = scen_7(data_pool$X_1)
    true_tailoring_vars = c("X_1")
  }
  
  # Nonlinear, nonmonotone continuous
  if (scenario == 8) {
    data_pool$truth = scen_8(data_pool$X_1)
    true_tailoring_vars = c("X_1")
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
  registerDoParallel(cores=ncores)
}

# Run the simulation study
results = foreach(i=1:sim_iters, .combine='rbind.fill') %dopar% {
  print(paste("iteration", i))
  datas = generate_data(scenario)
  family = gaussian()

  print("Park")
  start_time_park = Sys.time()
  results_park = runTrial(datas$data_pool,data_test,candsplinevars,
                          candbinaryvars,parseMHPark,getCutoffPark,
                          alpha,
                          B_1,
                          B_2,
                          B = 5000,
                          burnin = 10000,
                          thin = 1,
                          family=family,
                          true_tailoring_vars = datas$true_tailoring_vars,
                          interim_n = interim_n,
                          mu = mu,
                          tau = tau
  )
  end_time_park = Sys.time()
  results_park$time = difftime(end_time_park, start_time_park, units = "mins")
  colnames(results_park) = paste0("park_",colnames(results_park))
  
 
  print("Maleyeff")
  start_time_maley = Sys.time()
  results_maley = runTrial(datas$data_pool,data_test,candsplinevars,
                           candbinaryvars,parseMHMaley,getCutoffMaley,
                           alpha,
                           B_1,
                           B_2,
                           B = 2000,
                           burnin = 30000,
                           thin = 10,
                           family=family,
                           true_tailoring_vars = datas$true_tailoring_vars,
                           interim_n = interim_n,
                           lambda_1 = lambda_1,
                           lambda_2 = lambda_2
  )
  end_time_maley = Sys.time()
  results_maley$time = difftime(end_time_maley, start_time_maley, units = "mins")
  colnames(results_maley) = paste0("maley_",colnames(results_maley))
  return(cbind(results_park,
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
results$background = background
results$tau = tau
results$mu = mu
results$candbinaryvars = paste(candbinaryvars,collapse=",")
results$candsplinevars = paste(candsplinevars,collapse=",")

end_time = Sys.time()
results$time = difftime(end_time, start_time, units = "mins")
print(results$time)
setwd("results")
out.name = paste0("all_trials_", paste( "scen",    scenario,
                                        "B_1",     B_1,
                                        "B_2",     B_2,
                                        "no_cont", no_cont,
                                        "no_bin",  no_bin,
                                        "mu",      mu,
                                        "tau",     tau,
                                        "lambda_1", lambda_1,
                                        "lambda_2", lambda_2,
                                        "background",background,
                                        "alpha",    alpha,
                                        "index",    index, sep="_"),".csv")

write.csv(results, out.name, row.names=FALSE)
