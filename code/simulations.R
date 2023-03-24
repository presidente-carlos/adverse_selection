#-------------------------------------------------#
# Evaulation, Tests and Simulations of Algorithm 1
#------------------------------------------------#

## Description: 
# This code evaluates all DGP described in Gonzalez (2023)
# Additionally, it conducts further tests and simulations related to parameter
# sensibility

# 0. Load relevant files
#-----------------------

library(dplyr)        # For syntax
library(ggplot2)      # For plotting
library(future)       # General parallel computing tools
library(future.apply) # User friendly parallel computing functions
library(gridExtra)    # For general display
library(resample)     # For quick colVars()

# Load DGP
lambda = 0.7
source("algorithm.R")
source("dgp.R")


# 1. Benchmark Evaluation
#----------------------------------
K = 500 # Time periods
R = 100 # Replications


# Generate parameters
params = param_gen(K = K)


# 1.1. Uniform linear
# --------------------
final_display(R, data_function = unif_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = unif_lin_theory, label = "unif_lin")


# 1.2. Uniform non-linear
#--------------------------
# It seems to work better with smaller eta
final_display(R, data_function = unif_non_lin, 
              B = params$B, eta = 0.033, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = unif_non_lin_theory, label = "unif_non_lin")

# 1.3. Uniform - Uniform
#-----------------------
final_display(R, data_function = unif_unif, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = unif_unif_theory, label = "unif_unif")

# 1.4. Bernouilli linear
#-----------------------
# Probability p is imputed in dgp.R
final_display(R, data_function = bern_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = bern_lin_theory, label = "bern_lin")

# 1.5. Beta(2, 1) linear
# --------------------
final_display(R, data_function = beta21_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = beta21_lin_theory, label = "beta21_lin")

# 1.6. Beta(1/2, 1/2) linear
# --------------------
final_display(R, data_function = beta05_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = beta05_lin_theory, label = "beta05_lin")

# 1.7. Beta(2, 1) non-linear
# -------------------------
final_display(R, data_function = beta21_non_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = beta21_non_lin_theory, label = "beta21_non_lin")

# 2. Worst-case eta
#-------------------
# 2.1. Uniform linear
# --------------------
eta_vector = log_sequence(from = params$eta, to = 0.2, length.out = 2)
final_display_eta_loop(R = 100, data_function = unif_lin, B = params$B, 
                       gamma = params$gamma, K = 1000, lambda = 0.7,
                       theory_function = unif_lin_theory, label = "unif_lin_add",
                       K_store = 200, eta_vector = eta_vector)

# 2.2. Uniform non-linear
#------------------------
eta_vector = log_sequence(from = params$eta, to = 0.2, length.out = 10)
final_display_eta_loop(R = 100, data_function = unif_non_lin, B = params$B, 
                       gamma = params$gamma, K = 5000, lambda = 0.7,
                       theory_function = unif_non_lin_theory, label = "unif_non_lin_add",
                       K_store = 500, eta_vector = eta_vector)


# 2.Parameter Sensitivity
#------------------------

# Intuition, you 