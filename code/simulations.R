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
library(latex2exp)    # For LaTeX titles in graphs
library(pracma)       # For erf() access
library(forcats)  # For not-stacked area plots


# Load DGP
lambda = 0.7
source("dgp.R")
source("algorithm.R")
source("additional.R")


# 1. Benchmark Evaluation
#----------------------------------
K = 500 # Time periods
# R = 100 # Replications


# Generate parameters
params = param_gen(K = K)

# Set seed
set.seed(1234)

# 1.1. Uniform linear
# --------------------
final_display(R = 1000, data_function = unif_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = unif_lin_theory, label = "unif_lin")


# 1.2. Uniform non-linear (with small eta, larger K and K_store)
#--------------------------
final_display(R = 1000, data_function = unif_non_lin, 
              B = params$B, eta = 0.05, gamma = 0.029,
              K = 1000, lambda = lambda, K_store = 100,
              theory_function = unif_non_lin_theory, label = "unif_non_lin")

# 1.3. Uniform - Uniform
#-----------------------
final_display(R = 1000, data_function = unif_unif, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = unif_unif_theory, label = "unif_unif")

# 1.4. Bernouilli linear
#-----------------------
# Probability p is imputed in dgp.R
final_display(R = 1000, data_function = bern_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = K, lambda = lambda,
              theory_function = bern_lin_theory, label = "bern_lin")

# 1.5. Beta(2, 1) linear
# --------------------
final_display(R = 1000, data_function = beta21_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = 500, lambda = lambda,
              theory_function = beta21_lin_theory, label = "beta21_lin")

# 1.6. Beta(1/2, 1/2) linear
# --------------------
final_display(R = 1000, data_function = beta05_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = 500, lambda = lambda,
              theory_function = beta05_lin_theory, label = "beta05_lin")

# 1.7. Beta(2, 1) non-linear
# -------------------------
final_display(R = 1000, data_function = beta21_non_lin, 
              B = params$B, eta = 0.132, gamma = 0.029,
              K = 500, lambda = lambda,
              theory_function = beta21_non_lin_theory, label = "beta21_non_lin")


# 2. SD/eta Robustness Check
#----------------------------

set.seed(1234)
eta_values = log_sequence(0.0027, 2, 5)
sd_vector = c(0.4, 0.6, 0.8, 1)
for (eta_val in eta_values){
  multi_final_display_sd(R = 1000, data_function = unif_lin_random, B = params$B,
                         eta = eta_val, gamma = params$gamma, K = 500, lambda = lambda,
                         theory_function = unif_lin_random_theory, 
                         label = "unif_lin_random", K_store = 100, sd_vector = sd_vector)
}

# 3. Minimum Wage Analysis
#-------------------------
mw_vector = c(0, 0.2, 0.3, 0.4, 0.6)
multi_final_display_mw(R = 1000, B = params$B, data_function = unif_lin_high_u_low_v,
                       eta = 0.132, gamma = 0.029, K = 1000, lambda = lambda,
                       theory_function = unif_lin_theory_high_u_low_v, 
                       label = "unif_lin_mw_hl", K_store = 100, 
                       mw_vector = mw_vector) #x_optim = 1/4

multi_final_display_mw(R = 1000, B = params$B, data_function = unif_lin_high_u_high_v,
                       eta = 0.132, gamma = 0.029, K = 1000, lambda = lambda,
                       theory_function = unif_lin_theory_high_u_high_v, 
                       label = "unif_lin_mw_hh", K_store = 100, 
                       mw_vector = mw_vector) # x_optim = 1/2

multi_final_display_mw(R = 1000, B = params$B, data_function = unif_lin_low_u_low_v,
                       eta = 0.132, gamma = 0.029, K = 1000, lambda = lambda,
                       theory_function = unif_lin_theory_low_u_low_v, 
                       label = "unif_lin_mw_ll", K_store = 100, 
                       mw_vector = mw_vector) # x_optim = 3/16

multi_final_display_mw(R = 1000, B = params$B, data_function = unif_lin_low_u_high_v,
                       eta = 0.132, gamma = 0.029, K = 1000, lambda = lambda,
                       theory_function = unif_lin_theory_low_u_high_v, 
                       label = "unif_lin_mw_lh", K_store = 100, 
                       mw_vector = mw_vector) # x_optim = 3/8

# 3. LIPC
#-------------------------
p_theta_vector = c(1, 0.8, 0.6, 0.4, 0.2)
multi_final_display_lipc(R = 1000, B = params$B, data_function = unif_lin_high_u_low_v,
                       eta = 0.132, gamma = 0.029, K = 500, lambda = lambda,
                       theory_function = unif_lin_theory_high_u_low_v, 
                       label = "unif_lin_lipc_hl", K_store = 50, 
                       p_theta_vector = p_theta_vector) #x_optim = 1/4

multi_final_display_lipc(R = 1000, B = params$B, data_function = unif_lin_high_u_high_v,
                       eta = 0.132, gamma = 0.029, K = 500, lambda = lambda,
                       theory_function = unif_lin_theory_high_u_high_v, 
                       label = "unif_lin_lipc_hh", K_store = 50, 
                       p_theta_vector = p_theta_vector) # x_optim = 1/2

multi_final_display_lipc(R = 1000, B = params$B, data_function = unif_lin_low_u_low_v,
                       eta = 0.132, gamma = 0.029, K = 500, lambda = lambda,
                       theory_function = unif_lin_theory_low_u_low_v, 
                       label = "unif_lin_lipc_ll", K_store = 50, 
                       p_theta_vector = p_theta_vector) # x_optim = 3/16

multi_final_display_lipc(R = 1000, B = params$B, data_function = unif_lin_low_u_high_v,
                       eta = 0.132, gamma = 0.029, K = 500, lambda = lambda,
                       theory_function = unif_lin_theory_low_u_high_v, 
                       label = "unif_lin_lipc_lh", K_store = 50, 
                       p_theta_vector = p_theta_vector) # x_optim = 3/8

# 4. Productivity Shocks
#----------------------------
final_display_shocks(R = 1000, data_function = unif_lin_unchanged_fair, 
                     B = params$B, eta = 0.00267, 
                     gamma = 0.029, K = 10000, lambda = 0.1,
                     theory_function_1 = unif_lin_theory_shock_first,
                     theory_function_2 = unif_lin_theory_shock_second_unch,
                     label = "unif_lin_un_fair",
                     K_store = 500, sd = 0, mw = 0)

final_display_shocks(R = 1000, data_function = unif_lin_unchanged_unfair, 
                     B = params$B, eta = 0.00267, 
                     gamma = 0.029, K = 10000, lambda = 0.1,
                     theory_function_1 = unif_lin_theory_shock_first,
                     theory_function_2 = unif_lin_theory_shock_second_unch,
                     label = "unif_lin_un_unfair",
                     K_store = 500, sd = 0, mw = 0)

final_display_shocks(R = 1000, data_function = unif_lin_changed_fair, 
                     B = params$B, eta = 0.00267, 
                     gamma = 0.029, K = 10000, lambda = 0.1,
                     theory_function_1 = unif_lin_theory_shock_first,
                     theory_function_2 = unif_lin_theory_shock_second_ch,
                     label = "unif_lin_changed_fair",
                     K_store = 500, sd = 0, mw = 0)

final_display_shocks(R = 1000, data_function = unif_lin_changed_unfair, 
                     B = params$B, eta = 0.00267, 
                     gamma = 0.029, K = 10000, lambda = 0.1,
                     theory_function_1 = unif_lin_theory_shock_first,
                     theory_function_2 = unif_lin_theory_shock_second_ch,
                     label = "unif_lin_changed_unfair",
                     K_store = 500, sd = 0, mw = 0)

# To do list:
# - Understand theory welfare plot
# - Even more periods??
# - Lower eta (lowering eta does a bit the trick). Y eso que no lo estamos comparando con el unfair. Estamos jodidos!
# - Do the rest of them
# - Slides
# Study slides
