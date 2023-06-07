#-------------------------------------------------#
# Policy Simulations of Algorithm 1
#------------------------------------------------#

## Description: 
# This code conducts simulations to evaluate the theoretical policy
# prescriptions described in Gonzalez (2023)

# 0a. Load relevant files
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


# 0b. Set relevant parameters
#-----------------------------
K = 500 # Time periods
R = 100 # Replications


# Generate parameters
params = param_gen(K = K)

# 1. Minimum Wage Policies
#----------------------------

# 1a. High U, Low V
#------------------
mw_vector_hl = c(0, 1/8, 1/4, 1/2)
final_display_mw_loop(R = 100, data_function = unif_lin_high_u_low_v,
                      B = params$B, eta = 0.132, gamma = 0.029, K = K, lambda = lambda,
                      theory_function = unif_lin_theory_high_u_low_v,
                      label = "unif_lin_mw_hl", K_store = 50, 
                      mw_vector = mw_vector_hl)

# 1b. High U, High V
#------------------
mw_vector_hh = c(0, 1/4, 1/2, 1)
final_display_mw_loop(R = 100, data_function = unif_lin_high_u_high_v,
                      B = params$B, eta = 0.132, gamma = 0.029, K = K, lambda = lambda,
                      theory_function = unif_lin_theory_high_u_high_v,
                      label = "unif_lin_mw_hh", K_store = 50, 
                      mw_vector = mw_vector_hh)


# 1a. Low U, Low V
#------------------
mw_vector_ll = c(0, 3/32, 3/16, 3/8)
final_display_mw_loop(R = 100, data_function = unif_lin_low_u_low_v,
                      B = params$B, eta = 0.132, gamma = 0.029, K = K, lambda = lambda,
                      theory_function = unif_lin_theory_low_u_low_v,
                      label = "unif_lin_mw_ll", K_store = 50, 
                      mw_vector = mw_vector_ll)


# 1a. Low U, High V
#------------------
mw_vector_lh = c(0, 3/16, 3/8, 3/4)
final_display_mw_loop(R = 100, data_function = unif_lin_low_u_high_v,
                      B = params$B, eta = 0.132, gamma = 0.029, K = K, lambda = lambda,
                      theory_function = unif_lin_theory_low_u_high_v,
                      label = "unif_lin_mw_lh", K_store = 50, 
                      mw_vector = mw_vector_lh)

# 2a. Unchanged productivity, unfair
#-----------------------------------
final_display_shocks(R = 20, data_function = unif_lin_unchanged_fair, B, eta, 
                     gamma, K = 1500, lambda = 0.35,
                     theory_function_1, theory_function_2, label = "unif_lin_un_fair",
                     K_store = 100, sd = 0, mw = 0)

