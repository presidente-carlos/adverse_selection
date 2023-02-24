#--------------------------------------#
# Algorithm 1: Functions and evaluation
#--------------------------------------#

## Description: 
# This code prepares a series of functions including:
# - Tempered Exp3 (Algorithm 1)
# - Auxiliary (optimal) parameter generating functions
# - Data preparation and plotting functions

# Additionally, it evaluates all DGP described in Gonzalez (2023)

# 0. Libraries and load files
#----------------------------

library(dplyr)        # For syntax
library(ggplot2)      # For plotting
library(future)       # General parallel computing tools
library(future.apply) # User friendly parallel computing functions
library(gridExtra)    # For general display

# Load DGP
source("dgp.R")

# 1. Parameters Generator Function
# --------------------------------

param_gen = function(K){
  c_1 = 1/50 * ( K / log(K) )^(1/3)
  gamma = c_1 * (log(K) / K)^(1/3)
  c_3 = 10 * gamma 
  c_2 =  c_1 / (c_1 * c_3 + (K / log(K))^(1/3))
  eta = c_2 * gamma^2
  B = c_3 / gamma
  print(paste0("eta = ", eta,
               " B = ", B,
               " gamma = ", gamma))
  tibble("gamma" = gamma, "eta" = eta, "B" = B)
}

# 2. Algorithm 1
#---------------

# Tempered Exp3 for Monopolist Wage Setting

exp3_monop = function(B, eta, gamma, K, lambda, u, v, x_optim,
                      K_store = 50){
  
  # Inputs
  # B, eta, gamma, K, lambda = parameters
  # u, v = sequence of u and v realizations
  # x_optim = optimal policy in hindsight (given F_{U,v})
  
  # Initialize Inputs
  b_seq = 1:(B+1)
  x_disc = (b_seq - 1)/B
  
  # Initialize estimates
  g_hat = rep(0, B+1)
  u_hat = rep(0, B+1)
  
  # Initialize store vectors
  s_hat = rep(NA, B+1)
  p_i = rep(gamma / (B+1), B+1)

  store = tibble(b_i = rep(NA, K),
                 x_disc = NA,
                 v_i = NA,
                 u_i = NA, 
                 indic = NA,
                 s_i = NA,
                 s_optim = NA,
                 regret = 0,
                 cum_regret = NA,
                 avg_cum_regret = NA)
  
  # Store p_ib every 50 periods
  store_pib = tibble(x_disc = rep(x_disc, round(K / K_store)), 
                     p_ib = NA,
                     round = NA)
  
  for (i in 1:K){
    
    # Probability updates
    
    for (b in 1:(B+1)){
      if (b == 1){
            s_hat[b] = u_hat[b] - x_disc[b] * g_hat[b]
            }else{
            s_hat[b] = u_hat[b] - x_disc[b] * g_hat[b] + 
                          lambda / B * sum(g_hat[1:max((b-1), 1)])
            }
    }

    for (b in 1:(B+1)){
      # Importance Weighted estimator
      p_i[b] = (1 - gamma) * (exp(eta * s_hat[b]) / sum(exp(eta * s_hat))) +
               gamma / (B+1)
    }

    # Sample arm according to probabilities
    b_i = sample(b_seq, 1, prob = p_i)
    
    # Updating estimates
    x_disc_i = x_disc[b_i]
    indic = 1*(x_disc_i >= v[i])
    g_hat[b_i] = g_hat[b_i] + indic / p_i[b_i]
    u_hat[b_i] = u_hat[b_i] + u[i] * indic / p_i[b_i]
    
    # Regret
    s_i = indic * (u[i] - x_disc_i + lambda * (x_disc_i - v[i]))
    s_optim = 1*(x_optim[i] >= v[i])*(u[i] - x_optim[i] + 
                                        lambda * (x_optim[i] - v[i]))
    regret = s_optim - s_i

    store$b_i[i] = b_i          # Arm selected in period i
    store$v_i[i] = v[i]         # u in period i
    store$u_i[i] = u[i]         # v in period i
    store$indic[i] = indic      # x >= v in period i
    store$x_disc[i] = x_disc_i  # Policy selected in period i
    store$s_i[i] = s_i          # Policy welfare in period i
    store$s_optim[i] = s_optim  # Optimal policy welfare in period i
    store$regret[i] = regret    # Regret in period i
    
    # Cumulative regret
    store$cum_regret[i] = sum(store$regret)
    store$avg_cum_regret[i] = store$cum_regret[i] / i
    
    # Store probabilities of selected periods
    if ((i-1)%%50 == 0){
      t = round(i / 50)
      store_pib[(t*length(x_disc) + 1):(t*length(x_disc) +
                                          length(x_disc)), ]$p_ib = p_i
      store_pib[(t*length(x_disc) + 1):(t*length(x_disc) +
                                          length(x_disc)), ]$round = i
    }
    
  }
  
  # Return outcomes
  list("regret" = store$avg_cum_regret, "store_pib" = store_pib$p_ib)

}

# 3. Data preparation and Plotting
#---------------------------------

plotting = function(results, K, expression, K_store = 50, label){
  cum_regret = do.call(rbind, results[1,])
  
  # Avg cum regret
  cum_regret_tibble = tibble(x = seq(1:K),
                             y = colMeans(cum_regret))
  avg_cum_regret = cum_regret_tibble |> ggplot() +
                   geom_smooth(aes(x = x, y = y), se = FALSE,
                               color = "black") +
                   theme_minimal(base_size = 18) +
                   xlab("Period K") +
                   ylab("Average Cummulative Regret")

  # p_ib across K
  p_ib_across_K = do.call(rbind, results[2,])
  avg_p_ib_across_K = colMeans(p_ib_across_K)
  
  b_seq = 1:(params$B+1)
  x_disc = (b_seq - 1)/(params$B)
  times_arg = K/K_store
  
  p_ib_plotting = tibble("p_ib" = avg_p_ib_across_K, 
                         "K_round" = rep(seq(from = 1, to = K, by = K_store), 
                                   each = length(x_disc)),
                         "x_disc" = rep(x_disc, times_arg))
  
  # Empirical S_i(x)
  # Set sizes of titles. Ideally just one x title and one y title
  max_pib = max(p_ib_plotting$p_ib)
  P = list()
  t = 1
  k_minus = 0
  for (k in seq(K_store*3 + 1, K, by = K_store*3)){
    P[[t]] = p_ib_plotting |> filter(K_round <= k, K_round >= k_minus) |>
                              ggplot() +
                              geom_area(aes(x = x_disc,
                                            y = p_ib), color = "black", 
                                            fill = "coral") +
                              theme_minimal() +
                              xlab("Action Space x") +
                              ylab("Sampling probability") +
                              ylim(0, max_pib) + 
                              facet_grid(~K_round)
    t = t+1
    k_minus = k
  }

  # Theoretical S_i(x)
  x_seq_plot = seq(0, 1, by = 0.01)
  plot_tib = tibble("x_disc" = x_seq_plot, 
                    "s_ib" = expression(x_seq_plot))
  
  actual = plot_tib |> ggplot() +
                       geom_area(aes(x = x_disc,
                                     y = s_ib), 
                                     color = "black",
                                     fill = "coral2", size = 1) +
                       theme_minimal(base_size = 18) +
                       xlab("Action Space x") +
                       ylab("Expected Monopolist Profit S(x)")
  
  plot_list = list("avg_cum_regret" = plot(avg_cum_regret), 
                   "theory_welfare" = plot(actual), 
                   "empirical_probs" = do.call("grid.arrange", c(P, nrow = 3)))
  
  lapply(names(plot_list), 
         function(x)ggsave(filename=paste("../plots/", label, "_", x,
                                          ".jpeg",sep=""), 
                           plot = plot_list[[x]], bg = "white"))
  
}

final_display = function(R, data_function, B = params$B, eta = params$eta, 
                         gamma = params$gamma, K = K, lambda = lambda,
                         theory_function, label){
  
  simulation_data = data_function(K)
  
  # Prepare session for parallel computing
  plan(multisession(workers = 4))
  
  results = future_replicate(R, exp3_monop(B = params$B, eta = 0.132, 
                                           gamma = 0.029,
                                           K = K, lambda = lambda,
                                           u = simulation_data$u,
                                           v = simulation_data$v, 
                                           x_optim = simulation_data$x_optim))
  
  plotting(results = results, K = K, 
           expression = theory_function, label = label)
  
}

# 4. Evaluation
#----------------------------------
K = 500 # Time periods
R = 1000 # Replications

# Generate parameters
params = param_gen(K = K)



# 4.1. Uniform linear
# --------------------
lambda = 0.7
final_display(R, data_function = unif_lin, 
              B = params$B, eta = params$eta, gamma = params$gamma,
              K = K, lambda = lambda,
              theory_function = unif_lin_theory, label = "unif_lin")


# 4.2. Uniform non-linear
#--------------------------
lambda = 0.8
final_display(R, data_function = unif_non_lin, 
              B = params$B, eta = params$eta, gamma = params$gamma,
              K = K, lambda = lambda,
              theory_function = unif_non_lin_theory, label = "unif_non_lin")

# 4.3. Uniform - Uniform
#-----------------------
lambda = 0.7
final_display(R, data_function = unif_unif, 
              B = params$B, eta = params$eta, gamma = params$gamma,
              K = K, lambda = lambda,
              theory_function = unif_unif_theory, label = "unif_unif")

# 4.4. Bernouilli linear
#-----------------------
# Probability p is imputed in dgp.R
final_display(R, data_function = bern_lin, 
              B = params$B, eta = params$eta, gamma = params$gamma,
              K = K, lambda = lambda,
              theory_function = bern_lin_theory, label = "bern_lin")

# 4.5. Beta(2, 1) linear
# --------------------
final_display(R, data_function = beta21_lin, 
              B = params$B, eta = params$eta, gamma = params$gamma,
              K = K, lambda = lambda,
              theory_function = beta21_lin_theory, label = "beta21_lin")

# 4.6. Beta(1/2, 1/2) linear
# --------------------
final_display(R, data_function = beta05_lin, 
              B = params$B, eta = params$eta, gamma = params$gamma,
              K = K, lambda = lambda,
              theory_function = beta05_lin_theory, label = "beta05_lin")

# 4.7. Beta(2, 1) non-linear
# -------------------------
final_display(R, data_function = beta21_non_lin, 
              B = params$B, eta = params$eta, gamma = params$gamma,
              K = K, lambda = lambda,
              theory_function = beta21_non_lin_theory, label = "beta21_non_lin")
