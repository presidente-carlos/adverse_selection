#--------------------------------------#
# Algorithm 1: Functions
#--------------------------------------#

## Description: 
# This code prepares a series of functions including:
# - Tempered Exp3 (Algorithm 1)
# - Auxiliary (optimal) parameter generating functions
# - Data preparation and plotting functions

# 1. Parameters Generator Function
# --------------------------------

param_gen = function(K){
  c_1 = 1/34 * ( K / log(K) )^(1/3)
  gamma = c_1 * (log(K) / K)^(1/3)
  c_3 = 10 * gamma 
  c_2 =  c_1 / (c_1 * c_3 + (K / log(K))^(1/3))
  B = c_3 / gamma
  eta = gamma / (B+1)
  print(paste0("eta = ", eta,
               " B = ", B,
               " gamma = ", gamma))
  tibble("gamma" = gamma, "eta" = eta, "B" = B)
}

# 2. Algorithm 1
#---------------

# Tempered Exp3 for Monopolist Wage Setting

exp3_monop = function(B, eta, gamma, K, lambda, u, v, x_optim,
                      K_store){
  
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
    if ((i-1)%%K_store == 0){
      t = round(i / K_store)
      store_pib[(t*length(x_disc) + 1):(t*length(x_disc) +
                                          length(x_disc)), ]$p_ib = p_i
      store_pib[(t*length(x_disc) + 1):(t*length(x_disc) +
                                          length(x_disc)), ]$round = i
    }
    
  }
  
  # Return outcomes
  list("regret" = store$avg_cum_regret, "store_pib" = store_pib$p_ib,
       "indicator" = store$indic)

}

# 3. Data preparation and Plotting
#---------------------------------

plotting = function(results, K, expression, prob_work, label, 
                    K_store, lambda, eta){
  
  cum_regret = do.call(rbind, results[1,])
  
  # Avg cum regret (and additional statistics)
  mod_quantile = function(matrix){
    quantile(matrix, probs = c(.05, .25, .5, .75, .95))
  }
  cum_regret_quantiles = apply(cum_regret, 2, 
                               FUN = mod_quantile) |> t() |> as_tibble()
  colnames(cum_regret_quantiles) = c("p5", "p25", "p50", "p75", "p95")
  cum_regret_tibble = tibble(x = seq(1:K),
                             exp_regret = colMeans(cum_regret),
                             var_regret = colVars(cum_regret)) |> 
                      cbind(cum_regret_quantiles)
  
  # Plot average regret
  avg_cum_regret = cum_regret_tibble |> ggplot() +
                   geom_smooth(aes(x = x, y = exp_regret), se = FALSE,
                               color = "black") +
                   theme_minimal(base_size = 18) +
                   xlab("Period K") +
                   ylab("Average Cummulative Regret")
  
  # Plot regret additional percentiles
  extra_cum_regret = cum_regret_tibble |> ggplot() +
                   geom_smooth(aes(x = x, y = p50), se = FALSE,
                               color = "black") +
                   geom_smooth(aes(x = x, y = p25), se = FALSE, color = "blue") +
                   geom_smooth(aes(x = x, y = p75), se = FALSE, color = "blue") +
                   geom_smooth(aes(x = x, y = p5), se = FALSE, color = "red") +
                   geom_smooth(aes(x = x, y = p95), se = FALSE, color = "red") +
                   theme_minimal(base_size = 18) +
                   xlab("Period K") +
                   ylab("Cummulative Regret Percentiles")
  
  # Plot regret variance
  var_cum_regret = cum_regret_tibble |> ggplot() +
                   geom_smooth(aes(x = x, y = var_regret), se = FALSE,
                               color = "black") +
                   theme_minimal(base_size = 18) +
                   xlab("Period K") +
                   ylab("Cummulative Regret Variance")
  
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
                              ylab("Log Probability") +
                              scale_y_continuous(labels = NULL, breaks = NULL) +
                              facet_grid(~K_round)
    t = t+1
    k_minus = k
  }

  # Theoretical S_i(x)
  x_seq_plot = seq(0, 1, by = 0.01)
  plot_tib = tibble("x_disc" = x_seq_plot, 
                    "s_ib" = exp(expression(x_seq_plot, lambda = lambda)))
  # Downward shift for plotting considerations
  
  actual = plot_tib |> ggplot() +
                       geom_area(aes(x = x_disc,
                                     y = s_ib), 
                                     color = "black",
                                     fill = "coral2", size = 1) +
                       theme_minimal(base_size = 18) +
                       xlab("Action Space x") +
                       ylab("Expected Welfare S(x)")
  
  # Empirical probability of work
  indicator = do.call(rbind, results[3,])
  store_indicator_MA = matrix(NA, nrow = nrow(indicator), ncol = ncol(indicator))

  # Create Moving Averages (50 periods)
  for (i in 1:nrow(indicator)){
    for(j in 1:ncol(indicator)){
      store_indicator_MA[i, j] = mean(indicator[i:max((i-50),0), j])
    }
  }
  
  # Avg empirical probability of work
  indicator_tibble = tibble(x = seq(1:K),
                             y = colMeans(store_indicator_MA))
  indicator_plot = indicator_tibble |> ggplot() +
    geom_smooth(aes(x = x, y = y), se = FALSE,
                color = "black") +
    geom_hline(aes(yintercept = prob_work), color = "red", size = 1) +
    theme_minimal(base_size = 14) +
    xlab("Period K") +
    ylab("Moving Average (50 lags) - Probability of working")
  
  plot_list = list("avg_cum_regret" = plot(avg_cum_regret), 
                   "theory_welfare" = plot(actual), 
                   "empirical_probs" = do.call(grid.arrange, c(P, nrow = 3)),
                   "prob_work" = plot(indicator_plot),
                   "var_cum_regret" = plot(var_cum_regret),
                   "extra_cum_regret" = plot(extra_cum_regret))
  
  eta_label = round(eta, 3)
  lapply(names(plot_list), 
         function(x)ggsave(filename=paste("../plots/", label, "/",
                                           label, "_", x, eta_label,
                                          ".jpeg",sep=""), 
                           plot = plot_list[[x]], bg = "white"))
  
}

final_display = function(R, data_function, B, eta, 
                         gamma, K, lambda,
                         theory_function, label, K_store = 50){
  
  simulation_data = data_function(K)
  
  # Prepare session for parallel computing
  plan(multisession(workers = 4))
  
  results = future_replicate(R, exp3_monop(B = B, eta = eta, 
                                           gamma = gamma,
                                           K = K, lambda = lambda,
                                           u = simulation_data$u,
                                           v = simulation_data$v, 
                                           x_optim = simulation_data$x_optim,
                                           K_store = K_store))
  
  plotting(results = results, K = K, 
           expression = theory_function, prob_work = simulation_data$prob_work,
           label = label, K_store = K_store, lambda = lambda, eta = eta)
  
}

final_display_eta_loop = function(R, data_function, B = params$B, 
                                  gamma = params$gamma, K = K, lambda = lambda,
                                  theory_function, label, K_store = 50, eta_vector){
  
  future_sapply(X = eta_vector, FUN = final_display, 
                R = R, data_function = data_function, B = B, 
                gamma = gamma, K = K, lambda = lambda,
                theory_function = theory_function, label = label, 
                K_store = K_store, future.seed = T)
}

log_sequence = function(from, to, length.out){
  steps = c()
  for (i in 1:length.out){
    steps[i] = (to - from) / 2^i
  }
  rev(steps)
}