#--------------------------------------#
# Robustness Checks and Policy Functions
#--------------------------------------#

## Description: 
# This code contains additional functions following algorithm.R:
# - For variations across parameter values (eta)
# - For alternative sd in DGP
# - For variations in the policy space (MW policies)
# - For variations in the Inf Processing Capacity of firms
# - For productivity shocks

# 1. SD/eta robustness
#---------------------

# 1a) Multi-plotting sd/eta
#--------------------------
multi_plotting_sd = function(multi_results, K, expression, prob_work, label, 
                             K_store, lambda, eta, sd_vector){
  
  # Output multi_results recovers a tibble / list with length(sd_vector) columns
  # and K * 3 rows where rows include vectors of regrets, emp probs and prob of work 
  
  # Initialize store vectors
  cum_regret_store = c()
  indicator_store = c()
  
  for (iter in 1:ncol(multi_results)){
    
    # Average Regret
    index_regret = ((1:nrow(multi_results)) - 1)%%3 == 0
    cum_regret = do.call(rbind, multi_results[index_regret, iter]) |> colMeans()
    cum_regret_store = c(cum_regret_store, cum_regret)
    
    # Probability of work
    index_work = (1:nrow(multi_results))%%3 == 0
    indicator = do.call(rbind, multi_results[index_work, iter])
    
    # Create Moving Averages (100 periods) for prob of work
    store_indicator_MA = matrix(NA, nrow = nrow(indicator), ncol = ncol(indicator))
    for(i in 1:nrow(indicator)){
      for (j in 1:ncol(indicator)){
        store_indicator_MA[i, j] = mean(indicator[i, max((j-100), 1):j])
      }
    }
    
    indicator_store = c(indicator_store, colMeans(store_indicator_MA))
  }
  
  # Tibble for cum_regret
  cum_regret_tibble = tibble(x = rep(seq(1:K), ncol(multi_results)),
                             exp_regret = cum_regret_store,
                             sd = rep(sd_vector, each = K))
  
  # Plot average regret
  avg_cum_regret = cum_regret_tibble |> ggplot() +
                   geom_smooth(aes(x = x, y = exp_regret, 
                                   color = as.factor(sd)), se = FALSE) +
                    theme_minimal(base_size = 18) +
                    xlab("Period K") +
                    ylab("Average Cummulative Regret") +
                    theme(axis.title=element_text(size=16)) +
                    labs(color = TeX("$$\\sigma$$"))
  
  
  # Tibble for work_probability
  indicator_tibble = tibble(x = rep(seq(1:K), ncol(multi_results)),
                            y = indicator_store,
                            sd = rep(sd_vector, each = K))
  # Plot prob of work
  indicator_plot = indicator_tibble |> ggplot() +
                   geom_smooth(aes(x = x, y = y, 
                                   color = as.factor(sd)), se = FALSE) +
                    # geom_hline(aes(yintercept = prob_work), 
                    # color = "red", linewidth = 1) +
                    theme_minimal(base_size = 14) +
                    xlab("Period K") +
                    ylab("Probability of working - MA 100 periods") +
                    theme(axis.title=element_text(size=16)) +
                    labs(color = TeX("$$\\sigma$$"))
  
  # Plot cum_regret and prob_work
  plot_list = list("avg_cum_regret" = plot(avg_cum_regret),
                   "prob_work" = plot(indicator_plot))
  
  eta_label = round(eta, 2)
  # sd_label = paste0("sd", sd)
  # mw_label = paste0("mw", mw)
  
  lapply(names(plot_list), 
         function(x)ggsave(filename=paste("../plots/", label, "/",
                                          label, "_", x, eta_label,
                                          ".jpeg",sep=""), 
                           plot = plot_list[[x]], bg = "white",
                           width = 8, height = 6))
  
  # p_ib across K
  P = list() # Storing graphs along different sets of periods
  probability_plot_list = list() # Storing along different sd values
  
  for (iter in 1:ncol(multi_results)){
    
    # Initialize auxiliary objects for plotting
    P = list()
    t = 1
    k_minus = 0
    
    # Recover empirical probabilities
    index_prob = ((1:nrow(multi_results)) - 2)%%3 == 0
    avg_p_ib_across_K = do.call(rbind, 
                                multi_results[index_prob, iter]) |> colMeans()
    
    # Tibble for emp probs
    b_seq = 1:(params$B+1)
    x_disc = (b_seq - 1)/(params$B)
    times_arg = K/K_store
    p_ib_plotting = tibble("p_ib" = avg_p_ib_across_K,
                           "K_round" = rep(seq(from = 1, to = K, by = K_store), 
                                           each = length(x_disc)),
                           "x_disc" = rep(x_disc, times_arg))
    
    # Plot in groups of three (for display considerations)
    for (k in seq(K_store*3 + 1, K, by = K_store*3)){
      P[[t]] = p_ib_plotting |> filter(K_round <= k, K_round >= k_minus) |>
        ggplot() + geom_area(aes(x = x_disc, y = p_ib), color = "black", 
                             fill = "coral") +
        theme_minimal() +
        xlab("Action Space x") +
        ylab("Probability") +
        scale_y_continuous(labels = NULL, breaks = NULL) +
        facet_grid(~K_round) +
        theme(axis.title=element_text(size=16))
      
    # Update auxiliary objects  
      t = t+1
      k_minus = k
    }
    
    # Plot and save output
    probability_plot_list = list(do.call(grid.arrange, c(P, nrow = 3)))
    label_emp_probs = paste0("empirical_probs", sd_vector[iter])
    ggsave(filename = paste("../plots/", label, "/",
                            label, "_", label_emp_probs, ".jpeg",sep=""),
           plot = plot(probability_plot_list[[1]]), 
           bg = "white", width = 8, height = 6)
  }  
  
}

# Theoretical S_i(x)
multi_plotting_theory_sd = function(expression, lambda, sd, label){
  
  library(pracma) # For erf() access
  x_seq_plot = seq(0, 1, by = 0.01)
  s_ib = expression(x_seq_plot, lambda = lambda, sd = sd) |> exp()
  plot_tib = tibble("x_disc" = x_seq_plot, 
                    "s_ib" = s_ib)

  # Downward shift for plotting considerations
  min_sib = min(plot_tib$s_ib)
  plot_tib$s_ib = plot_tib$s_ib - min_sib
  
  actual = plot_tib |> ggplot() + geom_area(aes(x = x_disc, y = s_ib), 
                                            color = "black", fill = "coral2", linewidth = 1) +
                       theme_minimal() +
                       xlab("Action Space x") +
                       ylab("Expected Welfare S(x)") +
                       theme(axis.title=element_text(size=16),
                             axis.text.y=element_blank())
  
  ggsave(filename = paste("../plots/", label, "/",
                          label, "_", "theory_welfare_sd", sd, ".jpeg",sep=""),
         plot = plot(actual), bg = "white", width = 8, height = 6)
}

# 1b) Replication function (+ data generation)
#-----------------------------------------
rep_function_sd = function(R, B, data_function, 
                           eta, gamma, K, lambda, K_store, sd){
  
  # Simulate data
  simulation_data = data_function(K, sd)
  
  # Prepare session for parallel computing
  plan(multisession(workers = 6))
  
  future_replicate(R, exp3_monop(B = B, eta = eta, 
                                 gamma = gamma,
                                 K = K, lambda = lambda,
                                 u = simulation_data$u,
                                 v = simulation_data$v, 
                                 x_optim = simulation_data$x_optim,
                                 K_store = K_store))
}

# 1c) Final Display sd/eta
#--------------------------

multi_final_display_sd = function(R, B, data_function, eta, 
                                  gamma, K, lambda, theory_function, label, K_store = 50, 
                                  sd_vector){
  
  # Algorithm 1 (replicated)
  multi_results = future_sapply(X = sd_vector, FUN = rep_function_sd,
                                R = R, B = B, data_function = data_function,
                                eta = eta, gamma = gamma, K = K, lambda = lambda,
                                K_store = K_store)
  
  # Plotting and saving
  multi_plotting_sd(multi_results = multi_results, K = K,
                    expression = theory_function, prob_work = 0,
                    label = label, K_store = K_store, lambda = lambda, eta = eta,
                    sd_vector = sd_vector)
  
  # Plotting and saving theory function
  future_sapply(sd_vector, FUN = multi_plotting_theory_sd,
                expression = theory_function, lambda = lambda, label = label)
  
}

# 2. Minimum Wage Policy Analysis
#---------------------------------

multi_plotting_mw = function(multi_results, K, expression, prob_work, label, 
                             K_store, lambda, eta, mw_vector, x_optim){
  
  # Output multi_results recovers a tibble / list with length(mw_vector) columns
  # and K * 3 rows where rows include vectors of regrets, emp probs and prob of work 
  
  # Initialize store vectors
  cum_regret_store = c()
  indicator_store = c()
  
  for (iter in 1:ncol(multi_results)){
    
    # Average Regret
    index_regret = ((1:nrow(multi_results)) - 1)%%3 == 0
    cum_regret = do.call(rbind, multi_results[index_regret, iter]) |> colMeans()
    cum_regret_store = c(cum_regret_store, cum_regret)
    
    # Probability of work
    index_work = (1:nrow(multi_results))%%3 == 0
    indicator = do.call(rbind, multi_results[index_work, iter])
    
    # Create Moving Averages (100 periods) for prob of work
    store_indicator_MA = matrix(NA, nrow = nrow(indicator), ncol = ncol(indicator))
    for(i in 1:nrow(indicator)){
      for (j in 1:ncol(indicator)){
        store_indicator_MA[i, j] = mean(indicator[i, max((j-100), 1):j])
      }
    }
    
    indicator_store = c(indicator_store, colMeans(store_indicator_MA))
  }
  
  # Tibble for cum_regret
  cum_regret_tibble = tibble(x = rep(seq(1:K), ncol(multi_results)),
                             exp_regret = cum_regret_store,
                             mw = rep(round(mw_vector, 2), each = K),
                             binding = ifelse(mw < x_optim, "No", "Yes"))
  
  # Plot average regret
  avg_cum_regret = cum_regret_tibble |> ggplot() +
    geom_smooth(aes(x = x, y = exp_regret, 
                    color = as.factor(mw), linetype = binding), se = FALSE) +
    theme_minimal(base_size = 18) +
    xlab("Period K") +
    ylab("Average Cummulative Regret") +
    theme(axis.title=element_text(size=16)) +
    labs(color = "Minimum Wage", linetype = "Binding?")
  
  
  # Tibble for work_probability
  indicator_tibble = tibble(x = rep(seq(1:K), ncol(multi_results)),
                            y = indicator_store,
                            mw = rep(round(mw_vector, 2), each = K),
                            binding = ifelse(mw < x_optim, "No", "Yes"))
  # Plot prob of work
  indicator_plot = indicator_tibble |> ggplot() +
    geom_smooth(aes(x = x, y = y, 
                    color = as.factor(mw), linetype = binding), se = FALSE) +
    # geom_hline(aes(yintercept = prob_work), 
    # color = "red", linewidth = 1) +
    theme_minimal(base_size = 14) +
    xlab("Period K") +
    ylab("Probability of working - MA 100 periods") +
    theme(axis.title=element_text(size=16)) +
    labs(color = "Minimum Wage", linetype = "Binding?")
  
  # Plot cum_regret and prob_work
  plot_list = list("avg_cum_regret" = plot(avg_cum_regret),
                   "prob_work" = plot(indicator_plot))
  
  # eta_label = round(eta, 2)
  # sd_label = paste0("sd", sd)
  # mw_label = paste0("mw", mw)
  
  lapply(names(plot_list), 
         function(x)ggsave(filename=paste("../plots/", label, "/",
                                          label, "_", x,
                                          ".jpeg",sep=""), 
                           plot = plot_list[[x]], bg = "white",
                           width = 8, height = 6))
  
  # p_ib across K
  P = list() # Storing graphs along different sets of periods
  probability_plot_list = list() # Storing along different sd values
  
  for (iter in 1:ncol(multi_results)){
    
    # Initialize auxiliary objects for plotting
    P = list()
    t = 1
    k_minus = 0
    
    # Recover empirical probabilities
    index_prob = ((1:nrow(multi_results)) - 2)%%3 == 0
    avg_p_ib_across_K = do.call(rbind, 
                                multi_results[index_prob, iter]) |> colMeans()
    
    # Tibble for emp probs
    b_seq = 1:(params$B+1)
    x_disc = (b_seq - 1)/(params$B)
    times_arg = K/K_store
    p_ib_plotting = tibble("p_ib" = avg_p_ib_across_K,
                           "K_round" = rep(seq(from = 1, to = K, by = K_store), 
                                           each = length(x_disc)),
                           "x_disc" = rep(x_disc, times_arg))
    
    # Plot in groups of three (for display considerations)
    for (k in seq(K_store*3 + 1, K, by = K_store*3)){
      P[[t]] = p_ib_plotting |> filter(K_round <= k, K_round >= k_minus) |>
        ggplot() + geom_area(aes(x = x_disc, y = p_ib), color = "black", 
                             fill = "coral") +
        theme_minimal() +
        xlab("Action Space x") +
        ylab("Probability") +
        scale_y_continuous(labels = NULL, breaks = NULL) +
        facet_grid(~K_round) +
        theme(axis.title=element_text(size=16))
      
      # Update auxiliary objects  
      t = t+1
      k_minus = k
    }
    
    # Plot and save output
    probability_plot_list = list(do.call(grid.arrange, c(P, nrow = 3)))
    label_emp_probs = paste0("empirical_probs", round(mw_vector[iter], 2))
    ggsave(filename = paste("../plots/", label, "/",
                            label, "_", label_emp_probs, ".jpeg",sep=""),
           plot = plot(probability_plot_list[[1]]), 
           bg = "white", width = 8, height = 6)
  }
  
  # Theoretical S_i(x)
  x_seq_plot = seq(0, 1, by = 0.01)
  s_ib = expression(x_seq_plot, lambda = lambda) |> exp()
  plot_tib = tibble("x_disc" = x_seq_plot, 
                    "s_ib" = s_ib)
  
  # Downward shift for plotting considerations
  min_sib = min(plot_tib$s_ib)
  plot_tib$s_ib = plot_tib$s_ib - min_sib
  
  actual = plot_tib |> ggplot() + geom_area(aes(x = x_disc, y = s_ib), 
                                            color = "black", fill = "coral2", linewidth = 1) +
    theme_minimal() +
    xlab("Action Space x") +
    ylab("Expected Welfare S(x)") +
    theme(axis.title=element_text(size=16),
          axis.text.y=element_blank())
  
  ggsave(filename = paste("../plots/", label, "/",
                          label, "_", "theory_welfare", ".jpeg",sep=""),
         plot = plot(actual), bg = "white", width = 8, height = 6)
  
}

# 2b) Replication function (+ data generation)
#-----------------------------------------
rep_function_mw = function(R, B, data_function, 
                           eta, gamma, K, lambda, K_store, u, v, x_optim, mw){
  
  # Prepare session for parallel computing
  plan(multisession(workers = 6))
  future_replicate(R, exp3_monop(B = B, eta = eta, 
                                 gamma = gamma,
                                 K = K, lambda = lambda,
                                 u = u,
                                 v = simulation_data$v, 
                                 x_optim = simulation_data$x_optim,
                                 K_store = K_store,
                                 mw = mw))
}

multi_final_display_mw = function(R, B, data_function, eta, 
                                  gamma, K, lambda, theory_function, label, K_store = 50, 
                                  mw_vector){
  
  # Simulate data
  simulation_data = data_function(K)
  
  # Algorithm 1 (replicated)
  multi_results = future_sapply(X = mw_vector, FUN = rep_function_mw,
                                R = R, B = B, data_function = data_function,
                                eta = eta, gamma = gamma, K = K, lambda = lambda,
                                K_store = K_store, u = simulation_data$u,
                                v = simulation_data$v, 
                                x_optim = simulation_data$x_optim)
  
  # Plotting and saving
  multi_plotting_mw(multi_results = multi_results, K = K,
                    expression = theory_function, prob_work = 0,
                    label = label, K_store = K_store, lambda = lambda, eta = eta,
                    mw_vector = mw_vector, x_optim = simulation_data$x_optim)
  
}

# 3. Limited Information Processing Capacity (LIPC)
#--------------------------------------------------

multi_plotting_lipc = function(multi_results, K, expression, prob_work, label, 
                             K_store, lambda, eta, p_theta_vector){
  
  # Output multi_results recovers a tibble / list with length(p_theta_vector) columns
  # and K * 3 rows where rows include vectors of regrets, emp probs and prob of work 
  
  # Initialize store vectors
  cum_regret_store = c()
  indicator_store = c()
  
  for (iter in 1:ncol(multi_results)){
    
    # Average Regret
    index_regret = ((1:nrow(multi_results)) - 1)%%3 == 0
    cum_regret = do.call(rbind, multi_results[index_regret, iter]) |> colMeans()
    cum_regret_store = c(cum_regret_store, cum_regret)
    
    # Probability of work
    index_work = (1:nrow(multi_results))%%3 == 0
    indicator = do.call(rbind, multi_results[index_work, iter])
    
    # Create Moving Averages (100 periods) for prob of work
    store_indicator_MA = matrix(NA, nrow = nrow(indicator), ncol = ncol(indicator))
    for(i in 1:nrow(indicator)){
      for (j in 1:ncol(indicator)){
        store_indicator_MA[i, j] = mean(indicator[i, max((j-100), 1):j])
      }
    }
    
    indicator_store = c(indicator_store, colMeans(store_indicator_MA))
  }
  
  # Tibble for cum_regret
  cum_regret_tibble = tibble(x = rep(seq(1:K), ncol(multi_results)),
                             exp_regret = cum_regret_store,
                             p_theta = rep(round(p_theta_vector, 2), each = K))
  
  # Plot average regret
  avg_cum_regret = cum_regret_tibble |> ggplot() +
    geom_smooth(aes(x = x, y = exp_regret, 
                    color = as.factor(p_theta)), se = FALSE) +
    theme_minimal(base_size = 18) +
    xlab("Period K") +
    ylab("Average Cummulative Regret") +
    theme(axis.title=element_text(size=16)) +
    labs(color = TeX("$$p_\\theta"))
  
  
  # Tibble for work_probability
  indicator_tibble = tibble(x = rep(seq(1:K), ncol(multi_results)),
                            y = indicator_store,
                            p_theta = rep(round(p_theta_vector, 2), each = K))
  # Plot prob of work
  indicator_plot = indicator_tibble |> ggplot() +
    geom_smooth(aes(x = x, y = y, 
                    color = as.factor(p_theta)), se = FALSE) +
    # geom_hline(aes(yintercept = prob_work), 
    # color = "red", linewidth = 1) +
    theme_minimal(base_size = 14) +
    xlab("Period K") +
    ylab("Probability of working - MA 100 periods") +
    theme(axis.title=element_text(size=16)) +
    labs(color = TeX("$$p_\\theta"))
  
  # Plot cum_regret and prob_work
  plot_list = list("avg_cum_regret" = plot(avg_cum_regret),
                   "prob_work" = plot(indicator_plot))
  
  lapply(names(plot_list), 
         function(x)ggsave(filename=paste("../plots/", label, "/",
                                          label, "_", x,
                                          ".jpeg",sep=""), 
                           plot = plot_list[[x]], bg = "white",
                           width = 8, height = 6))
  
  # p_ib across K
  P = list() # Storing graphs along different sets of periods
  probability_plot_list = list() # Storing along different sd values
  
  for (iter in 1:ncol(multi_results)){
    
    # Initialize auxiliary objects for plotting
    P = list()
    t = 1
    k_minus = 0
    
    # Recover empirical probabilities
    index_prob = ((1:nrow(multi_results)) - 2)%%3 == 0
    avg_p_ib_across_K = do.call(rbind, 
                                multi_results[index_prob, iter]) |> colMeans()
    
    # Tibble for emp probs
    b_seq = 1:(params$B+1)
    x_disc = (b_seq - 1)/(params$B)
    times_arg = K/K_store
    p_ib_plotting = tibble("p_ib" = avg_p_ib_across_K,
                           "K_round" = rep(seq(from = 1, to = K, by = K_store), 
                                           each = length(x_disc)),
                           "x_disc" = rep(x_disc, times_arg))
    
    # Plot in groups of three (for display considerations)
    for (k in seq(K_store*3 + 1, K, by = K_store*3)){
      P[[t]] = p_ib_plotting |> filter(K_round <= k, K_round >= k_minus) |>
        ggplot() + geom_area(aes(x = x_disc, y = p_ib), color = "black", 
                             fill = "coral") +
        theme_minimal() +
        xlab("Action Space x") +
        ylab("Probability") +
        scale_y_continuous(labels = NULL, breaks = NULL) +
        facet_grid(~K_round) +
        theme(axis.title=element_text(size=16))
      
      # Update auxiliary objects  
      t = t+1
      k_minus = k
    }
    
    # Plot and save output
    probability_plot_list = list(do.call(grid.arrange, c(P, nrow = 3)))
    label_emp_probs = paste0("empirical_probs", round(p_theta_vector[iter], 2))
    ggsave(filename = paste("../plots/", label, "/",
                            label, "_", label_emp_probs, ".jpeg",sep=""),
           plot = plot(probability_plot_list[[1]]), 
           bg = "white", width = 8, height = 6)
  }
  
  # Theoretical S_i(x)
  x_seq_plot = seq(0, 1, by = 0.01)
  s_ib = expression(x_seq_plot, lambda = lambda) |> exp()
  plot_tib = tibble("x_disc" = x_seq_plot, 
                    "s_ib" = s_ib)
  
  # Downward shift for plotting considerations
  min_sib = min(plot_tib$s_ib)
  plot_tib$s_ib = plot_tib$s_ib - min_sib
  
  actual = plot_tib |> ggplot() + geom_area(aes(x = x_disc, y = s_ib), 
                                            color = "black", fill = "coral2", linewidth = 1) +
    theme_minimal() +
    xlab("Action Space x") +
    ylab("Expected Welfare S(x)") +
    theme(axis.title=element_text(size=16),
          axis.text.y=element_blank())
  
  ggsave(filename = paste("../plots/", label, "/",
                          label, "_", "theory_welfare", ".jpeg",sep=""),
         plot = plot(actual), bg = "white", width = 8, height = 6)
  
}

# 3b) Replication function (+ data generation)
#-----------------------------------------
rep_function_lipc = function(R, B, data_function, 
                           eta, gamma, K, lambda, K_store, u, v, x_optim, p_theta){
  
  # Prepare session for parallel computing
  plan(multisession(workers = 6))
  future_replicate(R, exp3_monop(B = B, eta = eta, 
                                 gamma = gamma,
                                 K = K, lambda = lambda,
                                 u = u,
                                 v = simulation_data$v, 
                                 x_optim = simulation_data$x_optim,
                                 K_store = K_store,
                                 p_theta = p_theta))
}

multi_final_display_lipc = function(R, B, data_function, eta, 
                                  gamma, K, lambda, theory_function, label, K_store = 50, 
                                  p_theta_vector, mw = 0){
  
  # Simulate data
  simulation_data = data_function(K)
  
  # Algorithm 2 (replicated)
  multi_results = future_sapply(X = p_theta_vector, FUN = rep_function_lipc,
                                R = R, B = B, data_function = data_function,
                                eta = eta, gamma = gamma, K = K, lambda = lambda,
                                K_store = K_store, u = simulation_data$u,
                                v = simulation_data$v, 
                                x_optim = simulation_data$x_optim)
  
  # Plotting and saving
  multi_plotting_lipc(multi_results = multi_results, K = K,
                    expression = theory_function, prob_work = 0,
                    label = label, K_store = K_store, lambda = lambda, eta = eta,
                    p_theta_vector = p_theta_vector)
  
}

# 5. Productivity Shocks (preliminary)
#--------------------------------------

plotting_prod_shock = function(results, K, expression_1, expression_2,
                               label, K_store, lambda, eta, x_optim){
  
  cum_regret = do.call(rbind, results[1,])
  cum_regret_tibble = tibble(x = seq(1:K),
                             exp_regret = colMeans(cum_regret))
  
  # Plot average regret
  avg_cum_regret = cum_regret_tibble |> ggplot() +
    geom_smooth(aes(x = x, y = exp_regret), se = FALSE,
                color = "black") +
    theme_minimal(base_size = 18) +
    geom_vline(aes(xintercept = K/10), linetype = "dashed", color = "red") +
    geom_hline(aes(yintercept = mean(exp_regret)), linetype = "dashed", color = "black") +
    xlab("Period K") +
    ylab("Average Cummulative Regret") +
    theme(axis.title=element_text(size=16))
  
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
      geom_vline(aes(xintercept = x_optim[K_round]), color = "black") +
      theme_minimal() +
      xlab(TeX("Action Space $\\Omega_X$")) +
      ylab("Probability") +
      theme(axis.title=element_text(size=16)) +
      scale_y_continuous(labels = NULL, breaks = NULL) +
      facet_grid(~K_round)
    t = t+1
    k_minus = k
  }
  
  x_seq_plot = seq(0, 1, by = 0.01)
  
  # Theoretical \exp(S)_i(x)
  s_ib_1 = expression_1(x_seq_plot, lambda = lambda) |> exp()
  s_ib_2 = expression_2(x_seq_plot, lambda = lambda) |> exp()
  plot_tib = tibble("x_disc" = rep(x_seq_plot, 2), 
                    "s_ib" = c(s_ib_1, s_ib_2),
                    "first_second" = c(rep("Before Shock", length(x_seq_plot)),
                                        rep("After Shock", length(x_seq_plot))))
  
  # Downward shift for plotting considerations
  min_sib = min(plot_tib$s_ib)
  plot_tib$s_ib = plot_tib$s_ib - min_sib
  
  actual = plot_tib |> ggplot() +
    geom_area(aes(x = x_disc,
                  y = s_ib, fill = fct_reorder(first_second, s_ib, .desc = TRUE)), alpha = 0.5,
              position = "identity",
              color = "black", linewidth = 1) +
    theme_minimal(base_size = 18) +
    xlab(TeX("Action Space $\\Omega_X$")) +
    ylab(TeX("Expected Welfare $$\\exp(S(x))$$")) +
    theme(axis.title=element_text(size=16),
          axis.text.y=element_blank()) +
    labs(fill = "Shock Timing")
    
  
  # Final plotting and saving
  plot_list = list("avg_cum_regret" = plot(avg_cum_regret), 
                   "empirical_probs" = do.call(grid.arrange, c(P, nrow = 3)))
  
  lapply(names(plot_list), 
         function(x)ggsave(filename=paste("../plots/", label, "/",
                                          label, "_", x, ".jpeg",sep=""), 
                           plot = plot_list[[x]], bg = "white",
                           width = 8, height = 6))
  
  # Theory functions
  ggsave(filename = paste("../plots/", label, "/",
                          label, "_", "theory_welfare", ".jpeg",sep=""),
         plot = plot(actual),
         bg = "white", width = 8, height = 6)
  
}

# 5b. Final display productivity shocks
final_display_shocks = function(R, data_function, B, eta, 
                         gamma, K, lambda,
                         theory_function_1, theory_function_2, label, K_store = 100, 
                         sd = 0, mw = 0){
  
  simulation_data = data_function(K)
  
  # Prepare session for parallel computing
  plan(multisession(workers = 6))
  
  results = future_replicate(R, exp3_monop(B = B, eta = eta, 
                                           gamma = gamma,
                                           K = K, lambda = lambda,
                                           u = simulation_data$u,
                                           v = simulation_data$v, 
                                           x_optim = simulation_data$x_optim,
                                           K_store = K_store))
  
  plotting_prod_shock(results = results, K = K, 
           expression_1 = theory_function_1,
           expression_2 = theory_function_2,
           label = label, K_store = K_store, lambda = lambda, eta = eta,
           x_optim = simulation_data$x_optim)
  
}
