#--------------------------------------#
# Sample from the middle: Graphical Heuristics
#--------------------------------------#

## Description: 
# This code creates a visual graph for showing the necessity of sampling
# from "the middle", a pool of (clearly) suboptimal policies
# to identify the optimal policy

# 0. Libraries
#--------------
library(dplyr)
library(ggplot2)

# 1. Formula
#-------------
# U = 1
# f_v^eps = (a, b*(1+epsilon), b*(1-epsilon), 1 - a - 2b)
# Supp_v = {0, 1/4, 1/2, 3/4}

s_x = function(x, a, b, epsilon, lambda){
  a*(1*(x>=0)*(1 - x + lambda*(x-0))) +
    b*(1+epsilon)*(1*(x>= 1/4)*(1-x + lambda*(x - 1/4))) +
      b*(1-epsilon)*(1*(x>= 1/2)*(1-x + lambda*(x - 1/2))) +
          (1-a-2*b)*(1*(x>= 3/4)*(1-x + lambda*(x - 3/4)))
}

# 2. Set parameters
#------------------
lower_bound_params = function(lambda){
  b = (1 - lambda) / (48 - 34*lambda) # 1/2 of b needed to set c_2 = 0
  a = (3*b*lambda + 1) / (4 - 3*lambda)
  list(a = a, b = b)
}

# 3. Display functions
#--------------------
wrap_lower_bound = function(lambda, x, epsilon){
  params = lower_bound_params(lambda)
  welfare = s_x(x, a = params$a, b = params$b, epsilon, lambda = lambda)
  list("welfare" = welfare, "epsilon" = epsilon, "x" = x)
}

grid = expand.grid(x = seq(0, 1, by = 0.05),
                   epsilon = c(1, -1))
lower_bound_tibble = mapply(wrap_lower_bound,
                          lambda = 0.7,
                          x = grid$x,
                          epsilon = grid$epsilon) |> t() |>
                     apply(2, unlist) |> as_tibble() |>
                     mutate(epsilon = as.factor(epsilon))

# 4. Plotting
#---------------
# Identify max by epsilon
maxs = lower_bound_tibble |> group_by(epsilon) |> summarize(max = max(welfare))

lb_plot = lower_bound_tibble |> ggplot() +
          geom_line(aes(x = x, y = welfare, color = epsilon), size = 1) +
          scale_color_manual(values = alpha(c("black", "coral1"))) +
          theme_minimal() +
          geom_hline(aes(yintercept = maxs$max[1]), color = "black",
                     size = 0.75, linetype = "dashed") +
          geom_hline(aes(yintercept = maxs$max[2]), color = "coral1",
                     size = 0.75, linetype = "dashed") +
          ylab("Welfare S(x)") +
          xlab("Action Space x")
  
  ggsave(lb_plot, filename="../plots/lower_bound.jpeg", bg = "white")

