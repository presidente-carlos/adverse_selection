#--------------------------------------#
# DGP for Algorithm 1
#--------------------------------------#

## Description: 
# This code generates data according to the different DGP in Gonzalez (2023)
# Additionally, it generates a theory prescription for the Expected Welfare
# according to the particular DGP
# These functions are then evaluated in "algorithm.R"

set.seed(1234)
library(pracma) # For erf() access

# 1. Uniform U with linear V = 1/2 * u
#-------------------------------------

unif_lin_theory = function(x, lambda){
  ifelse(x<1/2, lambda*x^2, 1/2 - x + lambda*x - lambda/4)
}  

unif_lin = function(K){
  u = runif(K)
  v = 1/2 * u
  tibble("u" = u, "v" = v,
         "x_optim" = 1/2,
         "prob_work" = 1)
}

unif_lin_random_theory = function(x, lambda, sd){
  -1/20 *(1-lambda)*(2*x-1)*erf(x/(sqrt(2)*sd)) + 
    lambda*sd*(exp(-x^2/(2*sd^2))-1) / (5*sqrt(2*pi))
} 

unif_lin_random = function(K, sd){
  # if(sd == 0){x_optim = 0}
  # if(sd == 0.2){x_optim = 0.1125}
  if(sd == 0.4){x_optim = 0.1146}
  if(sd == 0.6){x_optim = 0.1151}
  if(sd == 0.8){x_optim = 0.1152}
  if(sd == 1){x_optim = 0.1153}
  
  print(x_optim)
  
  u = runif(K, min = 0.4, max = 0.6)
  v = u + rnorm(K, mean = 0, sd = sd)
  v = apply(t(v), 2, function(x){max(x, 0)})
  v = apply(t(v), 2, function(x){min(x, 1)})
  tibble("u" = u, "v" = v,
         "x_optim" = x_optim,
         "prob_work" = 0)
}

# 2. Uniform U with non-linear V = * u^2
#---------------------------------------

unif_non_lin_theory = function(x, lambda){
 x/2 - x^(3/2) + 2/3 * lambda * x^(3/2)
} 

unif_non_lin = function(K){
  u = runif(K)
  v = u^2
  tibble("u" = u, "v" = v,
         "x_optim" = 0.3906,
         "prob_work" = sqrt(0.3906))
}

# 3. Uniform U, Uniform V
#-------------------------

unif_unif_theory = function(x, lambda){
  -x^2*(1- lambda / 2) + x/2
} 

unif_unif = function(K){
  u = runif(K)
  v = runif(K)
  tibble("u" = u, "v" = v,
         "x_optim" = 1 / (4-2*lambda),
         "prob_work" = 1 / (4-2*lambda))
}

# 4. Bern(p) U with linear V = 1/2*u
#-----------------------------------
p = 0.7 # Set probability
bern_lin_theory = function(x, prob = p, lambda){
  ifelse(x > 1/2, prob*(1-lambda/2) - (1-lambda)*x, -(1-prob)*(1-lambda)*x)
} 

bern_lin = function(K, prob = p){
  u = rbinom(K, size = 1, prob)
  v = u/2
  tibble("u" = u, "v" = v,
         "x_optim" = ifelse(p>=(1-2*lambda), 1/2, 0),
         "prob_work" = ifelse(p>=(1-2*lambda), 1, 1-p))
}

# 5. Beta(2,1) U with linear V = 1/2 * u
#---------------------------------------

beta21_lin_theory = function(x, alpha = 2, beta = 1, lambda){
  I_2x = 4*x^2
  I_2x_u = 16*(x^3) / 3
  ifelse(x<=1/2, I_2x*(I_2x_u - x) + lambda*I_2x*(x - (1/2) * I_2x_u), 
         (alpha / (alpha + beta)) - x + 
           lambda*(x-(1/2)* (alpha / (alpha + beta))))
}  

beta21_lin = function(K, alpha = 2, beta = 1){
  u = rbeta(K, alpha, beta)
  v = 1/2 * u
  tibble("u" = u, "v" = v,
         "x_optim" = 1/2,
         "prob_work" = 1)
}

# 6. Beta(1/2,1/2) U with linear V = 1/2 * u
#-------------------------------------------

beta05_lin_theory = function(x, alpha = 1/2, beta = 1/2, lambda){
  I_2x = 1 - 2/pi *acos(sqrt(2*x))
  I_2x_u = -1/pi * sqrt(2-4*x)*sqrt(x) + 1/pi*asin(sqrt(2*x)) 
  ifelse(x<1/2, I_2x*(I_2x_u - x) + lambda*I_2x*(x - (1/2) * I_2x_u), 
         alpha / (alpha + beta) - x + lambda*(x-(1/2)* alpha / (alpha + beta)))
}  

beta05_lin = function(K, alpha = 1/2, beta = 1/2){
  u = rbeta(K, alpha, beta)
  v = 1/2 * u
  tibble("u" = u, "v" = v,
         "x_optim" = 1/2,
         "prob_work" = 1)
}

# 7. Beta(2,1) U with non-linear V = = u^2
#-----------------------------------------

beta21_non_lin_theory = function(x, alpha = 2, beta = 1, lambda){
  2/3 * x^(3/2) -x^2 + lambda*x^2 / 2
}  

beta21_non_lin = function(K, alpha = 2, beta = 1){
  u = rbeta(K, alpha, beta)
  v = u^2
  tibble("u" = u, "v" = v,
         "x_optim" = 0.592,
         "prob_work" = 0.592)
}

# 8. MW Policy DGP
#------------------

# 1. Uniform U with linear V = 1/2 * u OR V = 1/4*u
#--------------------------------------------------

#Checked
unif_lin_theory_high_u_low_v = function(x, lambda){
  ifelse(x<1/8, 0, ifelse(x>=1/8 & x<1/4, -(1-lambda)*(8*x^2 - x) + 
                            (16*x^2 - 1/4)*(1-lambda/4), 
                          -(1-lambda)*x+3/4*(1-lambda/4)))
}

# Checked
unif_lin_theory_high_u_high_v = function(x, lambda){
  ifelse(x<1/4, 0, ifelse(x>=1/4 & x<1/2, 2*lambda*x^2 - (1-lambda)*x +
                            lambda/8 - 1/4,
                          -(1-lambda)*x + 3/4 * (1 - lambda/2)))
}

# Checked
unif_lin_theory_low_u_low_v = function(x, lambda){
  ifelse(x<1/16, 0, ifelse(x>=1/16 & x < 3/16, (16*x^2 - 1/16)*(1-lambda/4) -
                             2*x*(1-lambda)*(4*x-1/4),
                           -(1-lambda)*x + 1/2 - lambda/8))
}

# Checked
unif_lin_theory_low_u_high_v = function(x, lambda){
  ifelse(x<1/8, 0, ifelse(x>=1/8 & x < 3/8, (1-lambda/2)*(4*x^2-1/16) - 
                            (1-lambda)*(4*x^2-x/2),
                          -(1-lambda)*x + 1/2 - lambda/4))
}

# Checked
unif_lin_high_u_low_v = function(K){
  u = runif(K, 1/2, 1)
  v = 1/4 * u
  tibble("u" = u, "v" = v,
         "x_optim" = 1/4,
         "prob_work" = 1)
}

# Checked
unif_lin_high_u_high_v = function(K){
  u = runif(K, 1/2, 1)
  v = 1/2 * u
  tibble("u" = u, "v" = v,
         "x_optim" = 1/2,
         "prob_work" = 1)
}

#Checked
unif_lin_low_u_low_v = function(K){
  u = runif(K, 1/4, 3/4)
  v = 1/4 * u
  tibble("u" = u, "v" = v,
         "x_optim" = 3/16,
         "prob_work" = 1)
}

# Checked
unif_lin_low_u_high_v = function(K){
  u = runif(K, 1/4, 3/4)
  v = 1/2 * u
  tibble("u" = u, "v" = v,
         "x_optim" = 3/8,
         "prob_work" = 1)
}

#9. Productivity shock
#---------------------
# Setting lambda = 0.35
# Using 1,500 periods

# The thing here is that the optimal policy from the beginning coincides
# with the optimal policy after the change
# it would be great if the optimal policy ex-post would be different from the optimal policy ex-ante
# and the optimal policy ONLY ex-post

# Presentation idea: Present both overlaying with translucent colors
# Both change and no-change have the same theory function

unif_lin_theory_shock_first = function(x, lambda){
  ifelse(x<1/4, 0, ifelse(x>=1/4 & x<1/2, (2- lambda)/2 * (2*x^2 - 1/8)
                          - (1-lambda)*(2*x^2 - x/2), 
                          (2 - lambda)/2 * 3/4 - (1-lambda)*x))
}

unif_lin_theory_shock_second_unch = function(x, lambda){
  ifelse(x<1/4, 0, ifelse(x>=1/4 & x<1/2, (2 - lambda)/2 * (2*x^2 - 1/8)
                          - (1-lambda)*(2*x^2 - x/2) - 1/2 * (2*x - 1/2), 
                          (2 - lambda)/2 * 3/4 - (1-lambda)*x - 1/2)) 
}

unif_lin_unchanged_fair = function(K){
  u_orig = runif(K, 1/2, 1)
  u_first = u_orig[1:(K/10)]
  u_second = u_orig[(K/10 + 1):K] - 1/4
  v = 1/2 * u_orig
  tibble("u" = c(u_first, u_second),
         "v" = v, x_optim = 0.25)
}

unif_lin_unchanged_unfair = function(K){
  u_orig = runif(K, 1/2, 1)
  u_first = u_orig[1:(K/10)]
  u_second = u_orig[(K/10 + 1):K] - 1/4
  v = 1/2 * u_orig
  tibble("u" = c(u_first, u_second),
         "v" = v, 
         x_optim = c(rep(0.5, K/10), rep(0.25, 9*K/10)))
}

unif_lin_theory_shock_second_ch = function(x, lambda){
  ifelse(x<1/4, lambda * x^2, (2 - lambda) / 2 * 1/4 - (1- lambda)*x)
}

unif_lin_changed_fair = function(K){
  u_first = runif(K/10, 1/2, 1)
  u_second = runif(9*K/10, 0, 1/2)
  v = c(1/2 * u_first, 1/2 * u_second)
  tibble("u" = c(u_first, u_second),
         "v" = v, x_optim = 0.25)
}

unif_lin_changed_unfair = function(K){
  u_first = runif(K/10, 1/2, 1)
  u_second = runif(9*K/10, 0, 1/2)
  v = c(1/2 * u_first, 1/2 * u_second)
  tibble("u" = c(u_first, u_second),
         "v" = v, x_optim = c(rep(0.5, K/10), rep(0.25, 9*K/10)))
}
