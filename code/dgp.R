#--------------------------------------#
# DGP for Algorithm 1
#--------------------------------------#

## Description: 
# This code generates data according to the different DGP in Gonzalez (2023)
# Additionally, it generates a theory prescription for the Expected Welfare
# according to the particular DGP
# These functions are then evaluated in "algorithm.R"

set.seed(1234)

# 1. Uniform U with linear V = 1/2 * u
#-------------------------------------

unif_lin_theory = function(x){
  ifelse(x<1/2, lambda*x^2, 1/2 - x + lambda*x - lambda/4)
}  

unif_lin = function(K){
  u = runif(K)
  v = 1/2 * u
  tibble("u" = u, "v" = v,
         "x_optim" = 1/2)
}

# 2. Uniform U with non-linear V = * u^2
#---------------------------------------

unif_non_lin_theory = function(x){
 x^(3/2)*(2*lambda - 1)/2 - lambda*(x^2)/3
} 

unif_non_lin = function(K){
  u = runif(K)
  v = u^2
  tibble("u" = u, "v" = v,
         "x_optim" = 0.316)
}

# 3. Uniform U, Uniform V
#-------------------------

unif_unif_theory = function(x){
  -x^2*(1- lambda / 2) + x/2
} 

unif_unif = function(K){
  u = runif(K)
  v = runif(K)
  tibble("u" = u, "v" = v,
         "x_optim" = 1 / (4-2*lambda))
}

# 4. Bern(p) U with linear V = 1/2*u
#-----------------------------------
p = 0.7 # Set probability
bern_lin_theory = function(x, prob = p){
  ifelse(x > 1/2, (prob + lambda)/2 - (1-lambda)*x, -(1-prob)*(1-lambda)*x)
} 

bern_lin = function(K, prob = p){
  u = rbinom(K, size = 1, prob)
  v = u/2
  tibble("u" = u, "v" = v,
         "x_optim" = ifelse(p>=(1-2*lambda), 1/2, 0))
}

# 5. Beta(2,1) U with linear V = 1/2 * u
#---------------------------------------

beta21_lin_theory = function(x, alpha = 2, beta = 1){
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
         "x_optim" = 1/2)
}

# 6. Beta(1/2,1/2) U with linear V = 1/2 * u
#-------------------------------------------

beta05_lin_theory = function(x, alpha = 1/2, beta = 1/2){
  I_2x = 1 - 2/pi *acos(sqrt(2*x))
  I_2x_u = -1/pi * sqrt(2-4*x)*sqrt(x) + 1/pi*asin(sqrt(2*x)) 
  ifelse(x<1/2, I_2x*(I_2x_u - x) + lambda*I_2x*(x - (1/2) * I_2x_u), 
         alpha / (alpha + beta) - x + lambda*(x-(1/2)* alpha / (alpha + beta)))
}  

beta05_lin = function(K, alpha = 1/2, beta = 1/2){
  u = rbeta(K, alpha, beta)
  v = 1/2 * u
  tibble("u" = u, "v" = v,
         "x_optim" = 1/2)
}

# 7. Beta(2,1) U with non-linear V = = u^2
#-----------------------------------------

beta21_non_lin_theory = function(x, alpha = 2, beta = 1){

  x*((2/3) * x^(3/2) - x) + lambda*x*(x - x^2 / 2)
}  

beta21_non_lin = function(K, alpha = 2, beta = 1){
  u = rbeta(K, alpha, beta)
  v = u^2
  tibble("u" = u, "v" = v,
         "x_optim" = 1)
}