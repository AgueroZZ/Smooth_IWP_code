#### Functions for order selection section:

### Computation of reference marginal SD prior and prior conversion

### Input order and x vector, output covariance matrix of the 
### corresponding OS model:
compute_Sigma_OS <- function(order, x){
  designD <- compute_weights_design(x, p = order)[,-c(1:order)]
  Weight_cov <- solve(compute_weights_precision(svec = x))
  designD %*% Weight_cov %*% t(designD)
}

### Input covariance matrix, output reference marginal SD 
Compute_Ref_Margin <- function(Sigma){
  exp(mean(0.5*log(diag(Sigma))))
}

### Input marginal SD prior and reference marginal SD, output 
### corresponding prior for smoothing parameter
prior_conversion <- function(marginal, ref){
  alpha <- 1 - marginal$alpha 
  u <- marginal$u/ref
  list(u = u, alpha = alpha)
}


### Convert PC prior for order p to PC prior for order q: 
### by looking at the geometric average of marginal variance
prior_order_conversion_rue <- function(prior_p, x, p, q){
  Sigma_p <- compute_Sigma_OS(order = p, x = x)
  ref_p <- Compute_Ref_Margin(Sigma = Sigma_p)
  Sigma_q <- compute_Sigma_OS(order = q, x = x)
  ref_q <- Compute_Ref_Margin(Sigma = Sigma_q)
  prior_q <- list(alpha = prior_p$alpha, u = (prior_p$u * (ref_p /ref_q)))
  prior_q
}


### Convert PC prior for d-step-prediction variance to 
### smoothing variance at order p:
prior_order_conversion_predictive <- function(d, prior, p){
  Cp <- (d^((2*p) - 1))/(((2*p)-1)*(factorial(p-1)^2))
  prior_q <- list(alpha = prior$alpha, u = (prior$u * (1 /sqrt(Cp))))
  prior_q
}














