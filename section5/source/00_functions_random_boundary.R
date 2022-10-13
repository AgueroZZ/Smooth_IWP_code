### Functions:

### Implement WP with random boundary:
Imple_WP_rand_bound <- function(betaprec, sigmaE, theta, x, y, p = 3, optimizer = F){
  x_norm <- x - x[1]
  svec <- x_norm[-1] ## assume x[1] = 0
  B <- Compute_design_Aug(svec = svec, p = p)
  Q <- Compute_Aug_Wp_Prec(svec = svec,p = p)
  D_proposed <- compute_weights_design(x_norm, p)
  ## Design matrix for the global polynomial
  X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
  B <- rbind(rep(0, ncol(B)), B)
  tmbdat <- list(
    # Design matrix
    X = X,
    B = B,
    P = as(Q,'dgTMatrix'),
    logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
    # Response
    y = y,
    sigmaE = sigmaE,
    theta1 = theta,
    betaprec = betaprec
  )
  tmbparams <- list(
    W = c(rep(0, (ncol(Q) + ncol(X))))
  )
  
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    # random = "W",
    DLL = "00_Gaussian_Smoothing_unknown_bound",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  ## posterior mode and posterior precision matrix
  opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                control = list(eval.max = 200, iter.max = 200))
  prec_matrix <- forceSymmetric(ff$he(opt$par))
  
  if (optimizer == T){
    return(list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt))
  }
  
  ## Return posterior mean and posterior precision matrix
  return(list(mean = opt$par, prec = as.matrix(prec_matrix)))
}

### Implement WP with diffuse boundary:
Imple_WP_diff_bound <- function(sigmaE, theta, x, y, p = 3, optimizer = F){
  x_norm <- x - x[1]
  svec <- x_norm[-1] ## assume x[1] = 0
  B <- Compute_design_Aug(svec = svec, p = p)
  Q <- Compute_Aug_Wp_Prec(svec = svec,p = p)
  D_proposed <- compute_weights_design(x_norm, p)
  ## Design matrix for the global polynomial
  X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
  B <- rbind(rep(0, ncol(B)), B)
  tmbdat <- list(
    # Design matrix
    X = X,
    B = B,
    P = as(Q,'dgTMatrix'),
    logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
    # Response
    y = y,
    sigmaE = sigmaE,
    theta1 = theta
  )
  tmbparams <- list(
    W = c(rep(0, (ncol(Q) + ncol(X))))
  )
  
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    # random = "W",
    DLL = "00_Gaussian_Smoothing_diffuse_bound",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  ## posterior mode and posterior precision matrix
  opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                control = list(eval.max = 200, iter.max = 200))
  prec_matrix <- forceSymmetric(ff$he(opt$par))
  
  if (optimizer == T){
    return(list(mean = opt$par, prec = as.matrix(prec_matrix), opt = opt))
  }
  
  ## Return posterior mean and posterior precision matrix
  return(list(mean = opt$par, prec = as.matrix(prec_matrix)))
}

### Extract samples from fitted WP with random (or diffuse) boundary:
extract_samples_WP_rand_bound <- function(result, x, nsam = 3000, p = 3, target = "function"){
  samps <- rmvnp(n = nsam, mu = as.numeric(result$mean), Omega = as.matrix(result$prec))
  samps <- cbind(matrix(0,nrow = nsam, ncol = p), samps)
  xnorm <- x - x[1]
  if(target == "function"){
    samps_func <- data.frame(x = x)
    indx <- seq(1,(p*length(x)), by = p)
    wp_part <- t(samps[,indx])
    D_proposed <- compute_weights_design(xnorm,p)
    X = rbind(global_poly(xnorm[1], p = p), D_proposed[,1:p])
    poly_part <- as.matrix(X %*% t(samps[,((p*length(x)) + 1):ncol(samps)]))
    overall_part <- poly_part + wp_part
    samps_func <- cbind(samps_func, (overall_part))
    samps_func
  }
  else if(target == "derivative"){
    samps_deriv <- data.frame(x = x)
    indx <- seq(1,(p*length(x)), by = p) + 1
    D_proposed <- compute_weights_design(xnorm,p)
    degree = 1
    X = rbind(global_poly(xnorm[1], p = p), D_proposed[,1:p])
    X <- as.matrix(X[,1:(p-degree)])
    for (i in 1:ncol(X)) {
      X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
    }
    wp_part <- t(samps[,indx])
    poly_part <- as.matrix(X %*% t(samps[,((p*length(x)) + 2):ncol(samps)]))
    overall_part <- poly_part + wp_part
    samps_deriv <- cbind(samps_deriv, (overall_part))
    samps_deriv
  }
  else{
    stop(errorCondition("No selected target identified."))
  }
}

### Implement OS with random (or diffuse) boundary:
Imple_OS_rand_bound <- function(betaprec, sigmaE, theta, x, k, y, p = 3, optimizer = F, nsamps = 10000){
  if(k == length(x)){
    knots <- x
  }
  else{
    knots <- seq(0, max(x), length.out = k)
  }
  P_proposed <- compute_weights_precision(knots)
  D_proposed <- compute_weights_design(x,p)
  ## Design matrix for the global polynomial
  X = rbind(global_poly(x[1], p), D_proposed[,1:p])
  ## Design matrix for the spline basis weights
  B = as(local_poly(x = knots, refined_x = x, p = p),"dgTMatrix")
  tmbdat <- list(
    # Design matrix
    X = X,
    B = B,
    P = as(P_proposed,'dgTMatrix'),
    logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
    # Response
    y = y,
    betaprec = betaprec,
    sigmaE = sigmaE,
    theta1 = theta
  )
  tmbparams <- list(
    W = c(rep(0, (ncol(P_proposed) + ncol(X))))
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    # random = "W",
    DLL = "00_Gaussian_Smoothing_unknown_bound",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  
  ## posterior mode and posterior precision matrix
  opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                control = list(eval.max = 200, iter.max = 200))
  mean = opt$par
  prec_matrix <- forceSymmetric(ff$he(opt$par))
  sample_weights <- rmvnp(n = nsamps, mu = mean, Omega = as.matrix(prec_matrix))
  samples_fun <- data.frame(x = x)
  wp_part <- B %*% t(sample_weights[,1:ncol(P_proposed)])
  poly_part <- X %*% t(sample_weights[,(ncol(P_proposed) + 1):(ncol(sample_weights))])
  overall_part <- wp_part + poly_part
  samples_fun <- cbind(samples_fun, as.matrix(overall_part))
  
  
  Blower <- as(local_poly(x = knots, refined_x = x, p = (p-1)),"dgTMatrix")
  sample_deriv <- data.frame(x = x)
  wp_part <- as.matrix(Blower %*% t(sample_weights[,1:ncol(P_proposed)]))
  degree = 1
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  poly_part <- X %*% t(sample_weights[,(ncol(P_proposed) + 2):(ncol(sample_weights))])
  overall_part <- wp_part + poly_part
  sample_deriv <- cbind(sample_deriv, as.matrix(overall_part))
  if (optimizer == T){
    return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, opt = opt))
  }
  
  ## Return posterior mean and posterior precision matrix
  return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, mean = mean, prec = as.matrix(prec_matrix)))
}
Imple_OS_diffuse_bound <- function(sigmaE, theta, x, k, y, p = 3, optimizer = F, nsamps = 10000){
  if(k == length(x)){
    knots <- x
  }
  else{
    knots <- seq(0, max(x), length.out = k)
  }
  P_proposed <- compute_weights_precision(knots)
  D_proposed <- compute_weights_design(x,p)
  ## Design matrix for the global polynomial
  X = rbind(global_poly(x[1], p), D_proposed[,1:p])
  ## Design matrix for the spline basis weights
  B = as(local_poly(x = knots, refined_x = x, p = p),"dgTMatrix")
  tmbdat <- list(
    # Design matrix
    X = X,
    B = B,
    P = as(P_proposed,'dgTMatrix'),
    logPdet = as.numeric(determinant(P_proposed,logarithm = T)$modulus),
    # Response
    y = y,
    sigmaE = sigmaE,
    theta1 = theta
  )
  tmbparams <- list(
    W = c(rep(0, (ncol(P_proposed) + ncol(X))))
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    # random = "W",
    DLL = "00_Gaussian_Smoothing_diffuse_bound",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  
  ## posterior mode and posterior precision matrix
  opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                control = list(eval.max = 200, iter.max = 200))
  mean = opt$par
  prec_matrix <- forceSymmetric(ff$he(opt$par))
  sample_weights <- rmvnp(n = nsamps, mu = mean, Omega = as.matrix(prec_matrix))
  samples_fun <- data.frame(x = x)
  wp_part <- B %*% t(sample_weights[,1:ncol(P_proposed)])
  poly_part <- X %*% t(sample_weights[,(ncol(P_proposed) + 1):(ncol(sample_weights))])
  overall_part <- wp_part + poly_part
  samples_fun <- cbind(samples_fun, as.matrix(overall_part))
  
  
  Blower <- as(local_poly(x = knots, refined_x = x, p = (p-1)),"dgTMatrix")
  sample_deriv <- data.frame(x = x)
  wp_part <- as.matrix(Blower %*% t(sample_weights[,1:ncol(P_proposed)]))
  degree = 1
  X <- as.matrix(X[,1:(p-degree)])
  for (i in 1:ncol(X)) {
    X[,i] <- (factorial(i + degree - 1)/factorial(i-1))* X[,i]
  }
  poly_part <- X %*% t(sample_weights[,(ncol(P_proposed) + 2):(ncol(sample_weights))])
  overall_part <- wp_part + poly_part
  sample_deriv <- cbind(sample_deriv, as.matrix(overall_part))
  if (optimizer == T){
    return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, opt = opt))
  }
  
  ## Return posterior mean and posterior precision matrix
  return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, mean = mean, prec = as.matrix(prec_matrix)))
}

### Implement RW2 with random (or diffuse) boundary:
Imple_RW_rand_bound <- function(betaprec, sigmaE, theta, x, k, y, p = 3, optimizer = F, nsamps = 10000, diagonal = 0){
  if(k < length(x)){
    knots <- seq(0, max(x), length.out = k)
  }
  else {
    knots <- x
  }
  designB_rw2_reduced <- knots_RW2(x = knots, refined_x = x)
  d <- diff(knots)
  H <- compute_H_rue(d,n = length(knots))
  A <- compute_A(d,n = length(knots))
  Q_rw2_reduced <- as(t(H) %*% solve(A) %*% H, 'dgTMatrix')
  Q_rw2_reduced <- as(Q_rw2_reduced + Diagonal(k, x = diagonal), "dgTMatrix")
  D_proposed <- compute_weights_design(x,p)
  X = rbind(global_poly(x[1], p), D_proposed[,1:p])
  
  ### constraint for the two improper directions
  eigenQ <- eigen(Q_rw2_reduced, symmetric = T)
  log_det_P <-  sum(log(as.numeric(eigenQ$values[1:(ncol(Q_rw2_reduced)-2)])))
  # correction_matrix <- cbind(eigenQ$vectors[,c((ncol(Q_rw2_reduced)-1),ncol(Q_rw2_reduced))])
  correction_matrix <- cbind(c(1,rep(0,(k-1))),knots)
  correction_matrix <- rbind(correction_matrix, matrix(0, nrow = 2, ncol = p))
  
  tmbdat <- list(
    X = as(X,"dgTMatrix"),
    B = as(designB_rw2_reduced,"dgTMatrix"),
    # Penalty(Precision) matrix
    P = as(as.matrix(Q_rw2_reduced),"dgTMatrix"),
    # Log determinant of penalty matrix (without the sigma part)
    # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
    logPdet = log_det_P,
    # Response
    y = y,
    betaprec = betaprec,
    sigmaE = sigmaE,
    theta1 = theta
  )
  tmbparams <- list(
    W = c(rep(0, k + p))
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    # random = "W",
    DLL = "00_Gaussian_Smoothing_unknown_bound",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  ## posterior mode and posterior precision matrix
  opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                control = list(eval.max = 200, iter.max = 200))
  mean = opt$par
  prec_matrix <- forceSymmetric(ff$he(opt$par))
  
  sample_weights_uncorrected <- rmvnp(n = nsamps, mu = mean, Omega = as.matrix(prec_matrix))
  cov_mat <- solve(prec_matrix)
  sample_weights_uncorrected <- t(sample_weights_uncorrected)
  A <- t(correction_matrix)
  sample_weights_corrected <- sample_weights_uncorrected - cov_mat %*% t(A) %*% solve(A %*% cov_mat %*% t(A)) %*% (A %*% sample_weights_uncorrected)
  
  samples_fun <- data.frame(x = x)
  rw_part <- as.matrix(designB_rw2_reduced %*% sample_weights_corrected[1:k,])
  poly_part <- as.matrix(X %*% sample_weights_corrected[(k+1):(k+p),])
  samples_fun <- cbind(samples_fun, (poly_part + rw_part))
  
  sample_deriv_raw <- samples_fun[,-1] %>% apply(MARGIN = 2, compute_numeric_deriv, h = mean(diff(x)), degree = 1)
  sample_deriv_raw <- rbind(sample_weights_corrected[(k+p),], sample_deriv_raw)
  sample_deriv <- data.frame(x = x)
  sample_deriv <- cbind(sample_deriv, sample_deriv_raw)
  if (optimizer == T){
    return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, opt = opt))
  }
  
  ## Return posterior mean and posterior precision matrix
  return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, mean = mean, prec = as.matrix(prec_matrix)))
  
}
Imple_RW_diff_bound <- function(sigmaE, theta, x, k, y, p = 3, optimizer = F, nsamps = 10000, diagonal = 0){
  if(k < length(x)){
    knots <- seq(0, max(x), length.out = k)
  }
  else {
    knots <- x
  }
  designB_rw2_reduced <- knots_RW2(x = knots, refined_x = x)
  d <- diff(knots)
  H <- compute_H_rue(d,n = length(knots))
  A <- compute_A(d,n = length(knots))
  Q_rw2_reduced <- as(t(H) %*% solve(A) %*% H, 'dgTMatrix')
  Q_rw2_reduced <- as(Q_rw2_reduced + Diagonal(k, x = diagonal), "dgTMatrix")
  D_proposed <- compute_weights_design(x,p)

  ### constraint for the two improper directions
  eigenQ <- eigen(Q_rw2_reduced, symmetric = T)
  log_det_P <-  sum(log(as.numeric(eigenQ$values[1:(ncol(Q_rw2_reduced)-2)])))

  tmbdat <- list(
    B = as(designB_rw2_reduced,"dgTMatrix"),
    # Penalty(Precision) matrix
    P = as(as.matrix(Q_rw2_reduced),"dgTMatrix"),
    # Log determinant of penalty matrix (without the sigma part)
    # logPdet = as.numeric(sum(log(sort(all_values)[-(1:2)]))),
    logPdet = log_det_P,
    # Response
    y = y,
    sigmaE = sigmaE,
    theta1 = theta
  )
  tmbparams <- list(
    W = c(rep(0, k))
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    # random = "W",
    DLL = "00_RW2_Smoothing_diff_bound",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  ## posterior mode and posterior precision matrix
  opt <- nlminb(start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he, 
                control = list(eval.max = 200, iter.max = 200))
  mean = opt$par
  prec_matrix <- forceSymmetric(ff$he(opt$par))
  
  sample_weights <- t(rmvnp(n = nsamps, mu = mean, Omega = as.matrix(prec_matrix)))

  samples_fun <- data.frame(x = x)
  rw_part <- as.matrix(designB_rw2_reduced %*% sample_weights)
  samples_fun <- cbind(samples_fun, rw_part)
  
  sample_deriv_raw <- samples_fun[,-1] %>% apply(MARGIN = 2, compute_numeric_deriv, h = mean(diff(x)), degree = 1)
  sample_deriv_raw <- rbind(sample_deriv_raw[1,], sample_deriv_raw)
  sample_deriv <- data.frame(x = x)
  sample_deriv <- cbind(sample_deriv, sample_deriv_raw)
  if (optimizer == T){
    return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, opt = opt))
  }
  
  ## Return posterior mean and posterior precision matrix
  return(list(samples_fun = samples_fun, sample_deriv = sample_deriv, mean = mean, prec = as.matrix(prec_matrix)))
  
}

### Implement MGCV: B spline with diffuse boundary:
Imple_MGCV_diffuse_bound <- function(sigmaE, theta, x, k, y, p = 3, optimizer = F){
  data <- data.frame(x = x, y = y)
  sigmaS <- exp(-0.5*theta)
  ##### Fit P-Spline:
  ### cubic basis and second order penalty
  lambda <- (sigmaE^2)/((sigmaS^2))
  control_list <- list(scalePenalty = FALSE)
  BSmod <- gam(formula = y~ s(x,bs="bs",m=c(3,p), k = k, sp = lambda), data = data, family = gaussian(), scale = (sigmaE^2), control = control_list)
  BSmod
}

### Extract result (samples, intervals) from fitted MGCV:
extract_mgcv_samps <- function(BSmod, x, nsamps = 10000){
  newdata <- data.frame(x = x)
  ### Samples from the fitted GAM:
  gam_coef_samples <- t(rmvn(nsamps,coef(BSmod),vcov(BSmod)))
  gam_design <- predict.gam(BSmod, newdata = newdata, type = 'lpmatrix')
  gam_samples <- gam_design %*% gam_coef_samples
  gam_samples <- cbind(x, gam_samples)
  gam_samples <- as.data.frame(gam_samples)
  gam_samples
}
extract_mgcv_intervals <- function(BSmod, x, nsamps = 10000, control = list(method = "GET", GET = "rank", level = 0.8)){
  method = control$method
  GET = control$GET
  level = control$level
  newdata <- data.frame(x = x)
  ### Samples from the fitted GAM:
  gam_coef_samples <- t(rmvn(nsamps,coef(BSmod),vcov(BSmod)))
  gam_design <- predict.gam(BSmod, newdata = newdata, type = 'lpmatrix')
  gam_samples <- gam_design %*% gam_coef_samples
  gam_samples <- cbind(x, gam_samples)
  gam_samples <- as.data.frame(gam_samples)
  BS_GCR_g <- extract_GCR(gam_samples, level = level, pointwise = T, method = method, GET = GET)
  gam_conf <- confint(BSmod, parm = "s(x)", type = "confidence", shift = T, level = level, newdata = newdata)
  BS_GCR_g$mean <- gam_conf$est
  gam_deriv <- derivatives(
    BSmod,
    term = "s(x)",
    newdata = newdata,
    order = 1L,
    type = c("backward"),
    eps = 1e-10,
    interval = c("simultaneous"),
    n_sim = nsamps,
    level = level,
    unconditional = F,
    ncores = 1
  )
  gam_deriv_pt <- derivatives(
    BSmod,
    term = "s(x)",
    newdata = newdata,
    order = 1L,
    type = c("backward"),
    eps = 1e-10,
    interval = c("confidence"),
    n_sim = nsamps,
    level = level,
    unconditional = F,
    ncores = 1
  )
  
  BS_GCR_g1st <- data.frame(x = x, upper = gam_deriv$upper, lower = gam_deriv$lower, 
                            pupper = gam_deriv_pt$upper, plower = gam_deriv_pt$lower,
                            mean = gam_deriv_pt$derivative)
  list(BS_GCR_g = BS_GCR_g, BS_GCR_g1st = BS_GCR_g1st)
}



#### Summarize the result from the simulation


#### In this version, implement the new idea of coverage of posterior sample
#### pathes from the exact method (with diffuse boundary):

### Input:
# x: The location vector on which the data will be simulated
# x_compare: The location vector on which to evaluate the comparison,
#            should be choose to be far away from 0 to avoid boundary effect.
#            This have to be a subset of x vector.
# p: What is the smoothing degree to compare
# kvec: How many knots will be used in splines-based approximation (can be a vector)
# sigmaE: The residual SD in the simulation
# theta: The log precision of IWP to consider
# B: The number of replications
# nsamps: The number of samples to be drawn from each method

### Output: A list of two elements: overall_result_g and overall_result_g1st
### each element is a dataframe, with variable: Method, k, COV
compare_all_three_methods_B_times_sample_cov_diffuse <- function(x, x_compare = NULL, p, kvec, sigmaE, theta, seed = 12345,
                                                                 control = list(method = "GET", GET = "rank", level = 0.8), 
                                                                 nsamps = 8000, parallel = FALSE){
  if(is.null(x_compare)){
    x_compare <- x
  }
  method = control$method
  if(is.null(kvec)){
    kvec <- c(length(x))
  }
  GET = control$GET
  level = control$level
  indx <- which(x %in% x_compare)
  sigmaS <- sqrt(1/exp(theta))
  g_result <- data.frame()
  g1st_result <- data.frame()
  set.seed(seed)
  if(parallel == T){
    n <- length(x)
    PD_condition <- FALSE
    while(!PD_condition){
      RandomFunc <- sim_Wp_Var(mesh_size = 0.001, max_t = 10, p = p)
      RandomFunc[,-1] <- sigmaS * RandomFunc[,-1]
      
      # The function
      RandomFunc_IWPp <- RandomFunc[,1:2]
      # The first derivative, so q = p - 1
      RandomFunc_IWPq <- RandomFunc[,c(1,3)]
      
      # Plus the polynomial parts:
      beta <- rnorm(p, mean = 0, sd = 1/(factorial((1:p)-1)))
      poly_parts <- matrix(x, nrow = length(x), ncol = p)
      for (i in 1:p) {
        poly_parts[,i] <- poly_parts[,i]^(i-1)
      }
      poly_parts_gx <- poly_parts %*% beta
      poly_parts_1st <- matrix(poly_parts[,-1],ncol = (p-1))
      for (i in 2:p) {
        poly_parts_1st[,(i-1)] <- (i-1)*poly_parts[,(i-1)]
      }
      poly_parts_gx1st <- poly_parts_1st %*% beta[-1]
      gx <- evaluate_GP_Var(x, RandomFunc)[,1:2]
      names(gx) <- c("t", "IWp")
      gx$IWp <- gx$IWp + poly_parts_gx
      gx1st <- evaluate_GP_Var(x, RandomFunc)[,c(1,3)]
      names(gx1st) <- c("t", "IWq")
      gx1st$IWq <- gx1st$IWq + poly_parts_gx1st
      
      y <- gx$IWp + rnorm(n = n, mean = 0, sd = sigmaE)
      ## Fit the exact method
      exact_result <- Imple_WP_diff_bound(sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
      PD_condition <- is.positive.definite(exact_result$prec)
    }
    
    ### Convert to samples of functions
    samps_func <- extract_samples_WP_rand_bound(result = exact_result, x = x, target = "function", p = p, nsam = nsamps)
    samps_deriv <- extract_samples_WP_rand_bound(result = exact_result, x = x, target = "derivative", p = p, nsam = nsamps)
    exact_g_mean <- samps_func[,c(-1)] %>% apply(1, mean)
    exact_g1st_mean <- samps_deriv[,c(-1)] %>% apply(1, mean)
    
    
    do_once <- function(k){
      OSmod <- Imple_OS_diffuse_bound(sigmaE = sigmaE, theta = theta , x = x, k = (k+1), y = y, p = p, optimizer = T, nsamps = nsamps)
      OS_GCR_g <- extract_GCR(OSmod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      OS_GCR_g1st <- extract_GCR(OSmod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      diagonal <- 0
      RW2mod <- Imple_RW_diff_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = 2, optimizer = T, diagonal = diagonal, nsamps = nsamps)
      names(RW2mod$samples_fun)[1] <- "x"
      RW2_GCR_g <- extract_GCR(RW2mod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      RW2_GCR_g1st <- extract_GCR(RW2mod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      ##### Fit P-Spline:
      ### cubic basis and second order penalty
      BSmod <- Imple_MGCV_diffuse_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = p)
      
      ### Samples from the fitted GAM:
      BS_intervals <- extract_mgcv_intervals(BSmod = BSmod, x = x, nsamps = nsamps, control = control)
      BS_GCR_g <- BS_intervals$BS_GCR_g
      BS_GCR_g1st <- BS_intervals$BS_GCR_g1st
      
      
      OS_g_mean <- apply(OSmod$samples_fun[,-1], 1, mean)
      RW2_g_mean <- apply(RW2mod$samples_fun[,-1], 1, mean)
      OS_g1st_mean <- apply(OSmod$sample_deriv[,-1], 1, mean)
      RW2_g1st_mean <- apply(RW2mod$sample_deriv[,-1], 1, mean)
      
      g_RMSE_OS_compared <- sqrt(mean(((OS_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_OS_compared <- sqrt(mean(((OS_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_RW_compared <- sqrt(mean(((RW2_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_RW_compared <- sqrt(mean(((RW2_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_BS_compared <- sqrt(mean(((BS_GCR_g$mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_BS_compared <- sqrt(mean(((BS_GCR_g1st$mean-exact_g1st_mean)^2)[indx])) # by definition
      
      
      g_OS_cov <- compute_COV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_cov <- compute_COV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      g_OS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      
      g_result_OS <- data.frame(Method = "OS", k = k, COV = g_OS_cov, PCOV = g_OS_pcov, Type = "g", rMSE = g_RMSE_OS_compared)
      g1st_result_OS <- data.frame(Method = "OS", k = k, COV = g1st_OS_cov, PCOV = g1st_OS_pcov, Type = "g1st", rMSE = g1st_RMSE_OS_compared)
      result_OS <- rbind(g_result_OS, g1st_result_OS)
      
      g_RW2_cov <- compute_COV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_cov <- compute_COV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      g_RW2_pcov <- compute_pCOV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      
      g_result_RW <- data.frame(Method = "RW", k = k, COV = g_RW2_cov, PCOV = g_RW2_pcov, Type = "g", rMSE = g_RMSE_RW_compared)
      g1st_result_RW <- data.frame(Method = "RW", k = k, COV = g1st_RW2_cov, PCOV = g1st_RW2_pcov, Type = "g1st", rMSE = g1st_RMSE_RW_compared)
      result_RW <- rbind(g_result_RW, g1st_result_RW)
      
      g_BS_cov <- compute_COV_samps(samps = samps_func, intervals = BS_GCR_g, indx = indx)
      g1st_BS_cov <- compute_COV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      g_BS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = BS_GCR_g, indx = indx)
      g1st_BS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      
      g_result_BS <- data.frame(Method = "BS", k = k, COV = g_BS_cov, PCOV = g_BS_pcov, Type = "g", rMSE = g_RMSE_BS_compared)
      g1st_result_BS <- data.frame(Method = "BS", k = k, COV = g1st_BS_cov, PCOV = g1st_BS_pcov, Type = "g1st", rMSE = g1st_RMSE_BS_compared)
      result_BS <- rbind(g_result_BS, g1st_result_BS)
      
      overall_result <- rbind(result_OS, result_RW, result_BS)
      return(overall_result)
    }
    result_all <- foreach(i = 1:length(kvec), .combine = "rbind", .packages = c("mgcv", "gratia", "TMB", "aghq", "Matrix", "tidyverse", "GET", "tsDyn", "LaplacesDemon")) %dopar% {
      source("14_functions_defined.R") ## Only on Windows
      dyn.load(dynlib("00_Gaussian_Smoothing_known")) ## Only on Windows
      do_once(k = kvec[i])
    }
    g_result <- result_all %>% filter(Type == "g") %>% select(-Type)
    g1st_result <- result_all %>% filter(Type == "g1st") %>% select(-Type)
  }
  else{
    pb <- progress_bar$new(total = length(kvec))
    overall_result <- data.frame()
    n <- length(x)
    PD_condition <- FALSE
    while(!PD_condition){
      RandomFunc <- sim_Wp_Var(mesh_size = 0.001, max_t = 10, p = p)
      RandomFunc[,-1] <- sigmaS * RandomFunc[,-1]
      
      # The function
      RandomFunc_IWPp <- RandomFunc[,1:2]
      RandomFunc_IWPp$IWp <- RandomFunc_IWPp$IWp
      # The first derivative, so q = p - 1
      RandomFunc_IWPq <- RandomFunc[,c(1,3)]
      RandomFunc_IWPq$IWq <- RandomFunc_IWPq$IWq
      
      gx <- evaluate_GP_Var(x, RandomFunc)[,1:2]
      names(gx) <- c("t", "IWp")
      gx1st <- evaluate_GP_Var(x, RandomFunc)[,c(1,3)]
      names(gx1st) <- c("t", "IWq")
      y <- gx$IWp + rnorm(n = n, mean = 0, sd = sigmaE)
      ## Fit the exact method
      exact_result <- Imple_WP_diff_bound(sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
      PD_condition <- is.positive.definite(exact_result$prec)
    }
    
    ### Convert to samples of functions
    samps_func <- extract_samples_WP_rand_bound(result = exact_result, x = x, target = "function", p = p, nsam = nsamps)
    samps_deriv <- extract_samples_WP_rand_bound(result = exact_result, x = x, target = "derivative", p = p, nsam = nsamps)
    exact_g_mean <- samps_func[,c(-1)] %>% apply(1, mean)
    exact_g1st_mean <- samps_deriv[,c(-1)] %>% apply(1, mean)
    
    overall_result <- data.frame()
    for (i in 1:length(kvec)) {
      k <- kvec[i]
      pb$tick()
      OSmod <- Imple_OS_diffuse_bound(sigmaE = sigmaE, theta = theta , x = x, k = (k+1), y = y, p = p, optimizer = T)
      OS_GCR_g <- extract_GCR(OSmod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      OS_GCR_g1st <- extract_GCR(OSmod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      diagonal <- 0
      RW2mod <- Imple_RW_diff_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = 2, optimizer = T, diagonal = diagonal)
      names(RW2mod$samples_fun)[1] <- "x"
      RW2_GCR_g <- extract_GCR(RW2mod$samples_fun, level = level, pointwise = T, method = method, GET = GET)
      RW2_GCR_g1st <- extract_GCR(RW2mod$sample_deriv, level = level, pointwise = T, method = method, GET = GET)
      
      
      
      ##### Fit P-Spline:
      ### cubic basis and second order penalty
      
      BSmod <- Imple_MGCV_diffuse_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = p)
      
      ### Samples from the fitted GAM:
      BS_intervals <- extract_mgcv_intervals(BSmod = BSmod, x = x, nsamps = nsamps, control = control)
      BS_GCR_g <- BS_intervals$BS_GCR_g
      BS_GCR_g1st <- BS_intervals$BS_GCR_g1st
      
      OS_g_mean <- apply(OSmod$samples_fun[,-1], 1, mean)
      RW2_g_mean <- apply(RW2mod$samples_fun[,-1], 1, mean)
      OS_g1st_mean <- apply(OSmod$sample_deriv[,-1], 1, mean)
      RW2_g1st_mean <- apply(RW2mod$sample_deriv[,-1], 1, mean)
      
      g_RMSE_OS_compared <- sqrt(mean(((OS_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_OS_compared <- sqrt(mean(((OS_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_RW_compared <- sqrt(mean(((RW2_g_mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_RW_compared <- sqrt(mean(((RW2_g1st_mean-exact_g1st_mean)^2)[indx])) # by definition
      
      g_RMSE_BS_compared <- sqrt(mean(((BS_GCR_g$mean-exact_g_mean)^2)[indx])) # by definition
      g1st_RMSE_BS_compared <- sqrt(mean(((BS_GCR_g1st$mean-exact_g1st_mean)^2)[indx])) # by definition
      
      
      g_OS_cov <- compute_COV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_cov <- compute_COV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      g_OS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = OS_GCR_g, indx = indx)
      g1st_OS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = OS_GCR_g1st, indx = indx)
      
      g_result_OS <- data.frame(Method = "OS", k = k, COV = g_OS_cov, PCOV = g_OS_pcov, Type = "g", rMSE = g_RMSE_OS_compared)
      g1st_result_OS <- data.frame(Method = "OS", k = k, COV = g1st_OS_cov, PCOV = g1st_OS_pcov, Type = "g1st", rMSE = g1st_RMSE_OS_compared)
      result_OS <- rbind(g_result_OS, g1st_result_OS)
      
      g_RW2_cov <- compute_COV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_cov <- compute_COV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      g_RW2_pcov <- compute_pCOV_samps(samps = samps_func, intervals = RW2_GCR_g, indx = indx)
      g1st_RW2_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = RW2_GCR_g1st, indx = indx)
      
      g_result_RW <- data.frame(Method = "RW", k = k, COV = g_RW2_cov, PCOV = g_RW2_pcov, Type = "g", rMSE = g_RMSE_RW_compared)
      g1st_result_RW <- data.frame(Method = "RW", k = k, COV = g1st_RW2_cov, PCOV = g1st_RW2_pcov, Type = "g1st", rMSE = g1st_RMSE_RW_compared)
      result_RW <- rbind(g_result_RW, g1st_result_RW)
      
      g_BS_cov <- compute_COV_samps(samps = samps_func, intervals = BS_GCR_g, indx = indx)
      g1st_BS_cov <- compute_COV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      g_BS_pcov <- compute_pCOV_samps(samps = samps_func, intervals = BS_GCR_g, indx = indx)
      g1st_BS_pcov <- compute_pCOV_samps(samps = samps_deriv, intervals = BS_GCR_g1st, indx = indx)
      
      g_result_BS <- data.frame(Method = "BS", k = k, COV = g_BS_cov, PCOV = g_BS_pcov, Type = "g", rMSE = g_RMSE_BS_compared)
      g1st_result_BS <- data.frame(Method = "BS", k = k, COV = g1st_BS_cov, PCOV = g1st_BS_pcov, Type = "g1st", rMSE = g1st_RMSE_BS_compared)
      result_BS <- rbind(g_result_BS, g1st_result_BS)
      
      overall_result <- rbind(overall_result, result_OS, result_RW, result_BS)
      
    }
    g_result <- overall_result %>% filter(Type == "g") %>% select(-Type)
    g1st_result <- overall_result %>% filter(Type == "g1st") %>% select(-Type)
  }
  return(list(g = g_result, g1st = g1st_result))
}
