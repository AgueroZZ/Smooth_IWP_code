#### Fit IWP model with augmented method:
Imple_IWP_augment <- function(y, x, z, p, prior = NULL, starting_val = 0, aghq_k = 10){
  if(is.null(prior)){
    u1 = 1
    alpha1 = 0.5
    sigma2 = 0.2
    betaprec = 10^(-6)
  }
  else{
    u1 = prior$u1
    alpha1 = prior$alpha1
    sigma2 = prior$sigma2
    betaprec = prior$betaprec
  }
  svec <- sort(unique(c(x,z))) - starting_val ## bind x and z together
  if(min(svec) < 0){
    stop("Starting_val has to be smaller than all the locations.")
  }
  x_ind <- which(svec %in% x)
  if(min(svec) == 0){
    svec <- svec[-1]
    B <- Compute_design_Aug(svec = svec, p = p) 
    B <- rbind(rep(0, ncol(B)), B)
    B <- B[x_ind, ] ## only those in x vector affects likelihood
  }
  else{
    B <- Compute_design_Aug(svec = svec, p = p)[x_ind, ] ## only those in x vector affects likelihood
  }
  Q <- Compute_Aug_Wp_Prec(svec = svec,p = p)
  D_proposed <- compute_weights_design(x,p)
  ## Design matrix for the global polynomial
  X = rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
  tmbdat <- list(
    # Design matrix
    X = X,
    B = B,
    P = as(Q,'dgTMatrix'),
    logPdet = as.numeric(determinant(Q,logarithm = T)$modulus),
    # Response
    y = y,
    # PC Prior params
    u1 = u1,
    alpha1 = alpha1,
    sigma2 = sigma2,
    betaprec = betaprec
  )
  tmbparams <- list(
    W = c(rep(0, (ncol(Q) + ncol(X)))), 
    theta1 = 0
  )
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = tmbparams,
    random = "W",
    DLL = "00_Gaussian_Smoothing_PC_known",
    silent = TRUE
  )
  # Hessian not implemented for RE models
  ff$he <- function(w) numDeriv::jacobian(ff$gr,w)
  aghq::marginal_laplace_tmb(ff,aghq_k,c(0))
  
}


#### Extract samples from exact IWP:
extract_samples_IWP_exact <- function(mod, n_samps, x, p = 3){
  pos_samps_WP <- sample_marginal(mod, n_samps)
  if(min(x) == 0){
    x <- x[-1]
  }
  x_standardized <- c(0,x)
  B <- rbind(rep(0,length(x)), diag(1, nrow = length(x)))
  D_proposed <- compute_weights_design(x_standardized,p)
  ## Design matrix for the global polynomial
  designX <- rbind(c(1,rep(0,p-1)), D_proposed[,1:p])
  fitted_samps <- cbind(B,designX) %*% pos_samps_WP$samps 
  result <- cbind(x = x_standardized, data.frame(as.matrix(fitted_samps)))
  result
}








