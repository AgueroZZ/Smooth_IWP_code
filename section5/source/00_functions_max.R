identify_arg_max <- function(samps){
  indx <- apply(samps[,-1], MARGIN = 2, which.max)
  samps[,1][indx]
}

identify_max <- function(samps){
  max_var <- apply(samps[,-1], MARGIN = 2, max)
  max_var
}

compare_ks <- function(x, p, kvec, SNR, theta, seed = 12345, 
                       nsamps = 8000, g_type = "IWP"){
  set.seed(seed)
  n <- length(x)
  if(g_type == "IWP"){
    sigmaS <- sqrt(1/exp(theta))
    PD_condition <- FALSE
    while(!PD_condition){
      sigmaE <- sqrt((sigmaS^2)/SNR)
      RandomFunc <- sim_Wp_Var(mesh_size = 0.001, max_t = max(x), p = p, sigmaS = sigmaS)
      # The function
      RandomFunc_IWPp <- RandomFunc[,1:2]
      # The first derivative, so q = p - 1
      RandomFunc_IWPq <- RandomFunc[,c(1,3)]
      
      # Plus the polynomial parts:
      beta <- rnorm(p, mean = 0, sd = 1/(2*factorial((1:p)-1)))
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
  }
  else{
    sigmaE <- sqrt(1/SNR)
    PD_condition <- FALSE
    while(!PD_condition){
      g <- simulate_g_function(sigmaS = 1, x)
      y <- g(x) + rnorm(n = n, mean = 0, sd = sigmaE)
      ## Fit the exact method
      exact_result <- Imple_WP_diff_bound(sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
      PD_condition <- is.positive.definite(exact_result$prec)
    }
  }
  
  extracted_samples <- extract_samples_WP_rand_bound(exact_result, x, p = p, nsam = nsamps)
  arg_max_g_exact <- identify_arg_max(extracted_samples)
  max_g_exact <- identify_max(extracted_samples)
  extracted_samples_deriv <- extract_samples_WP_rand_bound(exact_result, x, p = p, target = "derivative", nsam = nsamps)
  arg_max_g1st_exact <- identify_arg_max(extracted_samples_deriv)
  max_g1st_exact <- identify_max(extracted_samples_deriv)
  
  
  g_result_arg <- data.frame()
  g1st_result_arg <- data.frame()
  
  g_result_max <- data.frame()
  g1st_result_max <- data.frame()
  
  pb <- progress_bar$new(total = length(kvec))
  
  for (i in 1:length(kvec)) {
    k <- kvec[i]
    pb$tick()
    OSmod <- Imple_OS_diffuse_bound(sigmaE = sigmaE, theta = theta , x = x, k = (k + 1), y = y, p = p, optimizer = T, nsamps = nsamps)
    arg_max_g_OS <- identify_arg_max(OSmod$samples_fun)
    arg_max_g1st_OS <- identify_arg_max(OSmod$sample_deriv)
    max_g_OS <- identify_max(OSmod$samples_fun)
    max_g1st_OS <- identify_max(OSmod$sample_deriv)
    
    diagonal <- 0
    RW2mod <- Imple_RW_diff_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = 2, optimizer = T, diagonal = diagonal, nsamps = nsamps)
    names(RW2mod$samples_fun)[1] <- "x"
    
    arg_max_g_RW <- identify_arg_max(RW2mod$samples_fun)
    arg_max_g1st_RW <- identify_arg_max(RW2mod$sample_deriv)
    max_g_RW <- identify_max(RW2mod$samples_fun)
    max_g1st_RW <- identify_max(RW2mod$sample_deriv)
    
    ##### Fit P-Spline:
    ### cubic basis and second order penalty
    BSmod <- Imple_MGCV_diffuse_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = p)
    fitted_mgcv_g_samps <- extract_mgcv_samps(BSmod, x = x, nsamps = nsamps)
    fitted_mgcv_g1st_samps <- fitted_mgcv_g_samps[,-1] %>% apply(MARGIN = 2, FUN = compute_numeric_deriv, h = diff(x)) 
    fitted_mgcv_g1st_samps <- rbind(fitted_mgcv_g1st_samps[1,],fitted_mgcv_g1st_samps)
    fitted_mgcv_g1st_samps <- cbind(x,fitted_mgcv_g1st_samps)
    
    
    arg_max_g_BS <- identify_arg_max(fitted_mgcv_g_samps)
    arg_max_g1st_BS <- identify_arg_max(fitted_mgcv_g1st_samps)
    max_g_BS <- identify_max(fitted_mgcv_g_samps)
    max_g1st_BS <- identify_max(fitted_mgcv_g1st_samps)
    
    ks_g_OS_arg <- suppressWarnings(ks.test(x = arg_max_g_exact, y = arg_max_g_OS))
    ks_g_RW_arg <- suppressWarnings(ks.test(x = arg_max_g_exact, y = arg_max_g_RW))
    ks_g_BS_arg <- suppressWarnings(ks.test(x = arg_max_g_exact, y = arg_max_g_BS))
    
    ks_g_OS_max <- suppressWarnings(ks.test(x = max_g_exact, y = max_g_OS))
    ks_g_RW_max <- suppressWarnings(ks.test(x = max_g_exact, y = max_g_RW))
    ks_g_BS_max <- suppressWarnings(ks.test(x = max_g_exact, y = max_g_BS))
    
    
    ks_g1st_OS_arg <- suppressWarnings(ks.test(x = arg_max_g1st_exact, y = arg_max_g1st_OS)) 
    ks_g1st_RW_arg <- suppressWarnings(ks.test(x = arg_max_g1st_exact, y = arg_max_g1st_RW))
    ks_g1st_BS_arg <- suppressWarnings(ks.test(x = arg_max_g1st_exact, y = arg_max_g1st_BS))
    
    ks_g1st_OS_max <- suppressWarnings(ks.test(x = max_g1st_exact, y = max_g1st_OS))
    ks_g1st_RW_max <- suppressWarnings(ks.test(x = max_g1st_exact, y = max_g1st_RW))
    ks_g1st_BS_max <- suppressWarnings(ks.test(x = max_g1st_exact, y = max_g1st_BS))
    
    
    g_result_arg <- rbind(g_result_arg, 
                          data.frame(method = "OS", KS = ks_g_OS_arg$statistic, P = ks_g_OS_arg$p.value, k = k))
    
    g_result_arg <- rbind(g_result_arg, 
                          data.frame(method = "RW", KS = ks_g_RW_arg$statistic, P = ks_g_RW_arg$p.value, k = k))
    
    g_result_arg <- rbind(g_result_arg, 
                          data.frame(method = "BS", KS = ks_g_BS_arg$statistic, P = ks_g_BS_arg$p.value, k = k))
    
    g1st_result_arg <- rbind(g1st_result_arg, 
                             data.frame(method = "OS", KS = ks_g1st_OS_arg$statistic, P = ks_g1st_OS_arg$p.value, k = k))
    
    g1st_result_arg <- rbind(g1st_result_arg, 
                             data.frame(method = "RW", KS = ks_g1st_RW_arg$statistic, P = ks_g1st_RW_arg$p.value, k = k))
    
    g1st_result_arg <- rbind(g1st_result_arg, 
                             data.frame(method = "BS", KS = ks_g1st_BS_arg$statistic, P = ks_g1st_BS_arg$p.value, k = k))
    
    
    g_result_max <- rbind(g_result_max, 
                          data.frame(method = "OS", KS = ks_g_OS_max$statistic, P = ks_g_OS_max$p.value, k = k))
    
    g_result_max <- rbind(g_result_max, 
                          data.frame(method = "RW", KS = ks_g_RW_max$statistic, P = ks_g_RW_max$p.value, k = k))
    
    g_result_max <- rbind(g_result_max, 
                          data.frame(method = "BS", KS = ks_g_BS_max$statistic, P = ks_g_BS_max$p.value, k = k))
    
    g1st_result_max <- rbind(g1st_result_max, 
                             data.frame(method = "OS", KS = ks_g1st_OS_max$statistic, P = ks_g1st_OS_max$p.value, k = k))
    
    g1st_result_max <- rbind(g1st_result_max, 
                             data.frame(method = "RW", KS = ks_g1st_RW_max$statistic, P = ks_g1st_RW_max$p.value, k = k))
    
    g1st_result_max <- rbind(g1st_result_max, 
                             data.frame(method = "BS", KS = ks_g1st_BS_max$statistic, P = ks_g1st_BS_max$p.value, k = k))
    
  }
  
  list(g_result_arg = g_result_arg, g1st_result_arg = g1st_result_arg, g_result_max = g_result_max, g1st_result_max = g1st_result_max)
  
}

simulate_g_function <- function(sigmaS, x){
  beta <- rnorm(n = 5, sd = 1)
  g <- function(x) {beta[1] * sin(x) + beta[2]*sin(beta[3] * x) + beta[4]*cos(beta[5] * x)}
  sdS <- sd(g(x))
  f <- function(x){sigmaS * g(x)/sdS}
  f
}


#### Aggregate the above simulation through B replications
agg_compare_diff_ks <- function(B = 300, x, p, kvec, SNR, theta, seed = 12345, 
                                nsamps = 8000, g_type = "IWP", parallel = FALSE){
  ### Aggregate the result B times:
  
  if(parallel == FALSE){
    result_cov_boundary_g_arg <- data.frame()
    result_cov_boundary_g1st_arg <- data.frame()
    result_cov_boundary_g_max <- data.frame()
    result_cov_boundary_g1st_max <- data.frame()
    
    for (i in 1:B) {
      result_cov_boundary_new <- compare_ks(x, p, kvec, SNR, theta, seed = i, 
                                            nsamps = nsamps, g_type = g_type)
      result_cov_boundary_new$g_result_arg$Rep = (i+1)
      result_cov_boundary_new$g1st_result_arg$Rep = (i+1)
      result_cov_boundary_new$g_result_max$Rep = (i+1)
      result_cov_boundary_new$g1st_result_max$Rep = (i+1)
      
      result_cov_boundary_g_arg <- rbind(result_cov_boundary_g_arg, result_cov_boundary_new$g_result_arg)
      result_cov_boundary_g1st_arg <- rbind(result_cov_boundary_g1st_arg, result_cov_boundary_new$g1st_result_arg)
      result_cov_boundary_g_max <- rbind(result_cov_boundary_g_max, result_cov_boundary_new$g_result_max)
      result_cov_boundary_g1st_max <- rbind(result_cov_boundary_g1st_max, result_cov_boundary_new$g1st_result_max)
    }
  }
  else{
    do_once <- function(seed){
      result_cov_boundary_new <- compare_ks(x, p, kvec, SNR, theta, seed = seed, 
                                            nsamps = nsamps, g_type = g_type)
      
      result_cov_boundary_new$g_result_arg$Type = "g_arg"
      result_cov_boundary_new$g1st_result_arg$Type = "g1st_arg"
      result_cov_boundary_new$g_result_max$Type = "g_max"
      result_cov_boundary_new$g1st_result_max$Type = "g1st_max"
      
      rbind(result_cov_boundary_new$g_result_arg, result_cov_boundary_new$g1st_result_arg, result_cov_boundary_new$g_result_max, result_cov_boundary_new$g1st_result_max)
    }
    result_all <- foreach(i = 1:B, .combine = "rbind", .packages = c("progress","mgcv", "gratia", "TMB", "aghq", "Matrix", "tidyverse", "GET", "tsDyn", "LaplacesDemon")) %dopar% {
      # ## Only on Windows
      # source("14_functions_defined.R")
      # source("00_functions_random_boundary.R")
      # dyn.load(dynlib("00_Gaussian_Smoothing_unknown_bound"))
      # dyn.load(dynlib("00_Gaussian_Smoothing_diffuse_bound"))
      # dyn.load(dynlib("00_RW2_Smoothing_diff_bound"))
      do_once(seed = i)
    }
    result_cov_boundary_g_arg <- result_all %>% filter(Type == "g_arg") %>% select(-Type)
    result_cov_boundary_g1st_arg <- result_all %>% filter(Type == "g1st_arg") %>% select(-Type)
    
    result_cov_boundary_g_max <- result_all %>% filter(Type == "g_max") %>% select(-Type)
    result_cov_boundary_g1st_max <- result_all %>% filter(Type == "g1st_max") %>% select(-Type)
  }
  
  return(list(g_result_arg = result_cov_boundary_g_arg, g1st_result_arg = result_cov_boundary_g1st_arg, 
              g_result_max = result_cov_boundary_g_max, g1st_result_max = result_cov_boundary_g1st_max))
  
}


### Consider the case where location vector x is random generated as well?

## call this version v2:
compare_ks_v2 <- function(n, xmax, p, kvec, SNR, theta, seed = 12345, 
                          nsamps = 8000, g_type = "IWP"){
  set.seed(seed)
  if(g_type == "IWP"){
    sigmaS <- sqrt(1/exp(theta))
    PD_condition <- FALSE
    while(!PD_condition){
      x_unique <- FALSE
      while (!x_unique) {
        x <- c(0,sort(runif((n-1), min = 0, max = xmax)))
        x <- round(x, digits = 3)
        x_unique <- (length(x) == length(unique(x)))
      }
      sigmaE <- sqrt((sigmaS^2)/SNR)
      RandomFunc <- sim_Wp_Var(mesh_size = 0.001, max_t = max(x), p = p, sigmaS = sigmaS)
      # The function
      RandomFunc_IWPp <- RandomFunc[,1:2]
      # The first derivative, so q = p - 1
      RandomFunc_IWPq <- RandomFunc[,c(1,3)]
      
      # Plus the polynomial parts:
      beta <- rnorm(p, mean = 0, sd = 1/(2*factorial((1:p)-1)))
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
  }
  else{
    sigmaE <- sqrt(1/SNR)
    PD_condition <- FALSE
    while(!PD_condition){
      x_unique <- FALSE
      while (!x_unique) {
        x <- c(0,sort(runif((n-1), min = 0, max = xmax)))
        x <- round(x, digits = 3)
        x_unique <- (length(x) == length(unique(x)))
      }
      g <- simulate_g_function(sigmaS = 1, x)
      y <- g(x) + rnorm(n = n, mean = 0, sd = sigmaE)
      ## Fit the exact method
      exact_result <- Imple_WP_diff_bound(sigmaE = sigmaE, theta = theta, p = p, x = x, y = y, optimizer = T)
      PD_condition <- is.positive.definite(exact_result$prec)
    }
  }
  
  extracted_samples <- extract_samples_WP_rand_bound(exact_result, x, p = p, nsam = nsamps)
  arg_max_g_exact <- identify_arg_max(extracted_samples)
  max_g_exact <- identify_max(extracted_samples)
  extracted_samples_deriv <- extract_samples_WP_rand_bound(exact_result, x, p = p, target = "derivative", nsam = nsamps)
  arg_max_g1st_exact <- identify_arg_max(extracted_samples_deriv)
  max_g1st_exact <- identify_max(extracted_samples_deriv)
  
  
  g_result_arg <- data.frame()
  g1st_result_arg <- data.frame()
  
  g_result_max <- data.frame()
  g1st_result_max <- data.frame()
  
  pb <- progress_bar$new(total = length(kvec))
  
  for (i in 1:length(kvec)) {
    k <- kvec[i]
    pb$tick()
    OSmod <- Imple_OS_diffuse_bound(sigmaE = sigmaE, theta = theta , x = x, k = (k + 1), y = y, p = p, optimizer = T, nsamps = nsamps)
    arg_max_g_OS <- identify_arg_max(OSmod$samples_fun)
    arg_max_g1st_OS <- identify_arg_max(OSmod$sample_deriv)
    max_g_OS <- identify_max(OSmod$samples_fun)
    max_g1st_OS <- identify_max(OSmod$sample_deriv)
    
    diagonal <- 0
    RW2mod <- Imple_RW_diff_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = 2, optimizer = T, diagonal = diagonal, nsamps = nsamps)
    names(RW2mod$samples_fun)[1] <- "x"
    
    arg_max_g_RW <- identify_arg_max(RW2mod$samples_fun)
    arg_max_g1st_RW <- identify_arg_max(RW2mod$sample_deriv)
    max_g_RW <- identify_max(RW2mod$samples_fun)
    max_g1st_RW <- identify_max(RW2mod$sample_deriv)
    
    ##### Fit P-Spline:
    ### cubic basis and second order penalty
    BSmod <- Imple_MGCV_diffuse_bound(sigmaE = sigmaE, theta = theta, x = x, k = k, y = y, p = p)
    fitted_mgcv_g_samps <- extract_mgcv_samps(BSmod, x = x, nsamps = nsamps)
    fitted_mgcv_g1st_samps <- fitted_mgcv_g_samps[,-1] %>% apply(MARGIN = 2, FUN = compute_numeric_deriv, h = diff(x)) 
    fitted_mgcv_g1st_samps <- rbind(fitted_mgcv_g1st_samps[1,],fitted_mgcv_g1st_samps)
    fitted_mgcv_g1st_samps <- cbind(x,fitted_mgcv_g1st_samps)
    
    
    arg_max_g_BS <- identify_arg_max(fitted_mgcv_g_samps)
    arg_max_g1st_BS <- identify_arg_max(fitted_mgcv_g1st_samps)
    max_g_BS <- identify_max(fitted_mgcv_g_samps)
    max_g1st_BS <- identify_max(fitted_mgcv_g1st_samps)
    
    ks_g_OS_arg <- suppressWarnings(ks.test(x = arg_max_g_exact, y = arg_max_g_OS))
    ks_g_RW_arg <- suppressWarnings(ks.test(x = arg_max_g_exact, y = arg_max_g_RW))
    ks_g_BS_arg <- suppressWarnings(ks.test(x = arg_max_g_exact, y = arg_max_g_BS))
    
    ks_g_OS_max <- suppressWarnings(ks.test(x = max_g_exact, y = max_g_OS))
    ks_g_RW_max <- suppressWarnings(ks.test(x = max_g_exact, y = max_g_RW))
    ks_g_BS_max <- suppressWarnings(ks.test(x = max_g_exact, y = max_g_BS))
    
    
    ks_g1st_OS_arg <- suppressWarnings(ks.test(x = arg_max_g1st_exact, y = arg_max_g1st_OS)) 
    ks_g1st_RW_arg <- suppressWarnings(ks.test(x = arg_max_g1st_exact, y = arg_max_g1st_RW))
    ks_g1st_BS_arg <- suppressWarnings(ks.test(x = arg_max_g1st_exact, y = arg_max_g1st_BS))
    
    ks_g1st_OS_max <- suppressWarnings(ks.test(x = max_g1st_exact, y = max_g1st_OS))
    ks_g1st_RW_max <- suppressWarnings(ks.test(x = max_g1st_exact, y = max_g1st_RW))
    ks_g1st_BS_max <- suppressWarnings(ks.test(x = max_g1st_exact, y = max_g1st_BS))
    
    
    g_result_arg <- rbind(g_result_arg, 
                          data.frame(method = "OS", KS = ks_g_OS_arg$statistic, P = ks_g_OS_arg$p.value, k = k))
    
    g_result_arg <- rbind(g_result_arg, 
                          data.frame(method = "RW", KS = ks_g_RW_arg$statistic, P = ks_g_RW_arg$p.value, k = k))
    
    g_result_arg <- rbind(g_result_arg, 
                          data.frame(method = "BS", KS = ks_g_BS_arg$statistic, P = ks_g_BS_arg$p.value, k = k))
    
    g1st_result_arg <- rbind(g1st_result_arg, 
                             data.frame(method = "OS", KS = ks_g1st_OS_arg$statistic, P = ks_g1st_OS_arg$p.value, k = k))
    
    g1st_result_arg <- rbind(g1st_result_arg, 
                             data.frame(method = "RW", KS = ks_g1st_RW_arg$statistic, P = ks_g1st_RW_arg$p.value, k = k))
    
    g1st_result_arg <- rbind(g1st_result_arg, 
                             data.frame(method = "BS", KS = ks_g1st_BS_arg$statistic, P = ks_g1st_BS_arg$p.value, k = k))
    
    
    g_result_max <- rbind(g_result_max, 
                          data.frame(method = "OS", KS = ks_g_OS_max$statistic, P = ks_g_OS_max$p.value, k = k))
    
    g_result_max <- rbind(g_result_max, 
                          data.frame(method = "RW", KS = ks_g_RW_max$statistic, P = ks_g_RW_max$p.value, k = k))
    
    g_result_max <- rbind(g_result_max, 
                          data.frame(method = "BS", KS = ks_g_BS_max$statistic, P = ks_g_BS_max$p.value, k = k))
    
    g1st_result_max <- rbind(g1st_result_max, 
                             data.frame(method = "OS", KS = ks_g1st_OS_max$statistic, P = ks_g1st_OS_max$p.value, k = k))
    
    g1st_result_max <- rbind(g1st_result_max, 
                             data.frame(method = "RW", KS = ks_g1st_RW_max$statistic, P = ks_g1st_RW_max$p.value, k = k))
    
    g1st_result_max <- rbind(g1st_result_max, 
                             data.frame(method = "BS", KS = ks_g1st_BS_max$statistic, P = ks_g1st_BS_max$p.value, k = k))
    
  }
  
  list(g_result_arg = g_result_arg, g1st_result_arg = g1st_result_arg, g_result_max = g_result_max, g1st_result_max = g1st_result_max)
  
}


#### Aggregate the above simulation through B replications
agg_compare_diff_ks_v2 <- function(B = 300, n, xmax, p, kvec, SNR, theta, seed = 12345, 
                                   nsamps = 8000, g_type = "IWP", parallel = FALSE){
  ### Aggregate the result B times:
  
  if(parallel == FALSE){
    result_cov_boundary_g_arg <- data.frame()
    result_cov_boundary_g1st_arg <- data.frame()
    result_cov_boundary_g_max <- data.frame()
    result_cov_boundary_g1st_max <- data.frame()
    
    for (i in 1:B) {
      result_cov_boundary_new <- compare_ks_v2(n, xmax, p, kvec, SNR, theta, seed = i, 
                                               nsamps = nsamps, g_type = g_type)
      result_cov_boundary_new$g_result_arg$Rep = (i+1)
      result_cov_boundary_new$g1st_result_arg$Rep = (i+1)
      result_cov_boundary_new$g_result_max$Rep = (i+1)
      result_cov_boundary_new$g1st_result_max$Rep = (i+1)
      
      result_cov_boundary_g_arg <- rbind(result_cov_boundary_g_arg, result_cov_boundary_new$g_result_arg)
      result_cov_boundary_g1st_arg <- rbind(result_cov_boundary_g1st_arg, result_cov_boundary_new$g1st_result_arg)
      result_cov_boundary_g_max <- rbind(result_cov_boundary_g_max, result_cov_boundary_new$g_result_max)
      result_cov_boundary_g1st_max <- rbind(result_cov_boundary_g1st_max, result_cov_boundary_new$g1st_result_max)
    }
  }
  else{
    do_once <- function(seed){
      result_cov_boundary_new <- compare_ks_v2(n, xmax, p, kvec, SNR, theta, seed = seed, 
                                               nsamps = nsamps, g_type = g_type)
      
      result_cov_boundary_new$g_result_arg$Type = "g_arg"
      result_cov_boundary_new$g1st_result_arg$Type = "g1st_arg"
      result_cov_boundary_new$g_result_max$Type = "g_max"
      result_cov_boundary_new$g1st_result_max$Type = "g1st_max"
      
      rbind(result_cov_boundary_new$g_result_arg, result_cov_boundary_new$g1st_result_arg, result_cov_boundary_new$g_result_max, result_cov_boundary_new$g1st_result_max)
    }
    result_all <- foreach(i = 1:B, .combine = "rbind", .packages = c("progress","mgcv", "gratia", "TMB", "aghq", "Matrix", "tidyverse", "GET", "tsDyn", "LaplacesDemon")) %dopar% {
      # ## Only on Windows
      # source("14_functions_defined.R")
      # source("00_functions_random_boundary.R")
      # dyn.load(dynlib("00_Gaussian_Smoothing_unknown_bound"))
      # dyn.load(dynlib("00_Gaussian_Smoothing_diffuse_bound"))
      # dyn.load(dynlib("00_RW2_Smoothing_diff_bound"))
      do_once(seed = i)
    }
    result_cov_boundary_g_arg <- result_all %>% filter(Type == "g_arg") %>% select(-Type)
    result_cov_boundary_g1st_arg <- result_all %>% filter(Type == "g1st_arg") %>% select(-Type)
    
    result_cov_boundary_g_max <- result_all %>% filter(Type == "g_max") %>% select(-Type)
    result_cov_boundary_g1st_max <- result_all %>% filter(Type == "g1st_max") %>% select(-Type)
  }
  
  return(list(g_result_arg = result_cov_boundary_g_arg, g1st_result_arg = result_cov_boundary_g1st_arg, 
              g_result_max = result_cov_boundary_g_max, g1st_result_max = result_cov_boundary_g1st_max))
  
}
