WORK_PATH <- getwd()
figure_path <- paste0(WORK_PATH, "/figure/")
result_path <- paste0(WORK_PATH, "/result/")
cpp_path <- paste0(WORK_PATH, "/cpp/")
source_path <- paste0(WORK_PATH, "/source/")
source(paste0(source_path,"00_load.R"))
source(paste0(source_path,"01_Ospline.R"))
source(paste0(source_path,"02_prior.R"))

compile(paste0(cpp_path,"00_Gaussian_Smoothing_PC.cpp"))
compile(paste0(cpp_path,"02_RW2Comparison.cpp"))

dyn.load(dynlib(paste0(cpp_path,"00_Gaussian_Smoothing_PC")))
dyn.load(dynlib(paste0(cpp_path,"02_RW2Comparison")))

GMM_simulation_agg <- function(B = 300, SNR, x, k, parallel = T, SYSTEM = "MAC"){
  do_once <- function(seed, SNR, x, k){
    set.seed(seed)
    pvec <- c(0.6,0.3,0.1)
    muvec <- rnorm(n = length(pvec), mean = 5, sd = 2)
    SDvec <- rep(1, length(pvec))
    gx <- dnormm(x = x, p = pvec, mu = muvec, sigma = SDvec)
    gx <- gx/sd(gx)
    gx1st <- diff(gx)/diff(x)
    gx1st <- c(gx1st[1], gx1st)
    gx2nd <- diff(gx1st)/diff(x)
    gx2nd <- c(gx2nd[1], gx2nd)
    y <- gx + rnorm(n = length(x), sd = sqrt(var(gx)/SNR))
    predictive_prior <- list(u = 1, alpha = 0.5)
    prior_set_1 <- prior_order_conversion_predictive(d = 1, predictive_prior, p = 3)
    prior_set_2 <- prior_order_conversion_predictive(d = 1, predictive_prior, p = 2)
    prior_set <- list(
      u1 = prior_set_1$u,
      alpha1 = prior_set_1$alpha,
      u2 = sqrt(var(gx)/SNR),
      alpha2 = 0.5,
      betaprec = 10^(-6)
    )
    prior_set_2 <- list(
      u1 = prior_set_2$u,
      alpha1 = prior_set_2$alpha,
      u2 = sqrt(var(gx)/SNR),
      alpha2 = 0.5,
      betaprec = 10^(-6)
    )
    if(is.null(k)){
      k <- length(x)
      OS3_fit <- Imple_BayesRegression(y = y, x = x, p = 3, prior.type = "PC", prior = prior_set)
      RW2_fit <- Imple_RW2(y = y, x = x, prior = prior_set_2)
      OS_samps <- extract_samples_Ospline_refined(OS3_fit, x = x, refined_x = x, p = 3)
      RW_samps <- extract_RW_samples(RW2_fit, x = x, refined_x = x)
      OS_deriv_samps <- extract_deriv_samples_OSpline(OS3_fit, p = 3, x = x, refined_x = x, degree = 1)
      RW_deriv_samps <- RW_samps[,-1] %>% apply(MARGIN = 2, FUN = compute_numeric_deriv, h = diff(x)) 
      RW_deriv_samps <- rbind(RW_deriv_samps[1,],RW_deriv_samps)
      RW_deriv_samps <- cbind(x,RW_deriv_samps)
      
      OS_2nd_deriv_samps <- extract_deriv_samples_OSpline(OS3_fit, p = 3, x = x, refined_x = x, degree = 2)
      RW_2nd_deriv_samps <- RW_deriv_samps[,-1] %>% apply(MARGIN = 2, FUN = compute_numeric_deriv, h = diff(x), degree = 1) 
      RW_2nd_deriv_samps <- rbind(RW_2nd_deriv_samps[1,],RW_2nd_deriv_samps)
      RW_2nd_deriv_samps <- cbind(x,RW_2nd_deriv_samps)
      
    }
    else{
      knots <- seq(from = 0, to = max(x), length.out = k)
      OS3_fit <- Imple_BayesRegression(y = y, x = x, p = 3, prior.type = "PC", prior = prior_set, knots = knots)
      RW2_fit <- Imple_RW2(y = y, x = x, prior = prior_set_2, knots = knots)
      data <- data.frame(x = x, y = y)
      OS_samps <- extract_samples_Ospline_refined(OS3_fit, x = knots, refined_x = x, p = 3)
      RW_samps <- extract_RW_samples(RW2_fit, x = knots, refined_x = x)
      
      OS_deriv_samps <- extract_deriv_samples_OSpline(OS3_fit, p = 3, x = knots, refined_x = x, degree = 1)
      RW_deriv_samps <- RW_samps[,-1] %>% apply(MARGIN = 2, FUN = compute_numeric_deriv, h = diff(x)) 
      RW_deriv_samps <- rbind(RW_deriv_samps[1,],RW_deriv_samps)
      RW_deriv_samps <- cbind(x,RW_deriv_samps)
      
      OS_2nd_deriv_samps <- extract_deriv_samples_OSpline(OS3_fit, p = 3, x = knots, refined_x = x, degree = 2)
      OS_2nd_deriv_samps <- extract_deriv_samples_OSpline(OS3_fit, p = 3, x = x, refined_x = x, degree = 2)
      RW_2nd_deriv_samps <- RW_deriv_samps[,-1] %>% apply(MARGIN = 2, FUN = compute_numeric_deriv, h = diff(x), degree = 1) 
      RW_2nd_deriv_samps <- rbind(RW_2nd_deriv_samps[1,],RW_2nd_deriv_samps)
      RW_2nd_deriv_samps <- cbind(x,RW_2nd_deriv_samps)
      
      
    }
    OS_g_mean <- OS_samps[,-1] %>% apply(MARGIN = 1, mean)
    RW_g_mean <- RW_samps[,-1] %>% apply(MARGIN = 1, mean)
    
    OS_g1st_mean <- OS_deriv_samps[,-1] %>% apply(MARGIN = 1, mean)
    RW_g1st_mean <- RW_deriv_samps[,-1] %>% apply(MARGIN = 1, mean)
    
    OS_g2nd_mean <- OS_2nd_deriv_samps[,-1] %>% apply(MARGIN = 1, mean)
    RW_g2nd_mean <- RW_2nd_deriv_samps[,-1] %>% apply(MARGIN = 1, mean)
    
    

    rMSE_g <- c(
    sqrt(mean((OS_g_mean - gx)^2)),
    sqrt(mean((RW_g_mean - gx)^2))
    )
    
    rMSE_g1st <- c(
    sqrt(mean((OS_g1st_mean - gx1st)^2)),
    sqrt(mean((RW_g1st_mean - gx1st)^2))
    )
    
    rMSE_g2nd <- c(
      sqrt(mean((OS_g2nd_mean - gx2nd)^2)),
      sqrt(mean((RW_g2nd_mean - gx2nd)^2))
    )
    
    
    result <- data.frame(method = c("OS", "RW"), k = k, seed = i, rMSE_g = rMSE_g,
                         rMSE_g1st = rMSE_g1st, rMSE_g2nd = rMSE_g2nd)
  }
  if(parallel == T){
    if(SYSTEM == "MAC"){
      result_all <- foreach(i = 1:B, .combine = "rbind", .packages = c("progress", "TMB", "aghq", "Matrix", "tidyverse", "GET", "tsDyn", "LaplacesDemon")) %dopar% {
        do_once(seed = i, SNR = SNR, x = x, k = k)
      } 
    }
    else{
      result_all <- foreach(i = 1:B, .combine = "rbind", .packages = c("progress", "TMB", "aghq", "Matrix", "tidyverse", "GET", "tsDyn", "LaplacesDemon")) %dopar% {
        WORK_PATH <- getwd()
        figure_path <- paste0(WORK_PATH, "/figure/")
        result_path <- paste0(WORK_PATH, "/result/")
        cpp_path <- paste0(WORK_PATH, "/cpp/")
        source_path <- paste0(WORK_PATH, "/source/")
        source(paste0(source_path,"00_load.R"))
        source(paste0(source_path,"01_Ospline.R"))
        source(paste0(source_path,"02_prior.R"))
        dyn.load(dynlib(paste0(cpp_path,"00_Gaussian_Smoothing_PC")))
        dyn.load(dynlib(paste0(cpp_path,"02_RW2Comparison")))
        do_once(seed = i, SNR = SNR, x = x, k = k)
      } 
    }

    result_all
  }
  else{
    result_all <- data.frame()
    for (i in 1:B) {
      result_all <- rbind(result_all, do_once(seed = i, SNR = SNR, x = x, k = k))
    }
    result_all
  }
}

### Parallel setup:
n.cores <- parallel::detectCores() - 2
my_cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK" # on mac/linux, use FORK; on windows, use PSOCK
)
doParallel::registerDoParallel(cl = my_cluster)

#########################################
#########################################
#########################################

## High SNR case: SNR = 100
x <- seq(from = 0, to = 10, length.out = 100)
start_time <- Sys.time()
result_high <- GMM_simulation_agg(SYSTEM = "WIN", B = 300, SNR = 100, x = x, k = 100, parallel = T)
end_time <- Sys.time()
end_time - start_time
save(file = paste0(result_path, "GMM_result_high.rda"), result_high)

## Med SNR case: SNR = 10
result_med <- GMM_simulation_agg(SYSTEM = "WIN", B = 300, SNR = 10, x = x, k = 100, parallel = T)
save(file = paste0(result_path, "GMM_result_med.rda"), result_med)

## Low SNR case: SNR = 1
result_low <- GMM_simulation_agg(SYSTEM = "WIN", B = 300, SNR = 1, x = x, k = 100, parallel = T)
save(file = paste0(result_path, "GMM_result_low.rda"), result_low)










