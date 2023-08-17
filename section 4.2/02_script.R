work_directory <- getwd()
source_location <- paste0(work_directory, "/source/")
cpp_location <- paste0(work_directory, "/cpp/")
figures_location <- paste0(work_directory, "/figures/")
data_location <- paste0(work_directory, "/data/")

source(paste0(source_location, "loads.R"))
source(paste0(source_location, "01_Ospline.R"))
source(paste0(source_location, "02_Prior.R"))
source(paste0(source_location, "03_IWP.R"))

compile(paste0(cpp_location, "00_Gaussian_Smoothing_PC_known.cpp"))
dyn.load(dynlib(paste0(cpp_location, "00_Gaussian_Smoothing_PC_known")))
Rcpp::sourceCpp(paste0(cpp_location, "OSpline_Design.cpp"), showOutput = FALSE)
aghq_k = 10
p = 3

############################################
############################################
### Time Comparison:  ######################
############################################
############################################

fit_once_time <- function(method = "IWP", y, x, z = NULL, p, prior, k = NULL, aghq_k = 10){
  locations <- sort(unique(c(x,z)))
  if(method == "IWP"){
    time_start <- Sys.time()
    ### Fit IWP3 with exact method:
    model_IWP_exact <- Imple_IWP_augment(y = y, x = x, z = NULL, prior = prior, starting_val = 0, p = p, aghq_k = aghq_k)
    ### Obtain IWP3 samples
    samps <- extract_samples_WP(model_IWP_exact, n_samp = 3000, p = p, x = locations)
    time_finish <- Sys.time()
    mylist <- model_IWP_exact$modesandhessians$H
    list(runtime = difftime(time1 = time_finish, time2 = time_start, units = "secs"), memory = pryr::object_size(model_IWP_exact),
         condition_number =  log10(unlist(lapply(mylist, kappa))))
  }
  else{
    knots <- seq(min(locations), max(locations), length.out = k)
    time_start <- Sys.time()
    model_OS3 <- Imple_BayesRegression(y = y, x = x, knots = knots, aghq_k = aghq_k, prior = prior, p = p)
    samps <- extract_samples_Ospline_refined(model_OS3, n_samp = 3000, p = p, x = knots, refined_x = locations)
    time_finish <- Sys.time()
    mylist <- model_OS3$modesandhessians$H
    list(runtime = difftime(time1 = time_finish, time2 = time_start, units = "secs"), memory = pryr::object_size(model_OS3),
         condition_number = log10(unlist(lapply(mylist, kappa))))
  }
}

### When N = 50:
N = 50
set.seed(12345)
x <- seq(0,20, length.out = N)
gx <- sqrt(3)*sin(x/2)
y <- gx + rnorm(n = length(x), sd = 1)
z <- NULL

N50_IWP3 <- c()
N50_OS3_small <- c()
N50_OS3_med <- c()
N50_OS3_large <- c()
N50_OS3_very_large <- c()

N50_IWP3_condition <- c()
N50_OS3_small_condition <- c()
N50_OS3_med_condition <- c()
N50_OS3_large_condition <- c()
N50_OS3_very_large_condition <- c()
for (i in 1:10) {
  predictive_prior <- list(u = 3, alpha = 0.01)
  prior_set_1 <- prior_order_conversion_predictive(d = 5, predictive_prior, p = p)
  prior_set <- list(
    u1 = prior_set_1$u,
    alpha1 = prior_set_1$alpha,
    sigma2 = 1,
    betaprec = 10^(-3)
  )
  IWP_fit <- fit_once_time(method = "IWP", y = y, x = x, z = z, prior = prior_set, p = p)
  OS3_small <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 10, p = p)
  OS3_med <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 30, p = p)
  OS3_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 50, p = p)
  OS3_very_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 100, p = p)
  
  N50_IWP3[i] <- IWP_fit$runtime
  N50_OS3_small[i] <- OS3_small$runtime
  N50_OS3_med[i] <- OS3_med$runtime
  N50_OS3_large[i] <- OS3_large$runtime
  N50_OS3_very_large[i] <- OS3_very_large$runtime
  
  N50_IWP3_condition[i] <- mean(IWP_fit$condition_number)
  N50_OS3_small_condition[i] <- mean(OS3_small$condition_number)
  N50_OS3_med_condition[i] <- mean(OS3_med$condition_number)
  N50_OS3_large_condition[i] <- mean(OS3_large$condition_number)
  N50_OS3_very_large_condition[i] <- mean(OS3_very_large$condition_number)
}

### When N = 100:
N = 100
set.seed(123456)
x <- seq(0,20, length.out = N)
gx <- sqrt(3)*sin(x/2)
y <- gx + rnorm(n = length(x), sd = 1)
z <- NULL

N100_IWP3 <- c()
N100_OS3_small <- c()
N100_OS3_med <- c()
N100_OS3_large <- c()
N100_OS3_very_large <- c()

N100_IWP3_condition <- c()
N100_OS3_small_condition <- c()
N100_OS3_med_condition <- c()
N100_OS3_large_condition <- c()
N100_OS3_very_large_condition <- c()
for (i in 1:10) {
  predictive_prior <- list(u = 3, alpha = 0.01)
  prior_set_1 <- prior_order_conversion_predictive(d = 5, predictive_prior, p = p)
  prior_set <- list(
    u1 = prior_set_1$u,
    alpha1 = prior_set_1$alpha,
    sigma2 = 1,
    betaprec = 10^(-3)
  )
  IWP_fit <- fit_once_time(method = "IWP", y = y, x = x, z = z, prior = prior_set, p = p)
  OS3_small <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 10, p = p)
  OS3_med <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 30, p = p)
  OS3_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 50, p = p)
  OS3_very_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 100, p = p)
  
  N100_IWP3[i] <- IWP_fit$runtime
  N100_OS3_small[i] <- OS3_small$runtime
  N100_OS3_med[i] <- OS3_med$runtime
  N100_OS3_large[i] <- OS3_large$runtime
  N100_OS3_very_large[i] <- OS3_very_large$runtime
  
  N100_IWP3_condition[i] <- mean(IWP_fit$condition_number)
  N100_OS3_small_condition[i] <- mean(OS3_small$condition_number)
  N100_OS3_med_condition[i] <- mean(OS3_med$condition_number)
  N100_OS3_large_condition[i] <- mean(OS3_large$condition_number)
  N100_OS3_very_large_condition[i] <- mean(OS3_very_large$condition_number)
}

### When N = 200:
N = 200
set.seed(123456)
x <- seq(0,20, length.out = N)
gx <- sqrt(3)*sin(x/2)
y <- gx + rnorm(n = length(x), sd = 1)
z <- NULL

N200_IWP3 <- c()
N200_OS3_small <- c()
N200_OS3_med <- c()
N200_OS3_large <- c()
N200_OS3_very_large <- c()

N200_IWP3_condition <- c()
N200_OS3_small_condition <- c()
N200_OS3_med_condition <- c()
N200_OS3_large_condition <- c()
N200_OS3_very_large_condition <- c()
for (i in 1:10) {
  predictive_prior <- list(u = 3, alpha = 0.01)
  prior_set_1 <- prior_order_conversion_predictive(d = 5, predictive_prior, p = p)
  prior_set <- list(
    u1 = prior_set_1$u,
    alpha1 = prior_set_1$alpha,
    sigma2 = 1,
    betaprec = 10^(-3)
  )
  IWP_fit <- fit_once_time(method = "IWP", y = y, x = x, z = z, prior = prior_set, p = p)
  OS3_small <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 10, p = p)
  OS3_med <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 30, p = p)
  OS3_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 50, p = p)
  OS3_very_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 100, p = p)
  
  N200_IWP3[i] <- IWP_fit$runtime
  N200_OS3_small[i] <- OS3_small$runtime
  N200_OS3_med[i] <- OS3_med$runtime
  N200_OS3_large[i] <- OS3_large$runtime
  N200_OS3_very_large[i] <- OS3_very_large$runtime
  
  N200_IWP3_condition[i] <- max(IWP_fit$condition_number)
  N200_OS3_small_condition[i] <- max(OS3_small$condition_number)
  N200_OS3_med_condition[i] <- max(OS3_med$condition_number)
  N200_OS3_large_condition[i] <- max(OS3_large$condition_number)
  N200_OS3_very_large_condition[i] <- max(OS3_very_large$condition_number)
  
  
}

### When N = 500:
N = 500
set.seed(12345)
x <- seq(0,20, length.out = N)
gx <- sqrt(3)*sin(x/2)
y <- gx + rnorm(n = length(x), sd = 1)
z <- NULL

N500_IWP3 <- c()
N500_OS3_small <- c()
N500_OS3_med <- c()
N500_OS3_large <- c()
N500_OS3_very_large <- c()

N500_IWP3_condition <- c()
N500_OS3_small_condition <- c()
N500_OS3_med_condition <- c()
N500_OS3_large_condition <- c()
N500_OS3_very_large_condition <- c()
for (i in 1:10) {
  predictive_prior <- list(u = 3, alpha = 0.01)
  prior_set_1 <- prior_order_conversion_predictive(d = 5, predictive_prior, p = p)
  prior_set <- list(
    u1 = prior_set_1$u,
    alpha1 = prior_set_1$alpha,
    sigma2 = 1,
    betaprec = 10^(-3)
  )
  IWP_fit <- fit_once_time(method = "IWP", y = y, x = x, z = z, prior = prior_set, p = p)
  OS3_small <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 10, p = p)
  OS3_med <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 30, p = p)
  OS3_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 50, p = p)
  OS3_very_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 100, p = p)
  
  N500_IWP3[i] <- IWP_fit$runtime
  N500_OS3_small[i] <- OS3_small$runtime
  N500_OS3_med[i] <- OS3_med$runtime
  N500_OS3_large[i] <- OS3_large$runtime
  N500_OS3_very_large[i] <- OS3_very_large$runtime
  
  N500_IWP3_condition[i] <- max(IWP_fit$condition_number)
  N500_OS3_small_condition[i] <- max(OS3_small$condition_number)
  N500_OS3_med_condition[i] <- max(OS3_med$condition_number)
  N500_OS3_large_condition[i] <- max(OS3_large$condition_number)
  N500_OS3_very_large_condition[i] <- max(OS3_very_large$condition_number)
}

### When N = 800:
N = 800
set.seed(1234)
x <- seq(0,20, length.out = N)
gx <- sqrt(3)*sin(x/2)
y <- gx + rnorm(n = length(x), sd = 1)
z <- NULL

N800_IWP3 <- c()
N800_OS3_small <- c()
N800_OS3_med <- c()
N800_OS3_large <- c()
N800_OS3_very_large <- c()

N800_IWP3_condition <- c()
N800_OS3_small_condition <- c()
N800_OS3_med_condition <- c()
N800_OS3_large_condition <- c()
N800_OS3_very_large_condition <- c()
for (i in 1:10) {
  predictive_prior <- list(u = 3, alpha = 0.01)
  prior_set_1 <- prior_order_conversion_predictive(d = 5, predictive_prior, p = p)
  prior_set <- list(
    u1 = prior_set_1$u,
    alpha1 = prior_set_1$alpha,
    sigma2 = 1,
    betaprec = 10^(-3)
  )
  IWP_fit <- fit_once_time(method = "IWP", y = y, x = x, z = z, prior = prior_set, p = p)
  # IWP_fit <- NULL
  OS3_small <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 10, p = p)
  OS3_med <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 30, p = p)
  OS3_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 50, p = p)
  OS3_very_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 100, p = p)
  N800_IWP3[i] <- IWP_fit$runtime
  # N800_IWP3[i] <- NULL
  N800_OS3_small[i] <- OS3_small$runtime
  N800_OS3_med[i] <- OS3_med$runtime
  N800_OS3_large[i] <- OS3_large$runtime
  N800_OS3_very_large[i] <- OS3_very_large$runtime
  
  N800_IWP3_condition[i] <- max(IWP_fit$condition_number)
  # N800_IWP3_condition[i] <- NULL
  N800_OS3_small_condition[i] <- max(OS3_small$condition_number)
  N800_OS3_med_condition[i] <- max(OS3_med$condition_number)
  N800_OS3_large_condition[i] <- max(OS3_large$condition_number)
  N800_OS3_very_large_condition[i] <- max(OS3_very_large$condition_number)
}

### When N = 2000:
N = 2000
set.seed(12345)
x <- seq(0,20, length.out = N)
gx <- sqrt(3)*sin(x/2)
y <- gx + rnorm(n = length(x), sd = 1)
z <- NULL

N2000_IWP3 <- c()
N2000_OS3_small <- c()
N2000_OS3_med <- c()
N2000_OS3_large <- c()
N2000_OS3_very_large <- c()

N2000_IWP3_condition <- c()
N2000_OS3_small_condition <- c()
N2000_OS3_med_condition <- c()
N2000_OS3_large_condition <- c()
N2000_OS3_very_large_condition <- c()
for (i in 1:10) {
  predictive_prior <- list(u = 3, alpha = 0.01)
  prior_set_1 <- prior_order_conversion_predictive(d = 5, predictive_prior, p = p)
  prior_set <- list(
    u1 = prior_set_1$u,
    alpha1 = prior_set_1$alpha,
    sigma2 = 1,
    betaprec = 10^(-3)
  )
  # IWP_fit <- fit_once_time(method = "IWP", y = y, x = x, z = z, prior = prior_set, p = p)
  IWP_fit <- NULL
  OS3_small <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 10, p = p)
  OS3_med <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 30, p = p)
  OS3_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 50, p = p)
  OS3_very_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 100, p = p)
  # N2000_IWP3[i] <- IWP_fit$runtime
  N2000_IWP3[i] <- NULL
  N2000_OS3_small[i] <- OS3_small$runtime
  N2000_OS3_med[i] <- OS3_med$runtime
  N2000_OS3_large[i] <- OS3_large$runtime
  N2000_OS3_very_large[i] <- OS3_very_large$runtime
  
  # N2000_IWP3_condition[i] <- max(IWP_fit$condition_number)
  N2000_IWP3_condition[i] <- NULL
  N2000_OS3_small_condition[i] <- max(OS3_small$condition_number)
  N2000_OS3_med_condition[i] <- max(OS3_med$condition_number)
  N2000_OS3_large_condition[i] <- max(OS3_large$condition_number)
  N2000_OS3_very_large_condition[i] <- max(OS3_very_large$condition_number)
}

### When N = 5000:
N = 5000
set.seed(12345)
x <- seq(0,20, length.out = N)
gx <- sqrt(3)*sin(x/2)
y <- gx + rnorm(n = length(x), sd = 1)
z <- NULL

N5000_IWP3 <- c()
N5000_OS3_small <- c()
N5000_OS3_med <- c()
N5000_OS3_large <- c()
N5000_OS3_very_large <- c()

N5000_IWP3_condition <- c()
N5000_OS3_small_condition <- c()
N5000_OS3_med_condition <- c()
N5000_OS3_large_condition <- c()
N5000_OS3_very_large_condition <- c()
for (i in 1:10) {
  predictive_prior <- list(u = 3, alpha = 0.01)
  prior_set_1 <- prior_order_conversion_predictive(d = 5, predictive_prior, p = p)
  prior_set <- list(
    u1 = prior_set_1$u,
    alpha1 = prior_set_1$alpha,
    sigma2 = 1,
    betaprec = 10^(-3)
  )
  # IWP_fit <- fit_once_time(method = "IWP", y = y, x = x, z = z, prior = prior_set, p = p)
  IWP_fit <- NULL
  OS3_small <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 10, p = p)
  OS3_med <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 30, p = p)
  OS3_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 50, p = p)
  OS3_very_large <- fit_once_time(method = "OS", y = y, x = x, z = z, prior = prior_set, k = 100, p = p)
  # N5000_IWP3[i] <- IWP_fit$runtime
  N5000_IWP3[i] <- NULL
  N5000_OS3_small[i] <- OS3_small$runtime
  N5000_OS3_med[i] <- OS3_med$runtime
  N5000_OS3_large[i] <- OS3_large$runtime
  N5000_OS3_very_large[i] <- OS3_very_large$runtime
  
  # N5000_IWP3_condition[i] <- max(IWP_fit$condition_number)
  N5000_IWP3_condition[i] <- NULL
  N5000_OS3_small_condition[i] <- max(OS3_small$condition_number)
  N5000_OS3_med_condition[i] <- max(OS3_med$condition_number)
  N5000_OS3_large_condition[i] <- max(OS3_large$condition_number)
  N5000_OS3_very_large_condition[i] <- max(OS3_very_large$condition_number)
}


### summary:

## N = 50
time_all_results_N50 <- cbind(N50_IWP3, N50_OS3_small, N50_OS3_med, N50_OS3_large, N50_OS3_very_large)
save(time_all_results_N50, file = paste0(data_location, "time_all_results_N50.Rda"))
# round(apply(time_all_results_N50, 2, mean),2)
# round(apply(time_all_results_N50, 2, sd),2)

condition_all_results_N50 <- cbind(N50_IWP3_condition, N50_OS3_small_condition, N50_OS3_med_condition, N50_OS3_large_condition, N50_OS3_very_large_condition)
# round(apply(condition_all_results_N50, 2, mean),2)
save(condition_all_results_N50, file = paste0(data_location, "condition_all_results_N50.Rda"))


## N = 100
time_all_results_N100 <- cbind(N100_IWP3, N100_OS3_small, N100_OS3_med, N100_OS3_large, N100_OS3_very_large)
save(time_all_results_N100, file = paste0(data_location, "time_all_results_N100.Rda"))
# round(apply(time_all_results_N100, 2, mean),2)
# round(apply(time_all_results_N100, 2, sd),2)

condition_all_results_N100 <- cbind(N100_IWP3_condition, N100_OS3_small_condition, N100_OS3_med_condition, N100_OS3_large_condition, N100_OS3_very_large_condition)
# round(apply(condition_all_results_N100, 2, mean),2)
save(condition_all_results_N100, file = paste0(data_location, "condition_all_results_N100.Rda"))

## N = 200
time_all_results_N200 <- cbind(N200_IWP3, N200_OS3_small, N200_OS3_med, N200_OS3_large, N200_OS3_very_large)
# round(apply(time_all_results_N200, 2, mean),2)
# round(apply(time_all_results_N200, 2, sd),2)
save(time_all_results_N200, file = paste0(data_location, "time_all_results_N200.Rda"))

condition_all_results_N200 <- cbind(N200_IWP3_condition, N200_OS3_small_condition, N200_OS3_med_condition, N200_OS3_large_condition, N200_OS3_very_large_condition)
# round(apply(condition_all_results_N200, 2, mean),2)
save(condition_all_results_N200, file = paste0(data_location, "condition_all_results_N200.Rda"))



## N = 500
time_all_results_N500 <- cbind(N500_IWP3, N500_OS3_small, N500_OS3_med, N500_OS3_large, N500_OS3_very_large)
# round(apply(time_all_results_N500, 2, mean),2)
# round(apply(time_all_results_N500, 2, sd),2)
save(time_all_results_N500, file = paste0(data_location, "time_all_results_N500.Rda"))

condition_all_results_N500 <- cbind(N500_IWP3_condition, N500_OS3_small_condition, N500_OS3_med_condition, N500_OS3_large_condition, N500_OS3_very_large_condition)
# round(apply(condition_all_results_N500, 2, mean),2)
save(condition_all_results_N500, file = paste0(data_location, "condition_all_results_N500.Rda"))



## N = 800
time_all_results_N800 <- cbind(N800_IWP3, N800_OS3_small, N800_OS3_med, N800_OS3_large, N800_OS3_very_large)
# round(apply(time_all_results_N800, 2, mean),2)
# round(apply(time_all_results_N800, 2, sd),2)
save(time_all_results_N800, file = paste0(data_location, "time_all_results_N800.Rda"))


condition_all_results_N800 <- cbind(N800_IWP3_condition, N800_OS3_small_condition, N800_OS3_med_condition, N800_OS3_large_condition, N800_OS3_very_large_condition)
# round(apply(condition_all_results_N800, 2, mean),2)
save(condition_all_results_N800, file = paste0(data_location, "condition_all_results_N800.Rda"))


## N = 2000
time_all_results_N2000 <- cbind(N2000_IWP3, N2000_OS3_small, N2000_OS3_med, N2000_OS3_large, N2000_OS3_very_large)
# round(apply(time_all_results_N2000, 2, mean),2)
# round(apply(time_all_results_N2000, 2, sd),2)
save(time_all_results_N2000, file = paste0(data_location, "time_all_results_N2000.Rda"))


condition_all_results_N2000 <- cbind(N2000_IWP3_condition, N2000_OS3_small_condition, N2000_OS3_med_condition, N2000_OS3_large_condition, N2000_OS3_very_large_condition)
# round(apply(condition_all_results_N2000, 2, mean),2)
save(condition_all_results_N2000, file = paste0(data_location, "condition_all_results_N2000.Rda"))



## N = 5000
time_all_results_N5000 <- cbind(N5000_IWP3, N5000_OS3_small, N5000_OS3_med, N5000_OS3_large, N5000_OS3_very_large)
# round(apply(time_all_results_N5000, 2, mean),2)
# round(apply(time_all_results_N5000, 2, sd),2)
save(time_all_results_N5000, file = paste0(data_location, "time_all_results_N5000.Rda"))


condition_all_results_N5000 <- cbind(N5000_IWP3_condition, N5000_OS3_small_condition, N5000_OS3_med_condition, N5000_OS3_large_condition, N5000_OS3_very_large_condition)
# round(apply(condition_all_results_N5000, 2, mean),2)
save(condition_all_results_N5000, file = paste0(data_location, "condition_all_results_N5000.Rda"))


###### Analysis:

## N = 50

load(paste0(data_location, "time_all_results_N50.Rda"))

reference = mean(time_all_results_N50[,2])

round(apply(time_all_results_N50/reference, 2, mean),2)
round(apply(time_all_results_N50/reference, 2, sd),2)

load(paste0(data_location, "condition_all_results_N50.Rda"))
round(apply(condition_all_results_N50, 2, mean),2)


## N = 100
load(paste0(data_location, "time_all_results_N100.Rda"))
round(apply(time_all_results_N100/reference, 2, mean),2)
round(apply(time_all_results_N100/reference, 2, sd),2)

load(paste0(data_location, "condition_all_results_N100.Rda"))
round(apply(condition_all_results_N100, 2, mean),2)


## N = 200
load(paste0(data_location, "time_all_results_N200.Rda"))
round(apply(time_all_results_N200/reference, 2, mean),2)
round(apply(time_all_results_N200/reference, 2, sd),2)


load(paste0(data_location, "condition_all_results_N200.Rda"))
round(apply(condition_all_results_N200, 2, mean),2)



## N = 500
load(paste0(data_location, "time_all_results_N500.Rda"))
round(apply(time_all_results_N500/reference, 2, mean),2)
round(apply(time_all_results_N500/reference, 2, sd),2)


load(paste0(data_location, "condition_all_results_N500.Rda"))
round(apply(condition_all_results_N500, 2, mean),2)



## N = 800
load(paste0(data_location, "time_all_results_N800.Rda"))
round(apply(time_all_results_N800/reference, 2, mean),2)
round(apply(time_all_results_N800/reference, 2, sd),2)

load(paste0(data_location, "condition_all_results_N800.Rda"))
round(apply(condition_all_results_N800, 2, mean),2)


## N = 2000
load(paste0(data_location, "time_all_results_N2000.Rda"))
round(apply(time_all_results_N2000/reference, 2, mean),2)
round(apply(time_all_results_N2000/reference, 2, sd),2)


load(paste0(data_location, "condition_all_results_N2000.Rda"))
round(apply(condition_all_results_N2000, 2, mean),2)



## N = 5000
load(paste0(data_location, "time_all_results_N5000.Rda"))
round(apply(time_all_results_N5000/reference, 2, mean),2)
round(apply(time_all_results_N5000/reference, 2, sd),2)

load(paste0(data_location, "condition_all_results_N5000.Rda"))
round(apply(condition_all_results_N5000, 2, mean),2)



