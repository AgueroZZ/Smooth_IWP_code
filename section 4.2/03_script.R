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

font_size <- 1.3
axis_size <- 1.3

############################################
############################################
### Single Comparison:(KL) #################
############################################
############################################


##########################
###### N = 100 ###########
##########################
N = 100
set.seed(12345)
x <- seq(0,20, length.out = N)

gx <- sqrt(3)*sin(x/2)
y <- gx + rnorm(n = length(x), sd = 1)
plot(y ~ x)
lines(gx ~ x, lty = 'dashed')

predictive_prior <- list(u = 3, alpha = 0.01)
prior_set_1 <- prior_order_conversion_predictive(d = 5, predictive_prior, p = p)
prior_set <- list(
  u1 = prior_set_1$u,
  alpha1 = prior_set_1$alpha,
  sigma2 = 1,
  betaprec = 10^(-3)
)

### Fit IWP3 with augmented method:
model_IWP_exact <- Imple_IWP_augment(y = y, x = x, z = NULL, prior = prior_set, starting_val = 0, p = p, aghq_k = aghq_k)

############################## 
############################## 
### KL comparison: ###########
############################## 
############################## 

### For PSD:
logpost_exact <- compute_pdf_and_cdf(model_IWP_exact$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
summ_exact <- logpost_exact %>% select(transparam, pdf_transparam)
colnames(summ_exact) <- c("PSD", "density")
summ_exact$density <- summ_exact$density/sum(summ_exact$density)
k_vec <- seq(from = 2, to = 50, by = 2)
KL_vec <- c()
for (i in 1:length(k_vec)) {
  knots <- seq(0,20, length.out = (k_vec[i]+1))
  model_OS3_k <- Imple_BayesRegression(y = y, x = x, knots = knots, aghq_k = aghq_k, prior = prior_set, p = p)
  logpost_k <- compute_pdf_and_cdf(model_OS3_k$marginals[[1]],list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline', finegrid = logpost_exact$theta)
  summ_k <- logpost_k %>% select(transparam, pdf_transparam)
  colnames(summ_k) <- c("PSD", "density")
  summ_k$density <- summ_k$density/sum(summ_k$density)
  KL_vec[i] <- suppressMessages(philentropy::KL(x = rbind(summ_exact$density, summ_k$density)))
}
pdf(paste0(figures_location, "/KL_PSD.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.6, 0.5, 0))
plot(KL_vec~k_vec, type = "o", xlab = "K", ylab = "KL", ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size)
dev.off()

### For g(x50):
### Obtain IWP3 samples
locations <- sort(unique(c(x)))
samps_IWP_f <- extract_samples_WP(model_IWP_exact, n_samp = 5000, p = p, x = locations)
samps_IWP_f50 <- as.numeric(samps_IWP_f %>% filter(x == locations[50]))[-1]
ecdf_IWP <- ecdf(samps_IWP_f50)
emp_dist_IWP <- data.frame(x = seq(from = -2.5, to = -0.5, length.out = 20))
emp_dist_IWP$cdf <- ecdf_IWP(emp_dist_IWP$x)
emp_dist_IWP$pmf <- c(0, diff(emp_dist_IWP$cdf))

k_vec <- seq(from = 2, to = 50, by = 2)
KL_vec <- c()
for (i in 1:length(k_vec)) {
  knots <- unique(sort(c(seq(0,20, length.out = k_vec[i]),locations[50])))
  model_OS3_k <- Imple_BayesRegression(y = y, x = x, knots = knots, aghq_k = aghq_k, prior = prior_set, p = p)
  samps_k_f <- extract_samples_Ospline_refined(model_OS3_k, n_samp = 5000, p = p, x = knots, refined_x = locations)
  samps_k_f50 <- as.numeric(samps_k_f %>% filter(x == locations[50]))[-1]
  ecdf_k <- ecdf(samps_k_f50)
  emp_dist_k <- data.frame(x = seq(from = -2.5, to = -0.5, length.out = 20))
  emp_dist_k$cdf <- ecdf_k(emp_dist_k$x)
  emp_dist_k$pmf <- c(0, diff(emp_dist_k$cdf))
  KL_vec[i] <- suppressMessages(philentropy::KL(x = rbind(emp_dist_IWP$pmf, emp_dist_k$pmf)))
}

pdf(paste0(figures_location, "/KL_g.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.6, 0.5, 0))
plot(KL_vec~k_vec, type = "o", xlab = "K", ylab = "KL", ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size)
dev.off()

### For g'(x50):
### Obtain IWP3 samples
samps_IWP_f1st <- extract_deriv_samples_WP(model_IWP_exact, n_samp = 5000, p = p, x = locations, degree = 1)
samps_IWP_f1st50 <- as.numeric(samps_IWP_f1st %>% filter(x == locations[50]))[-1]
ecdf_IWP <- ecdf(samps_IWP_f1st50)
emp_dist_IWP <- data.frame(x = seq(from = -1, to = 1, length.out = 20))
emp_dist_IWP$cdf <- ecdf_IWP(emp_dist_IWP$x)
emp_dist_IWP$pmf <- c(0, diff(emp_dist_IWP$cdf))

k_vec <- seq(from = 2, to = 50, by = 2)
KL_vec <- c()
for (i in 1:length(k_vec)) {
  knots <- unique(sort(c(seq(0,20, length.out = k_vec[i]),locations[50])))
  model_OS3_k <- Imple_BayesRegression(y = y, x = x, knots = knots, aghq_k = aghq_k, prior = prior_set, p = p)
  samps_k_f1st <- extract_deriv_samples_OSpline(quad = model_OS3_k, n_samp = 5000, p = p, x = knots, refined_x = locations, degree = 1)
  samps_k_f1st50 <- as.numeric(samps_k_f1st %>% filter(x == locations[50]))[-1]
  ecdf_k <- ecdf(samps_k_f1st50)
  emp_dist_k <- data.frame(x = seq(from = -1, to = 1, length.out = 20))
  emp_dist_k$cdf <- ecdf_k(emp_dist_k$x)
  emp_dist_k$pmf <- c(0, diff(emp_dist_k$cdf))
  KL_vec[i] <- suppressMessages(philentropy::KL(x = rbind(emp_dist_IWP$pmf, emp_dist_k$pmf)))
}
pdf(paste0(figures_location, "/KL_g1st.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.6, 0.5, 0))
plot(KL_vec~k_vec, type = "o", xlab = "K", ylab = "KL", ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size)
dev.off()

### For g''(x50):
### Obtain IWP3 samples
samps_IWP_f2nd <- extract_deriv_samples_WP(model_IWP_exact, n_samp = 5000, p = p, x = locations, degree = 2)
samps_IWP_f2nd50 <- as.numeric(samps_IWP_f2nd %>% filter(x == locations[50]))[-1]
ecdf_IWP <- ecdf(samps_IWP_f2nd50)
emp_dist_IWP <- data.frame(x = seq(from = -0.1, to = 1.0, length.out = 20))
emp_dist_IWP$cdf <- ecdf_IWP(emp_dist_IWP$x)
emp_dist_IWP$pmf <- c(0, diff(emp_dist_IWP$cdf))

k_vec <- seq(from = 2, to = 50, by = 2)
KL_vec <- c()
for (i in 1:length(k_vec)) {
  knots <- unique(sort(c(seq(0,20, length.out = k_vec[i]),locations[50])))
  model_OS3_k <- Imple_BayesRegression(y = y, x = x, knots = knots, aghq_k = aghq_k, prior = prior_set, p = p)
  samps_k_f2nd <- extract_deriv_samples_OSpline(quad = model_OS3_k, n_samp = 5000, p = p, x = knots, refined_x = locations, degree = 2)
  samps_k_f2nd50 <- as.numeric(samps_k_f2nd %>% filter(x == locations[50]))[-1]
  ecdf_k <- ecdf(samps_k_f2nd50)
  emp_dist_k <- data.frame(x = seq(from = -0.1, to = 1.0,length.out = 20))
  emp_dist_k$cdf <- ecdf_k(emp_dist_k$x)
  emp_dist_k$pmf <- c(0, diff(emp_dist_k$cdf))
  KL_vec[i] <- suppressMessages(philentropy::KL(x = rbind(emp_dist_IWP$pmf, emp_dist_k$pmf)))
}
pdf(paste0(figures_location, "/KL_g2nd.pdf"), width = 5, height = 5)
par(mar = c(2.5, 3, 0.1, 0.1), mgp = c(1.6, 0.5, 0))
plot(KL_vec~k_vec, type = "o", xlab = "K", ylab = "KL", ylim = c(0,1), cex.lab=font_size, cex.axis=axis_size)
dev.off()