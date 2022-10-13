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
### Single Comparison:  ####################
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

### Obtain IWP3 samples
locations <- sort(unique(c(x)))
samps <- extract_samples_WP(model_IWP_exact, n_samp = 3000, p = p, x = locations)
IWP3_summary_samps <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                 upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                 lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                 sd = samps[,-1] %>% apply(1, sd))
IWP3_summary_samps$bias <- IWP3_summary_samps$mean - IWP3_summary_samps$mean

samps <- extract_deriv_samples_WP(model_IWP_exact, n_samp = 3000, p = p, x = locations, degree = 1)
IWP3_summary_samps_deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                       upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                       lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                       sd = samps[,-1] %>% apply(1, sd))
IWP3_summary_samps_deriv$bias <- IWP3_summary_samps_deriv$mean - IWP3_summary_samps_deriv$mean

samps <- extract_deriv_samples_WP(model_IWP_exact, n_samp = 3000, p = p, x = locations, degree = 2)
IWP3_summary_samps_2deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                        upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                        lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                        sd = samps[,-1] %>% apply(1, sd))
IWP3_summary_samps_2deriv$bias <- IWP3_summary_samps_2deriv$mean - IWP3_summary_samps_2deriv$mean


### Fit OS3 method with small number of knots:
set.seed(12345)
knots <- seq(0,20, length.out = 10)
model_OS3_small <- Imple_BayesRegression(y = y, x = x, knots = knots, aghq_k = aghq_k, prior = prior_set, p = p)

### Obtain OS3 samples
samps <- extract_samples_Ospline_refined(model_OS3_small, n_samp = 3000, p = p, x = knots, refined_x = locations)

OS3_small_summary_samps <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                      upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                      lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                      sd = samps[,-1] %>% apply(1, sd))
OS3_small_summary_samps$bias <- OS3_small_summary_samps$mean - IWP3_summary_samps$mean

samps <- extract_deriv_samples_OSpline(quad = model_OS3_small, n_samp = 3000, p = p, x = knots, refined_x = locations, degree = 1)

OS3_small_summary_samps_deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                            upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                            lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                            sd = samps[,-1] %>% apply(1, sd))
OS3_small_summary_samps_deriv$bias <- OS3_small_summary_samps_deriv$mean - IWP3_summary_samps_deriv$mean


samps <- extract_deriv_samples_OSpline(quad = model_OS3_small, n_samp = 3000, p = p, x = knots, refined_x = locations, degree = 2)

OS3_small_summary_samps_2deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                             upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                             lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                             sd = samps[,-1] %>% apply(1, sd))
OS3_small_summary_samps_2deriv$bias <- OS3_small_summary_samps_2deriv$mean - IWP3_summary_samps_2deriv$mean



### Fit OS3 method with medium number of knots:
set.seed(12345)
knots <- seq(0,20, length.out = 30)
model_OS3_med <- Imple_BayesRegression(y = y, x = x, knots = knots, aghq_k = aghq_k, prior = prior_set, p = p)


### Obtain OS3 samples
samps <- extract_samples_Ospline_refined(model_OS3_med, n_samp = 3000, p = p, x = knots, refined_x = locations)

OS3_med_summary_samps <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                    upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                    lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                    sd = samps[,-1] %>% apply(1, sd))
OS3_med_summary_samps$bias <- OS3_med_summary_samps$mean - IWP3_summary_samps$mean

samps <- extract_deriv_samples_OSpline(quad = model_OS3_med, n_samp = 3000, p = p, x = knots, refined_x = locations, degree = 1)

OS3_med_summary_samps_deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                          upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                          lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                          sd = samps[,-1] %>% apply(1, sd))
OS3_med_summary_samps_deriv$bias <- OS3_med_summary_samps_deriv$mean - IWP3_summary_samps_deriv$mean

samps <- extract_deriv_samples_OSpline(quad = model_OS3_med, n_samp = 3000, p = p, x = knots, refined_x = locations, degree = 2)

OS3_med_summary_samps_2deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                           upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                           lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                           sd = samps[,-1] %>% apply(1, sd))
OS3_med_summary_samps_2deriv$bias <- OS3_med_summary_samps_2deriv$mean - IWP3_summary_samps_2deriv$mean



### Fit OS3 method with large number of knots:
set.seed(12345)
knots <- seq(0,20, length.out = 50)
model_OS3_large <- Imple_BayesRegression(y = y, x = x, knots = knots, aghq_k = aghq_k, prior = prior_set, p = p)


### Obtain OS3 samples
samps <- extract_samples_Ospline_refined(model_OS3_large, n_samp = 3000, p = p, x = knots, refined_x = locations)

OS3_large_summary_samps <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                      upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                      lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                      sd = samps[,-1] %>% apply(1, sd))
OS3_large_summary_samps$bias <- OS3_large_summary_samps$mean - IWP3_summary_samps$mean


samps <- extract_deriv_samples_OSpline(quad = model_OS3_large, n_samp = 3000, p = p, x = knots, refined_x = locations, degree = 1)

OS3_large_summary_samps_deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                            upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                            lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                            sd = samps[,-1] %>% apply(1, sd))
OS3_large_summary_samps_deriv$bias <- OS3_large_summary_samps_deriv$mean - IWP3_summary_samps_deriv$mean

samps <- extract_deriv_samples_OSpline(quad = model_OS3_large, n_samp = 3000, p = p, x = knots, refined_x = locations, degree = 2)

OS3_large_summary_samps_2deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                             upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                             lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                             sd = samps[,-1] %>% apply(1, sd))
OS3_large_summary_samps_2deriv$bias <- OS3_large_summary_samps_2deriv$mean - IWP3_summary_samps_2deriv$mean



### Fit OS3 method with very large number of knots:
set.seed(12345)
knots <- seq(0,20, length.out = 100)
model_OS3_very_large <- Imple_BayesRegression(y = y, x = x, knots = knots, aghq_k = aghq_k, prior = prior_set, p = p)


### Obtain OS3 samples
samps <- extract_samples_Ospline_refined(model_OS3_very_large, n_samp = 3000, p = p, x = knots, refined_x = locations)

OS3_very_large_summary_samps <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                           upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                           lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                           sd = samps[,-1] %>% apply(1, sd))
OS3_very_large_summary_samps$bias <- OS3_very_large_summary_samps$mean - IWP3_summary_samps$mean

samps <- extract_deriv_samples_OSpline(quad = model_OS3_very_large, n_samp = 3000, p = p, x = knots, refined_x = locations, degree = 1)

OS3_very_large_summary_samps_deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                                 upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                                 lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                                 sd = samps[,-1] %>% apply(1, sd))
OS3_very_large_summary_samps_deriv$bias <- OS3_very_large_summary_samps_deriv$mean - IWP3_summary_samps_deriv$mean

samps <- extract_deriv_samples_OSpline(quad = model_OS3_very_large, n_samp = 3000, p = p, x = knots, refined_x = locations, degree = 2)

OS3_very_large_summary_samps_2deriv <- data.frame(x = samps[,1], mean = samps[,-1] %>% apply(1, mean),
                                                  upper = samps[,-1] %>% apply(1, quantile, 0.975),
                                                  lower = samps[,-1] %>% apply(1, quantile, 0.025),
                                                  sd = samps[,-1] %>% apply(1, sd))
OS3_very_large_summary_samps_2deriv$bias <- OS3_very_large_summary_samps_2deriv$mean - IWP3_summary_samps_2deriv$mean


### Plotting:
line_size1 <- 0.8
line_size2 <- 0.8
line_size3 <- 0.8
line_size4 <- 0.8
line_size5 <- 1.5


color <- c("exact" = "black", "k = 10" = "green", "k = 30" = "brown", "k = 50" = "red", "k = 100" = "purple")
lintype <- c("exact" = "solid", "k = 10" = "dotted","k = 30" = "dashed", "k = 50" = "longdash", "k = 100" = "dotdash")
IWP3_summary_samps %>% ggplot(aes(x = x, y = sd)) + geom_line(size=line_size1, aes(color = "exact", lty = "exact"), alpha = 1) +
  geom_line(size=line_size2, aes(x = x, y = sd), data = OS3_med_summary_samps, color = "brown", lty = "dashed", alpha = 1) +
  geom_line(size=line_size3, aes(x = x, y = sd), data = OS3_large_summary_samps, color = "red", lty = "longdash", alpha = 1) +
  geom_line(size=line_size4, aes(x = x, y = sd), data = OS3_very_large_summary_samps, color = "purple", lty = "dotdash", alpha = 1) +
  geom_line(size=line_size5, aes(x = x, y = sd), data = OS3_small_summary_samps, color = "green", lty = "dotted", alpha = 1) +
  ylab("SD(g(x)|y)") +
  scale_color_manual(values = color, name = "Method") + 
  scale_linetype_manual(values = lintype, name = "Method") + 
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.35, 0.75), legend.key.width = unit(2, 'cm')) +
  ylim(c(0.0,0.6))

ggsave(paste0(figures_location,"compare_pos_sd_p3.pdf"), device = "pdf", width = 5, height = 5)



color <- c("exact" = "black", "k = 10" = "green", "k = 30" = "brown", "k = 50" = "red", "k = 100" = "purple")
lintype <- c("exact" = "solid", "k = 10" = "dotted","k = 30" = "dashed", "k = 50" = "longdash", "k = 100" = "dotdash")
IWP3_summary_samps %>% ggplot(aes(x = x, y = bias)) + geom_line(size=line_size1, aes(color = "exact", lty = "exact")) +
  geom_line(size=line_size2, aes(x = x, y = bias), data = OS3_med_summary_samps, color = "brown", lty = "dashed") +
  geom_line(size=line_size3, aes(x = x, y = bias), data = OS3_large_summary_samps, color = "red", lty = "longdash") +
  geom_line(size=line_size4, aes(x = x, y = bias), data = OS3_very_large_summary_samps, color = "purple", lty = "dotdash") +
  geom_line(size=line_size5, aes(x = x, y = bias), data = OS3_small_summary_samps, color = "green", lty = "dotted") +
  ylab("Bias of E(g(x)|y)") +
  scale_color_manual(values = color, name = "Method") + 
  scale_linetype_manual(values = lintype, name = "Method") + 
  theme_bw(base_size = 15) + ylim(c(-0.25,0.25)) +
  theme(legend.position = "none")

ggsave(paste0(figures_location,"compare_pos_bias_p3.pdf"), device = "pdf", width = 5, height = 5)



color <- c("exact" = "black", "k = 10" = "green", "k = 30" = "brown", "k = 50" = "red", "k = 100" = "purple")
lintype <- c("exact" = "solid", "k = 10" = "dotted","k = 30" = "dashed", "k = 50" = "longdash", "k = 100" = "dotdash")
IWP3_summary_samps %>% ggplot(aes(x = x, y = mean)) + geom_line(size=line_size1, aes(color = "exact", lty = "exact")) +
  geom_line(size=line_size2, aes(x = x, y = mean), data = OS3_med_summary_samps, color = "brown", lty = "dashed") +
  geom_line(size=line_size3, aes(x = x, y = mean), data = OS3_large_summary_samps, color = "red", lty = "longdash") +
  geom_line(size=line_size4, aes(x = x, y = mean), data = OS3_very_large_summary_samps, color = "purple", lty = "dotdash") +
  geom_line(size=line_size5, aes(x = x, y = mean), data = OS3_small_summary_samps, color = "green", lty = "dotted") +
  ylab("E(g(x)|y)") +
  scale_color_manual(values = color, name = "Method") + 
  scale_linetype_manual(values = lintype, name = "Method") + 
  theme_bw(base_size = 15) + ylim(c(-2,2.3)) +
  theme(legend.position = "none")

ggsave(paste0(figures_location,"compare_pos_mean_p3.pdf"), device = "pdf", width = 5, height = 5)


color <- c("exact" = "black", "k = 10" = "green", "k = 30" = "brown", "k = 50" = "red", "k = 100" = "purple")
lintype <- c("exact" = "solid", "k = 10" = "dotted","k = 30" = "dashed", "k = 50" = "longdash", "k = 100" = "dotdash")
IWP3_summary_samps_deriv %>% ggplot(aes(x = x, y = sd)) + geom_line(size=line_size1, aes(color = "exact", lty = "exact")) +
  theme_bw(base_size = 15) +
  geom_line(size=line_size2, aes(x = x, y = sd), data = OS3_med_summary_samps_deriv, color = "brown", lty = "dashed") +
  geom_line(size=line_size3, aes(x = x, y = sd), data = OS3_large_summary_samps_deriv, color = "red", lty = "longdash") +
  geom_line(size=line_size4, aes(x = x, y = sd), data = OS3_very_large_summary_samps_deriv, color = "purple", lty = "dotdash") +
  geom_line(size=line_size5, aes(x = x, y = sd), data = OS3_small_summary_samps_deriv, color = "green", lty = "dotted") +
  ylab("SD(g'(x)|y)") +
  scale_color_manual(values = color, name = "Method") + 
  scale_linetype_manual(values = lintype, name = "Method") + 
  theme(legend.position = c(0.25, 0.75), legend.key.width = unit(2, 'cm')) +
  ylim(c(0.0,0.6))

ggsave(paste0(figures_location,"compare_pos_sd_p3_deriv.pdf"), device = "pdf", width = 5, height = 5)


color <- c("exact" = "black", "k = 10" = "green", "k = 30" = "brown", "k = 50" = "red", "k = 100" = "purple")
lintype <- c("exact" = "solid", "k = 10" = "dotted","k = 30" = "dashed", "k = 50" = "longdash", "k = 100" = "dotdash")
IWP3_summary_samps_deriv %>% ggplot(aes(x = x, y = bias)) + geom_line(size=line_size1, aes(color = "exact", lty = "exact")) +
  geom_line(size=line_size2, aes(x = x, y = bias), data = OS3_med_summary_samps_deriv, color = "brown", lty = "dashed") +
  geom_line(size=line_size3, aes(x = x, y = bias), data = OS3_large_summary_samps_deriv, color = "red", lty = "longdash") +
  geom_line(size=line_size4, aes(x = x, y = bias), data = OS3_very_large_summary_samps_deriv, color = "purple", lty = "dotdash") +
  geom_line(size=line_size5, aes(x = x, y = bias), data = OS3_small_summary_samps_deriv, color = "green", lty = "dotted") +
  ylab("Bias of E(g'(x)|y)") +
  scale_color_manual(values = color, name = "Method") + 
  scale_linetype_manual(values = lintype, name = "Method") + 
  theme_bw(base_size = 15) + ylim(c(-0.5,0.5)) +
  theme(legend.position = "none")

ggsave(paste0(figures_location,"compare_pos_bias_p3_deriv.pdf"), device = "pdf", width = 5, height = 5)


color <- c("exact" = "black", "k = 10" = "green", "k = 30" = "brown", "k = 50" = "red", "k = 100" = "purple")
lintype <- c("exact" = "solid", "k = 10" = "dotted","k = 30" = "dashed", "k = 50" = "longdash", "k = 100" = "dotdash")
IWP3_summary_samps_deriv %>% ggplot(aes(x = x, y = mean)) + geom_line(size=line_size1, aes(color = "exact", lty = "exact")) +
  geom_line(size=line_size2, aes(x = x, y = mean), data = OS3_med_summary_samps_deriv, color = "brown", lty = "dashed") +
  geom_line(size=line_size3, aes(x = x, y = mean), data = OS3_large_summary_samps_deriv, color = "red", lty = "longdash") +
  geom_line(size=line_size4, aes(x = x, y = mean), data = OS3_very_large_summary_samps_deriv, color = "purple", lty = "dotdash") +
  geom_line(size=line_size5, aes(x = x, y = mean), data = OS3_small_summary_samps_deriv, color = "green", lty = "dotted") +
  ylab("E(g'(x)|y)") +
  scale_color_manual(values = color, name = "Method") + 
  scale_linetype_manual(values = lintype, name = "Method") + 
  theme_bw(base_size = 15) + ylim(c(-2,1)) +
  theme(legend.position = "none")

ggsave(paste0(figures_location,"compare_pos_mean_p3_deriv.pdf"), device = "pdf", width = 5, height = 5)



color <- c("exact" = "black", "k = 10" = "green", "k = 30" = "brown", "k = 50" = "red", "k = 100" = "purple")
lintype <- c("exact" = "solid", "k = 10" = "dotted","k = 30" = "dashed", "k = 50" = "longdash", "k = 100" = "dotdash")
IWP3_summary_samps_2deriv %>% ggplot(aes(x = x, y = sd)) + geom_line(size=line_size1, aes(color = "exact", lty = "exact")) +
  geom_line(size=line_size2, aes(x = x, y = sd), data = OS3_med_summary_samps_2deriv, color = "brown", lty = "dashed") +
  geom_line(size=line_size3, aes(x = x, y = sd), data = OS3_large_summary_samps_2deriv, color = "red", lty = "longdash") +
  geom_line(size=line_size4, aes(x = x, y = sd), data = OS3_very_large_summary_samps_2deriv, color = "purple", lty = "dotdash") +
  geom_line(size=line_size5, aes(x = x, y = sd), data = OS3_small_summary_samps_2deriv, color = "green", lty = "dotted") +
  ylab("SD(g''(x)|y)") +
  scale_color_manual(values = color, name = "Method") + 
  scale_linetype_manual(values = lintype, name = "Method") + 
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.25, 0.75), legend.key.width = unit(2, 'cm')) +
  ylim(c(0.0,0.6))

ggsave(paste0(figures_location,"compare_pos_sd_p3_2deriv.pdf"), device = "pdf", width = 5, height = 5)


color <- c("exact" = "black", "k = 10" = "green", "k = 30" = "brown", "k = 50" = "red", "k = 100" = "purple")
lintype <- c("exact" = "solid", "k = 10" = "dotted","k = 30" = "dashed", "k = 50" = "longdash", "k = 100" = "dotdash")
IWP3_summary_samps_2deriv %>% ggplot(aes(x = x, y = bias)) + geom_line(size=line_size1, aes(color = "exact", lty = "exact")) +
  geom_line(size=line_size2, aes(x = x, y = bias), data = OS3_med_summary_samps_2deriv, color = "brown", lty = "dashed") +
  geom_line(size=line_size3, aes(x = x, y = bias), data = OS3_large_summary_samps_2deriv, color = "red", lty = "longdash") +
  geom_line(size=line_size4, aes(x = x, y = bias), data = OS3_very_large_summary_samps_2deriv, color = "purple", lty = "dotdash") +
  geom_line(size=line_size5, aes(x = x, y = bias), data = OS3_small_summary_samps_2deriv, color = "green", lty = "dotted") +
  ylab("Bias of E(g''(x)|y)") +
  scale_color_manual(values = color, name = "Method") + 
  scale_linetype_manual(values = lintype, name = "Method") + 
  theme_bw(base_size = 15) + ylim(c(-0.1,0.1)) +
  theme(legend.position = "none")

ggsave(paste0(figures_location,"compare_pos_bias_p3_2deriv.pdf"), device = "pdf", width = 5, height = 5)


color <- c("exact" = "black", "k = 10" = "green", "k = 30" = "brown", "k = 50" = "red", "k = 100" = "purple")
lintype <- c("exact" = "solid", "k = 10" = "dotted","k = 30" = "dashed", "k = 50" = "longdash", "k = 100" = "dotdash")
IWP3_summary_samps_2deriv %>% ggplot(aes(x = x, y = mean)) + geom_line(size=line_size1, aes(color = "exact", lty = "exact")) +
  geom_line(size=line_size2, aes(x = x, y = mean), data = OS3_med_summary_samps_2deriv, color = "brown", lty = "dashed") +
  geom_line(size=line_size3, aes(x = x, y = mean), data = OS3_large_summary_samps_2deriv, color = "red", lty = "longdash") +
  geom_line(size=line_size4, aes(x = x, y = mean), data = OS3_very_large_summary_samps_2deriv, color = "purple", lty = "dotdash") +
  geom_line(size=line_size5, aes(x = x, y = mean), data = OS3_small_summary_samps_2deriv, color = "green", lty = "dotted") +
  ylab("E(g''(x)|y)") +
  scale_color_manual(values = color, name = "Method") + 
  scale_linetype_manual(values = lintype, name = "Method") + 
  theme_bw(base_size = 15) + ylim(c(-0.6,0.5)) +
  theme(legend.position = "none")

ggsave(paste0(figures_location,"compare_pos_mean_p3_2deriv.pdf"), device = "pdf", width = 5, height = 5)
