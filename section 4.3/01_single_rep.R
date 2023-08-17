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

### Try to simulate mixture normal density:
set.seed(1234) # set.seed(10)
pvec <- c(0.6,0.3,0.1)
muvec <- rnorm(n = length(pvec), mean = 5, sd = 2)
# SDvec <- abs(rnorm(n = length(pvec), mean = 1, sd = 1))
SDvec <- rep(1, length(pvec))
x <- seq(from = 0, to = 10, by = 0.1) 
gx <- dnormm(x = x, p = pvec, mu = muvec, sigma = SDvec)
gx <- gx/sd(gx)
plot(gx ~ x, type = 'l')
gx1st <- diff(gx)/diff(x)
gx1st <- c(gx1st[1], gx1st)
gx2nd <- diff(gx1st)/diff(x)
gx2nd <- c(gx2nd[1], gx2nd)
plot(gx1st ~ x, type = 'l')
SNR <- 100
y <- gx + rnorm(n = length(x), sd = sqrt(var(gx)/SNR))
plot(y ~ x)
lines(gx ~ x, lty = 'dashed')

### Fit OS3, RW2::
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

knots <- seq(from = 0, to = 10, length.out = 100) 
OS3_fit <- Imple_BayesRegression(y = y, x = x, knots = knots, p = 3, prior.type = "PC", prior = prior_set)
RW2_fit <- Imple_RW2(y = y, x = x, prior = prior_set_2, knots = knots)

### Compare performance: check rMSE for posterior mean of g,g'; check rMSE for argmax

OS_samps <- extract_samples_Ospline_refined(OS3_fit, x = knots, refined_x = x, p = 3)
RW_samps <- extract_RW_samples(RW2_fit, x = knots, refined_x = x)
OS_deriv_samps <- extract_deriv_samples_OSpline(OS3_fit, p = 3, x = knots, refined_x = x, degree = 1)
RW_deriv_samps <- RW_samps[,-1] %>% apply(MARGIN = 2, FUN = compute_numeric_deriv, h = diff(x)) 
RW_deriv_samps <- rbind(RW_deriv_samps[1,],RW_deriv_samps)
RW_deriv_samps <- cbind(x,RW_deriv_samps)

OS_2nd_deriv_samps <- extract_deriv_samples_OSpline(OS3_fit, p = 3, x = knots, refined_x = x, degree = 2)
RW_2nd_deriv_samps <- RW_deriv_samps[,-1] %>% apply(MARGIN = 2, FUN = compute_numeric_deriv, h = diff(x), degree = 1) 
RW_2nd_deriv_samps <- rbind(RW_2nd_deriv_samps[1,],RW_2nd_deriv_samps)
RW_2nd_deriv_samps <- cbind(x,RW_2nd_deriv_samps)

OS_g_mean <- OS_samps[,-1] %>% apply(MARGIN = 1, mean)
RW_g_mean <- RW_samps[,-1] %>% apply(MARGIN = 1, mean)
OS_g1st_mean <- OS_deriv_samps[,-1] %>% apply(MARGIN = 1, mean)
RW_g1st_mean <- RW_deriv_samps[,-1] %>% apply(MARGIN = 1, mean)
OS_g2nd_mean <- OS_2nd_deriv_samps[,-1] %>% apply(MARGIN = 1, mean)
RW_g2nd_mean <- RW_2nd_deriv_samps[,-1] %>% apply(MARGIN = 1, mean)

plot(OS_g_mean ~ x, type = 'l', col = 'red')
lines(RW_g_mean ~ x, type = 'l', col = 'blue')
points(y ~ x)
lines(gx ~ x, lty = 'dashed')

sqrt(mean((OS_g_mean - gx)^2))
sqrt(mean((RW_g_mean - gx)^2))


plot(OS_g1st_mean ~ x, type = 'l', col = 'red')
lines(RW_g1st_mean ~ x, type = 'l', col = 'blue')
lines(gx1st~x, lty = 'dashed')

sqrt(mean((OS_g1st_mean - gx1st)^2))
sqrt(mean((RW_g1st_mean - gx1st)^2))


plot(OS_g1st_mean ~ x, type = 'l', col = 'red')
lines(RW_g1st_mean ~ x, type = 'l', col = 'blue')
lines(gx1st~x, lty = 'dashed')

sqrt(mean((OS_g1st_mean - gx1st)^2))
sqrt(mean((RW_g1st_mean - gx1st)^2))


plot(OS_g2nd_mean ~ x, type = 'l', col = 'red', ylim = c(-5,3))
lines(RW_g2nd_mean ~ x, type = 'l', col = 'blue')
lines(gx2nd~x, lty = 'dashed')

sqrt(mean((OS_g2nd_mean - gx2nd)^2))
sqrt(mean((RW_g2nd_mean - gx2nd)^2))


#### Plotting of a single replication:
line_size <- 1

### for g:
g_plot <- data.frame(mean = c(OS_g_mean, RW_g_mean, gx), 
                     method = rep(c("OS", "RW", "Truth"), each = length(x)),
                     x = rep(x, 3))

g_plot$order <- c(rep(3, length(x)),rep(1, length(x)), rep(2, length(x)))

g_plot %>% ggplot(aes(x = x, y = mean, color = method)) + 
  geom_line(size=line_size, aes(linetype = method, alpha = method, order = order)) +
  theme_classic(base_size = 15) +
  scale_linetype_manual(values=c("solid", "dashed", "dotted")) + theme(legend.position="none") +
  scale_color_manual(values = c("red", "blue", "black")) +
  scale_alpha_manual(values = c(1,0.5,1))

ggsave(paste0(figure_path,"single_g_compare.pdf"), device = "pdf", width = 5, height = 5)

### for g':
g1st_plot <- data.frame(mean = c(OS_g1st_mean, RW_g1st_mean, gx1st), 
                     method = rep(c("OS", "RW", "Truth"), each = length(x)),
                     x = rep(x, 3))
g1st_plot$order <- c(rep(3, length(x)),rep(1, length(x)), rep(2, length(x)))

g1st_plot %>% ggplot(aes(x = x, y = mean, color = method)) + 
  geom_line(size=line_size, aes(linetype = method, alpha = method, order = order)) +
  theme_classic(base_size = 15) + theme(legend.position="none") +
  scale_linetype_manual(values=c("solid", "dashed", "dotted")) + 
  scale_color_manual(values = c("red", "blue", "black")) +
  scale_alpha_manual(values = c(1,0.5,1)) + 
  xlim(c(0.1,10))

ggsave(paste0(figure_path,"single_g1st_compare.pdf"), device = "pdf", width = 5, height = 5)

### for g'':
g2nd_plot <- data.frame(mean = c(OS_g2nd_mean, RW_g2nd_mean, gx2nd), 
                        method = rep(c("OS", "RW", "Truth"), each = length(x)),
                        x = rep(x, 3))
g2nd_plot$order <- c(rep(3, length(x)),rep(1, length(x)), rep(2, length(x)))
levels(as.factor(g2nd_plot$method))
g2nd_plot$method <- factor(g2nd_plot$method, levels = c("Truth","OS","RW"))

d2 <- g2nd_plot %>% ggplot(aes(x = x, y = mean, color = method)) +
  geom_line(size=line_size, aes(linetype = method, alpha = method, order = order)) +
  theme_classic(base_size = 15) +
  scale_linetype_manual(values=c("dotted", "solid", "dashed")) + 
  scale_color_manual(values = c("black","red", "blue")) +
  scale_alpha_manual(values = c(1,1,0.5)) + 
  xlim(c(0.2,10))

d2 <- lemon::reposition_legend(d2, "right bottom")
ggsave(paste0(figure_path,"single_g2nd_compare.pdf"), device = "pdf", width = 5, height = 5)



