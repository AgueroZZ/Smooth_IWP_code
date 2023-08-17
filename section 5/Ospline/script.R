working_path <- paste0(getwd(),"/")
figure_path <- paste0(working_path,"figures/")
library(OSplines)
library(tidyverse)
library(aghq)
library(TMB)

full_data_covid <- read.csv(file = paste0(working_path, "owid-covid-data.csv"), header = T)
full_data_covid <- full_data_covid %>% filter(date >= "2020-03-01")
full_data_covid$new_deaths[is.na(full_data_covid$new_deaths)] <- 0
full_data_covid$new_deaths <- round(full_data_covid$new_deaths)
full_data_covid$Date <- as.Date(full_data_covid$date)
full_data_covid$t <- (as.Date(full_data_covid$date) %>% as.numeric())/31
full_data_covid$weekdays <- weekdays(as.Date(full_data_covid$date))
full_data_covid$weekdays <- factor(full_data_covid$weekdays,
                                    levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"),
                                    ordered = F)
prior_set_conversion <- prior_conversion(d = (7/31), prior = list(u = log(2), alpha = 0.5), p = 3)
prior_set_overdis <- list(u1 = prior_set_conversion$u,
                          alpha1 = prior_set_conversion$a,
                          u2 = 0.1,
                          alpha2 = 0.5,
                          betaprec = 0.01,
                          Xfprec = 0.01)

###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
######################################################################     1. Canada Analysis:   ##############################################################################
##################################################a############################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################

### Processing data:
full_data_canada <- full_data_covid %>% filter(location == "Canada") %>% dplyr::select(Date, date, new_deaths, t, weekdays)
Xf <- model.matrix(~ weekdays, data = full_data_canada, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
final_data_canada <- cbind(full_data_canada, Xf)
final_data_canada$index <- 1:nrow(final_data_canada)

### Fitting model:
set.seed(123456)
OS3mod_Canada <- model_fit(formula = new_deaths ~ weekdays1 + weekdays2 + weekdays3 + weekdays4 + weekdays5 + weekdays6 + 
                             f(smoothing_var = t,
                               model = "IWP",
                               order = 3, k = 100, 
                               boundary.prior = list(prec = 0.01), 
                               sd.prior = list(prior = "exp", para = list(u = prior_set_conversion$u, alpha = prior_set_conversion$a))) + 
                             f(smoothing_var = index,
                               model = "IID",
                               sd.prior = list(prior = "exp", para = list(u = 0.1, alpha = 0.5))),
                           data = final_data_canada, family = "Poisson",
                           aghq_k = 10)

### Fitted Model:
IWP_part <- OS3mod_Canada$instances[[1]]
samps <- aghq::sample_marginal(OS3mod_Canada$mod, M = 3000)
global_samps <- samps$samps[OS3mod_Canada$boundary_samp_indexes$t, , drop = F]
coefsamps <- samps$samps[OS3mod_Canada$random_samp_indexes$t, , drop = F]
fixed_samps <- samps$samps[OS3mod_Canada$fixed_samp_indexes$Intercept, , drop = F]

f_can <- compute_post_fun(
  samps = coefsamps, global_samps = global_samps,
  knots = IWP_part@knots,
  refined_x = seq(min(IWP_part@observed_x), max(IWP_part@observed_x), length.out = 300),
  p = 3, degree = 0, 
  intercept_samps = fixed_samps
)
f_can_sum <- extract_mean_interval_given_samps(f_can)
f_can_sum$Date <- as.Date((f_can_sum$x)*31, origin = min(full_data_canada$Date))

f1st_can <- compute_post_fun(
  samps = coefsamps, global_samps = global_samps,
  knots = IWP_part@knots,
  refined_x = seq(min(IWP_part@observed_x), max(IWP_part@observed_x), length.out = 300),
  p = 3, degree = 1, 
  intercept_samps = fixed_samps
)
f1st_can_sum <- extract_mean_interval_given_samps(f1st_can)
f1st_can_sum$Date <- as.Date((f1st_can_sum$x)*31, origin = min(full_data_canada$Date))
g1st_can_sum <- extract_mean_interval_given_samps(cbind(f_can_sum$x, (exp(f_can[,-1]) * f1st_can[,-1])/31))
g1st_can_sum$Date <- as.Date((g1st_can_sum$x)*31, origin = min(full_data_canada$Date))

### Plotting:
full_data_canada %>% ggplot(aes(x = Date, y = log(new_deaths))) + geom_point(aes(color = weekdays), alpha = 0.4) + 
  theme_classic(base_size = 15) + ylab("Log New Deaths") + xlab("") + 
  theme(legend.position = "none"
  ) +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = (c(1,7)))
ggsave(paste0(figure_path, "canada_covid_death_raw.pdf"), device = "pdf", width = 5, height = 5)

CANADA_POP <- 38.01 ## millions
full_data_canada %>% ggplot(aes(x = Date)) + geom_point(aes(y = new_deaths/CANADA_POP), alpha = 0.3, size = 0.5) + 
  geom_line(aes(y = exp(mean)/CANADA_POP), data = f_can_sum, color = "blue") + 
  geom_ribbon(aes(ymax = exp(pupper)/CANADA_POP, ymin = exp(plower)/CANADA_POP), fill = "orange", alpha = 0.3, data = f_can_sum) +
  theme_classic(base_size = 15) + ylab("New Deaths (per million)") + xlab("") + guides(color= "none") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = c(0,9))
ggsave(paste0(figure_path, "canada_covid_death_adjusted.pdf"), device = "pdf", width = 5, height = 5)

full_data_canada %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = (mean)), data = f_can_sum, color = "blue") + 
  geom_ribbon(aes(ymax = (pupper), ymin = (plower)), fill = "orange", alpha = 0.3, data = f_can_sum) +
  theme_classic(base_size = 15) + ylab("g(t)") + xlab("") + guides(color= "none") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = c(-2,6))
ggsave(paste0(figure_path, "canada_covid_death_log.pdf"), device = "pdf", width = 5, height = 5)

g1st_can_sum %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean/CANADA_POP, color = "OS3"), linewidth = 1) + 
  geom_ribbon(aes(ymax = pupper/CANADA_POP, ymin = plower/CANADA_POP), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("Change in daily death (per million)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = c(-0.3,0.3))
ggsave(paste0(figure_path, "canada_covid_death_deriv_adjusted.pdf"), device = "pdf", width = 5, height = 5)

f1st_can_sum %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean, color = "OS3"), linewidth = 1) + 
  geom_ribbon(aes(ymax = pupper, ymin = plower), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("g'(t)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b") +
  coord_cartesian(ylim = c(-3,5))
ggsave(paste0(figure_path, "canada_covid_death_deriv_log.pdf"), device = "pdf", width = 5, height = 5)

week_day_can <- samps$samps[unlist(OS3mod_Canada$fixed_samp_indexes[-1]), , drop = F]
week_day_Sunday <- apply(-week_day_can, 2, sum); week_day_can <- rbind(week_day_can, week_day_Sunday)
rownames(week_day_can) <- c(colnames(Xf), "weekdays7")
week_mean <- week_day_can %>% apply(MARGIN = 1, mean); week_sd <- week_day_can %>% apply(MARGIN = 1, sd)
week_upper <- week_mean + 1.96 * week_sd; week_lower <- week_mean - 1.96 * week_sd
week_data_CA <- data.frame(mean = week_mean, sd = week_sd, upper = week_upper, lower = week_lower, 
                           days = substr(rownames(week_day_can), start = 9, stop = nchar(rownames(week_day_can))))


###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
######################################################################     2. Denmark Analysis:   #############################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################

### Processing data:
full_data_DK <- full_data_covid %>% filter(location == "Denmark") %>% dplyr::select(Date, date, new_deaths, t, weekdays)
Xf <- model.matrix(~ weekdays, data = full_data_DK, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
final_data_DK <- cbind(full_data_DK, Xf)
final_data_DK$index <- 1:nrow(final_data_DK)

### Fitting model:
set.seed(123456)
OS3mod_DK <- model_fit(formula = new_deaths ~ weekdays1 + weekdays2 + weekdays3 + weekdays4 + weekdays5 + weekdays6 + 
                             f(smoothing_var = t,
                               model = "IWP",
                               order = 3, k = 100, 
                               boundary.prior = list(prec = 0.01), 
                               sd.prior = list(prior = "exp", para = list(u = prior_set_conversion$u, alpha = prior_set_conversion$a))) + 
                             f(smoothing_var = index,
                               model = "IID",
                               sd.prior = list(prior = "exp", para = list(u = 0.1, alpha = 0.5))),
                           data = final_data_DK, family = "Poisson",
                           aghq_k = 10)

### Fitted Model:
IWP_part <- OS3mod_DK$instances[[1]]
samps <- aghq::sample_marginal(OS3mod_DK$mod, M = 3000)
global_samps <- samps$samps[OS3mod_DK$boundary_samp_indexes$t, , drop = F]
coefsamps <- samps$samps[OS3mod_DK$random_samp_indexes$t, , drop = F]
fixed_samps <- samps$samps[OS3mod_DK$fixed_samp_indexes$Intercept, , drop = F]

f_DK <- compute_post_fun(
  samps = coefsamps, global_samps = global_samps,
  knots = IWP_part@knots,
  refined_x = seq(min(IWP_part@observed_x), max(IWP_part@observed_x), length.out = 300),
  p = 3, degree = 0, 
  intercept_samps = fixed_samps
)
f_DK_sum <- extract_mean_interval_given_samps(f_DK)
f_DK_sum$Date <- as.Date((f_DK_sum$x)*31, origin = min(full_data_DK$Date))

f1st_DK <- compute_post_fun(
  samps = coefsamps, global_samps = global_samps,
  knots = IWP_part@knots,
  refined_x = seq(min(IWP_part@observed_x), max(IWP_part@observed_x), length.out = 300),
  p = 3, degree = 1, 
  intercept_samps = fixed_samps
)
f1st_DK_sum <- extract_mean_interval_given_samps(f1st_DK)
f1st_DK_sum$Date <- as.Date((f1st_DK_sum$x)*31, origin = min(full_data_DK$Date))
g1st_DK_sum <- extract_mean_interval_given_samps(cbind(f_DK_sum$x, (exp(f_DK[,-1]) * f1st_DK[,-1])/31))
g1st_DK_sum$Date <- as.Date((g1st_DK_sum$x)*31, origin = min(full_data_DK$Date))

### Plotting:
full_data_DK %>% ggplot(aes(x = Date, y = log(new_deaths))) + geom_point(aes(color = weekdays), alpha = 0.4) + 
  theme_classic(base_size = 15) + ylab("Log New Deaths") + xlab("") + 
  theme(legend.position = "none"
  ) +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = (c(1,7)))
ggsave(paste0(figure_path, "DK_covid_death_raw.pdf"), device = "pdf", width = 5, height = 5)

DK_POP <- 5.831 ## millions
full_data_DK %>% ggplot(aes(x = Date)) + geom_point(aes(y = new_deaths/DK_POP), alpha = 0.3, size = 0.5) + 
  geom_line(aes(y = exp(mean)/DK_POP), data = f_DK_sum, color = "blue") + 
  geom_ribbon(aes(ymax = exp(pupper)/DK_POP, ymin = exp(plower)/DK_POP), fill = "orange", alpha = 0.3, data = f_DK_sum) +
  theme_classic(base_size = 15) + ylab("New Deaths (per million)") + xlab("") + guides(color= "none") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = c(0,9))
ggsave(paste0(figure_path, "DK_covid_death_adjusted.pdf"), device = "pdf", width = 5, height = 5)

full_data_DK %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = (mean)), data = f_DK_sum, color = "blue") + 
  geom_ribbon(aes(ymax = (pupper), ymin = (plower)), fill = "orange", alpha = 0.3, data = f_DK_sum) +
  theme_classic(base_size = 15) + ylab("g(t)") + xlab("") + guides(color= "none") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = c(-2,6))
ggsave(paste0(figure_path, "DK_covid_death_log.pdf"), device = "pdf", width = 5, height = 5)

g1st_DK_sum %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean/DK_POP, color = "OS3"), linewidth = 1) + 
  geom_ribbon(aes(ymax = pupper/DK_POP, ymin = plower/DK_POP), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("Change in daily death (per million)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = c(-0.3,0.3))
ggsave(paste0(figure_path, "DK_covid_death_deriv_adjusted.pdf"), device = "pdf", width = 5, height = 5)

f1st_DK_sum %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean, color = "OS3"), linewidth = 1) + 
  geom_ribbon(aes(ymax = pupper, ymin = plower), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("g'(t)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b") +
  coord_cartesian(ylim = c(-3,5))
ggsave(paste0(figure_path, "DK_covid_death_deriv_log.pdf"), device = "pdf", width = 5, height = 5)

week_day_DK <- samps$samps[unlist(OS3mod_DK$fixed_samp_indexes[-1]), , drop = F]
week_day_Sunday <- apply(-week_day_DK, 2, sum); week_day_DK <- rbind(week_day_DK, week_day_Sunday)
rownames(week_day_DK) <- c(colnames(Xf), "weekdays7")
week_mean <- week_day_DK %>% apply(MARGIN = 1, mean); week_sd <- week_day_DK %>% apply(MARGIN = 1, sd)
week_upper <- week_mean + 1.96 * week_sd; week_lower <- week_mean - 1.96 * week_sd
week_data_DK <- data.frame(mean = week_mean, sd = week_sd, upper = week_upper, lower = week_lower, 
                           days = substr(rownames(week_day_DK), start = 9, stop = nchar(rownames(week_day_DK))))


###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
######################################################################   3. South Africa Analysis:   ########################################################################## 
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################

### Processing data:
full_data_SA <- full_data_covid %>% filter(location == "South Africa") %>% dplyr::select(Date, date, new_deaths, t, weekdays)
Xf <- model.matrix(~ weekdays, data = full_data_SA, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
final_data_SA <- cbind(full_data_SA, Xf)
final_data_SA$index <- 1:nrow(final_data_SA)

### Fitting model:
set.seed(123456)
OS3mod_SA <- model_fit(formula = new_deaths ~ weekdays1 + weekdays2 + weekdays3 + weekdays4 + weekdays5 + weekdays6 + 
                         f(smoothing_var = t,
                           model = "IWP",
                           order = 3, k = 100, 
                           boundary.prior = list(prec = 0.01), 
                           sd.prior = list(prior = "exp", para = list(u = prior_set_conversion$u, alpha = prior_set_conversion$a))) + 
                         f(smoothing_var = index,
                           model = "IID",
                           sd.prior = list(prior = "exp", para = list(u = 0.1, alpha = 0.5))),
                       data = final_data_SA, family = "Poisson",
                       aghq_k = 10)

### Fitted Model:
IWP_part <- OS3mod_SA$instances[[1]]
samps <- aghq::sample_marginal(OS3mod_SA$mod, M = 3000)
global_samps <- samps$samps[OS3mod_SA$boundary_samp_indexes$t, , drop = F]
coefsamps <- samps$samps[OS3mod_SA$random_samp_indexes$t, , drop = F]
fixed_samps <- samps$samps[OS3mod_SA$fixed_samp_indexes$Intercept, , drop = F]

f_SA <- compute_post_fun(
  samps = coefsamps, global_samps = global_samps,
  knots = IWP_part@knots,
  refined_x = seq(min(IWP_part@observed_x), max(IWP_part@observed_x), length.out = 300),
  p = 3, degree = 0, 
  intercept_samps = fixed_samps
)
f_SA_sum <- extract_mean_interval_given_samps(f_SA)
f_SA_sum$Date <- as.Date((f_SA_sum$x)*31, origin = min(full_data_SA$Date))

f1st_SA <- compute_post_fun(
  samps = coefsamps, global_samps = global_samps,
  knots = IWP_part@knots,
  refined_x = seq(min(IWP_part@observed_x), max(IWP_part@observed_x), length.out = 300),
  p = 3, degree = 1, 
  intercept_samps = fixed_samps
)
f1st_SA_sum <- extract_mean_interval_given_samps(f1st_SA)
f1st_SA_sum$Date <- as.Date((f1st_SA_sum$x)*31, origin = min(full_data_SA$Date))
g1st_SA_sum <- extract_mean_interval_given_samps(cbind(f_SA_sum$x, (exp(f_SA[,-1]) * f1st_SA[,-1])/31))
g1st_SA_sum$Date <- as.Date((g1st_SA_sum$x)*31, origin = min(full_data_SA$Date))

### Plotting:
full_data_SA %>% ggplot(aes(x = Date, y = log(new_deaths))) + geom_point(aes(color = weekdays), alpha = 0.4) + 
  theme_classic(base_size = 15) + ylab("Log New Deaths") + xlab("") + 
  theme(legend.position = "none"
  ) +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = (c(1,7)))
ggsave(paste0(figure_path, "SA_covid_death_raw.pdf"), device = "pdf", width = 5, height = 5)

SA_POP <- 59.31 ## millions
full_data_SA %>% ggplot(aes(x = Date)) + geom_point(aes(y = new_deaths/SA_POP), alpha = 0.3, size = 0.5) + 
  geom_line(aes(y = exp(mean)/SA_POP), data = f_SA_sum, color = "blue") + 
  geom_ribbon(aes(ymax = exp(pupper)/SA_POP, ymin = exp(plower)/SA_POP), fill = "orange", alpha = 0.3, data = f_SA_sum) +
  theme_classic(base_size = 15) + ylab("New Deaths (per million)") + xlab("") + guides(color= "none") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = c(0,9))
ggsave(paste0(figure_path, "SA_covid_death_adjusted.pdf"), device = "pdf", width = 5, height = 5)

full_data_SA %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = (mean)), data = f_SA_sum, color = "blue") + 
  geom_ribbon(aes(ymax = (pupper), ymin = (plower)), fill = "orange", alpha = 0.3, data = f_SA_sum) +
  theme_classic(base_size = 15) + ylab("g(t)") + xlab("") + guides(color= "none") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = c(-2,6))
ggsave(paste0(figure_path, "SA_covid_death_log.pdf"), device = "pdf", width = 5, height = 5)

g1st_SA_sum %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean/SA_POP, color = "OS3"), linewidth = 1) + 
  geom_ribbon(aes(ymax = pupper/SA_POP, ymin = plower/SA_POP), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("Change in daily death (per million)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = c(-0.3,0.3))
ggsave(paste0(figure_path, "SA_covid_death_deriv_adjusted.pdf"), device = "pdf", width = 5, height = 5)

f1st_SA_sum %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean, color = "OS3"), linewidth = 1) + 
  geom_ribbon(aes(ymax = pupper, ymin = plower), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("g'(t)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b") +
  coord_cartesian(ylim = c(-3,5))
ggsave(paste0(figure_path, "SA_covid_death_deriv_log.pdf"), device = "pdf", width = 5, height = 5)

week_day_SA <- samps$samps[unlist(OS3mod_SA$fixed_samp_indexes[-1]), , drop = F]
week_day_Sunday <- apply(-week_day_SA, 2, sum); week_day_SA <- rbind(week_day_SA, week_day_Sunday)
rownames(week_day_SA) <- c(colnames(Xf), "weekdays7")
week_mean <- week_day_SA %>% apply(MARGIN = 1, mean); week_sd <- week_day_SA %>% apply(MARGIN = 1, sd)
week_upper <- week_mean + 1.96 * week_sd; week_lower <- week_mean - 1.96 * week_sd
week_data_SA <- data.frame(mean = week_mean, sd = week_sd, upper = week_upper, lower = week_lower, 
                           days = substr(rownames(week_day_SA), start = 9, stop = nchar(rownames(week_day_SA))))



###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
######################################################################   4. South Korea Analysis:   ########################################################################### 
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################

### Processing data:
full_data_SK <- full_data_covid %>% filter(location == "South Korea") %>% dplyr::select(Date, date, new_deaths, t, weekdays)
Xf <- model.matrix(~ weekdays, data = full_data_SK, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
final_data_SK <- cbind(full_data_SK, Xf)
final_data_SK$index <- 1:nrow(final_data_SK)

### Fitting model:
set.seed(123456)
OS3mod_SK <- model_fit(formula = new_deaths ~ weekdays1 + weekdays2 + weekdays3 + weekdays4 + weekdays5 + weekdays6 + 
                         f(smoothing_var = t,
                           model = "IWP",
                           order = 3, k = 100, 
                           boundary.prior = list(prec = 0.01), 
                           sd.prior = list(prior = "exp", para = list(u = prior_set_conversion$u, alpha = prior_set_conversion$a))) + 
                         f(smoothing_var = index,
                           model = "IID",
                           sd.prior = list(prior = "exp", para = list(u = 0.1, alpha = 0.5))),
                       data = final_data_SK, family = "Poisson",
                       aghq_k = 10)

### Fitted Model:
IWP_part <- OS3mod_SK$instances[[1]]
samps <- aghq::sample_marginal(OS3mod_SK$mod, M = 3000)
global_samps <- samps$samps[OS3mod_SK$boundary_samp_indexes$t, , drop = F]
coefsamps <- samps$samps[OS3mod_SK$random_samp_indexes$t, , drop = F]
fixed_samps <- samps$samps[OS3mod_SK$fixed_samp_indexes$Intercept, , drop = F]

f_SK <- compute_post_fun(
  samps = coefsamps, global_samps = global_samps,
  knots = IWP_part@knots,
  refined_x = seq(min(IWP_part@observed_x), max(IWP_part@observed_x), length.out = 300),
  p = 3, degree = 0, 
  intercept_samps = fixed_samps
)
f_SK_sum <- extract_mean_interval_given_samps(f_SK)
f_SK_sum$Date <- as.Date((f_SK_sum$x)*31, origin = min(full_data_SK$Date))

f1st_SK <- compute_post_fun(
  samps = coefsamps, global_samps = global_samps,
  knots = IWP_part@knots,
  refined_x = seq(min(IWP_part@observed_x), max(IWP_part@observed_x), length.out = 300),
  p = 3, degree = 1, 
  intercept_samps = fixed_samps
)
f1st_SK_sum <- extract_mean_interval_given_samps(f1st_SK)
f1st_SK_sum$Date <- as.Date((f1st_SK_sum$x)*31, origin = min(full_data_SK$Date))
g1st_SK_sum <- extract_mean_interval_given_samps(cbind(f_SK_sum$x, (exp(f_SK[,-1]) * f1st_SK[,-1])/31))
g1st_SK_sum$Date <- as.Date((g1st_SK_sum$x)*31, origin = min(full_data_SK$Date))

### Plotting:
full_data_SK %>% ggplot(aes(x = Date, y = log(new_deaths))) + geom_point(aes(color = weekdays), alpha = 0.4) + 
  theme_classic(base_size = 15) + ylab("Log New Deaths") + xlab("") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = (c(1,7))) + 
  theme(legend.position = c(.2,.75), legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=12),
        legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm') #change legend key width
  )
ggsave(paste0(figure_path, "SK_covid_death_raw.pdf"), device = "pdf", width = 5, height = 5)

SK_POP <- 51.78 ## millions
full_data_SK %>% ggplot(aes(x = Date)) + geom_point(aes(y = new_deaths/SK_POP), alpha = 0.3, size = 0.5) + 
  geom_line(aes(y = exp(mean)/SK_POP), data = f_SK_sum, color = "blue") + 
  geom_ribbon(aes(ymax = exp(pupper)/SK_POP, ymin = exp(plower)/SK_POP), fill = "orange", alpha = 0.3, data = f_SK_sum) +
  theme_classic(base_size = 15) + ylab("New Deaths (per million)") + xlab("") + guides(color= "none") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = c(0,9))
ggsave(paste0(figure_path, "SK_covid_death_adjusted.pdf"), device = "pdf", width = 5, height = 5)

full_data_SK %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = (mean)), data = f_SK_sum, color = "blue") + 
  geom_ribbon(aes(ymax = (pupper), ymin = (plower)), fill = "orange", alpha = 0.3, data = f_SK_sum) +
  theme_classic(base_size = 15) + ylab("g(t)") + xlab("") + guides(color= "none") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = c(-2,6))
ggsave(paste0(figure_path, "SK_covid_death_log.pdf"), device = "pdf", width = 5, height = 5)

g1st_SK_sum %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean/SK_POP, color = "OS3"), linewidth = 1) + 
  geom_ribbon(aes(ymax = pupper/SK_POP, ymin = plower/SK_POP), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("Change in daily death (per million)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = c(-0.3,0.3))
ggsave(paste0(figure_path, "SK_covid_death_deriv_adjusted.pdf"), device = "pdf", width = 5, height = 5)

f1st_SK_sum %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean, color = "OS3"), linewidth = 1) + 
  geom_ribbon(aes(ymax = pupper, ymin = plower), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("g'(t)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b") +
  coord_cartesian(ylim = c(-3,5))
ggsave(paste0(figure_path, "SK_covid_death_deriv_log.pdf"), device = "pdf", width = 5, height = 5)

week_day_SK <- samps$samps[unlist(OS3mod_SK$fixed_samp_indexes[-1]), , drop = F]
week_day_Sunday <- apply(-week_day_SK, 2, sum); week_day_SK <- rbind(week_day_SK, week_day_Sunday)
rownames(week_day_SK) <- c(colnames(Xf), "weekdays7")
week_mean <- week_day_SK %>% apply(MARGIN = 1, mean); week_sd <- week_day_SK %>% apply(MARGIN = 1, sd)
week_upper <- week_mean + 1.96 * week_sd; week_lower <- week_mean - 1.96 * week_sd
week_data_SK <- data.frame(mean = week_mean, sd = week_sd, upper = week_upper, lower = week_lower, 
                           days = substr(rownames(week_day_SK), start = 9, stop = nchar(rownames(week_day_SK))))



###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
######################################################################   Additional Plots:   ##################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################
###############################################################################################################################################################################


### Overall Week Day Effects:
week_data_all <- rbind(week_data_CA, week_data_DK, week_data_SA, week_data_SK)
week_data_all$country <- rep(c("Canada", "Denmark", "South-Africa", "South-Korea"), each = 7)

week_data_all %>% ggplot(aes(x = as.numeric(days), y = mean, color = country, group = country)) + 
  scale_x_continuous(breaks = c(1:7), labels = c("Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun")) +
  theme_classic(base_size = 15) +
  xlab("Week Days") +
  ylab("Weekday effect in death rate") + 
  geom_point(position = position_dodge2(preserve = "total", width = 0.5)) + 
  geom_errorbar(aes(ymin=upper, ymax=lower), position = position_dodge2(preserve = "total"), width = 0.5) + ylim(c(-0.8,0.8)) +
  theme(legend.position="none")
ggsave(paste0(figure_path, "all_covid_death_weekdays.pdf"), device = "pdf", width = 5, height = 5)

### Hyper-parameters:
### Plotting the hyper-parameter of all countries:
library(lemon)

### For the smoothing parameter:
#### Figure for the posterior of smoothing parameter
theta_logprior <- function(theta,prior_alpha = c(.5),prior_u = c(prior_set_overdis$u1)) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2*log(x)))

### Inference for the smoothing parameter:
prec_marg_CA <- OS3mod_Canada$mod$marginals[[1]]
logpostsigmaCA <- compute_pdf_and_cdf(prec_marg_CA,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_CA <- data.frame(SD = logpostsigmaCA$transparam, 
                            density = logpostsigmaCA$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaCA$transparam))
smooth_var_CA <- rbind(smooth_var_CA, data.frame(SD = 0, density = 0, prior = priorfuncsigma(.Machine$double.eps)))

prec_marg_DK <- OS3mod_DK$mod$marginals[[1]]
logpostsigmaDK <- compute_pdf_and_cdf(prec_marg_DK,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_DK <- data.frame(SD = logpostsigmaDK$transparam, 
                            density = logpostsigmaDK$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaDK$transparam))


prec_marg_SA <- OS3mod_SA$mod$marginals[[1]]
logpostsigmaSA <- compute_pdf_and_cdf(prec_marg_SA,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_SA <- data.frame(SD = logpostsigmaSA$transparam, 
                            density = logpostsigmaSA$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaSA$transparam))


prec_marg_SK <- OS3mod_SK$mod$marginals[[1]]
logpostsigmaSK <- compute_pdf_and_cdf(prec_marg_SK,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_SK <- data.frame(SD = logpostsigmaSK$transparam, 
                            density = logpostsigmaSK$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaSK$transparam))

smooth_var_CA %>% ggplot(aes(x = SD, y = prior)) + geom_area(fill = "orange", alpha = 0.2) +
  theme_classic(base_size = 15) +
  geom_line(aes(y = density, color = "CA", linetype = "CA")) +
  geom_line(aes(x = SD, y = density, color = "DK", linetype = "DK"), data = smooth_var_DK) +
  geom_line(aes(x = SD, y = density, color = "SA", linetype = "SA"), data = smooth_var_SA) +
  geom_line(aes(x = SD, y = density, color = "SK", linetype = "SK"), data = smooth_var_SK) +
  ylab("Density") + xlim(c(0,15)) +
  labs(colour="Country") + 
  scale_linetype_manual(name = "Country", values = c("CA" = "solid", "DK" = "dashed", "SA" = "dotted", "SK" = "dotdash")) +
  theme(legend.position = c(0.8,0.5)) 
ggsave(paste0(figure_path, "smooth_para_covid.pdf"), device = "pdf", width = 5, height = 5)

#### Plot the prediction SD:

## Get the scaling number: sigma_s_h <- c * sigma_s
orig <- log(2)
c <- orig/prior_set_overdis$u1

smooth_var_CA %>% ggplot(aes(x = c*SD, y = prior/c)) + geom_area(fill = "orange", alpha = 0.2) +
  theme_classic(base_size = 15) +
  ylab("Density") + xlim(c(0,0.1)) +
  xlab("") +
  geom_line(aes(x = c*SD, y = density/c, color = "CA", linetype = "CA")) +
  geom_line(aes(x = c*SD, y = density/c, color = "DK", linetype = "DK"), data = smooth_var_DK) +
  geom_line(aes(x = c*SD, y = density/c, color = "SA", linetype = "SA"), data = smooth_var_SA) +
  geom_line(aes(x = c*SD, y = density/c, color = "SK", linetype = "SK"), data = smooth_var_SK) +
  labs(colour="Country") + 
  scale_linetype_manual(name = "Country", values = c("CA" = "solid", "DK" = "dashed", "SA" = "dotted", "SK" = "dotdash")) +
  theme(legend.position = c(0.8,0.5)) 
ggsave(paste0(figure_path, "pred_smooth_para_covid.pdf"), device = "pdf", width = 5, height = 5)

### Inference for the overdispersion parameter:
theta_logprior <- function(theta,prior_alpha = c(.5),prior_u = c(0.1)) {
  lambda <- -log(prior_alpha)/prior_u
  log(lambda/2) - lambda * exp(-theta/2) - theta/2
}
priorfunc <- function(x) exp(theta_logprior(x))
priorfuncsigma <- function(x) (2/x) * exp(theta_logprior(-2*log(x)))

prec_marg_CA <- OS3mod_Canada$mod$marginals[[2]]
logpostsigmaCA <- compute_pdf_and_cdf(prec_marg_CA,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_CA <- data.frame(SD = logpostsigmaCA$transparam, 
                            density = logpostsigmaCA$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaCA$transparam))

prec_marg_DK <- OS3mod_DK$mod$marginals[[2]]
logpostsigmaDK <- compute_pdf_and_cdf(prec_marg_DK,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_DK <- data.frame(SD = logpostsigmaDK$transparam, 
                            density = logpostsigmaDK$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaDK$transparam))

prec_marg_SA <- OS3mod_SA$mod$marginals[[2]]
logpostsigmaSA <- compute_pdf_and_cdf(prec_marg_SA,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_SA <- data.frame(SD = logpostsigmaSA$transparam, 
                            density = logpostsigmaSA$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaSA$transparam))

prec_marg_SK <- OS3mod_SK$mod$marginals[[2]]
logpostsigmaSK <- compute_pdf_and_cdf(prec_marg_SK,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_SK <- data.frame(SD = logpostsigmaSK$transparam, 
                            density = logpostsigmaSK$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaSK$transparam))

smooth_var_prior <- data.frame(SD = c(smooth_var_CA$SD, smooth_var_DK$SD, smooth_var_SA$SD, smooth_var_SK$SD), 
                               prior = c(smooth_var_CA$prior, smooth_var_DK$prior, smooth_var_SA$prior, smooth_var_SK$prior))

smooth_var_prior <- rbind(smooth_var_prior, data.frame(SD = 0, prior = priorfuncsigma(.Machine$double.eps)))

d2 <- smooth_var_prior %>% ggplot(aes(x = SD, y = prior)) + geom_area(fill = "orange", alpha = 0.2) +
  theme_classic(base_size = 15) +
  geom_line(aes(x = SD, y = density, color = "CA"), linetype="solid", data = smooth_var_CA) +
  geom_line(aes(x = SD, y = density, color = "DK"), linetype="dashed", data = smooth_var_DK) +
  geom_line(aes(x = SD, y = density, color = "SA"), linetype="dotted", data = smooth_var_SA) +
  geom_line(aes(x = SD, y = density, color = "SK"), linetype="dotdash", data = smooth_var_SK) +
  ylab("Density") +
  labs(colour="Country") +
  xlab("") + theme(legend.position = "none")
# d2 <- reposition_legend(d2, 'right')
d2
ggsave(plot = d2, paste0(figure_path, "overdis_para_covid.pdf"), device = "pdf", width = 5, height = 5)

