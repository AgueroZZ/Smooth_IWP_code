working_path <- "/Users/ziangzhang/Documents/GitHub/random-walks/Research_Progress/covid_data/v11/"
source_path <- "/Users/ziangzhang/Documents/GitHub/random-walks/Research_Progress/covid_data/v11/source/"
figure_path <- "/Users/ziangzhang/Documents/GitHub/random-walks/Research_Progress/covid_data/v11/figures/"
cpp_path <- "/Users/ziangzhang/Documents/GitHub/random-walks/Research_Progress/covid_data/v11/cpp/"

source(paste0(source_path,"08_loads.R"))
source(paste0(source_path,"18_functions_defined.R"))
source(paste0(source_path,"02_functions_model_selection.R"))
source(paste0(source_path,"bayes_regression_COVID.R"))
compile(paste0(cpp_path,"00_Poisson_Smoothing_PC_overdisp_covid.cpp"))
dyn.load(dynlib(paste0(cpp_path,"00_Poisson_Smoothing_PC_overdisp_covid")))
full_data_covid <- read.csv(file = paste0(working_path, "owid-covid-data.csv"), header = T)
full_data_covid <- full_data_covid %>% filter(date >= "2020-03-01")


###############################
###############################
###1.Canada Analysis:##########
###############################
###############################
full_data_canada <- full_data_covid %>% filter(location == "Canada") %>% dplyr::select(date, new_deaths, new_deaths_smoothed)
full_data_canada$new_deaths[is.na(full_data_canada$new_deaths)] <- 0
full_data_canada$t <- as.Date(full_data_canada$date) %>% as.numeric()
full_data_canada$t <- full_data_canada$t - full_data_canada$t[1]
full_data_canada$t <- full_data_canada$t/31
full_data_canada$weekdays <- weekdays(as.Date(full_data_canada$date))
full_data_canada$weekdays <- factor(full_data_canada$weekdays,
                                    levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"),
                                    ordered = F)
Xf <- model.matrix(~ weekdays, data = full_data_canada, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
p = 3
prior_set_conversion <- prior_order_conversion_predictive(d = (7/31), prior = list(u = log(2), alpha = 0.5), p = 3)
prior_set_overdis <- list(u1 = prior_set_conversion$u,
                          alpha1 = prior_set_conversion$alpha,
                          u2 = 0.1,
                          alpha2 = 0.5,
                          betaprec = 0.01,
                          Xfprec = 0.01)
knots <- seq(0,max(full_data_canada$t), length.out = 100)
OS3mod_Canada <- Imple_BayesRegression_COVID(y = round(full_data_canada$new_deaths), x = full_data_canada$t,
                                           knots = knots, p = p, prior = prior_set_overdis, aghq_k = 10, Xf = Xf)
OS3g_result_Over <- extract_mean_interval_refined_overd(OS3mod_Canada, x = knots, 
                                                        refined_x = full_data_canada$t,
                                                        p = p, type = "function", Xf = Xf)
full_data_canada$OS3smoothed_overd <- OS3g_result_Over$mean
full_data_canada$OS3smoothed_upper_overd <- OS3g_result_Over$upper
full_data_canada$OS3smoothed_lower_overd <- OS3g_result_Over$lower
full_data_canada$Date <- as.Date(full_data_canada$date)

full_data_canada %>% ggplot(aes(x = Date, y = log(new_deaths))) + geom_point(aes(color = weekdays), alpha = 0.4) + 
  theme_classic(base_size = 15) + ylab("Log New Deaths") + xlab("") + 
  theme(legend.position = "none"
        ) +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = (c(1,7)))
ggsave(paste0(figure_path, "canada_covid_death_raw.pdf"), device = "pdf", width = 5, height = 5)


### Also present the population adjusted figure:
CANADA_POP <- 38.01 ## millions
full_data_canada %>% ggplot(aes(x = Date)) +
  geom_line(aes(y = exp(OS3smoothed_overd)/CANADA_POP, color = "OS3")) + 
  geom_ribbon(aes(ymax = exp(OS3smoothed_upper_overd)/CANADA_POP, ymin = exp(OS3smoothed_lower_overd)/CANADA_POP), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + ylab("New Deaths (per million)") + xlab("") + 
  scale_color_manual(values = c("blue")) + guides(color= "none") +
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = c(0,9))
ggsave(paste0(figure_path, "canada_covid_death_adjusted.pdf"), device = "pdf", width = 5, height = 5)


full_data_canada %>% ggplot(aes(x = Date)) +
  geom_line(aes(y = (OS3smoothed_overd), color = "OS3")) + 
  geom_ribbon(aes(ymax = (OS3smoothed_upper_overd), ymin = (OS3smoothed_lower_overd)), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + ylab("g(t)") + xlab("") + 
  scale_color_manual(values = c("blue")) + guides(color= "none") +
  scale_x_date(date_labels = "%y %b")
ggsave(paste0(figure_path, "canada_covid_log_death.pdf"), device = "pdf", width = 5, height = 5)


OS3f_samps_over <- extract_samples_Ospline_refined_overd(OS3mod_Canada, x = knots, 
                                                         refined_x = full_data_canada$t,
                                                         p = 3, type = "function")
OS3f1st_samps_over <- extract_deriv_samples_OSpline_overd(OS3mod_Canada, x = knots, 
                                                          refined_x = full_data_canada$t,
                                                          p = 3, degree = 1)
OS3g1st_samps_over <- cbind(OS3f_samps_over[,1]*31, (exp(OS3f_samps_over[,-1]) * OS3f1st_samps_over[,-1])/31)


OS3f1st_result_over <- data.frame(Date = full_data_canada$Date, mean = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, mean), 
                                 upper = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.975),
                                 lower = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.025))


OS3g1st_result_over <- data.frame(x = OS3g1st_samps_over[,1], mean = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, mean), 
                                  upper = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.975),
                                  lower = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.025))

full_data_canada$derivOS3smoothed_over <- OS3g1st_result_over$mean
full_data_canada$derivOS3smoothed_upper_over <- OS3g1st_result_over$upper
full_data_canada$derivOS3smoothed_lower_over <- OS3g1st_result_over$lower


full_data_canada %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = derivOS3smoothed_over/CANADA_POP, color = "OS3"), size = 1) + 
  geom_ribbon(aes(ymax = derivOS3smoothed_upper_over/CANADA_POP, ymin = derivOS3smoothed_lower_over/CANADA_POP), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("Change in daily death (per million)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = c(-0.3,0.3))

ggsave(paste0(figure_path, "canada_covid_death_deriv_adjusted.pdf"), device = "pdf", width = 5, height = 5)


OS3f1st_result_over %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean, color = "OS3"), size = 1) + 
  geom_ribbon(aes(ymax = upper, ymin = lower), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("g'(t)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")

ggsave(paste0(figure_path, "canada_covid_death_deriv_log.pdf"), device = "pdf", width = 5, height = 5)

poster <- sample_marginal(OS3mod_Canada, M = 3000)
pos_samps <- poster$samps
pos_samps_fixed <- pos_samps[(length(knots) + p): (nrow(pos_samps) - nrow(full_data_canada)), ]
pos_samps_fixed_Sunday <- apply(-pos_samps_fixed, 2, sum)
pos_samps_fixed <- rbind(pos_samps_fixed, pos_samps_fixed_Sunday)
rownames(pos_samps_fixed) <- c(colnames(Xf), "weekdays7")

week_mean <- pos_samps_fixed %>% apply(MARGIN = 1, mean)
week_sd <- pos_samps_fixed %>% apply(MARGIN = 1, sd)
week_upper <- week_mean + 1.96 * week_sd
week_lower <- week_mean - 1.96 * week_sd

week_data_CA <- data.frame(mean = week_mean, sd = week_sd, upper = week_upper, lower = week_lower, 
                        days = substr(rownames(pos_samps_fixed), start = 9, stop = nchar(rownames(pos_samps_fixed))))


###############################
###############################
###2.Denmark Analysis:#########
###############################
###############################

full_data_DK <- full_data_covid %>% filter(location == "Denmark") %>% dplyr::select(date, new_deaths, new_deaths_smoothed)
full_data_DK$new_deaths[is.na(full_data_DK$new_deaths)] <- 0
full_data_DK$t <- as.Date(full_data_DK$date) %>% as.numeric()
full_data_DK$t <- full_data_DK$t - full_data_DK$t[1]
full_data_DK$t <- full_data_DK$t/31
full_data_DK$weekdays <- weekdays(as.Date(full_data_DK$date))
full_data_DK$weekdays <- factor(full_data_DK$weekdays,
                                levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"),
                                ordered = F)
Xf <- model.matrix(~ weekdays, data = full_data_DK, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
p = 3
prior_set_conversion <- prior_order_conversion_predictive(d = (7/31), prior = list(u = log(2), alpha = 0.5), p = 3)

OS3mod_Denmark <- Imple_BayesRegression_COVID(y = round(full_data_DK$new_deaths), x = full_data_DK$t,
                                              knots = knots, p = p, prior = prior_set_overdis, aghq_k = 10, Xf = Xf)
OS3g_result_Over <- extract_mean_interval_refined_overd(OS3mod_Denmark, x = knots, 
                                                        refined_x = full_data_DK$t,
                                                        p = p, type = "function", Xf = Xf)
full_data_DK$OS3smoothed_overd <- OS3g_result_Over$mean
full_data_DK$OS3smoothed_upper_overd <- OS3g_result_Over$upper
full_data_DK$OS3smoothed_lower_overd <- OS3g_result_Over$lower
full_data_DK$Date <- as.Date(full_data_DK$date)


full_data_DK %>% ggplot(aes(x = Date, y = log(new_deaths))) + geom_point(aes(color = weekdays), alpha = 0.4) + 
  theme_classic(base_size = 15) + ylab("Log New Deaths") + xlab("") + 
  theme(legend.position = "none") + 
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = (c(1,7)))

ggsave(paste0(figure_path, "DK_covid_death_raw.pdf"), device = "pdf", width = 5, height = 5)



DK_SIZE <- 5.831 # millions
full_data_DK %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = exp(OS3smoothed_overd)/DK_SIZE, color = "OS3")) + 
  geom_ribbon(aes(ymax = exp(OS3smoothed_upper_overd)/DK_SIZE, ymin = exp(OS3smoothed_lower_overd)/DK_SIZE), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + ylab("New Deaths (per million)") + xlab("") + 
  scale_color_manual(values = c("blue"))  +
  scale_x_date(date_labels = "%y %b") + guides(color= "none") + coord_cartesian(ylim = c(0,9))
ggsave(paste0(figure_path, "DK_covid_death_adjusted.pdf"), device = "pdf", width = 5, height = 5)


full_data_DK %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = (OS3smoothed_overd), color = "OS3")) + 
  geom_ribbon(aes(ymax = (OS3smoothed_upper_overd), ymin = (OS3smoothed_lower_overd)), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + ylab("g(t)") + xlab("") + 
  scale_color_manual(values = c("blue"))  +
  scale_x_date(date_labels = "%y %b") + guides(color= "none")
ggsave(paste0(figure_path, "DK_covid_log_death.pdf"), device = "pdf", width = 5, height = 5)


OS3g_result_Over <- extract_mean_interval_refined_overd(OS3mod_Denmark, x = knots, 
                                                        refined_x = full_data_DK$t,
                                                        p = p, type = "function", Xf = Xf)

OS3f_samps_over <- extract_samples_Ospline_refined_overd(OS3mod_Denmark, x = knots, 
                                                         refined_x = full_data_DK$t,
                                                         p = 3, type = "function")
OS3f1st_samps_over <- extract_deriv_samples_OSpline_overd(OS3mod_Denmark, x = knots, 
                                                          refined_x = full_data_DK$t,
                                                          p = 3, degree = 1)
OS3g1st_samps_over <- cbind(OS3f_samps_over[,1]*31, (exp(OS3f_samps_over[,-1]) * OS3f1st_samps_over[,-1])/31)
OS3g1st_result_over <- data.frame(x = OS3g1st_samps_over[,1], mean = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, mean), 
                                  upper = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.975),
                                  lower = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.025))
full_data_DK$derivOS3smoothed_over <- OS3g1st_result_over$mean
full_data_DK$derivOS3smoothed_upper_over <- OS3g1st_result_over$upper
full_data_DK$derivOS3smoothed_lower_over <- OS3g1st_result_over$lower

full_data_DK %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = derivOS3smoothed_over/DK_SIZE, color = "OS3"), size = 1) + 
  geom_ribbon(aes(ymax = derivOS3smoothed_upper_over/DK_SIZE, ymin = derivOS3smoothed_lower_over/DK_SIZE), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("Change in daily death (per million)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  scale_x_date(date_labels = "%y %b") + guides(color= "none")  + coord_cartesian(ylim = c(-0.3,0.3))
ggsave(paste0(figure_path,"DK_covid_death_deriv_adjusted.pdf"), device = "pdf", width = 5, height = 5)


OS3f1st_result_over <- data.frame(Date = full_data_DK$Date, mean = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, mean), 
                                  upper = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.975),
                                  lower = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.025))

OS3f1st_result_over %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean, color = "OS3"), size = 1) + 
  geom_ribbon(aes(ymax = upper, ymin = lower), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("g'(t)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")

ggsave(paste0(figure_path, "DK_covid_death_deriv_log.pdf"), device = "pdf", width = 5, height = 5)

poster <- sample_marginal(OS3mod_Denmark, M = 3000)
pos_samps <- poster$samps
pos_samps_fixed <- pos_samps[(length(knots) + p): (nrow(pos_samps) - nrow(full_data_DK)), ]
pos_samps_fixed_Sunday <- apply(-pos_samps_fixed, 2, sum)
pos_samps_fixed <- rbind(pos_samps_fixed, pos_samps_fixed_Sunday)
rownames(pos_samps_fixed) <- c(colnames(Xf), "weekdays7")

week_mean <- pos_samps_fixed %>% apply(MARGIN = 1, mean)
week_sd <- pos_samps_fixed %>% apply(MARGIN = 1, sd)
week_upper <- week_mean + 1.96 * week_sd
week_lower <- week_mean - 1.96 * week_sd

week_data_DK <- data.frame(mean = week_mean, sd = week_sd, upper = week_upper, lower = week_lower, 
                        days = substr(rownames(pos_samps_fixed), start = 9, stop = nchar(rownames(pos_samps_fixed))))


###############################
###############################
###3.South Africa Analysis:####
###############################
###############################

full_data_SA <- full_data_covid %>% filter(location == "South Africa") %>% dplyr::select(date, new_deaths, new_deaths_smoothed)
full_data_SA$new_deaths[is.na(full_data_SA$new_deaths)] <- 0
full_data_SA$t <- as.Date(full_data_SA$date) %>% as.numeric()
full_data_SA$t <- full_data_SA$t - full_data_SA$t[1]
full_data_SA$t <- full_data_SA$t/31
full_data_SA$weekdays <- weekdays(as.Date(full_data_SA$date))
full_data_SA$weekdays <- factor(full_data_SA$weekdays,
                                levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"),
                                ordered = F)
Xf <- model.matrix(~ weekdays, data = full_data_SA, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
p = 3
knots <- seq(0,max(full_data_SA$t), length.out = 100)
OS3mod_SA <- Imple_BayesRegression_COVID(y = round(full_data_SA$new_deaths), x = full_data_SA$t,
                                           knots = knots, p = p, prior = prior_set_overdis, aghq_k = 10, Xf = Xf)
OS3g_result_Over <- extract_mean_interval_refined_overd(OS3mod_SA, x = knots, 
                                                        refined_x = full_data_SA$t,
                                                        p = p, type = "function", Xf = Xf)

full_data_SA$OS3smoothed_overd <- OS3g_result_Over$mean
full_data_SA$OS3smoothed_upper_overd <- OS3g_result_Over$upper
full_data_SA$OS3smoothed_lower_overd <- OS3g_result_Over$lower
full_data_SA$Date <- as.Date(full_data_SA$date)


full_data_SA %>% ggplot(aes(x = Date, y = log(new_deaths))) + geom_point(aes(color = weekdays), alpha = 0.4) + 
  theme_classic(base_size = 15) + ylab("Log New Deaths") + xlab("") + 
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = (c(1,7))) +
  theme(legend.position = "none")
ggsave(paste0(figure_path, "SA_covid_death_raw.pdf"), device = "pdf", width = 5, height = 5)


SA_SIZE <- 59.31
full_data_SA %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = exp(OS3smoothed_overd)/SA_SIZE, color = "OS3")) + 
  geom_ribbon(aes(ymax = exp(OS3smoothed_upper_overd)/SA_SIZE, ymin = exp(OS3smoothed_lower_overd)/SA_SIZE), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + ylab("New Deaths (per million)") + xlab("") + 
  scale_color_manual(values = c("blue")) + guides(color= "none") +
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = c(0,9))
ggsave(paste0(figure_path, "SA_covid_death_adjusted.pdf"), device = "pdf", width = 5, height = 5)

full_data_SA %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = (OS3smoothed_overd), color = "OS3")) + 
  geom_ribbon(aes(ymax = (OS3smoothed_upper_overd), ymin = (OS3smoothed_lower_overd)), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + ylab("g(t)") + xlab("") + 
  scale_color_manual(values = c("blue")) + guides(color= "none") +
  scale_x_date(date_labels = "%y %b")

ggsave(paste0(figure_path, "SA_covid_log_death.pdf"), device = "pdf", width = 5, height = 5)


OS3f_samps_over <- extract_samples_Ospline_refined_overd(OS3mod_SA, x = knots, 
                                                         refined_x = full_data_SA$t,
                                                         p = 3, type = "function")
OS3f1st_samps_over <- extract_deriv_samples_OSpline_overd(OS3mod_SA, x = knots, 
                                                          refined_x = full_data_SA$t,
                                                          p = 3, degree = 1)
OS3g1st_samps_over <- cbind(OS3f_samps_over[,1]*31, (exp(OS3f_samps_over[,-1]) * OS3f1st_samps_over[,-1])/31)


OS3g1st_result_over <- data.frame(x = OS3g1st_samps_over[,1], mean = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, mean), 
                                  upper = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.975),
                                  lower = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.025))
full_data_SA$derivOS3smoothed_over <- OS3g1st_result_over$mean
full_data_SA$derivOS3smoothed_upper_over <- OS3g1st_result_over$upper
full_data_SA$derivOS3smoothed_lower_over <- OS3g1st_result_over$lower

full_data_SA %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = derivOS3smoothed_over/SA_SIZE, color = "OS3"), size = 1) + 
  geom_ribbon(aes(ymax = derivOS3smoothed_upper_over/SA_SIZE, ymin = derivOS3smoothed_lower_over/SA_SIZE), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("Change in daily death (per million)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = c(-0.3,0.3))
ggsave(paste0(figure_path, "SA_covid_death_deriv_adjusted.pdf"), device = "pdf", width = 5, height = 5)

OS3f1st_result_over <- data.frame(Date = full_data_SA$Date, mean = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, mean), 
                                  upper = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.975),
                                  lower = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.025))

OS3f1st_result_over %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean, color = "OS3"), size = 1) + 
  geom_ribbon(aes(ymax = upper, ymin = lower), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("g'(t)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")

ggsave(paste0(figure_path, "SA_covid_death_deriv_log.pdf"), device = "pdf", width = 5, height = 5)


poster <- sample_marginal(OS3mod_SA, M = 3000)
pos_samps <- poster$samps
pos_samps_fixed <- pos_samps[(length(knots) + p): (nrow(pos_samps) - nrow(full_data_SA)), ]
pos_samps_fixed_Sunday <- apply(-pos_samps_fixed, 2, sum)
pos_samps_fixed <- rbind(pos_samps_fixed, pos_samps_fixed_Sunday)
rownames(pos_samps_fixed) <- c(colnames(Xf), "weekdays7")

week_mean <- pos_samps_fixed %>% apply(MARGIN = 1, mean)
week_sd <- pos_samps_fixed %>% apply(MARGIN = 1, sd)
week_upper <- week_mean + 1.96 * week_sd
week_lower <- week_mean - 1.96 * week_sd

week_data_SA <- data.frame(mean = week_mean, sd = week_sd, upper = week_upper, lower = week_lower, 
                        days = substr(rownames(pos_samps_fixed), start = 9, stop = nchar(rownames(pos_samps_fixed))))

###############################
###############################
###4.South Korea Analysis:####
###############################
###############################
full_data_SK <- full_data_covid %>% filter(location == "South Korea") %>% dplyr::select(date, new_deaths, new_deaths_smoothed)
full_data_SK$new_deaths[is.na(full_data_SK$new_deaths)] <- 0
full_data_SK$t <- as.Date(full_data_SK$date) %>% as.numeric()
full_data_SK$t <- full_data_SK$t - full_data_SK$t[1]
full_data_SK$t <- full_data_SK$t/31
full_data_SK$weekdays <- weekdays(as.Date(full_data_SK$date))
full_data_SK$weekdays <- factor(full_data_SK$weekdays,
                                levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"),
                                ordered = F)
Xf <- model.matrix(~ weekdays, data = full_data_SK, contrasts.arg = list(weekdays = "contr.sum"))[,-1]
knots <- seq(0,max(full_data_SK$t), length.out = 100)
OS3mod_SK <- Imple_BayesRegression_COVID(y = round(full_data_SK$new_deaths), x = full_data_SK$t,
                                           knots = knots, p = p, prior = prior_set_overdis, aghq_k = 10, Xf = Xf)
OS3g_result_Over <- extract_mean_interval_refined_overd(OS3mod_SK, x = knots, 
                                                        refined_x = full_data_SK$t,
                                                        p = p, type = "function", Xf = Xf)
full_data_SK$OS3smoothed_overd <- OS3g_result_Over$mean
full_data_SK$OS3smoothed_upper_overd <- OS3g_result_Over$upper
full_data_SK$OS3smoothed_lower_overd <- OS3g_result_Over$lower
full_data_SK$Date <- as.Date(full_data_SK$date)

full_data_SK %>% ggplot(aes(x = Date, y = log(new_deaths))) + geom_point(aes(color = weekdays), alpha = 0.4) + 
  theme_classic(base_size = 15) + ylab("Log New Deaths") + xlab("") + 
  scale_x_date(date_labels = "%y %b") + coord_cartesian(ylim = (c(0.5,7))) + 
  theme(legend.position = c(.2,.75), legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=12),
        legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm') #change legend key width
  )
ggsave(paste0(figure_path, "SK_covid_death_raw.pdf"), device = "pdf", width = 5, height = 5)


SK_SIZE <- 51.78 ## millions
full_data_SK %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = exp(OS3smoothed_overd)/SK_SIZE, color = "OS3")) + 
  geom_ribbon(aes(ymax = exp(OS3smoothed_upper_overd)/SK_SIZE, ymin = exp(OS3smoothed_lower_overd)/SK_SIZE), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + ylab("New Deaths (per million)") + xlab("") + 
  scale_color_manual(values = c("blue")) + guides(color= "none") +
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = c(0,9))
ggsave(paste0(figure_path, "SK_covid_death_adjusted.pdf"), device = "pdf", width = 5, height = 5)

full_data_SK %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = (OS3smoothed_overd), color = "OS3")) + 
  geom_ribbon(aes(ymax = (OS3smoothed_upper_overd), ymin = (OS3smoothed_lower_overd)), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + ylab("g(t)") + xlab("") + 
  scale_color_manual(values = c("blue")) + guides(color= "none") +
  scale_x_date(date_labels = "%y %b")

ggsave(paste0(figure_path, "SK_covid_log_death.pdf"), device = "pdf", width = 5, height = 5)


OS3f_samps_over <- extract_samples_Ospline_refined_overd(OS3mod_SK, x = knots, 
                                                         refined_x = full_data_SK$t,
                                                         p = 3, type = "function")
OS3f1st_samps_over <- extract_deriv_samples_OSpline_overd(OS3mod_SK, x = knots, 
                                                          refined_x = full_data_SK$t,
                                                          p = 3, degree = 1)
OS3g1st_samps_over <- cbind(OS3f_samps_over[,1]*31, (exp(OS3f_samps_over[,-1]) * OS3f1st_samps_over[,-1])/31)


OS3f1st_result_over <- data.frame(Date = full_data_SK$Date, mean = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, mean), 
                                  upper = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.975),
                                  lower = OS3f1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.025))

OS3g1st_result_over <- data.frame(x = OS3g1st_samps_over[,1], mean = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, mean), 
                                  upper = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.975),
                                  lower = OS3g1st_samps_over[,-1] %>% apply(MARGIN = 1, quantile, p = 0.025))

full_data_SK$derivOS3smoothed_over <- OS3g1st_result_over$mean
full_data_SK$derivOS3smoothed_upper_over <- OS3g1st_result_over$upper
full_data_SK$derivOS3smoothed_lower_over <- OS3g1st_result_over$lower



full_data_SK %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = derivOS3smoothed_over/SK_SIZE, color = "OS3"), size = 1) + 
  geom_ribbon(aes(ymax = derivOS3smoothed_upper_over/SK_SIZE, ymin = derivOS3smoothed_lower_over/SK_SIZE), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("Change in daily death (per million)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")  + coord_cartesian(ylim = c(-0.3,0.3))

ggsave(paste0(figure_path, "SK_covid_death_deriv_adjusted.pdf"), device = "pdf", width = 5, height = 5)

OS3f1st_result_over %>% ggplot(aes(x = Date)) + 
  geom_line(aes(y = mean, color = "OS3"), size = 1) + 
  geom_ribbon(aes(ymax = upper, ymin = lower), fill = "orange", alpha = 0.3) +
  theme_classic(base_size = 15) + 
  ylab("g'(t)") + xlab("") +
  scale_color_manual(values = c("blue")) +
  guides(color= "none") +
  scale_x_date(date_labels = "%y %b")

ggsave(paste0(figure_path, "SK_covid_death_deriv_log.pdf"), device = "pdf", width = 5, height = 5)


poster <- sample_marginal(OS3mod_SK, M = 3000)
pos_samps <- poster$samps
pos_samps_fixed <- pos_samps[(length(knots) + p): (nrow(pos_samps) - nrow(full_data_SK)), ]
pos_samps_fixed_Sunday <- apply(-pos_samps_fixed, 2, sum)
pos_samps_fixed <- rbind(pos_samps_fixed, pos_samps_fixed_Sunday)
rownames(pos_samps_fixed) <- c(colnames(Xf), "weekdays7")

week_mean <- pos_samps_fixed %>% apply(MARGIN = 1, mean)
week_sd <- pos_samps_fixed %>% apply(MARGIN = 1, sd)
week_upper <- week_mean + 1.96 * week_sd
week_lower <- week_mean - 1.96 * week_sd

week_data_SK <- data.frame(mean = week_mean, sd = week_sd, upper = week_upper, lower = week_lower, 
                        days = substr(rownames(pos_samps_fixed), start = 9, stop = nchar(rownames(pos_samps_fixed))))




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
prec_marg_CA <- OS3mod_Canada$marginals[[1]]
logpostsigmaCA <- compute_pdf_and_cdf(prec_marg_CA,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_CA <- data.frame(SD = logpostsigmaCA$transparam, 
                            density = logpostsigmaCA$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaCA$transparam))
smooth_var_CA <- rbind(smooth_var_CA, data.frame(SD = 0, density = 0, prior = priorfuncsigma(.Machine$double.eps)))

prec_marg_DK <- OS3mod_Denmark$marginals[[1]]
logpostsigmaDK <- compute_pdf_and_cdf(prec_marg_DK,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_DK <- data.frame(SD = logpostsigmaDK$transparam, 
                            density = logpostsigmaDK$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaDK$transparam))


prec_marg_SA <- OS3mod_SA$marginals[[1]]
logpostsigmaSA <- compute_pdf_and_cdf(prec_marg_SA,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_SA <- data.frame(SD = logpostsigmaSA$transparam, 
                            density = logpostsigmaSA$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaSA$transparam))


prec_marg_SK <- OS3mod_SK$marginals[[1]]
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


prec_marg_CA <- OS3mod_Canada$marginals[[2]]
logpostsigmaCA <- compute_pdf_and_cdf(prec_marg_CA,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_CA <- data.frame(SD = logpostsigmaCA$transparam, 
                            density = logpostsigmaCA$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaCA$transparam))

prec_marg_DK <- OS3mod_Denmark$marginals[[2]]
logpostsigmaDK <- compute_pdf_and_cdf(prec_marg_DK,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_DK <- data.frame(SD = logpostsigmaDK$transparam, 
                            density = logpostsigmaDK$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaDK$transparam))


prec_marg_SA <- OS3mod_SA$marginals[[2]]
logpostsigmaSA <- compute_pdf_and_cdf(prec_marg_SA,list(totheta = function(x) -2*log(x),fromtheta = function(x) exp(-x/2)),interpolation = 'spline')
smooth_var_SA <- data.frame(SD = logpostsigmaSA$transparam, 
                            density = logpostsigmaSA$pdf_transparam,
                            prior = priorfuncsigma(logpostsigmaSA$transparam))


prec_marg_SK <- OS3mod_SK$marginals[[2]]
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



