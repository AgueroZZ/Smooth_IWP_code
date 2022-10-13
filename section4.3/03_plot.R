WORK_PATH <- getwd()
figure_path <- paste0(WORK_PATH, "/figure/")
result_path <- paste0(WORK_PATH, "/result/")
cpp_path <- paste0(WORK_PATH, "/cpp/")
source_path <- paste0(WORK_PATH, "/source/")
source(paste0(source_path,"00_load.R"))
source(paste0(source_path,"01_Ospline.R"))
source(paste0(source_path,"02_prior.R"))

### Plotting the result in the paper:
load(paste0(result_path, "GMM_result_high.rda"))

c <- median(result_high$rMSE_g[result_high$method == "OS"])
result_high$rMSE_g_rt <- (result_high$rMSE_g)/c
result_high %>% ggplot(aes(x = method, color = method, y = rMSE_g_rt)) + geom_boxplot() +
  theme_classic(base_size = 15) + theme(legend.position="none") + ylab("scaled-rMSE") + ylim(0.5,4)
ggsave(paste0(figure_path, "high_agg_g_rt.pdf"), width = 5, height = 5, device = "pdf")



c <- median(result_high$rMSE_g1st[result_high$method == "OS"])
result_high$rMSE_g1st_rt <- (result_high$rMSE_g1st)/c
result_high %>% ggplot(aes(x = method, color = method, y = rMSE_g1st_rt)) + geom_boxplot() +
  theme_classic(base_size = 15) + theme(legend.position="none") + ylab("scaled-rMSE") + ylim(0.5,4)
ggsave(paste0(figure_path, "high_agg_g1st_rt.pdf"), width = 5, height = 5, device = "pdf")



c <- median(result_high$rMSE_g2nd[result_high$method == "OS"])
result_high$rMSE_g2nd_rt <- (result_high$rMSE_g2nd)/c
result_high %>% ggplot(aes(x = method, color = method, y = rMSE_g2nd_rt)) + geom_boxplot() +
  theme_classic(base_size = 15) + theme(legend.position="none") + ylab("scaled-rMSE") + ylim(0.5,4)
ggsave(paste0(figure_path, "high_agg_g2nd_rt.pdf"), width = 5, height = 5, device = "pdf")






#### Appendix for other choices of SNR:

load(paste0(result_path, "GMM_result_med.rda"))

c <- median(result_med$rMSE_g[result_med$method == "OS"])
result_med$rMSE_g_rt <- (result_med$rMSE_g)/c
result_med %>% ggplot(aes(x = method, color = method, y = rMSE_g_rt)) + geom_boxplot() +
  theme_classic(base_size = 15) + theme(legend.position="none") + ylab("scaled-rMSE") + ylim(0.5,4)
ggsave(paste0(figure_path, "med_agg_g_rt.pdf"), width = 5, height = 5, device = "pdf")



c <- median(result_med$rMSE_g1st[result_med$method == "OS"])
result_med$rMSE_g1st_rt <- (result_med$rMSE_g1st)/c
result_med %>% ggplot(aes(x = method, color = method, y = rMSE_g1st_rt)) + geom_boxplot() +
  theme_classic(base_size = 15) + theme(legend.position="none") + ylab("scaled-rMSE") + ylim(0.5,4)
ggsave(paste0(figure_path, "med_agg_g1st_rt.pdf"), width = 5, height = 5, device = "pdf")



c <- median(result_med$rMSE_g2nd[result_med$method == "OS"])
result_med$rMSE_g2nd_rt <- (result_med$rMSE_g2nd)/c
result_med %>% ggplot(aes(x = method, color = method, y = rMSE_g2nd_rt)) + geom_boxplot() +
  theme_classic(base_size = 15) + theme(legend.position="none") + ylab("scaled-rMSE") + ylim(0.5,4)
ggsave(paste0(figure_path, "med_agg_g2nd_rt.pdf"), width = 5, height = 5, device = "pdf")



load(paste0(result_path, "GMM_result_low.rda"))

c <- median(result_low$rMSE_g[result_med$method == "OS"])
result_low$rMSE_g_rt <- (result_low$rMSE_g)/c
result_low %>% ggplot(aes(x = method, color = method, y = rMSE_g_rt)) + geom_boxplot() +
  theme_classic(base_size = 15) + theme(legend.position="none") + ylab("scaled-rMSE") + ylim(0.5,4)
ggsave(paste0(figure_path, "low_agg_g_rt.pdf"), width = 5, height = 5, device = "pdf")



c <- median(result_low$rMSE_g1st[result_med$method == "OS"])
result_low$rMSE_g1st_rt <- (result_low$rMSE_g1st)/c
result_low %>% ggplot(aes(x = method, color = method, y = rMSE_g1st_rt)) + geom_boxplot() +
  theme_classic(base_size = 15) + theme(legend.position="none") + ylab("scaled-rMSE") + ylim(0.5,4)
ggsave(paste0(figure_path, "low_agg_g1st_rt.pdf"), width = 5, height = 5, device = "pdf")


c <- median(result_low$rMSE_g2nd[result_med$method == "OS"])
result_low$rMSE_g2nd_rt <- (result_low$rMSE_g2nd)/c
result_low %>% ggplot(aes(x = method, color = method, y = rMSE_g2nd_rt)) + geom_boxplot() +
  theme_classic(base_size = 15) + theme(legend.position="none") + ylab("scaled-rMSE") + ylim(0.5,4)
ggsave(paste0(figure_path, "low_agg_g2nd_rt.pdf"), width = 5, height = 5, device = "pdf")







