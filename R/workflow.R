
# ============================================================================================================= -

rm(list = ls())

.libPaths("T:/R3UserLibs")

list.of.packages <- c("tidyverse", "mvoutlier", "data.table", "prospectr", "ggsci", "doParallel", "foreach")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos = 'https://stat.ethz.ch/CRAN/')

library(tidyverse)
library(data.table)
library(mvoutlier)
library(prospectr)
library(ggsci)
library(foreach)
library(doParallel)
library(foreach)

# set working directory
path_to_data <- "O:/Projects/KP0011/8/test/data/"
path_to_data <- "O:/Projects/KP0011/8/Spectral_Data/"
path_to_funcs <- "O:/Projects/KP0011/8/test/funcs/"

# load functions
source(paste0(path_to_funcs, "spectra_proc.R"))

# # load all spectral data and save as .rds
# data <- load_spectra(dir = path_to_data, format = "sed")
# saveRDS(data, paste0(path_to_data, "alldat.rds"))

# ============================================================================================================= -

# TRT COMPARISON

# load spectral data
data0 <- readRDS(paste0(path_to_data, "alldat.rds"))
data_2d <- data0 %>% filter(meas_date %in% c("20200615"))
data <- data0 %>% filter(meas_date == "20200622")

# pre-process spectra
spc_pp <- data %>% 
  preprocess_spc(m = 0, new_col = TRUE, 
                 trim = c(1350, 1475, 1781, 1990, 2400, 2500), 
                 binning = 3) %>% 
  preprocess_spc(m = 1, new_col = TRUE,
                 trim = c(1350, 1475, 1781, 1990, 2400, 2500))

meta <- read_csv(paste0(path_to_data, "spectral_files_fusarium_asign.csv"))
dd <- right_join(meta, spc_pp, by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")

# plot spectra
plots <- plot_spectra(dd, treatment = "Treatment", topdf = F)

plots[[2]] <- plots[[2]] +
  ylab("1st derivative of reflectance")

png("O:/Projects/KP0011/8/Figures/spc_trt_der1.png", width = 6, height = 4, units = 'in', res = 300)
plot(plots[[2]])
dev.off()

# ============================================================================================================= -

# SINGLE WVLT

spc_pp <- data_2d %>% 
  preprocess_spc(m = 0, new_col = TRUE, 
                 trim = c(1350, 1475, 1781, 1990, 2400, 2500), 
                 binning = 6) 

data_ref <- read.csv(paste0(path_to_data, "20200220_Fusarium_Ratings.csv"), sep = ";") %>% as_tibble() %>% 
  dplyr::select(9, 2, 4:6, 8)
dd <- right_join(meta, spc_pp, by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")
ddd <- right_join(data_ref, dd)
infected <- ddd %>% dplyr::filter(Treatment != "EG")

# get data
cordat <- infected[c("FUS20200616", "meas_date", "rflt_p3w21m0_trim_bin6")] %>% 
  unnest(rflt_p3w21m0_trim_bin6) %>% filter(FUS20200616 != 9) %>% 
  group_by(meas_date) %>% group_nest() %>% 
  mutate(coefs = purrr::map(data, get_cor_coefs)) %>% 
  dplyr::select(meas_date, coefs) %>% unnest(coefs)

p1 <- ggplot(cordat) + 
  geom_line(aes(x = wvlt, y = cor_coef, col = meas_date,  group = interaction(meas_date, spc_range)), size = 0.6) + 
  xlab("Wavelength (nm)") + ylab("Correlation coefficient") +
  scale_x_continuous(breaks = seq(0,2500,500), limits = c(350, 2550), expand = c(0.01, 0.01)) +
  geom_abline(slope = 0, intercept = 0) +
  ggsci::scale_color_npg() +
  theme(axis.line = element_line(colour = "black"),
        legend.position=c(0.8,0.2),
        panel.border = element_rect(fill = NA),
        panel.background = element_blank(),
        axis.title.y = element_text(size = 15, angle = 90),
        axis.title.x = element_text(size = 15, angle = 0),
        axis.text.x = element_text(size = 12.5),
        axis.text.y = element_text(size = 12.5),
        legend.text = element_text(size = 12.5),
        legend.title = element_text(size = 15),
        plot.title = element_text( size=20, face="bold"))

png("O:/Projects/KP0011/8/Figures/cor_1.png", width = 6, height = 4, units = 'in', res = 300)
plot(p1)
dev.off()

# trait data
get_cor_coefs <- function(data){
  data <- as.data.frame(data)
  traitvector <- as.matrix(as.numeric(data[,"FUS20200616"]))
  # variables with which to correlate
  variables <- as.matrix(data[, -which(names(data) == "FUS20200616")])
  # calculate correlation coefficient
  cor_coef <- -(as.vector(cor(traitvector, variables, method = "spearman")))
  wvlt <- as.numeric(gsub("_", "", stringr::str_sub(names(data)[-c(1)], -4, -1)))
  coefs <- as.data.frame(cbind(wvlt, cor_coef))
  # add spectral range for plotting
  coefs <- add_spc_range(coefs)
  return(coefs)
}

# ============================================================================================================= -

# SPECTRAL INDICES

ddd <- ddd %>% preprocess_spc()

# calculate spectral indices: use all data (including outliers...)
SVI0 <- calculate_SVI(data = ddd, col_in = "rflt_p3w21m0") %>% 
  dplyr::filter(Treatment != "EG") %>% 
  dplyr::filter(FUS20200623 != 9)

cors <- dplyr::select(SVI0, FUS20200623, SVI) %>% unnest(SVI) %>% 
  gather(Index, val, 2:length(.)) %>% group_by(Index) %>% group_nest() %>% 
  mutate(cor = purrr::map_dbl(data, ~cor(.$FUS20200623, .$val, use = "everything"))) %>% 
  arrange(desc(cor))

# ============================================================================================================= -
data <- readRDS(paste0(path_to_data, "alldat.rds"))
spc_pp <- preprocess_spc(data)
meta <- read_csv(paste0(path_to_data, "spectral_files_fusarium_asign.csv"))
data_ref <- read.csv(paste0(path_to_data, "20200220_Fusarium_Ratings.csv"), sep = ";") %>% as_tibble() %>% 
  dplyr::select(9, 2, 4:6, 8)
dd <- right_join(meta, spc_pp, by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")
ddd <- right_join(data_ref, dd)

SVI <- calculate_SVI(ddd, col_in = "rflt_p3w21m0") %>% scale_SVI(mean = "Treatment")
parameters <- get_svi_dynamics(SVI, timevar = "dafm")

ll <- parameters %>% ungroup() %>% 
  tidyr::gather(param, value, t80:dur2) %>% 
  group_by(Treatment, variable, param) %>% 
  summarise(mean = mean(value, na.rm = TRUE)) %>% 
  spread(Treatment, mean) %>% 
  mutate(diff = EG - FUS) %>% 
  arrange(desc(abs(diff)))
  
# ============================================================================================================= -

plot_SVI(SVI, col_in = "SVI_sc",
         x = "meas_date", x_is_date = TRUE, 
         svi = c("SIPI", "PSRI", "NDVI"), 
         groups = "Treatment")

# ============================================================================================================= -

# get data
moddat <- infected[c("FUS20200623", "meas_date", "rflt_p3w21m0_trim_bin6")] %>% 
  unnest(rflt_p3w21m0_trim_bin6) %>% filter(FUS20200623 != 9) %>% as.data.frame() %>% dplyr::select(-meas_date)

names(moddat)[2:length(moddat)] <- paste0("rflt_", names(moddat)[2:length(moddat)])

subsets = c(200, 150, 100, 75, 50, 40, 30, 20, 16, 14, 12, 10:1) 

# Set up cluster
y = detectCores()
cluster = makeCluster(y-2)
registerDoParallel(cluster)
clusterEvalQ(cluster, .libPaths("C:/Users/anjonas/R3UserLibs"))
clusterEvalQ(cluster, {
  library(ggplot2)
  library(MASS)
  library(caret)
  library(Cubist)
  library(tidyverse)
  library(foreach)
})
clusterExport(cluster, c("moddat", "subsets", "perform_rfe_par", 
                         "get_acc", "rmse"), envir=environment())

# Perform recursive feature elimination
rfe <- perform_rfe_par(response = "FUS20200623", base_learner = "ranger", type = "regression",
                       p = 0.75, times = 30, 
                       subsets = subsets, data = moddat,
                       importance = "permutation",
                       num.trees = 200)
saveRDS(rfe, paste0(path_to_data, "Output/rfe_spc.rds"))

tidy2 <- tidy[[2]] %>% slice(1:20)

tidy <- tidy_rfe_output(rfe, base_learner = "ranger")
p <- plot_perf_profile(tidy[[1]])
plot_feature_ranks(tidy2)

# ============================================================================================================= -

mtry <- c(1, 2, 5, 8, 12, 16, 20, 30)
min.node.size = c(1, 2, 3, 5, 8, 11)
tune_grid <- expand.grid(mtry = mtry,
                         splitrule = "variance", #default
                         min.node.size = min.node.size)
ctrl <- caret::trainControl(method = "repeatedcv",
                            number = 10,
                            rep = 5,
                            savePredictions = T,
                            verbose = TRUE)
rf_ranger <- caret::train(FUS20200623 ~ .,
                          data = moddat,
                          preProc = c("center", "scale"),
                          method = "ranger",
                          num.trees = 200,
                          tuneGrid = tune_grid,
                          importance = "permutation",
                          trControl = ctrl)

rf_ranger$results

perf <- assess_mod_perf(rf_ranger)
plot <- plot_predobs(perf)

png("O:/Projects/KP0011/8/Figures/fullspc.png", width = 6, height = 4, units = 'in', res = 300)
plot(plot)
dev.off()

#===================================================================================================


tune_grid <- expand.grid(committees = c(1, 2, 5, 10),
                         neighbors = c(0))

cubist <- caret::train(FUS20200623 ~ .,
                       data = cordat,
                       preProc = c("center", "scale"),
                       method = "cubist",
                       tuneGrid = tune_grid,
                       trControl = ctrl)

cubist$results

perf <- assess_mod_perf(cubist)
plot <- plot_predobs(perf)

#===================================================================================================

data_ref

pd <- data_ref %>% dplyr::select(-Heading_JulianDay) %>% 
  gather(date, score, starts_with("FUS")) %>% 
  mutate(date = gsub("FUS", "", date) %>% as.Date(format = "%Y%m%d")) %>% 
  dplyr::filter(score != 9)

p <- ggplot(pd) +
  geom_boxplot(aes(x = date, y = score, group = interaction(date, Treatment), col = Treatment)) +
  xlab("Date") + ylab("Visual FHB Score") +
  scale_y_continuous(breaks = c(1:8)) +
  scale_x_date(breaks = unique(pd$date)) +
  scale_color_npg() +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

png("O:/Projects/KP0011/8/Figures/FHB_dev.png", width = 6, height = 4, units = 'in', res = 300)
plot(p)
dev.off()

