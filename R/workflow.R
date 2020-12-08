
# ============================================================================================================= -

rm(list = ls())

.libPaths("T:/R3UserLibs")

list.of.packages <- c("tidyverse", "mvoutlier", "data.table", "prospectr", "ggsci")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos = 'https://stat.ethz.ch/CRAN/')

library(tidyverse)
library(data.table)
library(mvoutlier)
library(prospectr)
library(ggsci)
library(foreach)

# set working directory
path_to_data <- "O:/Projects/KP0011/8/test/data/"
path_to_funcs <- "O:/Projects/KP0011/8/test/funcs/"

# load functions
source(paste0(path_to_funcs, "spectra_proc.R"))

# # load all spectral data and save as .rds
# data <- load_spectra(dir = path_to_data, format = "sed")
# saveRDS(data, paste0(path_to_data, "alldat.rds"))

# ============================================================================================================= -

# load spectral data
data <- readRDS(paste0(path_to_data, "alldat.rds"))

plot_spectra(data, col_in = "rflt")

# pre-process spectra
spc_pp <- data %>% 
  # default
  preprocess_spc(new_col = TRUE) %>%
  # alternative parameters
  preprocess_spc(snv = T, new_col = TRUE) %>% 
  preprocess_spc(m = 1, new_col = TRUE) %>% 
  preprocess_spc(m = 1, snv = T, new_col = TRUE) %>% 
  preprocess_spc(m = 0, snv = F, new_col = TRUE, 
                 trim = c(1350, 1475, 1781, 1990, 2400, 2500), 
                 binning = 3)

# ============================================================================================================= -

# detect multivariate outliers, plot and exclude from db
spc_pp_o <- spc_pp %>% 
  detect_outlier_spectra(grouping = NULL, 
                         outliers_rm = "rflt_p3w21m0_trim_bin3", 
                         create_plot = F)

# add measurement meta data
meta <- read_csv(paste0(path_to_data, "spectral_files_fusarium_asign.csv"))
dd <- right_join(meta, spc_pp_o[[1]], by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")

# plot spectra
plot_spectra(dd, facets = "meas_date",
             treatment = "Treatment", mark_outliers = T)

# ============================================================================================================= -

# add measurement meta data
dd <- right_join(meta, spc_pp, by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")

# ============================================================================================================= -

parameters <- get_svi_dynamics(SVI, timevar = "dafm")

# ============================================================================================================= -

# calculate spectral indices: use all data (including outliers...)
SVI <- calculate_SVI(data = dd, col_in = "rflt_p3w21m0") %>% 
  scale_SVI(plot = F)
saveRDS(SVI, paste0(basedir, "Spectral_Data/spc_pp.rds"))

plot_SVI(SVI, col_in = "SVI_sc",
         x = "meas_date", x_is_date = TRUE, 
         svi = c("SIPI", "PSRI", "NDVI"), 
         groups = "Treatment")

# ============================================================================================================= -

dd <- data %>% 
  preprocess_spc(m = 0, snv = F, new_col = T, 
               trim = c(1350, 1475, 1781, 1990, 2400, 2500), 
               binning = 6)

dd <- right_join(meta, dd, by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")

data_ref <- read.csv(paste0(path_to_data, "20200220_Fusarium_Ratings.csv"), sep = ";") %>% as_tibble() %>% 
  dplyr::select(9, 2, 4:6, 8)

ddd <- right_join(data_ref, dd)
infected <- ddd %>% dplyr::filter(Treatment != "EG")

cors <- spectra_cor(infected, trait = "FUS20200623", col_in = "rflt_p3w21m0_trim_bin6", topdf = F)

# ============================================================================================================= -

cordat <- infected[c("FUS20200623", "rflt_p3w21m0_trim_bin6")] %>% unnest(rflt_p3w21m0_trim_bin6) %>% as.data.frame()

names(cordat)[2:length(cordat)] <- paste0("rflt_", names(cordat)[2:length(cordat)])

subsets = c(200, 150, 100, 50, 40, 30, 20, 15, 10, 7, 5, 2, 1) 

# Perform recursive feature elimination
rfe <- perform_rfe(response = "FUS20200623", base_learner = "ranger", type = "regression",
                   p = 0.75, times = 30, 
                   subsets = subsets, data = cordat,
                   importance = "permutation",
                   num.trees = 200)

tidy <- tidy_rfe_output(rfe, base_learner = "ranger")
plot_perf_profile(tidy[[1]])

