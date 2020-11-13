
# ============================================================================================================= -

rm(list = ls())

.libPaths("T:/R3UserLibs")

list.of.packages <- c("tidyverse", "mvoutlier", "data.table", "prospectr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos='https://stat.ethz.ch/CRAN/')

library(tidyverse)
library(data.table)
library(mvoutlier)
library(prospectr)

# set working directory
basedir <- "O:/Projects/KP0011/8/"

# load functions
source(paste0(basedir, "spectra-analysis-showcase/spectra_proc.R"))

# load all spectral data
data <- load_spectra(dir = paste0(basedir, "Spectral_Data/"), format = "sed")
saveRDS(data, "alldat_list.rds")

# ============================================================================================================= -

# load spectral data
data <- readRDS("alldat_list.rds")

# pre-process spectra
spc_pp <- data %>% 
  # default
  preprocess_spc() %>%
  # for calculation of SVI
  preprocess_spc(trim = NULL, bin = NULL) %>% 
  # alternative parameters
  preprocess_spc(snv = T) %>% 
  preprocess_spc(m = 1) %>% 
  preprocess_spc(m = 1, snv = T) 

# ============================================================================================================= -

# detect multivariate outliers, plot and exclude from db
spc_pp_o <- spc_pp %>% 
  detect_outlier_spectra(grouping = "meas_date", 
                         outliers_rm= c("rflt_p3w21m0_trim_bin3", "rflt_p3w21m1_trim_bin3"),
                         outliers_rm = NULL, # if outliers are to be plotted later
                         create_plot = T)

# add measurement meta data
meta <- read_csv("spectral_files_fusarium_asign.csv")
dd <- right_join(meta, spc_pp_o[[1]], by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")

# plot spectra
plot_spc(dd, treatment = "Treatment", mark_outliers = F)

# ============================================================================================================= -

# add measurement meta data
dd <- right_join(meta, spc_pp, by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")

# calculate spectral indices: use all data (including outliers...)
SVI <- calculate_SVI(data = dd, col_in = "rflt_p3w21m0") %>% 
  scale_SVI()

# ============================================================================================================= -
