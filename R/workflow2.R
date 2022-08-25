
# ============================================================================================================= -

rm(list = ls())

.libPaths("T:/R4UserLibs")

list.of.packages <- c("tidyverse", "mvoutlier", "data.table", "prospectr", 
                      "ggsci", "doParallel", "foreach", "furrr", "nls.mulstart",
                      "scam")
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
library(furrr)
library(nls.multstart)
library(scam)

# set working directory
path_to_data <- "Z:/Public/Jonas/Data/ESWW006/ASD/"
path_to_funcs <- "O:/Evaluation/Projects/KP0011/8/spectra-analysis-showcase/R/"

# load functions
source(paste0(path_to_funcs, "spectra_proc.R"))

# # load all spectral data and save as .rds
# data <- load_spectra(dir = path_to_data, subdir = "renamed", format = "asd")
# saveRDS(data, paste0(path_to_data, "alldat.rds"))

# load design (as measurement meta data)
design <- read_csv("Z:/Public/Jonas/003_ESWW/Design_20211001/final_design_gen_rev.csv") %>% 
  dplyr::select(1:8, gen_name, dis_treat_new)

# ============================================================================================================= -

data <- readRDS(paste0(path_to_data, "alldat.rds"))

spc_pp <- data %>%
  preprocess_spc(average = 5, 
                 bin = 12, 
                 w = 21, p = 3,
                 trim = c(350, 375, 1350, 1475, 1781, 1990, 2400, 2500),
                 new_col = F)

# p <- plot_spectra(spc_pp, col_in = "rflt_p3w21m0_trim_bin6")

# ============================================================================================================= -

# add design and drop the reference measurements
spc_pp$plot_UID <- strsplit(spc_pp$meas_id, "_") %>% lapply("[[", 3) %>% unlist()
spc_pp <- spc_pp %>% 
  right_join(design, .,by = "plot_UID") %>% 
  dplyr::filter(!grepl("Ref", plot_UID))

# plot spectra
plot_spectra(spc_pp, facets = "meas_date", col_in = "rflt_p3w21m0_trim_bin12",
             treatment = "dis_treat_new", mark_outliers = F,
             topdf = T)

# ============================================================================================================= -

data$plot_UID <- strsplit(data$meas_id, "_") %>% lapply("[[", 3) %>% unlist()

spc_ref <- data %>% 
  right_join(design, ., by = "plot_UID") %>% 
  dplyr::filter(!grepl("Ref", plot_UID))

SVI <- spc_ref %>%
  preprocess_spc(average = 5, 
                 w = 21, p = 3, m = 0, 
                 new_col = F) %>% 
  calculate_SVI(col_in = "rflt_p3w21m0") %>% 
  scale_SVI(data = SVI, 
            plotid = "plot_UID", 
            treatid = "dis_treat_new", 
            plot = F)

saveRDS(SVI, paste0(path_to_data, "SVIdat.rds"))

plot_SVI(SVI, 
         svi = c("MCARI2", "FII"),
         col_in = "SVI_sc",
         x = "meas_date",
         x_is_date = TRUE,
         groups = "Treatment")

# ============================================================================================================= -

SVI <- readRDS(paste0(path_to_data, "SVIdat.rds"))

svi <- c("NDVI_nb_ASD", "MCARI2", "FII", "NDWI1650", "mND705", "780_740",
         "PSRI", "CARG", "PRI$", "WI_NDVI")

parameters <- get_svi_dynamics(data = SVI, 
                               svi = svi,
                               method = c("pspl"),
                               timevar = "dafm", 
                               plot_dynamics = T)

# ============================================================================================================= -


