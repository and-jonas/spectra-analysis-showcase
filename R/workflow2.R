
# ============================================================================================================= -
# Prepare workspace ----
# ============================================================================================================= -

rm(list = ls())

.libPaths("T:/R4UserLibs")

list.of.packages <- c("tidyverse", "mvoutlier", "data.table", "prospectr", 
                      "ggsci", "doParallel", "foreach", "furrr", "nls.multstart",
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

# ============================================================================================================= -
# Pre-process and plot spectra ----
# ============================================================================================================= -

data <- readRDS(paste0(path_to_data, "alldat.rds"))

# remove data from plots with blue background 
p <- read_csv("Z:/Public/Jonas/Data/ESWW006/ImagesNadir/Meta/images_exclude.csv")
pdates <- strsplit(p$Exclude, "_") %>% lapply("[[", 1) %>% unlist() %>% unique()
pplots <- strsplit(p$Exclude, "_") %>% lapply("[[", 3) %>% unlist() %>% unique()
data <- data %>% 
  dplyr::filter(!(meas_date %in% pdates & Plot_ID %in% pplots))

spc_pp <- data %>%
  preprocess_spc(data = ., 
                 col_in = "rflt", 
                 average = 5, 
                 bin = 6, 
                 w = 21, p = 3,
                 trim = c(350, 375, 1350, 1475, 1781, 1990, 2400, 2500),
                 new_col = T) %>% 
  preprocess_spc(data = ., 
                 col_in = "rflt", 
                 average = 5, 
                 w = 21, p = 3, m = 1, 
                 trim = c(350, 375, 1350, 1475, 1781, 1990, 2400, 2500),
                 new_col = T)


p <- plot_spectra(data = spc_pp, 
                  facets = NULL,
                  treatment = NULL,
                  col_in = "all",
                  topdf = T)

# ============================================================================================================= -
# Add design ----
# ============================================================================================================= -

# load design (as measurement meta data)
design <- read_csv("Z:/Public/Jonas/003_ESWW/Design_20211001/final_design_gen_rev.csv") %>% 
  dplyr::select(1:8, gen_name, dis_treat_new)

# add design and drop the reference measurements
spc_pp$plot_UID <- strsplit(spc_pp$meas_id, "_") %>% lapply("[[", 3) %>% unlist()
spc_pp <- spc_pp %>% 
  right_join(design, .,by = "plot_UID") %>% 
  dplyr::filter(!grepl("Ref", plot_UID))

# plot spectra
plot_spectra(data = spc_pp, 
             facets = "meas_date", 
             col_in = "all",
             treatment = "dis_treat_new", 
             mark_outliers = F,
             topdf = T)

# ============================================================================================================= -
# Calculate SVI ----
# ============================================================================================================= -

data$plot_UID <- strsplit(data$meas_id, "_") %>% lapply("[[", 3) %>% unlist()

spc_ref <- data %>% 
  right_join(design, ., by = "plot_UID") %>% 
  dplyr::filter(!grepl("Ref", plot_UID))

SVI <- spc_ref %>%
  preprocess_spc(average = 5, 
                 w = 21, p = 3, m = 0, 
                 new_col = F) %>% 
  calculate_SVI(col_in = "rflt_p3w21m0")
  
SVI_ <- scale_SVI(data = SVI, 
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
# Process spectra ----
# ============================================================================================================= -

SVI <- readRDS(paste0(path_to_data, "SVIdat.rds"))

SVI_ <- scale_SVI(data = SVI, 
                  plotid = "plot_UID", 
                  treatid = "dis_treat_new", 
                  plot = F)

# svi <- c("NDVI_nb_ASD", "MCARI2", "FII", "NDWI1650", "mND705", "780_740",
#          "PSRI", "CARG", "PRI$", "WI_NDVI")

svi <- names(SVI_$SVI[[1]])

parameters <- get_svi_dynamics(data = SVI_, 
                               svi = svi,
                               method = c("cgom", "pspl"),
                               timevar = "dafm", 
                               plot_dynamics = F)

saveRDS(parameters, paste0(path_to_data, "PARdat.rds"))

fits <- parameters$fits
parameters <- parameters$parameters

indis <- extract_sensitivity_indicators(fits, parameters)

integrals <- indis[[1]]
differences <- indis[[2]]

# ============================================================================================================= -
# Process scorings ---- 
# ============================================================================================================= -

path_to_data = "Z:/Public/Jonas/Data/ESWW006/RefData/GLA/"
SEN <- readRDS(paste0(path_to_data, "all_sen_dat.rds"))
# load design (as measurement meta data)
design <- read_csv("Z:/Public/Jonas/003_ESWW/Design_20211001/final_design_gen_rev.csv") %>% 
  dplyr::select(1:8, gen_name, dis_treat_new)

SEN <- SEN %>% 
  right_join(design, ., by = "plot_UID")

names(SEN)[12] <- c("SVI")
 
SEN_ <- scale_SVI(data = SEN, 
                  plotid = "plot_UID", 
                  treatid = "dis_treat_new", 
                  plot = F)

svi <- c("Senescence_plot", "Senescence_Fl0")

# subset <- SEN_ %>% filter(Plot_ID %in% c("ESWW0060146"))

parameters <- get_svi_dynamics(data = SEN_, 
                               svi = svi,
                               method = c("linear", "cgom", "fgom", "pspl"),
                               timevar = "dafm", 
                               plot_dynamics = T)

plot_parameters(data = parameters$parameters)

# ============================================================================================================= -

