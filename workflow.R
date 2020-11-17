
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

# set working directory
path_to_data <- "O:/Projects/KP0011/8/test/data/"
path_to_funcs <- "O:/Projects/KP0011/8/test/funcs/"

# load functions
source(paste0(path_to_funcs, "spectra_proc.R"))

# load all spectral data and save as .rds
data <- load_spectra(dir = path_to_data, format = "sed")
saveRDS(data, paste0(path_to_data, "alldat.rds"))

# ============================================================================================================= -

# load spectral data
data <- readRDS(paste0(path_to_data, "alldat.rds"))

plot_spectra(data, col_in = "rflt")

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
  detect_outlier_spectra(grouping = NULL, 
                         outliers_rm = "rflt_p3w21m0_trim_bin3", 
                         create_plot = T)

# add measurement meta data
meta <- read_csv(paste0(basedir, "Spectral_Data/spectral_files_fusarium_asign.csv"))
dd <- right_join(meta, spc_pp_o[[1]], by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")

# plot spectra
plot_spc(dd, facets = "meas_date",
         treatment = "Treatment", mark_outliers = T)

# ============================================================================================================= -

# add measurement meta data
dd <- right_join(meta, spc_pp, by = c("sed_new" = "meas_id")) %>% 
  dplyr::rename("meas_id" = "sed_new")

# calculate spectral indices: use all data (including outliers...)
SVI <- calculate_SVI(data = dd, col_in = "rflt_p3w21m0") %>% 
  scale_SVI(plot = F)
saveRDS(SVI, paste0(basedir, "Spectral_Data/spc_pp.rds"))

# ============================================================================================================= -
# ============================================================================================================= -




# add model fits
dings <- SVI %>% dplyr::select(Plot_ID, Treatment, dafm, SVI_sc) %>% unnest(SVI_sc) %>% 
  tidyr::gather(SVI, value, starts_with("SI_"))

# add model fits
data_fits <- dings %>% 
  group_by(Plot_ID, Treatment, SVI) %>% 
  # remove plots with incomplete time series
  nest() %>% 
  # interpolate linearly between measurement time points 
  mutate(fit_lin = purrr::map(data, lin_approx, n_meas = 6) %>% purrr::map(cbind.data.frame))

# extract dynamics parameters from nls fits
dynpars_lin <- data_fits %>% 
  mutate(pars_lin = purrr::map(fit_lin, extract_pars))

bla <- dynpars_lin %>% dplyr::select(Plot_ID, Treatment, SVI, pars_lin) %>% 
  mutate(pars_lin = purrr::map(pars_lin, tibble::as_tibble)) %>% 
  unnest(pars_lin)

dur <- bla %>% dplyr::select(Plot_ID, Treatment, SVI, dur1, dur2) %>% 
  tidyr::gather(param, value, dur1:dur2)

tps <- bla %>% dplyr::select(Plot_ID, Treatment, SVI, t80:t20) %>% 
  tidyr::gather(param, value, t80:t20)

plot <- ggplot(dur) +
  geom_boxplot(aes(x = param, y = value, group = interaction(param, Treatment), fill = Treatment)) +
  facet_wrap(~SVI) +
  ggsci::scale_fill_npg() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank())

plot <- ggplot(tps) +
  geom_boxplot(aes(x = param, y = value, group = interaction(param, Treatment), fill = Treatment)) +
  facet_wrap(~SVI) +
  ggsci::scale_fill_npg() +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank())
