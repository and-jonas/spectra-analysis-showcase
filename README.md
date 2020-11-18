# Spectral Data Analysis Showcase

This repository contains functions for a few standard steps in spectral data processing, visualizaiton, and analysis. The aim is to enable people with little experience in spectral data analysis and R programming to perform basic quality control analysis of spectral datasets. 

## Author

> Jonas Anderegg  
> Crop Science Group  
> ETH Zürich  


## Dependencies

The R-packages `data.table`, `mvoutlier`, `prospectr`, `ggsci` and `tidyverse` are required. These are automatically loaded, if necessary, in the demo `workflow`. 

## Content

### Spectra processing functions

`load_spectra()`: Loads spectral data from .asd or .sed files into R
`preprocess_spc()`: Implements different forms of signal pre-processing, such as smoothing, calculation of derivatives, calculation of standard normal variate, signal binning.
`detect_outlier_spectra()`: Performs multivariate outlier detection based on robust PCA, using the `mvoutlier::pcout`function. 
`Calculate_SVI()`: Calculates a large number of published spectral vegetation indices. 
`scale_SVI()`: Scales spectral vegetation indices at the plot level to range from 10 to 0, representing the highest and lowest value recorded during a measurement campaign comprising several consecutive measurements over time. 
`plot_spectra()`: Creates plots of spectra or derived data (such as e.g., their derivatives). 

### Example workflow

The file workflow.R implements an example workflow. 

### Demo

The file demo.pdf demonstrates some of the capabilities of the functions and carries out an example analysis. 

