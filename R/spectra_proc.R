
# ============================================================================================================= -

# main functions ----

#' Loads spectra files from working directory
#' @param dir Directory to read files from
#' @param format file extension as a character string, "sed" or "asd"
#' @return A tibble with measurement_id, measurement date and spectra in a list column of data.tables
#' @details Files tp read must be stored in subfolders with measurement dates as folder names, e.g. "20200709"
#' File names must start with a four digit number indicating the year of measurement 
#' and finish with the specified file extension, i.e. ".sed" or ".asd". 
#' All files meeting these criteria are read. 
#' @export
load_spectra <- function(dir, subdir = NULL, format = "sed"){
  # print("loading spetra ...")
  # get all subdirectories
  subdirs <- dir(path = dir, full.names = TRUE, recursive = FALSE, pattern = "^[0-9]{8}")
  # if spectra are stored in subdirectories
  if(!is.null(subdir)){
    subdirs <- paste(subdirs, subdir, sep = "/")
  }
  # get all filenames
  dirs_spc_files <- list.files(subdirs, pattern = paste0("[0-9]{4}.*", ".", format), full.names = TRUE)
  # load spectral data
  data <- dirs_spc_files %>% 
    # read list of files
    lapply(., read_spectrum, format = format) %>% 
    # bind to tible
    data.table::rbindlist(.,) %>% tibble::as_tibble()
  # print("done")
  return(data)
}

#' Pre-processes spectra 
#' @param data A tibble with spectra in list columns of data.tables, as returned by "load_spectra()".  
#' @param p Order of the polynomial used for smoothing of spectra, defaults to p = 3
#' @param w Window size in wavebands, defaults to w = 21
#' @param m Order of the derivative calculated, defaults to m = 0, i.e. original spectra
#' @param col_in Column of the tibble to which the pre-processing is applied, defaults to "rflt", i.e. raw reflectance spectra
#' @param new_col Boolean, whether to add output as a new list column or overwrite the input column
#' @param snv Boolean, whether to calculate the standard normal variate or not
#' @param trim A numeric vector specifying the lower and upper bounds of spectral regions to exclude, or NULL 
#' @param binning numeric or NULL, specifying the bin width
#' @return The modified input tibble, with THE overwritten or added output in a list column
#' @export
preprocess_spc <- function(data,
                           average = NULL,
                           p = 3, w = 21, m = 0, 
                           col_in = "rflt", 
                           new_col = FALSE, 
                           snv = FALSE,
                           trim = NULL,
                           binning = NULL){
  
  # print("pre-processing spectra ...")
  
  # Convert list of data.tables to one data.table
  spc_raw <- data.table::rbindlist(data[[col_in]])
  
  if(!is.null(average)){
    spc_raw <- setDT(spc_raw)[, as.list(colMeans(.SD)), 
                              by = as.integer(c(0, rep(1:(nrow(spc_raw)-1)%/%average)))]
  }
  
  # apply SavitzkyGolay filter to smooth spectra
  if(p != 0){
    # set column names depending on selected parameters
    ptype <- paste0("p", p, "w", w, "m", m)
    spc_pp <- with(data, prospectr::savitzkyGolay(spc_raw, p = p, w = w, m = m))
  } else {
    ptype = gsub("rflt_", "", col_in)
    spc_pp <- spc_raw %>% as.matrix()
  }
  
  # calculate the standard normal variate (AFTER filtering; Fearn, 2008)
  if(snv){
    spc_pp <- with(data, prospectr::standardNormalVariate(spc_pp))
  }
  
  # remove regions with low signal-to-noise ratio
  if(!is.null(trim)){
    drop_vars <- get_regions_drop(trim) 
    spc_pp <- spc_pp[, !colnames(spc_pp) %in% drop_vars]
  }
  
  # bin spectra
  if(!is.null(binning)){
    spc_pp <- prospectr::binning(X = spc_pp, bin.size = binning)
  }
  
  # convert to list of data.tables
  spc_pp <- data.table::as.data.table(spc_pp)
  spc_pre_list <- map(purrr::transpose(spc_pp), data.table::as.data.table)
  
  # define column name to represent pre-processing
  colname_out <- paste0("rflt_", ptype)
  if(snv){
    colname_out <- paste0(colname_out, "_snv")
  }
  if(!is.null(trim)){
    colname_out <- paste0(colname_out, "_trim")
  }
  if(!is.null(binning)){
    colname_out <- paste0(colname_out, "_bin", binning)
  }
  
  # replace raw spectra with pre-processed spectra
  if(!is.null(average)){
    data <- data[seq(1, nrow(data), average), ]
  }
  if(new_col == F){
    spc_pp <- data[,-which(colnames(data) == col_in)]
  } else {
    spc_pp <- data
  }
  spc_pp[, colname_out] <- list(spc_pre_list)
  
  return(spc_pp)
  
}


#' Detects multivariate outliers in spectral datasets 
#' @param data A tibble with spectra in list columns of data.tables, as returned by "load_spectra()" and "preprocess_spc". 
#' @param col_in Column of the tibble to which the pre-processing is applied, defaults to "all", 
#' i.e. using all types of spectra contained in list columns of the tibble
#' @param grouping A character vector or NULL, specifying the structure of the spectral dataset. 
#' If a grouping variable is supplied, outliers are detected within each group
#' @param outliers_rm Character vector or NULL, specifying which type of outliers should be removed (i.e. which of the col_in).
#' If "all", all samples detected as an outlier in any of the col_in will be removed from the dataset
#' @param create_plot Boolean, whether the outliers should be plotted
#' @return A list, consisting of a tibble with list columns, and a list of measurement_ids detected as mvoutliers,
#' depending on the spectral input data type (col_in).
detect_outlier_spectra <- function(data, col_in = "all",
                                   grouping = NULL,
                                   outliers_rm = NULL,
                                   create_plot = T,
                                   topdf = F){
  
  # print("detecting multivariate outliers ...")
  
  # fix grouping structure if exists
  if(!is.null(grouping)){
    group_list <- split(data, data[grouping])
  } else{
    group_list <- list(data)
  }
  
  if(col_in == "all"){
    col_in <- grep("^rflt", names(data), value = TRUE)
  }
  
  # iterate over data types
  out_meas <- list()
  for (j in col_in){
    
    # detect outliers within each group of measurements
    v <- list()
    for (i in 1:length(group_list)){
      dat <- group_list[[i]]
      # Convert list of data.tables to one data.table
      spc <- data.table::rbindlist(dat[[j]])
      out <- mvoutlier::pcout(spc, makeplot = F, outbound = 0.05)
      outliers <- which(out$wfinal01 == 0)
      v[[i]] <- numeric(nrow(dat))
      v[[i]][outliers] <- 1
    }
    
    otype <- paste0("out_", j)
    outlier <- unlist(v)
    data[, otype] <- outlier
    
    if(create_plot){
      plot_outliers(data, j, mark_outliers = T, 
                    facets = grouping, topdf = topdf)
    }
    
    out_idx <- which(outlier == 1)
    out_meas[[j]] <- data$meas_id[out_idx]
    
  }
  
  # drop outliers
  if(!is.null(outliers_rm)){
    
    # print("removing outliers from dataset ...")
    
    # get index of columns to drop
    if("all" %in% outliers_rm){
      out_colidx <- grep("^out_", colnames(data))
    } else {
      out_colidx <- grep(paste0("out_", outliers_rm, collapse = "|"), colnames(data))
    }
    outids <- data[out_colidx]
    drop_idx <- which(rowSums(outids) > 0)
    data <- data[-drop_idx,]
    data <- data[,-grep("^out_", names(data))]
  }
  
  return(list(data, out_meas))
  
}

#' Calculates spectral vegetation indices
#' @param data A tibble with spectra in list columns of data.tables, as returned by "load_spectra()" and "preprocess_spc"
#' @param col_in Column name of the tibble containing the spectra to use as character string
#' @param new_col Boolean, whether to add output as a new list column or overwrite the input column
#' @return A tibble with list columns
#' @export
calculate_SVI <- function(data, col_in, new_col = T) {
  
  # print("calculating spectral indices ...")
  
  d <- data.table::rbindlist(data[[col_in]]) %>% tibble::as_tibble()
  
  # BANDS ########################################################################## -
  ALI5 <- rowMeans(d[775:805])
  ALI6 <- rowMeans(d[ ,match(845:890, colnames(d))])
  ALI7 <- rowMeans(d[ ,match(1200:1300, colnames(d))])
  ASTER4 <- rowMeans(d[ ,match(1600:1700, colnames(d))])
  ASTER5 <- rowMeans(d[ ,match(2145:2185, colnames(d))])
  ASTER6 <- rowMeans(d[ ,match(2185:2225, colnames(d))])
  ASTER7 <- rowMeans(d[ ,match(2235:2285, colnames(d))])
  ASTER8 <- rowMeans(d[ ,match(2295:2365, colnames(d))])
  MODIS5 <- rowMeans(d[ ,match(1230:1250, colnames(d))])
  MODIS6 <- rowMeans(d[ ,match(1628:1652, colnames(d))])
  MODIS7 <- rowMeans(d[ ,match(2105:2155, colnames(d))])
  R1050_ASD <- rowMeans(d[ ,match(1046:1054, colnames(d))])
  R1090to1110 <- rowMeans(d[ ,match(1090:1110, colnames(d))])
  R1094_HYP <- rowMeans(d[ ,match(1089:1099, colnames(d))]) #B095 1094.09 10.99
  R1100_5 <- rowMeans(d[ ,match(1098:1102, colnames(d))]) #interpolated ASD
  R1180to1220 <- rowMeans(d[ ,match(1180:1220, colnames(d))])
  R1200_5 <- rowMeans(d[ ,match(1198:1202, colnames(d))]) #interpolated ASD
  R1205_HYP <- rowMeans(d[ ,match(1201:1209, colnames(d))]) #B106 1205.07 10.79
  R1220 <- rowMeans(d[ ,match(1216:1224, colnames(d))])
  R1225 <- rowMeans(d[ ,match(1221:1229, colnames(d))])
  R1240_AVIRIS <- rowMeans(d[ ,match(1236:1244, colnames(d))])
  R1250_ASD <- rowMeans(d[ ,match(1246:1254, colnames(d))])
  R1265to1285 <- rowMeans(d[ ,match(1265:1285, colnames(d))])
  R1300_100 <- rowMeans(d[ ,match(1250:1350, colnames(d))])
  R1450_100 <- rowMeans(d[ ,match(1400:1500, colnames(d))])
  R1510_AVIRIS <- rowMeans(d[ ,match(1506:1514, colnames(d))])
  R1555_ASD <- rowMeans(d[ ,match(1551:1559, colnames(d))])
  R1650 <- rowMeans(d[ ,match(1649:1651, colnames(d))])
  R1650_CS <- rowMeans(d[ ,match(1550:1750, colnames(d))])
  R1659_HYP <- rowMeans(d[ ,match(1654:1664, colnames(d))]) #B151 1659 11.53
  R1680_AVIRIS <- rowMeans(d[ ,match(1676:1684, colnames(d))])
  R1754_AVIRIS <- rowMeans(d[ ,match(1750:1758, colnames(d))])
  R2020_IRIS <- rowMeans(d[ ,match(2019:2021, colnames(d))])
  R2100_IRIS <- rowMeans(d[ ,match(2099:2101, colnames(d))])
  R2220_IRIS <- rowMeans(d[ ,match(2219:2221, colnames(d))])
  R360to370 <- rowMeans(d[ ,match(360:370, colnames(d))])
  R430_6 <- rowMeans(d[ ,match(428:432, colnames(d))])
  R440_SD <- rowMeans(d[ ,match(438:442, colnames(d))]) #bandwidth 5nm
  R445_SIRIS <- rowMeans(d[ ,match(444:446, colnames(d))]) 
  R470_SIRIS <- rowMeans(d[ ,match(469:471, colnames(d))])
  R470to480 <- rowMeans(d[ ,match(470:480, colnames(d))])
  R480_ASD <- rowMeans(d[ ,match(476:484, colnames(d))])
  R500 <- rowMeans(d[ ,match(499:501, colnames(d))])
  R500_SIRIS <- rowMeans(d[ ,match(499:501, colnames(d))])
  R500to599 <- rowMeans(d[ ,match(500:599, colnames(d))])
  R510 <- rowMeans(d[ ,match(509:511, colnames(d))])
  R510_10 <- rowMeans(d[ ,match(506:514, colnames(d))])
  R510to520 <- rowMeans(d[ ,match(510:520, colnames(d))])
  R513_ASD <- rowMeans(d[ ,match(512:514, colnames(d))])
  R520_ASD <- rowMeans(d[ ,match(519:521, colnames(d))])
  R520to585 <- rowMeans(d[ ,match(520:585, colnames(d))])
  R531_10 <- rowMeans(d[ ,match(527:535, colnames(d))])
  R531_6 <- rowMeans(d[ ,match(529:533, colnames(d))])
  R534_ASD <- rowMeans(d[ ,match(533:535, colnames(d))])
  R540to560 <- rowMeans(d[ ,match(540:560, colnames(d))])
  R550 <- rowMeans(d[ ,match(549:551, colnames(d))])
  R550_30 <- rowMeans(d[ ,match(535:565, colnames(d))])
  R550_ASD <- rowMeans(d[ ,match(549:551, colnames(d))])
  R550_CASI <- rowMeans(d[ ,match(547:553, colnames(d))])
  R550_HYP <- rowMeans(d[ ,match(543:553, colnames(d))]) #B020 548.92 11.02
  R550_SE3 <- rowMeans(d[ ,match(549:551, colnames(d))]) #doublons!!
  R550_YARA <- rowMeans(d[ ,match(549:551, colnames(d))]) #doublons!!
  R560_10 <- rowMeans(d[ ,match(556:564, colnames(d))])
  R560_MSR <- rowMeans(d[ ,match(556:564, colnames(d))]) #bandwidth 8.7nm
  R570_10 <- rowMeans(d[ ,match(566:574, colnames(d))])
  R570_6 <- rowMeans(d[ ,match(568:572, colnames(d))])
  R570_ASD <- rowMeans(d[ ,match(569:571, colnames(d))])
  R573_SD <- rowMeans(d[ ,match(567:579, colnames(d))]) #bandwidth 14nm
  R584_ASD <- rowMeans(d[ ,match(583:585, colnames(d))])
  R600 <- rowMeans(d[ ,match(599:601, colnames(d))])
  R600to699 <- rowMeans(d[ ,match(600:699, colnames(d))])
  R630 <- rowMeans(d[ ,match(629:631, colnames(d))])
  R650_SIRIS <- rowMeans(d[ ,match(649:651, colnames(d))])
  R660_MSR <- rowMeans(d[ ,match(656:664, colnames(d))]) #bandwidth 9.4nm
  R670 <- rowMeans(d[ ,match(669:671, colnames(d))])
  R670_10 <- rowMeans(d[ ,match(666:674, colnames(d))])
  R670_5 <- rowMeans(d[ ,match(668:672, colnames(d))]) #interpolated ASD
  R670_3 <- rowMeans(d[ ,match(669:671, colnames(d))])
  R670_CASI <- rowMeans(d[ ,match(667:673, colnames(d))]) #see above Rr_CASI
  R670_SE3 <- rowMeans(d[ ,match(669:671, colnames(d))]) #doublons!!
  R670_YARA <- rowMeans(d[ ,match(669:671, colnames(d))]) #doublons!!
  R675_SIRIS <- rowMeans(d[ ,match(674:676, colnames(d))])
  R678 <- rowMeans(d[ ,match(677:679, colnames(d))])
  R680 <- rowMeans(d[ ,match(679:681, colnames(d))])
  R680_6 <- rowMeans(d[ ,match(678:682, colnames(d))])
  R680_MERIS <- rowMeans(d[ ,match(678:684, colnames(d))]) #B08 681.25 7.5
  R680_SIRIS <- rowMeans(d[ ,match(679:681, colnames(d))]) #doublons!!
  R680_UNI <- rowMeans(d[ ,match(676:684, colnames(d))])
  R681_HYP <- rowMeans(d[ ,match(676:686, colnames(d))]) #B033 681.2 10.33
  R683 <- rowMeans(d[ ,match(682:684, colnames(d))])
  R690to710 <- rowMeans(d[ ,match(690:710, colnames(d))])
  R690to720 <- rowMeans(d[ ,match(690:720, colnames(d))])
  R695to740 <- rowMeans(d[ ,match(695:740, colnames(d))])
  R698_ASD <- rowMeans(d[ ,match(697:699, colnames(d))])
  R700 <- rowMeans(d[ ,match(699:701, colnames(d))])
  R700_10 <- rowMeans(d[ ,match(696:704, colnames(d))])
  R700_15 <- rowMeans(d[ ,match(693:707, colnames(d))])
  R700_ASD <- rowMeans(d[ ,match(699:701, colnames(d))]) #doublons!!
  R700_CASI <- rowMeans(d[ ,match(697:703, colnames(d))])
  R700_SE3 <- rowMeans(d[ ,match(699:701, colnames(d))]) #doublons!!
  R700_YARA <- rowMeans(d[ ,match(699:701, colnames(d))]) #doublons!!
  R700_ASD <- rowMeans(d[ ,match(699:701, colnames(d))])
  R704_ASD <- rowMeans(d[ ,match(703:705, colnames(d))])
  R704 <- rowMeans(d[ ,match(703:705, colnames(d))])
  R705_HYP <- rowMeans(d[ ,match(698:706, colnames(d))]) #B035 701.55 10.46
  R705 <- rowMeans(d[ ,match(698:706, colnames(d))])
  R710 <- rowMeans(d[ ,match(704:706, colnames(d))])
  R710_MERIS <- rowMeans(d[ ,match(705:713, colnames(d))]) #B09 708.75 10
  R715_ASD <- rowMeans(d[ ,match(714:716, colnames(d))])
  R720_10 <- rowMeans(d[ ,match(716:724, colnames(d))])
  R720_ASD <- rowMeans(d[ ,match(719:721, colnames(d))])
  R720_CASI <- rowMeans(d[ ,match(717:723, colnames(d))])
  R724_ASD <- rowMeans(d[ ,match(723:725, colnames(d))])
  R726_ASD <- rowMeans(d[ ,match(725:727, colnames(d))])
  R730_YARA <- rowMeans(d[ ,match(729:731, colnames(d))])
  R734_ASD <- rowMeans(d[ ,match(733:735, colnames(d))])
  R737_ASD <- rowMeans(d[ ,match(736:738, colnames(d))])
  R740_ASD <- rowMeans(d[ ,match(739:741, colnames(d))]) #doublons!!
  R740_YARA <- rowMeans(d[ ,match(739:741, colnames(d))])
  R745_ASD <-  rowMeans(d[ ,match(744:746, colnames(d))])
  R747_ASD <- rowMeans(d[ ,match(746:748, colnames(d))])
  R750 <- rowMeans(d[ ,match(749:751, colnames(d))])
  R750_CASI <- rowMeans(d[ ,match(747:753, colnames(d))])
  R750_HYP <- rowMeans(d[ ,match(748:756, colnames(d))]) #B040 752.43 10.71
  R750_MERIS <- rowMeans(d[ ,match(751:757, colnames(d))]) #B10 753.75 7.5
  R750to800 <- rowMeans(d[ ,match(750:800, colnames(d))]) #doublons!! see Rn_1
  R760_YARA <- rowMeans(d[ ,match(759:761, colnames(d))])
  R760to800 <- rowMeans(d[ ,match(760:800, colnames(d))])
  R780_YARA <- rowMeans(d[ ,match(779:781, colnames(d))])
  R790_10 <- rowMeans(d[ ,match(786:794, colnames(d))])
  R800 <- rowMeans(d[ ,match(799:801, colnames(d))])
  R800_10 <- rowMeans(d[ ,match(796:804, colnames(d))])
  R800_ASD <- rowMeans(d[ ,match(799:801, colnames(d))])
  R800_CASI <- rowMeans(d[ ,match(797:803, colnames(d))])
  R800_HYP <- rowMeans(d[ ,match(795:805, colnames(d))]) #B045 803.3 11.1
  R800_SIRIS <- rowMeans(d[ ,match(799:801, colnames(d))])
  R800_UNI <- rowMeans(d[ ,match(796:804, colnames(d))])
  R810_10 <- rowMeans(d[ ,match(806:814, colnames(d))])
  R810_MSR <- rowMeans(d[ ,match(805:815, colnames(d))]) #bandwidth 11.2nm
  R830 <- rowMeans(d[ ,match(829:831, colnames(d))])
  R840_CS <- rowMeans(d[ ,match(770:910, colnames(d))])
  R850 <- rowMeans(d[ ,match(849:851, colnames(d))])
  R850_5 <- rowMeans(d[ ,match(848:852, colnames(d))]) #interpolated ASD
  R860_AVIRIS <- rowMeans(d[ ,match(856:864, colnames(d))])
  R874 <- rowMeans(d[ ,match(856:864, colnames(d))])
  R880 <- rowMeans(d[ ,match(879:881, colnames(d))])
  R893_HYP <- rowMeans(d[ ,match(888:898, colnames(d))]) #B075 892.28 11.05
  R900 <- rowMeans(d[ ,match(899:901, colnames(d))])
  R900_UNI <- rowMeans(d[ ,match(896:904, colnames(d))])
  R900_YARA <- rowMeans(d[ ,match(899:901, colnames(d))])
  R920 <- rowMeans(d[ ,match(919:921, colnames(d))])
  R920to940 <- rowMeans(d[ ,match(920:940, colnames(d))])
  R955 <-  rowMeans(d[ ,match(954:956, colnames(d))])
  R960to990 <- rowMeans(d[ ,match(960:990, colnames(d))])
  R970 <- rowMeans(d[ ,match(969:971, colnames(d))])
  R970_UNI <- rowMeans(d[ ,match(966:974, colnames(d))])
  R970_YARA <- rowMeans(d[ ,match(969:971, colnames(d))])
  Rb_1 <- rowMeans(d[ ,match(460:480, colnames(d))]) #Kaufman and Tanre 1992
  Rb_bb <- rowMeans(d[ ,match(450:520, colnames(d))])
  Rb_cb <- rowMeans(d[ ,match(400:520, colnames(d))])
  Rb_MODIS3 <- rowMeans(d[ ,match(459:479, colnames(d))])
  Re_MERIS <- rowMeans(d[ ,match(700:710, colnames(d))])
  Rg_1 <- rowMeans(d[ ,match(530:570, colnames(d))]) #MODIS band?
  Rg_bb <- rowMeans(d[ ,match(520:600, colnames(d))])
  Rg_cb <- rowMeans(d[ ,match(580:610, colnames(d))])
  Rg_MODIS12 <- rowMeans(d[ ,match(546:556, colnames(d))])
  Rg_MODIS4 <- rowMeans(d[ ,match(545:565, colnames(d))])
  Rg_RE2 <- rowMeans(d[ ,match(520:590, colnames(d))]) #RapidEye B2
  Rn_1 <- rowMeans(d[ ,match(845:885, colnames(d))]) #Kaufman and Tanre 1992
  Rn_2 <- rowMeans(d[ ,match(750:800, colnames(d))]) #doublons!!
  Rn_3 <- rowMeans(d[ ,match(800:900, colnames(d))])
  Rn_4 <- rowMeans(d[ ,match(750:900, colnames(d))])
  Rn_AVHRR <- rowMeans(d[ ,match(750:1000, colnames(d))])
  Rn_bb <- rowMeans(d[ ,match(760:900, colnames(d))]) #doublons!!
  Rn_CASI <- rowMeans(d[ ,match(857:863, colnames(d))])
  Rn_MERIS <- rowMeans(d[ ,match(750:757, colnames(d))])
  Rn_MODIS2 <- rowMeans(d[ ,match(841:876, colnames(d))])
  Rn_RE5 <- rowMeans(d[ ,match(760:850, colnames(d))]) #RapidEye B5
  Rn_SPOT <- rowMeans(d[ ,match(790:890, colnames(d))])
  Rn_TM4 <- rowMeans(d[ ,match(760:900, colnames(d))]) #doublons!!
  Rr_1 <- rowMeans(d[ ,match(635:685, colnames(d))]) #Kaufman and Tanre 1992 
  Rr_ALI4 <- rowMeans(d[ ,match(630:690, colnames(d))]) #doublons!!
  Rr_ASTER2 <- rowMeans(d[ ,match(630:690, colnames(d))]) #doublons!!
  Rr_AVHRR <- rowMeans(d[ ,match(580:680, colnames(d))])
  Rr_bb <- rowMeans(d[ ,match(630:690, colnames(d))]) #doublons!!
  Rr_CASI <- rowMeans(d[ ,match(667:673, colnames(d))])
  Rr_cb <- rowMeans(d[ ,match(580:660, colnames(d))])
  Rr_MODIS1 <- rowMeans(d[ ,match(620:670, colnames(d))])
  Rr_RE3 <- rowMeans(d[ ,match(630:685, colnames(d))]) #RapidEye B3
  Rr_SPOT <- rowMeans(d[ ,match(610:690, colnames(d))])
  Rr_TM3 <- rowMeans(d[ ,match(630:690, colnames(d))]) #doublons!!
  Rre_bb <- rowMeans(d[ ,match(700:730, colnames(d))])
  Rre_MERIS <- rowMeans(d[ ,match(704:714, colnames(d))])
  Rre_RE4 <- rowMeans(d[ ,match(690:730, colnames(d))]) #RapidEye B4
  TM5 <- rowMeans(d[ ,match(1550:1750, colnames(d))])
  TM7 <- rowMeans(d[ ,match(2080:2350, colnames(d))])
  # A_vis <- rowMeans(d[ ,match(350:700, colnames(d))])  #ALBEDO in respective region
  A_NIR <- rowMeans(d[ ,match(700:1000, colnames(d))])
  A_SWIR1 <- rowMeans(d[ ,match(1000:1350, colnames(d))]) 
  A_SWIR2 <- rowMeans(d[ ,match(1440:1810, colnames(d))])
  A_SWIR3 <- rowMeans(d[ ,match(1940:2400, colnames(d))])
  R <- rowMeans(d[ ,match(675:685, colnames(d))])
  G <- rowMeans(d[ ,match(545:555, colnames(d))])
  B <- rowMeans(d[ ,match(495:505, colnames(d))])
  
  # INDICES ########################################################################## -
  
  #ARVI
  SI_ARVI <- (Rn_1-(Rr_1-Rb_1+Rr_1))/(Rn_1+(Rr_1-Rb_1+Rr_1)) #with gamma equals 1 here
  
  #EVI/SARVI2
  SI_EVI <- 2.5*(Rn_MODIS2-Rr_MODIS1)/(Rn_MODIS2+6*Rr_MODIS1-7.5*Rb_MODIS3+1)
  
  #NDVI narrow band
  SI_NDVI_nb_ASD <- (R800_ASD-R670_3)/(R800_ASD+R670_3)
  SI_NDVI_nb_CASI <- (R800_CASI-R670_CASI)/(R800_CASI+R670_CASI)
  
  #NDVI MODIS
  SI_NDVI_MODIS <- (Rn_MODIS2-Rr_MODIS1)/(Rn_MODIS2+Rr_MODIS1)
  
  #NDVI broad band
  SI_NDVI_bb <- (Rn_2-Rr_bb)/(Rn_2+Rr_bb)
  SI_NDVI_bb2 <- (Rn_3-Rr_bb)/(Rn_3+Rr_bb)
  SI_NDVI_bb3 <- (Rn_4-Rr_bb)/(Rn_4+Rr_bb)
  
  #SAVI
  L <- 0.5
  SI_SAVI <- (1+L)*(Rn_CASI-Rr_CASI)/(Rn_CASI+Rr_CASI+L)
  
  #NGRDI 
  SI_NGRDI <- (Rg_bb-Rr_bb)/(Rg_bb+Rr_bb)
  
  #WDRVI
  a <- 0.1 #weighting coeff
  SI_WDRVI <- (a*Rn_AVHRR-Rr_AVHRR)/(a*Rn_AVHRR+Rr_AVHRR)
  
  #CAI
  f <- 100  #bands reflectance translation in %
  SI_CAI <-  0.5*(R2020_IRIS*f+R2220_IRIS*f)-R2100_IRIS*f
  
  #RRDI
  SI_RRDI <- (R745_ASD - R740_ASD)/(R740_ASD - R700_ASD)
  
  #mND705
  SI_mND705 <- (R750 - R705)/(R750 + R705 - 2*R445_SIRIS)
  
  #NDTI
  SI_NDTI <- (TM5-TM7)/(TM5+TM7)
  
  #LCA
  SI_LCA <- 100*((ASTER6-ASTER5)/(ASTER6-ASTER8))
  
  #NDRI
  SI_NDRI <- (Rr_TM3-TM7)/(Rr_TM3+TM7)
  
  #ANSI
  SI_ANSI <- (ASTER5-ASTER4)/(ASTER5+ASTER4)
  
  #SINDRI
  SI_SINDRI <- (ASTER6-ASTER7)/(ASTER6+ASTER7)
  
  #NDLI
  SI_NDLI = (log(1/R1754_AVIRIS)-log(1/R1680_AVIRIS))/(log(1/R1754_AVIRIS)+log(1/R1680_AVIRIS))
  
  #SR_dm
  SI_SR_dm <- R780_YARA/R670_YARA
  
  #REIP
  SI_REIP <- 700 + 40*((R670_YARA+R780_YARA)/2-R700_YARA)/(R740_YARA-R700_YARA)
  
  #760/730 ratio
  SI_760_730 <- R760_YARA/R730_YARA
  
  #780/740 ratio
  SI_780_740 <- R780_YARA/R740_YARA
  
  #780/700 ratio  
  SI_780_700	<-R780_YARA/R700_YARA
  
  #780/550 ratio
  SI_780_550 <-R780_YARA/R550_YARA
  
  #970/900 ratio
  SI_970_900 <- R970_YARA/R900_YARA
  
  #MCARI1
  SI_MCARI1 <- 1.2*(2.5*(R800_CASI-Rr_CASI)-1.3*(R800_CASI-R550_CASI))
  
  #MTVI1
  SI_MTVI1 <- 1.2*(1.2*(R800_CASI-R550_CASI)-2.5*(R670_CASI-R550_CASI))
  
  #MCARI2
  SI_MCARI2 <- (1.5*(2.5*(R800_CASI-R670_CASI)-1.3*(R800_CASI-R550_CASI)))/sqrt((2*R800_CASI+1)^2-(6*R800_CASI-5*sqrt(R670_CASI))-0.5)
  
  #MTVI2
  SI_MTVI2 <- (1.5*(1.2*(R800_CASI-R550_CASI)-2.5*(R670_CASI-R550_CASI)))/sqrt((2*R800_CASI+1)^2-(6*R800_CASI-5*sqrt(R670_CASI))-0.5)
  
  #GLI
  SI_GLI <- (2*Rg_cb-Rr_cb-Rb_cb)/(2*Rg_cb+Rr_cb+Rb_cb)
  
  #OSAVI
  SI_OSAVI <- (1+0.16)*(R800_CASI-R670_CASI)/(R800_CASI+R670_CASI+0.16)
  
  #MSAVI
  SI_MSAVI <- 0.5*(2*Rn_SPOT+1-sqrt((2*Rn_SPOT+1)^2-8*(Rn_SPOT-Rr_SPOT)))
  
  #SARVI
  L <- 0.5
  SI_SARVI <- (1+L)*(Rn_1-(Rr_1-Rb_1+Rr_1))/(Rn_1+(Rr_1-Rb_1+Rr_1+L)) #with gamma equals 1 here  #VARIgreen
  SI_VARIgreen <- (Rg_MODIS12-Rr_MODIS1)/(Rg_MODIS12+Rr_MODIS1-Rb_MODIS3)
  
  #VARI700
  SI_VARI700 <- (Re_MERIS-1.7*Rr_MODIS1+0.7*Rb_MODIS3)/(Re_MERIS+2.3*Rr_MODIS1-1.3*Rb_MODIS3)
  
  #VIgreen
  SI_VIgreen <- (Rg_MODIS12-Rr_MODIS1)/(Rg_MODIS12+Rr_MODIS1)
  
  #VI700
  SI_VI700 <- (Re_MERIS-Rr_MODIS1)/(Re_MERIS+Rr_MODIS1)
  
  #SGR
  SI_SGR <- R500to599
  
  #SR_lai
  SI_SR_lai_bb <- Rn_bb/Rr_bb
  SI_SR_lai_cb <- Rn_bb/Rr_cb #no NIR-band for camera
  
  #SLAIDI
  S <- 5 #scaling factor
  SI_SLAIDI1 <- S*(R1050_ASD-R1250_ASD)/(R1050_ASD+R1250_ASD)
  S <- 40 #scaling factor
  SI_SLAIDI2 <- S*R1555_ASD*(R1050_ASD-R1250_ASD)/(R1050_ASD+R1250_ASD)
  
  #SIPI
  SI_SIPI <- (R800_SIRIS-R445_SIRIS)/(R800_SIRIS-R680_SIRIS)
  
  #PSSR
  SI_PSSR1 <- R800_SIRIS/R675_SIRIS
  SI_PSSR2 <- R800_SIRIS/R650_SIRIS
  SI_PSSR3 <- R800_SIRIS/R500_SIRIS
  
  #PSND
  SI_PSND1 <- (R800_SIRIS-R675_SIRIS)/(R800_SIRIS+R675_SIRIS)
  SI_PSND2 <- (R800_SIRIS-R650_SIRIS)/(R800_SIRIS+R650_SIRIS)
  SI_PSND3 <- (R800_SIRIS-R500_SIRIS)/(R800_SIRIS+R500_SIRIS)
  SI_PSND4 <- (R800_SIRIS-R470_SIRIS)/(R800_SIRIS+R470_SIRIS)
  
  SI_PSRI <- (R678-R500)/R750
  
  #NPCI
  SI_NPCI <- (R680_6-R430_6)/(R680_6+R430_6)
  
  #PRI
  SI_PRI <- (R531_10-R570_10)/(R531_10+R570_10)
  
  #PRI570
  SI_PRI570_10 <- (R570_10-R531_10)/(R570_10+R531_10)
  SI_PRI570_6 <- (R570_6-R531_6)/(R570_6+R531_6)
  
  #PRInorm
  SI_PRInorm <- ((R570_10-R531_10)/(R570_10+R531_10))/(((R800_10-R670_10)/sqrt(R800_10+R670_10))*(R700_10/R670_10))
  
  #RGR
  SI_RGR <- R683/R510
  SI_RGR2 <- R600to699/R500to599
  
  #ANTH
  SI_ANTH <- R760to800*(1/R540to560-1/R690to710)
  
  #ARI
  SI_ARI <- 1/(R550*f)-1/(R700*f)
  
  #CARG
  SI_CARG <- R760to800*(1/R510to520-1/R540to560)
  
  #CARRE
  SI_CARRE <- R760to800*(1/R510to520-1/R690to710)
  
  #CRI1
  SI_CRI1 <- 1/R510_10-1/R550_30
  
  #CRI2
  SI_CRI2 <- 1/R510_10-1/R700_15
  
  #CHLG
  SI_CHLG <- R760to800/R540to560
  
  #CHLRE
  SI_CHLRE <- R760to800/R690to720-1
  
  #LCI
  SI_LCI <- (R850-R710)/(R850-R680)
  SI_LCI2 <- (R850-R710)/(R850+R680) #in Pu 2012*
  
  #GNDVI
  SI_GNDVI_HI <- (R750-Rg_1)/(R750+Rg_1) #in Gitelson, Kaufman and Merzlyak 1996
  SI_GNDVI_IRIS <- (R750-R550)/(R750+R550) #in datat 1998
  
  #NDREI
  SI_NDREI <- (R750-R705)/(R750+R705)
  #SI_NDREI2 <- (Rn_bb-Rre_bb)/(Rn_bb+Rre_bb) #in Hunt et al 2011
  
  #NDRE
  SI_NDRE <- (R790_10-R720_10)/(R790_10+R720_10)
  
  #CIG
  SI_CIG <- R750to800/R520to585-1
  #SI_CIG2 <- Rn_bb/Rg_bb-1 #in Hunt et al 2011
  
  #CIRE
  SI_CIRE <- R750to800/R695to740-1
  #SI_CIRE2 <- Rn_bb/Rre_bb-1 #in Hunt et al 2011
  
  #PBI
  SI_PBI <- R810_10/R560_10
  
  #TGI
  SI_TGI_cb <- -0.5*(190*(Rr_cb-Rg_cb)-120*(Rr_cb-Rb_cb))
  SI_TGI_bb <- -0.5*(190*(Rr_bb-Rg_bb)-120*(Rr_bb-Rb_bb))
  SI_TGI_nb <- -0.5*(190*(R670_10-R550_ASD)-120*(R670_10-R480_ASD))
  
  #TCARI
  SI_TCARI <- 3*((R700_CASI-R670_CASI)-0.2*(R700_CASI-R550_CASI)*(R700_CASI/R670_CASI))
  
  #MCARI revised
  SI_MCARI_rev <- ((R750_HYP-R705_HYP)-0.2*(R750_HYP-R550_HYP))/(R750_HYP/R705)
  
  #MSR revised
  SI_MSR_rev <- ((R750_HYP/R705_HYP)-1)/sqrt((R750_HYP/R705_HYP)+1)
  
  #TCARI/OSAVI
  SI_TCARItemp <- 3*((R700_CASI*f-R670_CASI*f)-0.2*(R700_CASI*f-R550_CASI*f)*(R700_CASI*f/(R670_CASI*f)))
  SI_OSAVItemp <- (1+0.16)*(R800_CASI*f-R670_CASI*f)/(R800_CASI*f+R670_CASI*f+0.16)
  SI_TCARI_OSAVI <- SI_TCARItemp/SI_OSAVItemp
  
  #TCARI/OSAVI revised
  TCARItemp <- 3*((R750_HYP-R705_HYP)-0.2*(R750_HYP-R550_HYP)*(R750_HYP/R705_HYP))
  OSAVItemp <- (1+0.16)*(R750_HYP-R705_HYP)/(R750_HYP+R705_HYP+0.16)
  SI_TCARI_OSAVI_rev <- TCARItemp/OSAVItemp
  
  #MCARI/OSAVI revised
  SI_MCARI_OSAVI_rev <- SI_MCARI_rev/OSAVItemp
  
  #MTCI
  SI_MTCI <- (R750_MERIS-R710_MERIS)/(R710_MERIS-R680_MERIS)
  
  #OCAR
  SI_OCAR <- R630/R680
  
  #YCAR
  SI_YCAR <- R600/R680
  
  #MCARI
  SI_MCARI <- ((R700_SE3-R670_SE3)-0.2*(R700_SE3-R550_SE3))*(R700_SE3/R670_SE3)
  
  #TVI
  SI_TVI <- 0.5*(120*(R750_CASI-R550_CASI)-200*(R670_CASI-R550_CASI))
  
  #REM
  SI_REM <- (Rn_MERIS/Rre_MERIS)-1
  #SI_REM2 <- R750/R720-1 #in Chen et al 2010
  
  #GM
  SI_GM <- Rn_MODIS2/Rg_MODIS4-1
  
  #VIopt
  SI_VIopt <- (1+0.45)*((Rn_bb)*2+1)/(Rr_bb+0.45)
  #SI_VIopt2 <- (1+0.45)*((R800)*2+1)/(R670+0.45) #in Chen et al 2010
  
  #RVI1
  SI_RVI1 <- R810_MSR/R660_MSR
  
  #RVI2
  SI_RVI2 <- R810_MSR/R660_MSR
  
  #MCARI/MTVI2
  MCARItemp1 <- (R700-R670-0.2*(R700-R550))*(R700-R670)
  MTVI2temp1 <- (1.5*(1.2*(R800-R550)-2.5*(R670-R550)))/sqrt((2*R800+1)^2-(6*R800-5*sqrt(R670))-0.5)
  MCARItemp2 <- (Rre_RE4-Rr_RE3-0.2*(Rre_RE4-Rg_RE2))*(Rre_RE4-Rr_RE3)
  MTVI2temp2 <- (1.5*(1.2*(Rn_RE5-Rg_RE2)-2.5*(Rr_RE3-Rg_RE2)))/sqrt((2*Rn_RE5+1)^2-(6*R800-5*sqrt(Rr_RE3))-0.5)
  SI_MCARI_MTVI2_ASD <- MCARItemp1/MTVI2temp1
  SI_MCARI_MTVI2_RE <- MCARItemp2/MTVI2temp2
  
  #NDNI
  SI_NDNI <- (log(1/R1510_AVIRIS)-log(1/R1680_AVIRIS))/(log(1/R1510_AVIRIS)+log(1/R1680_AVIRIS))
  
  #DCNI
  n <- 0.03 #soil constant
  SI_DCNI_CASI <- (R720_CASI-R700_CASI)/(R700_CASI-R670_CASI)/(R720_CASI-R670_CASI+n)
  SI_DCNI_ASD <- (R720_ASD-R700_ASD)/(R700_ASD-R670_3)/(R720_ASD-R670_3+n)
  
  #GBNDVI
  SI_GBNDVI <- (R573_SD-R440_SD)/(R573_SD+R440_SD)
  
  #SRWI
  SI_SRWI <- Rn_MODIS2/MODIS5
  
  #WI-WBI
  SI_WI <- R900_UNI/R970_UNI
  
  #WI/NDVI
  SI_WI_NDVI <- (R900_UNI/R970_UNI)/((R800_UNI-R680_UNI)/(R800_UNI+R680_UNI))
  
  #NDWI
  SI_NDWI1 <- (R860_AVIRIS-R1240_AVIRIS)/(R860_AVIRIS+R1240_AVIRIS)
  
  #Ratio at 975
  SI_975 <- (2*R960to990)/(R920to940+R1090to1110)
  
  #Ratio at 1200
  SI_1200 <- (2*R1180to1220)/(R1090to1110+R1265to1285)
  
  #NDMI
  SI_NDMI <- (Rn_TM4-TM5)/(Rn_TM4+TM5)
  
  #LWVI1
  SI_LWVI1 <- (R1094_HYP-R893_HYP)/(R1094_HYP+R893_HYP)
  
  #LWVI2
  SI_LWVI2 <- (R1094_HYP-R1205_HYP)/(R1094_HYP+R1205_HYP)
  
  #LWI
  SI_LWI <- R1300_100/R1450_100
  
  #SIWSI-NDWI1640
  SI_SIWSI <- (Rn_MODIS2-MODIS6)/(Rn_MODIS2+MODIS6)
  
  #NDWI1650
  SI_NDWI1650 <- (R840_CS-R1650_CS)/(R840_CS+R1650_CS)
  
  #NDWI2130
  SI_NDWI2130 <- (Rn_MODIS2-MODIS7)/(Rn_MODIS2+MODIS7)
  
  #PRI
  SI_PRI_10 <- (R531_10-R570_10)/(R531_10+R570_10) #SE590 (10nm)
  SI_PRI_6 <- (R531_6-R570_6)/(R531_6+R570_6) #SE590 (6nm)
  
  #DSWI
  SI_DSWI <- (R800_HYP+R550_HYP)/(R1659_HYP+R681_HYP)
  
  #HI
  SI_HI <- (R534_ASD-R698_ASD)/(R534_ASD+R698_ASD)-R704_ASD/2
  
  #CLSI
  SI_CLSI <- (R698_ASD-R570_ASD)/(R698_ASD+R570_ASD)-R734_ASD
  
  #SBRI
  SI_SBRI <- (R570_ASD-R513_ASD)/(R570_ASD+R513_ASD)-R704_ASD/2
  
  #PMI
  SI_PMI <- (R520_ASD-R584_ASD)/(R520_ASD+R584_ASD)+R724_ASD
  
  #NHI
  SI_NHI_ASD <- (R1100_5-R1200_5)/(R1100_5+R1200_5)
  SI_NHI_ALI <- (ALI6-ALI7)/(ALI6+ALI7)
  
  #Corrected NHI
  SI_CNHI_ASD <- ((R1100_5-R1200_5)/(R1100_5+R1200_5))/((R850_5-R670_5)/(R850_5+R670_5))
  SI_CNHI_ALI <- ((ALI6-ALI7)*(ALI5+Rr_ALI4))/((ALI6+ALI7)*(ALI5-Rr_ALI4))
  
  #FII
  SI_FII <- (R470to480-R360to370)/(R470to480+R360to370)
  
  #NDSVI
  SI_NDSVI <- (TM5-Rr_TM3)/(TM5+Rr_TM3)
  #SI_NDSVI2 <- (ASTER4-ASTER2/ASTER4+ASTER2) #in Pena-Barragan et al 2011
  
  #ALBEDO
  
  # ALBEDO <- (A_vis+A_NIR+A_SWIR1+A_SWIR2+A_SWIR3)/5
  ALBEDO_SWIR <- (A_SWIR1+A_SWIR2+A_SWIR3)/5
  # ALBEDO_VIS <- A_vis*1
  ALBEDO_NIR <- A_NIR*1
  
  ######## SENESCENCE RELEVANT INDICES
  
  #Simple Ratio Index 
  SI_SR <- R800/R670
  
  #Vogelmann Red Edge Index 1
  SI_VOG1 <- R740_ASD/R720_ASD
  
  #Vogelmann Red Edge Index 2
  SI_VOG2 <- (R734_ASD - R747_ASD)/(R715_ASD + R726_ASD)
  
  #Vogelmann Red Edge Index 3
  SI_VOG3 <- (R734_ASD - R747_ASD)/(R715_ASD + R720_ASD)
  
  #R550
  SI_R550 <- R550_ASD
  
  #Gnyli
  SI_GNYLI <- (R900*R1050_ASD-R955*R1220)/(R900*R1050_ASD+R955*R1220)
  
  #NRI
  SI_NRI <- (R874-R1225)/(R874+R1225)
  
  #chormatographic indices (provided by HA)
  #ExG
  SI_Exg_HA <- 2*G-R-B
  
  #NGRDI
  SI_NGRDI_HA <- (G-R)/(G+R)
  
  #GCC
  SI_GCC_HA <- G/(R + G + B)
  
  # CREATE TIBBLE ########################################################################## -
  
  DF <- do.call(cbind.data.frame, mget(ls(pattern = "SI_"))) %>% tibble::as_tibble()
  
  # convert to list of data.tables
  svi <- data.table::as.data.table(DF)
  svi_pre_list <- map(purrr::transpose(svi), data.table::as.data.table)
  
  colname_out <- "SVI"
  
  # replace raw spectra with pre-processed spectra
  if(new_col == F){
    SVI <- data[,-grep(col_in, colnames(data))]
  } else {
    SVI <- data
  }
  SVI[, colname_out] <- list(svi_pre_list)
  
  return(SVI)
  
}

#' Scales spectral indices to range from 0 to 10 
#' representing the minimum and maximum value observed in a time series of measurements
#' @param data A tibble with SVI values in a list column of data.tables, as returned by "calculate_SVI"
#' @param plotid The variable name of the plot identifier
#' @return A tibble with list columns
#' @details This function reverts the scale for SVI with an increase during the measurement period
scale_SVI <- function(data, plotid = "Plot_ID", plot = T, topdf = F) {
  
  # fix plot identifier 
  colid <- which(names(data)==plotid)
  names(data)[colid] <- "Plot_ID"
  
  # keep only plots with measurements covering the entire measurement period
  # get start and end date
  min_date <- min(data$meas_date)
  max_date <- max(data$meas_date)
  # transform measurement date in days after first measurement
  # this should be replaced by a measure of chronological or thermal time after heading
  data$dafm <- as.numeric(as.Date(data$meas_date, "%Y%m%d") - as.Date(min_date, "%Y%m%d"))
  
  # extract SVI from list column and add required metadata
  meta <- data[c(plotid, "meas_date", "dafm", "Treatment")]
  SVI_dat <- data.table::rbindlist(data[["SVI"]])
  d <- cbind(meta, SVI_dat)
  
  # filter dataset
  d_ <- d %>% group_by(Plot_ID) %>% nest() %>% 
    mutate(max = purrr::map_chr(data,  ~max(.$meas_date)),
           min = purrr::map_chr(data,  ~min(.$meas_date))) %>% 
    filter(max == max_date & min == min_date)
  # remove helper columns
  d_ <- d_[!(names(d_) %in% c("max", "min"))] 
  
  # scale SVI
  d_scaled <- d_ %>% 
    transmute(SVI_sc = purrr::map(data, col_scaling)) %>% 
    unnest(cols = c(SVI_sc))
  
  # revert scale where required
  ids <- d_scaled[c("Plot_ID", "meas_date", "dafm", "Treatment")] %>% ungroup()
  SVI_sc_dat <- d_scaled[grepl("^SI_", names(d_scaled))]
  SVI_sc_dat <- SVI_sc_dat %>% 
    mutate_all(funs(r = revert)) %>% 
    #select original or reversed values
    dplyr::select_if(function(col) col[1] > 5) %>% 
    data.table::as.data.table()
  
  if(plot){
    
    # transform measurement date in days after first measurement
    # this should be replaced by a measure of chronological or thermal time after heading
    pd <- cbind(ids, SVI_sc_dat)
    pd$dafm <- as.numeric(as.Date(pd$meas_date, "%Y%m%d") - as.Date(min_date, "%Y%m%d"))
    # reshape data for plotting
    pd <- pd %>% dplyr::select(1:Treatment, dafm, starts_with("SI_")) %>% 
      tidyr::gather(SVI, value, starts_with("SI_"))
  
    p <- ggplot(pd) +
      geom_line(aes(x = dafm, y = value, group = Plot_ID, col = Treatment), alpha  = 0.3) +
      facet_wrap(~SVI) +
      ggsci::scale_color_npg() +
      theme_bw(base_size = 7)
    
    # check if Output directory exists
    if(!file.exists(paste0(path_to_data, "Output"))){
      dir.create(paste0(path_to_data, "Output"))
    }
    
    if(topdf){
      # save plot to pdf
      pdf(paste0(path_to_data, "Output/SVI_ts.pdf"), width = 35, height = 35)
      plot(p)
      dev.off()
    } else {
      plot(p)
    }

  }
  
  # create data output (list column)
  spc_pre_list <- map(purrr::transpose(SVI_sc_dat), data.table::as.data.table)
  ids[, "SVI_sc"] <- list(spc_pre_list)
  
  # join with input spectral data
  data_out <- full_join(data, ids, by = c("Plot_ID", "meas_date", "dafm", "Treatment"))
  
} 

#' Gets dynamics parameters for each index and plot
#' @param data A tibble with scaled SVI values in a list column of data.tables, as returned by "scale_SVI"
#' @param timevar A character string specifying the name of the column containing the time variable
#' @param method A character string specifying the method to be used for modelling of the temporal index dynamics. 
#' So far, only "interpolate" is supported.
#' @param plot Boolean, indicating whether or not to create a plot
#' @param topdf Boolean, indicating whether or not to save the plot to pdf
#' @return A tibble containing the dynamics parameters for each Plot and SVI. 
get_svi_dynamics <- function(data, timevar, method = "interpolate",
                             plot = T, topdf = F){
  
  # unnest 
  if(!is.null(data[["SVI_sc"]])){
    dat_svi <- data.table::rbindlist(data[["SVI_sc"]])
  } else {
    stop("Dynamics parameters can only be extracted from scaled SVI values!")
  }
  
  # fix variable names
  names(data)[which(names(data)==timevar)] <- "timevar"
  
  # reshape
  meta <- data[, names(data) %in% c("Plot_ID", "Treatment", "timevar")]
  dat <- cbind(meta, dat_svi)
  dat_long <- melt(dat, 
                   id.vars = c("Plot_ID", "Treatment", "timevar"),
                   measure.vars = c(grep("^SI_", names(dat), value = T)))
  
  # fit parametric model or perform linear interpolation
  dat_mod <- dat_long %>%
    dplyr::group_by(Plot_ID, variable, Treatment) %>%
    tidyr::nest()
  if(method == "interpolate"){
    fits <- dat_mod %>%
      mutate(fit = purrr::map(data, lin_approx, n_meas = 6) %>% purrr::map(cbind.data.frame)) %>%
      transmute(pars = purrr::map(fit, extract_pars) %>% purrr::map(cbind.data.frame)) %>%
      unnest(pars) %>% 
      # to (full) long
      tidyr::gather(par, value, t80:dur2) %>% 
      mutate(param = paste(variable, par, sep = "-")) %>% ungroup()  %>%  dplyr::select(-variable, -par) %>% 
      tidyr::spread(param, value)
  } else {
    stop("Only linear interpolation between measurement timepoints is implemented so far. Specify by setting method = interpolate")
  }

  if(plot){

    dur <- fits %>% dplyr::select(Plot_ID, Treatment, variable, dur1, dur2) %>%
      tidyr::gather(param, value, dur1:dur2)

    tps <- fits %>% dplyr::select(Plot_ID, Treatment, variable, t80:t20) %>%
      tidyr::gather(param, value, t80:t20)

    plot1 <- ggplot(dur) +
      geom_boxplot(aes(x = param, y = value, group = interaction(param, Treatment), fill = Treatment)) +
      facet_wrap(~variable) +
      ggsci::scale_fill_npg() +
      theme_bw(base_size = 7) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank())

    plot2 <- ggplot(tps) +
      geom_boxplot(aes(x = param, y = value, group = interaction(param, Treatment), fill = Treatment)) +
      facet_wrap(~variable) +
      ggsci::scale_fill_npg() +
      theme_bw(base_size = 7) +
      theme(panel.grid = element_blank(),
            panel.background = element_blank())

    plot(plot1)
    plot(plot2)

  }

  return(fits)

}

#' Plots reflectance spectra and pre-processing products
#' @param data A tibble with spectra in list columns of data.tables, as returned by "load_spectra()" and "preprocess_spc" 
#' @param col_in Column name of the tibble containing the spectra to use as character string
#' @param treatment Character string specifying the name of the variable indicating the treatments
#' @param xlim numeric vector, specifying the range along the x-axis to plot, defaults to c(350, 2500)
#' @param facets Character vector specifying the variable to use as facetting variable, passed to ggplot::facet_wrap()
#' @param mark_outliers Boolean, whether or not to mark multivariate outliers
#' @param topdf Boolean, whether or not to save the plot as a pdf
#' @return A ggplot object
plot_spectra <- function(data, col_in = "all", treatment = NULL,
                         xlim = c(350, 2500),
                         facets = NULL,
                         mark_outliers = FALSE,
                         topdf = FALSE){
  
  # print("plotting spectra ...")
  
  if("all" %in% col_in){
    col_in <- grep("^rflt_", names(data), value = TRUE)
  }
  
  # factorise and rename structuring variables
  if(!is.null(treatment)){
    data[[treatment]] <- as.factor(data[[treatment]])
    data <- dplyr::rename(data, "Treatment" = treatment)
  }
  if(!is.null(facets)){
    data[[facets]] <- as.factor(data[[facets]])
    data <- dplyr::rename(data, "group" = facets)
  } else{
    group <- NULL
  }
  
  # iterate over data types
  for (column in col_in) {
    
    # reshape data for plotting
    pd <- data %>% dplyr::select(meas_id, group, Treatment, !!as.symbol(column)) %>% 
      unnest(c(!!as.symbol(column))) %>% 
      tidyr::gather(wvlt, rflt, grep("^[0-9]", names(.)), factor_key = TRUE) %>% 
      mutate(wvlt = as.numeric(as.character(wvlt))) %>% 
      add_spc_range()
    
    otype = paste0("out_", column)
    
    ylim1 <-max(pd[pd$wvlt %in% seq(xlim[1], xlim[2], 1),]$rflt)
    ylim2 <- min(pd[pd$wvlt %in% seq(xlim[1], xlim[2], 1),]$rflt)
    ylim <- c(ylim2, ylim1)
    
    # create empty plot
    plot <- ggplot() + 
      xlab("Wavelength (nm)") + ylab(col_in) +
      xlim(xlim) + ylim(ylim) +
      theme_bw() + 
      theme(panel.grid = element_blank())
    
    # add spectra
    if(!is.null(treatment)){
      
      plotf <- plot +
        geom_line(data = pd, aes(x = wvlt, y = rflt,  group = interaction(meas_id, spc_range),
                                 col = Treatment), size = 0.05, alpha = 0.25) +
        ggsci::scale_color_npg()
      
      # calculate mean and sd of reflectance
      if(!is.null(facets)){
        pdmeans <- pd %>% group_by(Treatment, group, wvlt, spc_range)
      } else {
        pdmeans <- pd %>% group_by(Treatment, wvlt, spc_range)
      }
      pdmeans <- pdmeans %>% 
        dplyr::summarise(mean_rflt = mean(rflt, na.rm = T),
                         sd_rflt = sd(rflt, na.rm = T))
      
      # plot mean and sd per treatment and group
      plotmeans <- plot +
        geom_line(data = pdmeans, aes(x = wvlt, y = mean_rflt,  group = interaction(Treatment, spc_range),
                                      col = Treatment), size = 1, alpha = 1) +
        geom_ribbon(data=pdmeans, aes(x = wvlt, ymin = mean_rflt - sd_rflt, 
                                      ymax = mean_rflt + sd_rflt, 
                                      fill = Treatment, 
                                      group = interaction(Treatment, spc_range)), alpha = 0.25) +
        ggsci::scale_color_npg() +
        ggsci::scale_fill_npg()
    } else {
      plotf <- plot + 
        geom_line(data = pd, aes(x = wvlt, y = rflt,  group = interaction(meas_id, spc_range)), size = 0.05, alpha = 0.25)
    }
    
    # add facets
    if(!is.null(facets)){
      plotf <- plotf +
        facet_wrap(~group, ncol = 1, scales = "free")
      if(!is.null(treatment)){
        plotmeans <- plotmeans +
          facet_wrap(~group, ncol = 1, scales = "free")
      }
    }
    
    # mark outliers
    if(mark_outliers){
      
      if(otype %in% names(data)){
        
        # print("marking outliers ...")
        
        # reshape data for plotting
        pd_out <- data %>% dplyr::select(meas_id, group, treatment, column, otype) %>% unnest(c(column)) %>% 
          tidyr::gather(wvlt, rflt, grep("^[0-9]", names(.)), factor_key = TRUE) %>% 
          mutate(wvlt = as.numeric(as.character(wvlt))) %>% 
          add_spc_range()
        
        # mark outliers
        plotf <- plotf +
          geom_line(data = pd_out[pd_out[otype] == 1,], aes(x = wvlt, y = rflt, group = interaction(meas_id, spc_range)),
                    size = 0.05, alpha = 1, col = "black", lty = "dashed")
        
      } else {
        print("no data for outliers found!")
      }
      
    }
    
    # check if Output directory exists
    if(!file.exists(paste0(path_to_data, "Output"))){
      dir.create(paste0(path_to_data, "Output"))
    }
    
    if(topdf){
      # save plot to pdf
      # individual spectra
      levels = nrow(unique(pd["group"]))
      width = 10
      pdf(paste0(path_to_data, "Output/spectra_trt_", column,".pdf"), width = width, height = width * (levels/2))
      # pdf(paste0("C:/Users/anjonas/output/spectra_", column,".pdf"), width = width, height = width * (levels/2))
      plot(plotf)
      dev.off()
      # spectra means
      if(exists("plotmeans")){
        pdf(paste0(path_to_data, "Output/spectra_means_", column,".pdf"), width = width, height = width * (levels/2))
        plot(plotmeans)
        dev.off()
      }
    } else{
      plot(plotf)
      if(exists("plotmeans")){
        plot(plotmeans)
      }
    }
    
  }
  
}

#' Creates plots of spectral index values
#' @param data A tibble with original and/or scaled SVI values in list column(s) of data.tables
#' @param svi A vector of character strings specifying the name(s) of indices to show
#' @param col_in A character string specifying which type of Index values to use
#' @param x A character string specifying the x-variable
#' @param x_is_date Boolean, whether or not the x variable is a date
#' @param groups A character string specifying the grouping structure of the data
#' @param topdf Boolean, whether or not to save the plot to pdf
#' @return A ggplot
plot_SVI <- function(data, 
                     svi, col_in,
                     x, x_is_date = FALSE, 
                     groups = NULL,
                     topdf = FALSE){
  
  # for variable selection
  svis <- paste0("_", svi)
  pattern <- paste0(svis, collapse = "|")
  # unnest SVI data
  dat_svi <- data.table::rbindlist(data[[col_in]])
  namevec <- grep(pattern = pattern, names(dat_svi), value = T)
  dat_svi <- dat_svi[, ..namevec]
  # extract meta data
  meta <- data[, !grepl("^rflt|^SVI", names(SVI))]
  dat <- cbind(meta, dat_svi)
  # reshape for plotting
  dat_long <- melt(dat, 
                   id.vars = c("Plot_ID", "meas_date", "Treatment"),
                   measure.vars = c(grep("^SI_", names(dat), value = T)))
  
  # fix variable names
  if(!is.null(groups)){
    names(dat_long)[which(names(dat_long) == groups)] <- "groups"
  }
  names(dat_long)[which(names(dat_long) == x)] <- "xvar"
  
  if(x_is_date){
    dat_long$xvar <- as.Date(dat_long$xvar, format = "%Y%m%d")
  }
  
  # create plot of index values
  p <- ggplot() +
    xlab(x) + ylab("Index value") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  if(is.null(x)){
    p <- p + geom_boxplot(data = dat_long, aes(y = value))
  } else {
    if(!is.null(groups)){
      p <- p + geom_boxplot(data = dat_long, aes(x = xvar, y = value, group = interaction(xvar, groups), fill = groups)) +
        ggsci::scale_fill_npg() +
        guides(fill = guide_legend(title = groups))
      
    } else {
      p <- p + geom_boxplot(data = dat_long, aes(x = xvar, y = value, group = interaction(xvar)))
    }
  }
  
  p <- p + facet_wrap(~dat_long$variable, scales = "free")
  
  
  if(topdf){
    # check if Output directory exists
    if(!file.exists(paste0(path_to_data, "Output"))){
      dir.create(paste0(path_to_data, "Output"))
    }
    pdf(paste0(path_to_data, "Output/SVI.pdf"))
    plot(p)
    dev.off()
  } else {
    plot(p)
  }
  
}

spectra_cor <- function(data, trait, col_in, topdf = F){
  
  # get data
  cordat <- data[c(trait, col_in)] %>% unnest(col_in) %>% as.data.frame()
  # trait data
  traitvector <- as.matrix(as.numeric(cordat[,trait]))
  # variables with which to correlate
  variables <- as.matrix(cordat[, -which(names(cordat) == trait)])
  # calculate correlation coefficient
  cor_coef <- -(as.vector(cor(traitvector, variables, method = "spearman")))
  wvlt <- as.numeric(gsub("_", "", stringr::str_sub(names(cordat)[-c(1)], -4, -1)))
  coefs <- as.data.frame(cbind(wvlt, cor_coef))
  # add spectral range for plotting
  coefs <- add_spc_range(coefs)
  # create plot
  p1 <- ggplot(coefs) + 
    geom_line(aes(x = wvlt, y = cor_coef, group = spc_range), size = 0.6) + 
    xlab(col_in) + ylab("Correlation coefficient") +
    scale_x_continuous(breaks = seq(0,2500,500), limits = c(350, 2550), expand = c(0.01, 0.01)) +
    geom_abline(slope = 0, intercept = 0) +
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
  
  plot(p1)
  
  if(topdf){
    pdf(paste0(path_to_data, "Output/spectra_corr_", col_in,".pdf"))
    plot(p1)
    dev.off()
  }
  
  return(coefs)
  
}


# ============================================================================================================= -

# helper functions ---- 

#' Read an individual spectrum from an .sed or an .asd file
#' @param dir full file name
#' @param format file extension as a character string, "sed" or "asd"
#' @return A tibble with measurement_id, measurement date and the spectrum in a data.table
read_spectrum <- function(dir, subdir = NULL, format = "sed"){
  if(format == "sed"){
    # read file 
    spc <- data.table::fread(dir, skip = 26)
    # transpose to wide
    spct <- dcast(melt(spc, id.vars = "Wvl"), variable ~ Wvl)
    spct[, "variable":=NULL]
    spct <- spct/100
    colnames(spct) <- as.character(spc$Wvl)
    # add measurement name and measurement date as identifiers 
    filename <- base::basename(dir)
    l <- str_split(dir, "/") %>% unlist()
    meas_date <- l[length(l)-1]
    # Return spectra as tibble
    t <- tibble::tibble(
      meas_id = filename,
      meas_date = meas_date,
      rflt = list(spct)
    )
  } else if (format == "asd"){
    spc <- data.table(prospectr::readASD(dir))
    # add measurement name and measurement date as identifiers 
    filename <- base::basename(dir)
    l <- str_split(dir, "/") %>% unlist()
    meas_date <- grep("[0-9]{8}$", l, value = TRUE)
    t <- tibble::tibble(
      meas_id = filename,
      meas_date = meas_date,
      rflt = list(spc)
    )
  } else{
    print("data type not recognized!")
  }
}

get_regions_drop <- function(limits){
  # bands to drop
  wb1 <- seq(limits[1], limits[2], 1)
  wb2 <- seq(limits[3], limits[4], 1)
  fb <- seq(limits[5], limits[6], 1)
  # corresponding variable names
  drop_vars <- as.character(c(wb1, wb2, fb))
  return(drop_vars)
}  

add_spc_range <- function(data){
  data <- data %>% 
    mutate(spc_range = ifelse(data$wvlt <= 1350, "visnir", 
                              ifelse(data$wvlt <= 1990, "swir1",
                                     "swir2")))
  
  return(data)
  
} 

plot_outliers <- function(data, col_in, 
                          mark_outliers = TRUE, 
                          facets = "meas_date",
                          topdf = FALSE){
  
  # print("plotting spectra ...")
  
  otype <- paste0("out_", col_in)
  
  # reshape data for plotting
  pd <- data %>% dplyr::select(meas_id, meas_date, col_in, otype) %>% unnest(c(col_in)) %>% 
    tidyr::gather(wvlt, rflt, grep("^[0-9]", names(.)), factor_key = TRUE) %>% 
    mutate(wvlt = as.numeric(as.character(wvlt))) %>% 
    add_spc_range()
  
  plot <- ggplot() + 
    geom_line(data = pd, aes(x = wvlt, y = rflt,  group = interaction(meas_id, spc_range)), size = 0.05, alpha = 0.25) + 
    # geom_abline(intercept = 0, slope = 0) +
    # facet_wrap(as.formula(paste0("~", facets)), ncol = 1, scales = "free") +
    theme_bw() + 
    theme(panel.grid = element_blank())
  
  if(!is.null(facets)){
    plot <- plot +
      facet_wrap(as.formula(paste0("~", facets)), ncol = 1, scales = "free")
  }
  
  if(mark_outliers){
    plot <- plot +
      geom_line(data = pd[pd[otype] == 1,], aes(x = wvlt, y = rflt, group = interaction(meas_id, spc_range)), 
                size = 0.05, alpha = 1, col = "darkred")
  }
  
  # check if Output directory exists
  if(!file.exists(paste0(path_to_data, "Output"))){
    dir.create(paste0(path_to_data, "Output"))
  }
  
  # save plot to pdf
  if(!is.null(facets)){
    levels <- nrow(unique(pd[facets]))
  } else {
    levels <- 1
  }
  
  if(topdf){
    width = 10
    pdf(paste0(path_to_data, "Output/spectra_outs_", col_in,".pdf"), width = width, height = width * (levels/2))
    plot(plot)
    dev.off()
  } else {
    plot(plot)
  }

}

revert <- function(x){10 - x}

col_scaling <- function(d) {
  ids <- d[!grepl("^SI_", names(d))]
  d <- d[grep("^SI_", names(d))]
  d <- purrr::map_df(d, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE))*10)
  out <- bind_cols(ids, d)
  return(out)
}

lin_approx <- function(data, n_meas){
  
  data <- as.data.frame(data)
  
  # linearly interpolate between measurement time points
  out <- approx(data[, "timevar"], data[,"value"],
                xout = seq(round(min(data[, "timevar"], na.rm = TRUE), 0),
                           round(max(data[, "timevar"], na.rm = TRUE), 0), 1))
  names(out) <- c("timevar", ".fitted")
  
  # check that time series is complete
  # and first measurement is larger than the critical value for onset
  n <- nrow(data)
  init <- data[1,2]
  if(n < n_meas || init <= 8){
    out$.fitted <- rep(NA, length(out["timevar"]))
  }
  return(out)
}

extract_pars <- function(data){
  t80 <- data[which(data[".fitted"] < 8)[1], "timevar"]
  t50 <- data[which(data[".fitted"] < 5)[1], "timevar"]
  t20 <- data[which(data[".fitted"] < 2)[1], "timevar"]
  dur1 <- t50 - t80
  dur2 <- t20 - t80
  pars <- cbind(t80, t50, t20, dur1, dur2)
  return(pars)
}

# ============================================================================================================= -

#' Perform recursive feature elimination
#' @param response A character string indicating the name of the variable to predict
#' @param base_learner A character string indicating the algorithm to use as base learner (so far either "ranger" or "cubist")
#' @param type A character string, either "regression" or "classification"
#' @param p A numeric in the range (0,1) indicating the split proportions (the larger fraction is used for feature elimination)
#' @param times A numeric, indicating the number of times the dataset is to be resampled
#' @param groups A numeric, indicating the number of groups for stratified sampling (ensures balanced evaluation datasets)
#' @param subsets A numeric vector, defining the number of variables to be used for modelling
#' @param data Dataset
#' @param ... Other arguments passed to caret::train, depending on the base learner used
#' @return A list of length length(subsets), hwere each list element holds the output for one resample of the data. 
#' Each list element is a list of 4, 
#' holding (i) feature ranks, (ii) cross-validate training RMSE, (iii) test RMSE, (iv) the corresponding number of features.
perform_rfe <- function(response, base_learner = "ranger", type = "regression",
                        p = 0.75, times = 30, groups = 9, 
                        subsets, data,
                        ...) {
  
  #create multifolds for repeated n-fold cross validation
  index <- caret::createDataPartition(pull(data[response]), p = p, times = times, groups = ifelse(is.numeric(groups), groups, 2))
  
  #outer resampling
  #CV of feature selection
  out <- list()
  for(i in 1:length(index)){
    
    #Verbose
    print(paste("resample ", i, "/", length(index), sep = ""))
    
    #use indices to create train and test data sets for the resample
    ind <- as.numeric(index[[i]])
    train <- data[ind,]
    test <- data[-ind, ]
    
    #for each subset of decreasing size
    #tune/train rf and select variables to retain
    keep_vars <- drop_vars <- test_perf <- train_perf <- npred <- NULL
    for(j in 1:length(subsets)){
      
      #define new training data
      #except for first iteration, where the full data set ist used
      if(exists("newtrain")) {train = newtrain}
      
      #Verbose iter
      print(paste("==> subset size = ", length(train)-1, sep = ""))
      
      #define tune grid
      if(base_learner == "ranger"){
        #adjust mtry parameter to decreasing predictor set
        #maximum mtry at 200
        mtry <- ceiling(seq(1, length(train[-1]), len = 7)) %>% unique()
        if(any(mtry > 250)){
          mtry <- mtry[-which(mtry >= 250)]
        }
        min.node.size <- c(5)
        tune_grid <- expand.grid(mtry = mtry,
                                 splitrule = ifelse(type == "regression", "variance", "gini"),
                                 min.node.size = ifelse(type == "regression", 5, 1)) 
      } else if(base_learner == "cubist"){
        tune_grid <- expand.grid(committees = c(1, 2, 5, 10),
                                 neighbors = c(0))
      }
      
      #define inner resampling procedure
      ctrl <- caret::trainControl(method = "repeatedcv",
                                  number = 10,
                                  rep = 1,
                                  verbose = FALSE,
                                  allowParallel = TRUE,
                                  savePredictions = TRUE,
                                  classProbs = ifelse(type == "classification", TRUE, FALSE))
      
      #define model to fit
      formula <- as.formula(paste(response, " ~ .", sep = ""))
      
      #tune/train random forest
      fit <- caret::train(formula,
                          data = train,
                          preProc = c("center", "scale"),
                          method = base_learner,
                          tuneGrid = tune_grid,
                          trControl = ctrl,
                          ...)
      
      if(type == "regression"){
        #extract predobs of each cv fold
        predobs_cv <- plyr::match_df(fit$pred, fit$bestTune, on = names(fit$bestTune))
        #Average predictions of the held out samples;
        predobs <- predobs_cv %>% 
          group_by(rowIndex) %>% 
          dplyr::summarize(obs = mean(obs),
                           mean_pred = mean(pred))
        #get train performance
        train_perf[j] <- caret::getTrainPerf(fit)$TrainRMSE
        #get test performance
        test_perf[j] <- rmse(test %>% pull(response), caret::predict.train(fit, test))
      } else if (type == "classification"){
        #get train accuracy
        train_perf[j] <- caret::getTrainPerf(fit)$TrainAccuracy
        #get test accuracy
        test_perf[j] <- get_acc(fit, test)
      }
      
      #number of preds used
      npred[[j]] <- length(train)-1
      
      #extract retained variables
      #assign ranks
      #define reduced training data set
      if(j < length(subsets)){
        #extract top variables to keep for next iteration
        keep_vars[[j]] <- varImp(fit)$importance %>% 
          tibble::rownames_to_column() %>% 
          as_tibble() %>% dplyr::rename(var = rowname) %>%
          arrange(desc(Overall)) %>% slice(1:subsets[j+1]) %>% pull(var)
        #extract variables dropped from dataset
        drop_vars[[j]] <- names(train)[!names(train) %in% c(keep_vars[[j]], response)] %>% 
          tibble::enframe() %>% mutate(rank = length(subsets)-j+1) %>% 
          dplyr::select(value, rank) %>% dplyr::rename(var = value)
        #define new training data
        newtrain <- dplyr::select(train, response, keep_vars[[j]])
        #last iteration
      } else {
        drop_vars[[j]] <- names(train)[names(train) != response] %>% 
          tibble::enframe() %>% mutate(rank = length(subsets)-j+1) %>% 
          dplyr::select(value, rank) %>% rename(var = value)
      }
    } #END OF FEATURE ELIMINATION ON RESAMPLE i
    #clean environment 
    rm("newtrain")
    #gather results for resample i
    ranks <- drop_vars %>% do.call("rbind", .)
    out[[i]] <- list(ranks, train_perf, test_perf, npred)
  } #END OF OUTER RESAMPLING
  return(out)
}

#Perform recursive feature elimination
perform_rfe_par <- function(response, base_learner = "ranger", type = "regression",
                            p = 0.75, times = 30, groups = 9, parallel = T, 
                            subsets, data,
                            ...) {
  
  #create multifolds for repeated n-fold cross validation
  index <- caret::createDataPartition(pull(data[response]), p = p, times = times, groups = ifelse(is.numeric(groups), groups, 2))
  
  #outer resampling
  #CV of feature selection
  `%infix%` <- ifelse(parallel, `%dopar%`, `%do%`)
  foreach(i=1:length(index)) %infix% {
    
    #Verbose
    print(paste("resample ", i, "/", length(index), sep = ""))
    
    #use indices to create train and test data sets for the resample
    ind <- as.numeric(index[[i]])
    train <- data[ind,]
    test <- data[-ind, ]
    
    #for each subset of decreasing size
    #tune/train rf and select variables to retain
    keep_vars <- drop_vars <- test_perf <- train_perf <- npred <- NULL
    for(j in 1:length(subsets)){
      
      #define new training data
      #except for first iteration, where the full data set ist used
      if(exists("newtrain")) {train = newtrain}
      
      #Verbose iter
      print(paste("==> subset size = ", length(train)-1, sep = ""))
      
      #define tune grid
      if(base_learner == "ranger"){
        #adjust mtry parameter to decreasing predictor set
        #maximum mtry at 200
        mtry <- ceiling(seq(1, length(train[-1]), len = 7)) %>% unique()
        if(any(mtry > 250)){
          mtry <- mtry[-which(mtry >= 250)]
        }
        min.node.size <- c(5)
        tune_grid <- expand.grid(mtry = mtry,
                                 splitrule = ifelse(type == "regression", "variance", "gini"),
                                 min.node.size = ifelse(type == "regression", 5, 1)) 
      } else if(base_learner == "cubist"){
        tune_grid <- expand.grid(committees = c(1, 2, 5, 10),
                                 neighbors = c(0))
      }
      
      #define inner resampling procedure
      ctrl <- caret::trainControl(method = "repeatedcv",
                                  number = 10,
                                  rep = 1,
                                  verbose = FALSE,
                                  allowParallel = TRUE,
                                  savePredictions = TRUE,
                                  classProbs = ifelse(type == "classification", TRUE, FALSE))
      
      #define model to fit
      formula <- as.formula(paste(response, " ~ .", sep = ""))
      
      #tune/train random forest
      fit <- caret::train(formula,
                          data = train,
                          preProc = c("center", "scale"),
                          method = base_learner,
                          tuneGrid = tune_grid,
                          trControl = ctrl,
                          ...)
      
      if(type == "regression"){
        #extract predobs of each cv fold
        predobs_cv <- plyr::match_df(fit$pred, fit$bestTune, on = names(fit$bestTune))
        #Average predictions of the held out samples;
        predobs <- predobs_cv %>% 
          group_by(rowIndex) %>% 
          dplyr::summarize(obs = mean(obs),
                           mean_pred = mean(pred))
        #get train performance
        train_perf[j] <- caret::getTrainPerf(fit)$TrainRMSE
        #get test performance
        test_perf[j] <- rmse(test %>% pull(response), caret::predict.train(fit, test))
      } else if (type == "classification"){
        #get train accuracy
        train_perf[j] <- caret::getTrainPerf(fit)$TrainAccuracy
        #get test accuracy
        test_perf[j] <- get_acc(fit, test)
      }
      
      #number of preds used
      npred[[j]] <- length(train)-1
      
      #extract retained variables
      #assign ranks
      #define reduced training data set
      if(j < length(subsets)){
        #extract top variables to keep for next iteration
        keep_vars[[j]] <- varImp(fit)$importance %>% 
          tibble::rownames_to_column() %>% 
          as_tibble() %>% dplyr::rename(var = rowname) %>%
          arrange(desc(Overall)) %>% slice(1:subsets[j+1]) %>% pull(var)
        #extract variables dropped from dataset
        drop_vars[[j]] <- names(train)[!names(train) %in% c(keep_vars[[j]], response)] %>% 
          tibble::enframe() %>% mutate(rank = length(subsets)-j+1) %>% 
          dplyr::select(value, rank) %>% dplyr::rename(var = value)
        #define new training data
        newtrain <- dplyr::select(train, response, keep_vars[[j]])
        #last iteration
      } else {
        drop_vars[[j]] <- names(train)[names(train) != response] %>% 
          tibble::enframe() %>% mutate(rank = length(subsets)-j+1) %>% 
          dplyr::select(value, rank) %>% rename(var = value)
      }
    } #END OF FEATURE ELIMINATION ON RESAMPLE i
    #clean environment 
    rm("newtrain")
    #gather results for resample i
    ranks <- drop_vars %>% do.call("rbind", .)
    return(list(ranks, train_perf, test_perf, npred))
  } #END OF OUTER RESAMPLING
}


#' Create a tidy output
#' @param data The output of perform_rfe, a list of length length(subsets), with each list element a list of 4.
#' @param base_learner A character string, indicating the base learner used.
#' @return A list holding the performance profile and the robust feature ranks 
#' (as averages and standard errors across resamples)
tidy_rfe_output <- function(data, base_learner){
  #tidy up list output
  subsets <- data[[1]][[length(data[[1]])]]
  ranks <- lapply(data, "[[", 1) %>% 
    Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = "var"), .) %>% 
    purrr::set_names(., c("var", paste("Resample", 1:length(data), sep = "")))
  RMSEtrain <- lapply(data, "[[", 2) %>% lapply(., cbind, subsets) %>%
    lapply(., as_tibble) %>% Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = "subsets"), .) %>% 
    dplyr::select(subsets, everything()) %>% 
    purrr::set_names(c("subset_size", paste("Resample", 1:length(data), sep = "")))
  RMSEtest <- lapply(data, "[[", 3) %>% lapply(., cbind, subsets) %>%
    lapply(., as_tibble) %>% Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = "subsets"), .) %>% 
    dplyr::select(subsets, everything()) %>% 
    purrr::set_names(c("subset_size", paste("Resample", 1:length(data), sep = "")))
  #average across resamples, get sd and means
  Trainperf <- RMSEtrain %>%
    gather(resample, RMSE, contains("Resample")) %>%
    group_by(subset_size) %>%
    arrange(subset_size) %>%
    summarise_at(vars(RMSE), funs(mean, sd), na.rm = TRUE) %>%
    mutate(set = "Train")
  Testperf <- RMSEtest %>% 
    gather(resample, RMSE, contains("Resample")) %>%
    group_by(subset_size) %>%
    arrange(subset_size) %>%
    summarise_at(vars(RMSE), funs(mean, sd), na.rm = TRUE) %>%
    mutate(set = "Test")
  Perf <- bind_rows(Trainperf, Testperf) %>% mutate(algorithm = base_learner)
  #average ranks
  robranks <- ranks %>% 
    gather(resample, rank, contains("Resample")) %>%
    group_by(var) %>%
    summarise_at(vars(rank), funs(mean, sd), na.rm = TRUE) %>%
    arrange(mean)
  tidy_out <- list(Perf, robranks)
  return(tidy_out)
}

#' Plot performance profiles
#' @param The first element of the tidy rfe output 
#' @return A ggplot
plot_perf_profile <- function(data){
  pd <- position_dodge(0.5) # move them .05 to the left and right
  #plot performance profiles
  ggplot(data, aes(x = subset_size, y = mean, group = set, colour = set)) +
    geom_point(position = pd) + geom_line() +
    geom_errorbar(position = pd, aes(ymin = mean - sd, ymax = mean + sd), width = 1, alpha = 0.5) + 
    xlab("#Features") + ylab("RMSE") +
    scale_x_continuous(limits = c(-2.5, 32.5)) +
    facet_wrap(~algorithm) +
    theme_bw() +
    theme(legend.title=element_blank(),
          plot.title = element_text(size = 15, face = "bold"),
          strip.text = element_text(face = "bold"))
}

plot_feature_ranks <- function(data, topdf = F){
  # order features
  ranks <- tidy[[2]] %>% tibble::rowid_to_column("order") %>% 
    mutate(order = as.numeric(1-order))
  # create plot
  ranks <- ggplot(ranks, aes(x=order, y=mean)) +
    geom_point() +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.5) +
    ylab("Feature rank") +
    xlab("Feature")+
    # Add categories to axis
    scale_x_continuous(
      breaks = ranks$order,
      labels = ranks$var,
      expand = c(0,0)
    ) +
    coord_flip() +
    theme_bw() %+replace%
    theme(axis.title.y =  element_blank(),
          plot.title = element_text(size=15, face="bold"),
          panel.grid.minor = element_blank())
  plot(ranks)
  # save to png
  if(topdf){
    png(paste0(path_to_data, "Output/feature_ranks.png"), width = 7, height = 8, units = 'in', res = 400)
    plot(ranks)
    dev.off()
  }
  
}


# ============================================================================================================= -

#Helper function to calculate accuracy
get_acc <- function(model, testdata) {
  preds_class <- caret::predict.train(model, newdata = testdata[ , names(testdata) != "trt"])
  true_class <- testdata$trt
  res <- cbind(preds_class, true_class) %>% data.frame()
  match <- ifelse(res$preds_class == res$true_class, 1, 0) %>% sum()
  acc <- match/nrow(testdata)
}

#Helper function to get RMSE
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}

# ============================================================================================================= -
# ============================================================================================================= -
# ============================================================================================================= -
