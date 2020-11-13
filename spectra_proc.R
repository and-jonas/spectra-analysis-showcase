
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
load_spectra <- function(dir, format = "sed"){
  print("loading spetra ...")
  # get all subdirectories
  subdirs <- dir(path = dir, full.names = TRUE, recursive = FALSE, pattern = "^[0-9]{8}")
  # get all filenames
  dirs_spc_files <- list.files(subdirs, pattern = paste0("^[0-9]{4}.*", ".", format), full.names = TRUE)
  # load spectral data
  data <- dirs_spc_files %>% 
    # read list of files
    lapply(., read_spectrum, format = format) %>% 
    # bind to tible
    data.table::rbindlist(.,) %>% tibble::as_tibble()
  print("done")
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
                           p = 3, w = 21, m = 0, 
                           col_in = "rflt", 
                           new_col = TRUE, 
                           snv = FALSE,
                           trim = c(1350, 1475, 1781, 1990, 2400, 2500),
                           binning = 3){
  
  print("pre-processing spectra ...")
  
  # Convert list of data.tables to one data.table
  spc_raw <- data.table::rbindlist(data[[col_in]])
  
  # set column names depending on selected parameters
  ptype <- paste0("p", p, "w", w, "m", m)
  
  # apply SavitzkyGolay filter to smooth spectra
  spc_pp <- with(data, prospectr::savitzkyGolay(spc_raw, p = p, w = w, m = m))
  
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
  if(new_col == F){
    spc_pp <- data[,-grep(col_in, colnames(data))]
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
#' @return A tibble with list columns, 
detect_outlier_spectra <- function(data, col_in = "all",
                                   grouping = NULL,
                                   outliers_rm = NULL,
                                   create_plot = T){
  
  print("detecting multivariate outliers ...")
  
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
                    facets = grouping)
    }
    
    out_idx <- which(outlier == 1)
    out_meas[[j]] <- data$meas_id[out_idx]
    
  }
  
  # drop outliers
  if(!is.null(outliers_rm)){
    
    print("removing outliers from dataset ...")
    
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
#' @param col_in Column of the tibble to containing the spectra to use
#' @param new_col Boolean, whether to add output as a new list column or overwrite the input column
#' @return A tibble with list columns
#' @export
calculate_SVI <- function(data, col_in, new_col = T) {
  
  print("calculating spectral indices ...")
  
  d <- data.table::rbindlist(dd[[col_in]]) %>% tibble::as_tibble()
  
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
  
  DF <- do.call(cbind.data.frame, mget(ls(pattern = "SI_"))) %>% tibble::as.tibble()
  
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

#' Sclaes spectral indices to range from 0 to 10 
#' representing the minimum and maximum value observed in a time series of measurements
#' @param data A tibble with SVI values in a list column of data.tables, as returned by "calculate_SVI"
#' @param plotid The variable name of the plot identifier
#' @return A tibble with list columns
#' @details This function reverts the scale for SVI with an increase during the measurement period
#' @export
scale_SVI <- function(data, plotid = "Plot_ID") {
  
  # extract SVI from list column and add required metadata
  meta <- data[c(plotid, "meas_date")]
  SVI_dat <- data.table::rbindlist(data[["SVI"]])
  d <- cbind(meta, SVI_dat)
  
  # fix plot identifier 
  colid <- which(names(d)==plotid)
  names(d)[colid] <- "Plot_ID"
  
  # keep only plots with measurements covering the entire measurement period
  # get start and end date
  min_date <- min(d$meas_date)
  max_date <- max(d$meas_date)
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
  ids <- d_scaled[c("Plot_ID", "meas_date")] %>% ungroup()
  SVI_sc_dat <- d_scaled[grepl("^SI_", names(d_scaled))]
  SVI_sc_dat <- SVI_sc_dat %>% 
    mutate_all(funs(r = revert)) %>% 
    #select original or reversed values
    dplyr::select_if(function(col) col[1] > 5) %>% 
    data.table::as.data.table()
  
  spc_pre_list <- map(purrr::transpose(SVI_sc_dat), data.table::as.data.table)
  ids[, "SVI_sc"] <- list(spc_pre_list)
  
  # join with spectral data
  data_out <- full_join(SVI, ids, by = c("Plot_ID", "meas_date"))
  
} 


# ============================================================================================================= -

# helper functions ---- 

#' Read an individual spectrum from an .sed or an .asd file
#' @param dir full file name
#' @param format file extension as a character string, "sed" or "asd"
#' @return A tibble with measurement_id, measurement date and the spectrum in a data.table
#' @export
read_spectrum <- function(dir, format = "sed"){
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
    tibble::tibble(
      meas_id = filename,
      meas_date = meas_date,
      rflt = list(spct)
    )
  } else {
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
                          mark_outliers = T, 
                          facets = "meas_date"){
  
  print("plotting spectra ...")
  
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
  if(!file.exists("Output")){
    dir.create("Output")
  }
  
  # save plot to pdf
  levels = nrow(unique(pd[facets]))
  width = 10
  pdf(paste0("Output/spectra_", col_in,".pdf"), width = width, height = width * (levels/2))
  # pdf(paste0("C:/Users/anjonas/output/spectra_", col_in,".pdf"), width = width, height = width * (levels/2))
  plot(plot)
  dev.off()
  
}

plot_spc <- function(data, col_in = "all", treatment = NULL,
                     facets = "meas_date",
                     mark_outliers = FALSE){
  
  print("plotting spectra ...")
  
  if("all" %in% col_in){
    col_in <- grep("^rflt_", names(data), value = TRUE)
  }
  
  # iterate over data types
  for (colunm in col_in) {
    
    # reshape data for plotting
    pd <- data %>% dplyr::select(meas_id, meas_date, treatment, colunm) %>% unnest(c(colunm)) %>% 
      tidyr::gather(wvlt, rflt, grep("^[0-9]", names(.)), factor_key = TRUE) %>% 
      mutate(wvlt = as.numeric(as.character(wvlt))) %>% 
      add_spc_range()
    
    otype = paste0("out_", colunm)
    
    # create empty plot
    plot <- ggplot() + 
      # geom_abline(intercept = 0, slope = 0) +
      theme_bw() + 
      theme(panel.grid = element_blank())
    
    # add spectra
    if(!is.null(treatment)){
      # factorise
      pd[[treatment]] <- as.factor(pd[[treatment]])
      # rename treatment column 
      colnames(pd[treatment]) <- "Treatment"
      plot <- plot +
        geom_line(data = pd, aes(x = wvlt, y = rflt,  group = interaction(meas_id, spc_range), col = Treatment), size = 0.05, alpha = 0.25) +
        scale_color_brewer(palette = "YlGnBu")
    } else {
      plot <- plot + 
        geom_line(data = pd, aes(x = wvlt, y = rflt,  group = interaction(meas_id, spc_range)), size = 0.05, alpha = 0.25)
    }
    
    # add facets
    if(!is.null(facets)){
      plot <- plot +
        facet_wrap(as.formula(paste0("~", facets)), ncol = 1, scales = "free")
    }
    
    # mark outliers
    if(mark_outliers){
      
      if(otype %in% names(data)){
        
        print("marking outliers ...")
        
        # reshape data for plotting
        pd_out <- data %>% dplyr::select(meas_id, meas_date, treatment, colunm, otype) %>% unnest(c(colunm)) %>% 
          tidyr::gather(wvlt, rflt, grep("^[0-9]", names(.)), factor_key = TRUE) %>% 
          mutate(wvlt = as.numeric(as.character(wvlt))) %>% 
          add_spc_range()
        
        # mark outliers
        plot <- plot +
          geom_line(data = pd_out[pd_out[otype] == 1,], aes(x = wvlt, y = rflt, group = interaction(meas_id, spc_range)),
                    size = 0.05, alpha = 1, col = "darkred", lty = "dashed")
        
      } else {
        print("no data for outliers found!")
      }
      
    }
    
    # check if Output directory exists
    if(!file.exists("Output")){
      dir.create("Output")
    }
    
    # save plot to pdf
    levels = nrow(unique(pd[facets]))
    width = 10
    pdf(paste0("Output/spectra_trt_", colunm,".pdf"), width = width, height = width * (levels/2))
    # pdf(paste0("C:/Users/anjonas/output/spectra_", colunm,".pdf"), width = width, height = width * (levels/2))
    plot(plot)
    dev.off()
    
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

# ============================================================================================================= -
  

