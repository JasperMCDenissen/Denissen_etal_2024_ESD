# By: Jasper Denissen
# 2024-01-06
# Script to read CMIP6 data
# 1981 - 2100, 2.0x2.0 grid cell resolution WITH THE TERRA PACKAGE!
# 1) Read in all data and average to monthly timescale
# 2) detrend
# 3) anomalies
# 4) keep only anomalies when T > threshold (10 deg C) and LAI > .2
# 5) calculate correlation difference three hottest months per decade
pdf(NULL)

#####################################################################################################################
#####################################################################################################################
######################## !!! Don't forget to reset the working directory to a directory of your choosing!!! #########
#####################################################################################################################
#####################################################################################################################
setwd('/Net/Groups/BGI/work_3/HydroBioClim/archive/Denissen_etal_2024_ESD/')

# packages
source('Scripts/to_be_loaded_packages.R')

# load table containing source_id and member_id of the regridded CMIP6
path_cmip6 <- "Data/cmip6_202312_climex_ESD_acccmip6/"
path_cmip6_2.0x2.0 <- "Data/cmip6_202312_climex_ESD_acccmip6_2.0x2.0/"
# RData
path_RData <- "RData/"
# # testdir
# path_RData <- "testdir/"
list_all <- list.files(path_cmip6)
list_tas <- list_all[grep("tas_*",list_all)]

source_id <- member_id <- c()
for(i in 1:length(list_tas)){
  source_id[i] <- strsplit(list_tas[i],"_")[[1]][3]
  member_id[i] <- strsplit(list_tas[i],"_")[[1]][5]
}

cmip6_data.df <- setNames(data.frame(matrix(ncol = 2, nrow = 0)),
                          c("source_id","member_id"))
for(source in unique(source_id)){
  cmip6_data.df <- rbind(cmip6_data.df,
                         data.frame("source_id" = source,
                                    "member_id" = member_id[which(source_id == source)[1]]))
}
cmip6_data.df$SM <- 'mrsol'
cmip6_data.df$SM[c(3,10:13)] <- 'mrso'
cmip6_data.df <- cmip6_data.df[which(cmip6_data.df$SM == 'mrsol'),1:2]

##########################################################################################
##########################################################################################
########################################## 2d00 ##########################################
##########################################################################################
##########################################################################################

# load functions
source('Scripts/parallel_raster_functions.R')

files <- list.files(path_cmip6_2.0x2.0,
                    pattern = "_2.0x2.0.nc$")
# make sure to fetch the right tas and tasmax: first half of grep files are tas and second half are tasmax
files_tas <- files[grepl("tas_*",files)][1:12]
files_tasmax <- files[grepl("tasmax_*",files)]
files_mrsol <- files[grepl("mrsol_*",files)]
files_hfls <- files[grepl("hfls_*",files)]
files_rlus <- files[grepl("rlus_*",files)]
files_rsus <- files[grepl("rsus_*",files)]
files_rlds <- files[grepl("rlds_*",files)]
files_rsds <- files[grepl("rsds_*",files)]
files_lai <- files[grepl("lai_*",files)]

max_tasmax_moy.list <- list()
moy <- rep(1:12,120) # all months consecutively for 120 years. This allows to pick the maximum temperature per month of year
for(i in 1:length(cmip6_data.df$source_id)){
  stack_tasmax <- rast(paste0(path_cmip6_2.0x2.0,files_tasmax[i]))
  max_tasmax_moy <- tapp(stack_tasmax, index = moy, fun = max, na.rm = T)
  max_tasmax_moy.list[[i]] <- as.array(max_tasmax_moy)
  save(max_tasmax_moy.list,
       file = paste0(path_RData, "202401_max_tasmax_moy_terra.RData"))
}

load(paste0(path_RData, "202401_max_tasmax_moy_terra.RData"))
hottest_3_moy <- array(NaN,c(90,180,12,12)) # lat,lon,#source,#moy
for(i in 1:length(cmip6_data.df$source_id)){
  for(y in 1:90){
    for(x in 1:180){
      hottest_3_moy[y,x,i,which(rank(max_tasmax_moy.list[[i]][y,x,],ties.method = 'random') > 9)] <- 1
    }
  }
}

# consistency check
check <- array(NaN,c(90,180,12))
for(i in 1:length(cmip6_data.df$source_id)){
  for(y in 1:90){
    for(x in 1:180){
      check[y,x,i] <- sum(hottest_3_moy[y,x,i,], na.rm=T)
    }
  }
}

save(hottest_3_moy, file = paste0(path_RData, "hottest_3_moy_acccmip6.RData"))
load(paste0(path_RData, "hottest_3_moy_acccmip6.RData"))

for(i in c(1:12)){
  if(i == 1){
    hottest_3_moy_stack <- rast(hottest_3_moy[,,i,])
  }else{
    hottest_3_moy_stack <- c(hottest_3_moy_stack, rast(hottest_3_moy[,,i,]))
  }
}
crs(hottest_3_moy_stack) <- crs(stack_tasmax)
xmin(hottest_3_moy_stack) <- 0
xmax(hottest_3_moy_stack) <- 360
ymin(hottest_3_moy_stack) <- -90
ymax(hottest_3_moy_stack) <- 90

moy <- rep(1:12,10) # all months consecutively for 10 years.
month_year <- rep(1:10,each=12)
hottest_3_moy_index <- c(1:12) # change this to whatever model you start at

# when script is run for the first time, allocate all lists
av_tas.list <-
  av_EF.list <-
  av_netrad.list <-
  av_rsds.list <-
  Tmaxbar.list <-
  dcorr.list <-
  corr_rgy_veg.list <-
  corr_wtr_veg.list <-
  dcorr_netrad.list <-
  corr_rgy_veg_netrad.list <-
  dcorr_rsds.list <-
  corr_rgy_veg_rsds.list <-
  tas_mask.list <-
  lai_mask.list <-
  combined_mask.list <-
  list()

for(i in 1:length(cmip6_data.df$source_id)){
  stack_tasmax <- rast(paste0(path_cmip6_2.0x2.0,files_tasmax[i]))
  stack_tas <- rast(paste0(path_cmip6_2.0x2.0,files_tas[i]))
  stack_mrsol <- rast(paste0(path_cmip6_2.0x2.0,files_mrsol[i]))
  stack_hfls <- rast(paste0(path_cmip6_2.0x2.0,files_hfls[i]))
  stack_rsus <- rast(paste0(path_cmip6_2.0x2.0,files_rsus[i]))
  stack_rlus <- rast(paste0(path_cmip6_2.0x2.0,files_rlus[i]))
  stack_rsds <- rast(paste0(path_cmip6_2.0x2.0,files_rsds[i]))
  stack_rlds <- rast(paste0(path_cmip6_2.0x2.0,files_rlds[i]))
  stack_lai <- rast(paste0(path_cmip6_2.0x2.0,files_lai[i]))
  
  # Mask all data when tas < 273.15+10 K
  tas_mask <- app(stack_tas, fun=temp_K_mask_fun)
  # Mask all data when stack_lai < 0.5
  lai_mask <- app(stack_lai, fun=lai_mask_fun02)
  # when either too cold or not enough vegetation, Na will appear in the combined_mask
  combined_mask <- tas_mask * lai_mask
  hottest_3_moy_rep <- rep(hottest_3_moy_stack[[hottest_3_moy_index]], 10)
  # analyze every 10 year block and shift by 10 years. starting at 1981-1990, 1991-2000, 2001-2010, 2011-2020, etc
  for(t in seq(1,1440,120)){ # so the loop moves 10 years and in the loop we read in 10 years of data.
    # tasmax
    tasmax_10yr <- stack_tasmax[[t:(t+119)]]
    tasmax_10yr_wrm_h3m <- tasmax_10yr * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
    month_year_max_tasmax_wrm_3hm <- tapp(tasmax_10yr_wrm_h3m, index = month_year, fun = max, na.rm = T)
    Tmaxbar_10yr <- app(month_year_max_tasmax_wrm_3hm, fun = mean, na.rm = T)
    # tas
    tas_10yr <- stack_tas[[t:(t+119)]]
    detrend_tas_10yr <- app(tas_10yr, fun = detrend_fun_10yr)
    moy_tas_10yr <- tapp(detrend_tas_10yr, index = moy, fun = mean, na.rm = T)
    moy_tas_10yr_rep <- rep(moy_tas_10yr, 10)
    anom_tas_wrm_10yr <- (detrend_tas_10yr - moy_tas_10yr_rep) * tas_mask[[t:(t+119)]]
    anom_tas_wrm_10yr_h3m <-  anom_tas_wrm_10yr * hottest_3_moy_rep * combined_mask[[t:(t+119)]]
    tas_10yr_wrm_h3m <- tas_10yr * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
    av_tas_10yr <- app(tas_10yr_wrm_h3m, fun = mean, na.rm = T)
    # hfls (conversion to mm/d: /26.15741)
    hfls_10yr <- stack_hfls[[t:(t+119)]]
    detrend_hfls_10yr <- app(hfls_10yr, fun = detrend_fun_10yr)
    moy_hfls_10yr <- tapp(detrend_hfls_10yr, index = moy, fun = mean, na.rm = T)
    moy_hfls_10yr_rep <- rep(moy_hfls_10yr, 10)
    anom_hfls_wrm_10yr <- (detrend_hfls_10yr - moy_hfls_10yr_rep) * tas_mask[[t:(t+119)]]
    anom_hfls_wrm_10yr_h3m <- anom_hfls_wrm_10yr * hottest_3_moy_rep * combined_mask[[t:(t+119)]]
    
    
    # mrsol
    mrsol_10yr <- stack_mrsol[[t:(t+119)]]
    detrend_mrsol_10yr <- app(mrsol_10yr, fun = detrend_fun_10yr)
    moy_mrsol_10yr <- tapp(detrend_mrsol_10yr, index = moy, fun = mean, na.rm = T)
    moy_mrsol_10yr_rep <- rep(moy_mrsol_10yr, 10)
    anom_mrsol_wrm_10yr <- (detrend_mrsol_10yr - moy_mrsol_10yr_rep) * tas_mask[[t:(t+119)]]
    anom_mrsol_wrm_10yr_h3m <- anom_mrsol_wrm_10yr * hottest_3_moy_rep * combined_mask[[t:(t+119)]]
    
    # rlds
    rlds_10yr <- stack_rlds[[t:(t+119)]]
    # rlus
    rlus_10yr <- stack_rlus[[t:(t+119)]]
    # rsds
    rsds_10yr <- stack_rsds[[t:(t+119)]]
    rsds_10yr_wrm_h3m <- rsds_10yr * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
    av_rsds_10yr <- app(rsds_10yr_wrm_h3m, fun = mean, na.rm = T)
    # rsus
    rsus_10yr <- stack_rsus[[t:(t+119)]]
    
    # calculate net radiation
    netrad_10yr <- (rsds_10yr - rsus_10yr) + (rlds_10yr - rlus_10yr)
    detrend_netrad_10yr <- app(netrad_10yr, fun = detrend_fun_10yr)
    moy_netrad_10yr <- tapp(detrend_netrad_10yr, index = moy, fun = mean, na.rm = T)
    moy_netrad_10yr_rep <- rep(moy_netrad_10yr, 10)
    anom_netrad_wrm_10yr <- (detrend_netrad_10yr - moy_netrad_10yr_rep) * tas_mask[[t:(t+119)]]
    anom_netrad_wrm_10yr_h3m <- anom_netrad_wrm_10yr * hottest_3_moy_rep * combined_mask[[t:(t+119)]]
    
    # calculations with rsds
    detrend_rsds_10yr <- app(rsds_10yr, fun = detrend_fun_10yr)
    moy_rsds_10yr <- tapp(detrend_rsds_10yr, index = moy, fun = mean, na.rm = T)
    moy_rsds_10yr_rep <- rep(moy_rsds_10yr, 10)
    anom_rsds_wrm_10yr <- (detrend_rsds_10yr - moy_rsds_10yr_rep) * tas_mask[[t:(t+119)]]
    anom_rsds_wrm_10yr_h3m <- anom_rsds_wrm_10yr * hottest_3_moy_rep * combined_mask[[t:(t+119)]]
    
    # calculate EF
    av_netrad_wrm_10yr_h3m <- app((netrad_10yr * combined_mask[[t:(t+119)]] * hottest_3_moy_rep), fun = mean, na.rm = T)
    av_hfls_wrm_10yr_h3m <- app((hfls_10yr * combined_mask[[t:(t+119)]] * hottest_3_moy_rep), fun = mean, na.rm = T)
    av_EF_10yr <- av_hfls_wrm_10yr_h3m/av_netrad_wrm_10yr_h3m
    
    # Ta > 10 deg C & LAI > 0.5 & warm season (h3m)
    # ELI
    stack_anom_tas_hfls_wrm <- c(anom_tas_wrm_10yr_h3m, anom_hfls_wrm_10yr_h3m)
    stack_anom_mrsol_hfls_wrm <- c(anom_mrsol_wrm_10yr_h3m, anom_hfls_wrm_10yr_h3m)
    corr_rgy_veg_10yr <- app(stack_anom_tas_hfls_wrm, fun = corr_fun_h3m)
    corr_wtr_veg_10yr <- app(stack_anom_mrsol_hfls_wrm, fun = corr_fun_h3m)
    dcorr_10yr <- corr_rgy_veg_10yr - corr_wtr_veg_10yr
    
    
    # ELI with netrad
    stack_anom_netrad_hfls_wrm <- c(anom_netrad_wrm_10yr_h3m, anom_hfls_wrm_10yr_h3m)
    corr_rgy_veg_10yr_netrad <- app(stack_anom_netrad_hfls_wrm, fun = corr_fun_h3m)
    dcorr_10yr_netrad <- corr_rgy_veg_10yr_netrad - corr_wtr_veg_10yr    
    
    # ELI with rsds
    stack_anom_rsds_hfls_wrm <- c(anom_rsds_wrm_10yr_h3m, anom_hfls_wrm_10yr_h3m)
    corr_rgy_veg_10yr_rsds <- app(stack_anom_rsds_hfls_wrm, fun = corr_fun_h3m)
    dcorr_10yr_rsds <- corr_rgy_veg_10yr_rsds - corr_wtr_veg_10yr
    
    if(t == 1){
      dcorr <- dcorr_10yr
      corr_rgy_veg <- corr_rgy_veg_10yr
      corr_wtr_veg <- corr_wtr_veg_10yr
      dcorr_netrad <- dcorr_10yr_netrad
      corr_rgy_veg_netrad <- corr_rgy_veg_10yr_netrad
      dcorr_rsds <- dcorr_10yr_rsds
      corr_rgy_veg_rsds <- corr_rgy_veg_10yr_rsds
      Tmaxbar <- Tmaxbar_10yr
      av_tas <- av_tas_10yr
      av_EF <- av_EF_10yr
      av_netrad <- av_netrad_wrm_10yr_h3m
      av_rsds <- av_rsds_10yr
      tas_mask_r <- tas_mask
      lai_mask_r <- lai_mask
      combined_mask_r <- combined_mask
    }else{
      dcorr <- c(dcorr, dcorr_10yr)
      corr_rgy_veg <- c(corr_rgy_veg, corr_rgy_veg_10yr)
      corr_wtr_veg <- c(corr_wtr_veg, corr_wtr_veg_10yr)
      dcorr_netrad <- c(dcorr_netrad, dcorr_10yr_netrad)
      corr_rgy_veg_netrad <- c(corr_rgy_veg_netrad, corr_rgy_veg_10yr_netrad)
      dcorr_rsds <- c(dcorr_rsds, dcorr_10yr_rsds)
      corr_rgy_veg_rsds <- c(corr_rgy_veg_rsds, corr_rgy_veg_10yr_rsds)
      Tmaxbar <- c(Tmaxbar, Tmaxbar_10yr)
      av_tas <- c(av_tas, av_tas_10yr)
      av_EF <- c(av_EF, av_EF_10yr)
      av_netrad <- c(av_netrad, av_netrad_wrm_10yr_h3m)
      av_rsds <- c(av_rsds, av_rsds_10yr)
      tas_mask_r <- c(tas_mask_r, tas_mask)
      lai_mask_r <- c(lai_mask_r, lai_mask)
      combined_mask_r <- c(combined_mask_r, combined_mask)
    }
    print(paste(cmip6_data.df$source_id[i]," and month ",t," are done...",sep=''))
  }
  dcorr.list[[i]] <- aperm(as.array(dcorr), c(2,1,3))[,90:1,]
  corr_rgy_veg.list[[i]] <- aperm(as.array(corr_rgy_veg), c(2,1,3))[,90:1,]
  corr_wtr_veg.list[[i]] <- aperm(as.array(corr_wtr_veg), c(2,1,3))[,90:1,]
  dcorr_netrad.list[[i]] <- aperm(as.array(dcorr_netrad), c(2,1,3))[,90:1,]
  corr_rgy_veg_netrad.list[[i]] <- aperm(as.array(corr_rgy_veg_netrad), c(2,1,3))[,90:1,]
  dcorr_rsds.list[[i]] <- aperm(as.array(dcorr_rsds), c(2,1,3))[,90:1,]
  corr_rgy_veg_rsds.list[[i]] <- aperm(as.array(corr_rgy_veg_rsds), c(2,1,3))[,90:1,]
  Tmaxbar.list[[i]] <- aperm(as.array(Tmaxbar), c(2,1,3))[,90:1,]
  av_tas.list[[i]] <- aperm(as.array(av_tas), c(2,1,3))[,90:1,]
  av_EF.list[[i]] <- aperm(as.array(av_EF), c(2,1,3))[,90:1,]
  av_netrad.list[[i]] <- aperm(as.array(av_netrad), c(2,1,3))[,90:1,]
  av_rsds.list[[i]] <- aperm(as.array(av_rsds), c(2,1,3))[,90:1,]
  tas_mask.list[[i]] <- aperm(as.array(tas_mask), c(2,1,3))[,90:1,]
  lai_mask.list[[i]] <- aperm(as.array(lai_mask), c(2,1,3))[,90:1,]
  combined_mask.list[[i]] <- aperm(as.array(combined_mask), c(2,1,3))[,90:1,]
  
  
  save(tas_mask.list,
       lai_mask.list,
       combined_mask.list,
       dcorr.list,
       corr_rgy_veg.list,
       corr_wtr_veg.list,
       dcorr_netrad.list,
       corr_rgy_veg_netrad.list,
       dcorr_rsds.list,
       corr_rgy_veg_rsds.list,
       Tmaxbar.list,
       av_tas.list,
       av_EF.list,
       av_netrad.list,
       av_rsds.list,
       cmip6_data.df,
       
       file = paste0(path_RData, "202401_dcorr_cmip6_10yr_h3m_combined_mask_acccmip6_rsds.RData"))
  hottest_3_moy_index <- hottest_3_moy_index + 12
}
