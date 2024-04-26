# By: Jasper Denissen
# 2021-10-04
# Script to read ERA5 data
# 1981 - 2020, 2.0x2.0 grid cell resolution WITH THE RASTER PACKAGE!
# 1) Read in all data and average to monthly timescale
# 2) detrend
# 3) anomalies
# 4) keep only anomalies when LAI_hv & LAI_lv > threshold (.2)
# 5) calculate correlation difference three hottest months per decade
pdf(NULL)

###############################################################################################################
###############################################################################################################
######################## !!! Don't forget to reset the working directory to a directory of your choosing!!! ###
###############################################################################################################
###############################################################################################################
setwd('/Net/Groups/BGI/work_3/HydroBioClim/archive/Denissen_etal_2024_ESD/')

# packages
source('Scripts/to_be_loaded_packages.R')

# RData
path_RData <- "RData/"
# # testdir
# path_RData <- "testdir/"

# load functions
source('Scripts/parallel_raster_functions.R')
swvlall_fun <- function(a){
  ifelse(sum(!is.na(a)) == 4,
         swvl1234 <- (a[1]*.07+a[2]*.21+a[3]*.72+a[4]*1.89)/2.89,
         swvl1234 <- NaN)
}

swvl123_fun <- function(a){
  ifelse(sum(!is.na(a)) == 3,
         swvl123 <- (a[1]*.07+a[2]*.21+a[3]*.72),
         swvl123 <- NaN)
}

path_ERA5_t2mmax_2d00 <- "Data/ERA5_t2m_2d00_daily/"
all_files_ERA5_t2mmax_2d00 <- list.files(path_ERA5_t2mmax_2d00,
                                         pattern = ".nc$")
all_files_ERA5_t2mmax_2d00[c(1,70)]

t2mmax_d <- rast(paste0(path_ERA5_t2mmax_2d00, all_files_ERA5_t2mmax_2d00))
day1 <- dmy("01-01-1951")
ndays <- day1 + c(0:(dim(t2mmax_d)[3]-1))
m_i <- c()
month_count <- 0
for(year in 1951:2020){
  m_i[(length(m_i)+1):(length(which(year(ndays) == year))+length(m_i))] <- month(ndays[which(year(ndays) == year)]) + month_count
  month_count <- month_count + 12
}
stack_t2mmax_raw <- tapp(t2mmax_d, index = m_i, fun = max, na.rm = T)
moy <- rep(1:12,70)
t2mmax_moy <- tapp(stack_t2mmax_raw, index = moy, fun = max, na.rm = T)
t2mmax_moy.array <- as.array(t2mmax_moy)

hottest_3_moy <- array(NaN,c(90,180,12))
for(y in 1:90){
  for(x in 1:180){
    hottest_3_moy[y,x,which(rank(t2mmax_moy.array[y,x,],ties.method = 'random') > 9)] <- 1
  }
}
# apparently the array needs to be flipped upside down...
hottest_3_moy <- hottest_3_moy

# consistency check
check <- array(NaN,c(90,180))
for(y in 1:90){
  for(x in 1:180){
    check[y,x] <- sum(hottest_3_moy[y,x,], na.rm=T)
  }
}

# convert this into a stack!
save(hottest_3_moy, file = paste0(path_RData, "ERA5_hottest_3_moy.RData"))
load(paste0(path_RData, "ERA5_hottest_3_moy.RData"))

hottest_3_moy_stack <- rast(hottest_3_moy)
crs(hottest_3_moy_stack) <- crs(stack_t2mmax_raw)
xmin(hottest_3_moy_stack) <- -180
xmax(hottest_3_moy_stack) <- 180
ymin(hottest_3_moy_stack) <- -90
ymax(hottest_3_moy_stack) <- 90

############################## READING THE REGRIDDED DATA ##########################################

# soil moisture
path_swvl1 <- "Data/era5land/swvl1/"
list_swvl1 <- list.files(path = path_swvl1, pattern = ".nc$")[2:71]
list_swvl1[c(1,70)]
stack_swvl1_raw <- rast(paste(path_swvl1,list_swvl1, sep=''))
path_swvl2 <- "Data/era5land/swvl2/"
list_swvl2 <- list.files(path = path_swvl2, pattern = ".nc$")[2:71]
list_swvl2[c(1,70)]
stack_swvl2_raw <- rast(paste(path_swvl2,list_swvl2, sep=''))
path_swvl3 <- "Data/era5land/swvl3/"
list_swvl3 <- list.files(path = path_swvl3, pattern = ".nc$")[2:71]
list_swvl3[c(1,70)]
stack_swvl3_raw <- rast(paste(path_swvl3,list_swvl3, sep=''))
path_swvl4 <- "Data/era5land/swvl4/"
list_swvl4 <- list.files(path = path_swvl4, pattern = ".nc$")[2:71]
list_swvl4[c(1,70)]
stack_swvl4_raw <- rast(paste(path_swvl4,list_swvl4, sep=''))

# evapotranspiration
path_slhf <- "Data/era5land/slhf/"
list_slhf <- list.files(path = path_slhf, pattern = ".nc$")[2:71]
list_slhf[c(1,70)]
stack_slhf_raw <- rast(paste(path_slhf,list_slhf,sep=''))

# lai_hv
path_lai_hv <- "Data/era5land/lai_hv/"
list_lai_hv <- list.files(path = path_lai_hv, pattern = ".nc$")[2:71]
list_lai_hv[c(1,70)]
stack_lai_hv_raw <- rast(paste(path_lai_hv,list_lai_hv,sep=''))

# lai_lv
path_lai_lv <- "Data/era5land/lai_lv/"
list_lai_lv <- list.files(path = path_lai_lv, pattern = ".nc$")[2:71]
list_lai_lv[c(1,70)]
stack_lai_lv_raw <- rast(paste(path_lai_lv,list_lai_lv,sep=''))

# net radiation
path_ssrd <- "Data/era5land/ssrd/"
list_ssrd <- list.files(path = path_ssrd, pattern = ".nc$")[2:71]
list_ssrd[c(1,70)]
stack_ssrd_raw <- rast(paste(path_ssrd,list_ssrd,sep=''))

# t2m
path_t2m <- "Data/era5land/t2m/"
list_t2m <- list.files(path = path_t2m, pattern = ".nc$")[2:71]
list_t2m[c(1,70)]
stack_t2m_raw <- rast(paste(path_t2m,list_t2m,sep=''))


for(i in 1:840){
  swvlall_pre <- c(stack_swvl1_raw[[i]],stack_swvl2_raw[[i]],stack_swvl3_raw[[i]],stack_swvl4_raw[[i]])
  swvl123_pre <- c(stack_swvl1_raw[[i]],stack_swvl2_raw[[i]],stack_swvl3_raw[[i]])
  if(i == 1){
    stack_swvlall_raw <- app(swvlall_pre, swvlall_fun)
    stack_swvl123_raw <- app(swvl123_pre, swvl123_fun)
  }else{
    stack_swvlall_raw <- c(stack_swvlall_raw, app(swvlall_pre, swvlall_fun))
    stack_swvl123_raw <- c(stack_swvl123_raw, app(swvl123_pre, swvl123_fun))
  }
  print(i)
}

# Mask all data when t2m_mon < 273.15+10
t2m_mask <- app(stack_t2m_raw, temp_K_mask_fun)

# First, calculate both masks for lai_hv and lai_lv, which gives a 1 when there is more than 0.5 fraction and 0 when not
lai_hv_mask <- app(stack_lai_hv_raw, lai_hlv_mask_fun02)
lai_lv_mask <- app(stack_lai_lv_raw, lai_hlv_mask_fun02)
lai_hlv_mask <- lai_hv_mask + lai_lv_mask
lai_hlv_mask_combined <- app(lai_hlv_mask, lai_hlv_mask_fun_combined)
combined_mask <- t2m_mask * lai_hlv_mask_combined

# # Because ERA5 LE is negative, multiply by -1
slhf_mon <- stack_slhf_raw * -1

moy <- rep(1:12,10)
month_year <- rep(1:10,each=12)

hottest_3_moy_rep <- rep(hottest_3_moy_stack, 10)
for(t in seq(1,840,120)){ # so the loop moves 10 years and in the loop we read in 5 years of data.
  # t2mmax
  t2mmax_10yr <- stack_t2mmax_raw[[t:(t+119)]]
  t2mmax_10yr_wrm_h3m <- t2mmax_10yr * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
  month_year_max_t2mmax_wrm_h3m <- tapp(t2mmax_10yr_wrm_h3m, index = month_year, fun = max, na.rm = T)
  T2maxbar_10yr <- app(month_year_max_t2mmax_wrm_h3m, fun = mean, na.rm = T)
  
  
  # t2m
  t2m_10yr <- stack_t2m_raw[[t:(t+119)]]
  detrend_t2m_10yr <- app(t2m_10yr, fun = detrend_fun_10yr)
  moy_t2m_10yr <- tapp(detrend_t2m_10yr, index = moy, fun = mean, na.rm = T)
  moy_t2m_10yr_rep <- rep(moy_t2m_10yr, 10)
  anom_t2m_wrm_10yr_h3m <- (detrend_t2m_10yr - moy_t2m_10yr_rep) * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
  t2m_10yr_wrm_h3m <- t2m_10yr * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
  av_t2m_10yr <- app(t2m_10yr_wrm_h3m, fun = mean, na.rm = T)
  
  
  # ssrd
  ssrd_10yr <- stack_ssrd_raw[[t:(t+119)]]
  detrend_ssrd_10yr <- app(ssrd_10yr, fun = detrend_fun_10yr)
  moy_ssrd_10yr <- tapp(detrend_ssrd_10yr, index = moy, fun = mean, na.rm = T)
  moy_ssrd_10yr_rep <- rep(moy_ssrd_10yr, 10)
  anom_ssrd_wrm_10yr_h3m <- (detrend_ssrd_10yr - moy_ssrd_10yr_rep) * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
  
  
  
  
  # slhf
  slhf_10yr <- slhf_mon[[t:(t+119)]]
  detrend_slhf_10yr <- app(slhf_10yr, fun = detrend_fun_10yr)
  moy_slhf_10yr <- tapp(detrend_slhf_10yr, index = moy, fun = mean, na.rm = T)
  moy_slhf_10yr_rep <- rep(moy_slhf_10yr, 10)
  anom_slhf_wrm_10yr_h3m <- (detrend_slhf_10yr - moy_slhf_10yr_rep) * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
  
  # swvlall
  swvlall_10yr <- stack_swvlall_raw[[t:(t+119)]]
  detrend_swvlall_10yr <- app(swvlall_10yr, fun = detrend_fun_10yr)
  moy_swvlall_10yr <- tapp(detrend_swvlall_10yr, index = moy, fun = mean, na.rm = T)
  moy_swvlall_10yr_rep <- rep(moy_swvlall_10yr, 10)
  anom_swvlall_wrm_10yr_h3m <- (detrend_swvlall_10yr - moy_swvlall_10yr_rep) * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
  
  # swvl123
  swvl123_10yr <- stack_swvl123_raw[[t:(t+119)]]
  detrend_swvl123_10yr <- app(swvl123_10yr, fun = detrend_fun_10yr)
  moy_swvl123_10yr <- tapp(detrend_swvl123_10yr, index = moy, fun = mean, na.rm = T)
  moy_swvl123_10yr_rep <- rep(moy_swvl123_10yr, 10)
  anom_swvl123_wrm_10yr_h3m <- (detrend_swvl123_10yr - moy_swvl123_10yr_rep) * combined_mask[[t:(t+119)]] * hottest_3_moy_rep
  
  
  
  # ELI swvlall
  stack_anom_t2m_slhf_wrm <- c(anom_t2m_wrm_10yr_h3m, anom_slhf_wrm_10yr_h3m)
  stack_anom_ssrd_slhf_wrm <- c(anom_ssrd_wrm_10yr_h3m, anom_slhf_wrm_10yr_h3m)
  stack_anom_swvlall_slhf_wrm <- c(anom_swvlall_wrm_10yr_h3m, anom_slhf_wrm_10yr_h3m)
  stack_anom_swvl123_slhf_wrm <- c(anom_swvl123_wrm_10yr_h3m, anom_slhf_wrm_10yr_h3m)
  corr_rgy_veg_10yr <- app(stack_anom_t2m_slhf_wrm, fun = corr_fun_h3m)
  corr_rgy_veg_10yr_ssrd <- app(stack_anom_ssrd_slhf_wrm, fun = corr_fun_h3m)
  corr_wtr_veg_10yr <- app(stack_anom_swvlall_slhf_wrm, fun = corr_fun_h3m)
  corr_wtr_veg_10yr_swvl123 <- app(stack_anom_swvl123_slhf_wrm, fun = corr_fun_h3m)
  
  
  if(t == 1){
    corr_rgy_veg <- corr_rgy_veg_10yr
    corr_rgy_veg_ssrd <- corr_rgy_veg_10yr_ssrd
    corr_wtr_veg <- corr_wtr_veg_10yr
    corr_wtr_veg_swvl123 <- corr_wtr_veg_10yr_swvl123
    T2maxbar <- T2maxbar_10yr
    av_t2m <- av_t2m_10yr
  }else{
    corr_rgy_veg <- c(corr_rgy_veg, corr_rgy_veg_10yr)
    corr_rgy_veg_ssrd <- c(corr_rgy_veg_ssrd, corr_rgy_veg_10yr_ssrd)
    corr_wtr_veg <- c(corr_wtr_veg, corr_wtr_veg_10yr)
    corr_wtr_veg_swvl123 <- c(corr_wtr_veg_swvl123, corr_wtr_veg_10yr_swvl123)
    T2maxbar <- c(T2maxbar, T2maxbar_10yr)
    av_t2m <- c(av_t2m, av_t2m_10yr)
  }
  print(paste("Month ",t," is done...",sep=''))
}

corr_rgy_veg.array <- aperm(as.array(corr_rgy_veg), c(2,1,3))[,90:1,]
corr_rgy_veg_ssrd.array <- aperm(as.array(corr_rgy_veg_ssrd), c(2,1,3))[,90:1,]
corr_wtr_veg.array <- aperm(as.array(corr_wtr_veg), c(2,1,3))[,90:1,]
corr_wtr_veg_swvl123.array <- aperm(as.array(corr_wtr_veg_swvl123), c(2,1,3))[,90:1,]
T2maxbar.array <- aperm(as.array(T2maxbar), c(2,1,3))[,90:1,]
av_t2m.array <- aperm(as.array(av_t2m), c(2,1,3))[,90:1,]
save(corr_rgy_veg.array,
     corr_rgy_veg_ssrd.array,
     corr_wtr_veg.array,
     corr_wtr_veg_swvl123.array,
     T2maxbar.array,
     av_t2m.array,
     
     
     file = paste0(path_RData, "202207_dcorr_ERA5_10yr_h3m_daily_combined_mask_terra_tas_ssrd.RData"))
