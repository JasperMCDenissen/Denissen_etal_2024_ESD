# # By: Jasper Denissen
# 2024-03-12
# Script to make paper figures
pdf(NULL)

############################################################################################################################
############################################################################################################################
######################## !!! Don't forget to reset the working directory to a directory of your choosing!!! ################
############################################################################################################################
###############################################################################################################
setwd('/Net/Groups/BGI/work_3/HydroBioClim/archive/Denissen_etal_2024_ESD/')

# Get the proper packages
source('Scripts/to_be_loaded_packages.R')
# source functions
source('Scripts/calc_boxes.R')
source('Scripts/plot_discrete_cbar.R')


##########################################################################################
##########################################################################################
########################################## 2d00 ##########################################
##########################################################################################
##########################################################################################

# RData
load("RData/202401_dcorr_cmip6_10yr_h3m_combined_mask_acccmip6_rsds.RData")
load('RData/total_land_area_CMIP6_climex_acccmip6.RData')
path_RData <- "RData/"
# # testdir
# load("testdir/202401_dcorr_cmip6_10yr_h3m_combined_mask_acccmip6_rsds.RData")
# load('testdir/total_land_area_CMIP6_climex_acccmip6.RData')
# path_RData <- "testdir/"

# maybe now per surface area instead of % of grid cells?
r <- raster()  # by default 1 by 1 degree
res(r) <- 2 # so change the resolution
a <- raster::area(r) # calculate the area of a 2x2 degree grid from N - S, as area varies only by latitude, not longitude
area <- a[,1]
area.array <- array(NaN,c(180,90))
for(x in 1:180){
  area.array[x,] <- area
}

# lon <- seq(1,359,2)
lon <- seq(-179,179,2)
lat <- seq(-89,89,2)

# reg <- c("SAM", "NAM", "CEU", "NAS")
hotspot_regs.df <- data.frame("x" = c(lon[53]-1, lon[69]+1, lon[53]-1, lon[53]-1,
                                      lon[28]-1, lon[45]+1, lon[28]-1, lon[28]-1,
                                      lon[91]-1, lon[114]+1, lon[91]-1, lon[91]-1,
                                      lon[118]-1, lon[147]+1, lon[118]-1, lon[118]-1),
                              "y" = c(lat[34]-1, lat[34]-1, lat[34]-1, lat[46]+1,
                                      lat[67]-1, lat[67]-1, lat[67]-1, lat[73]+1,
                                      lat[64]-1, lat[64]-1, lat[64]-1, lat[72]+1,
                                      lat[75]-1, lat[75]-1, lat[75]-1, lat[79]+1),
                              "xend" = c(lon[53]-1, lon[69]+1, lon[69]+1, lon[69]+1,
                                         lon[28]-1, lon[45]+1, lon[45]+1, lon[45]+1,
                                         lon[91]-1, lon[114]+1, lon[114]+1, lon[114]+1,
                                         lon[118]-1, lon[147]+1, lon[147]+1, lon[147]+1),
                              "yend" = c(lat[46]+1, lat[46]+1, lat[34]-1, lat[46]+1,
                                         lat[73]+1, lat[73]+1, lat[67]-1, lat[73]+1,
                                         lat[72]+1, lat[72]+1, lat[64]-1, lat[72]+1,
                                         lat[79]+1, lat[79]+1, lat[75]-1, lat[79]+1))





count <- 1
for(i in 1:length(cmip6_data.df$source_id)){
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- dcorr.list[[count]][91:180,,]; test[91:180,,] <- dcorr.list[[count]][1:90,,]; dcorr.list[[i]] <- test; dcorr.list[[count]] <- -1*dcorr.list[[count]]
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- corr_rgy_veg.list[[count]][91:180,,]; test[91:180,,] <- corr_rgy_veg.list[[count]][1:90,,]; corr_rgy_veg.list[[i]] <- test
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- corr_wtr_veg.list[[count]][91:180,,]; test[91:180,,] <- corr_wtr_veg.list[[count]][1:90,,]; corr_wtr_veg.list[[i]] <- test
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- dcorr_rsds.list[[count]][91:180,,]; test[91:180,,] <- dcorr_rsds.list[[count]][1:90,,]; dcorr_rsds.list[[i]] <- test; dcorr_rsds.list[[count]] <- -1*dcorr_rsds.list[[count]]
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- corr_rgy_veg_rsds.list[[count]][91:180,,]; test[91:180,,] <- corr_rgy_veg_rsds.list[[count]][1:90,,]; corr_rgy_veg_rsds.list[[i]] <- test
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- Tmaxbar.list[[count]][91:180,,]; test[91:180,,] <- Tmaxbar.list[[count]][1:90,,]; Tmaxbar.list[[i]] <- test
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- av_tas.list[[count]][91:180,,]; test[91:180,,] <- av_tas.list[[count]][1:90,,]; av_tas.list[[i]] <- test
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- av_EF.list[[count]][91:180,,]; test[91:180,,] <- av_EF.list[[count]][1:90,,]; av_EF.list[[i]] <- test
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- av_rsds.list[[count]][91:180,,]; test[91:180,,] <- av_rsds.list[[count]][1:90,,]; av_rsds.list[[i]] <- test
  count <- count + 1
}

# now the strongest average energy proxy
dcorr_max.list <- list()
corr_rgy_veg_max.list <- corr_rgy_veg_rsds.list
tas_or_rsds.array <- array(NaN,c(180,90,12))
tas_cum.array <- array(0,c(180,90))
for(i in 1:length(cmip6_data.df$source_id)){
  for(x in 1:180){
    for(y in 1:90){
      if(sum(!is.na(corr_rgy_veg.list[[i]][x,y,])) > 0){
        if(mean(corr_rgy_veg.list[[i]][x,y,], na.rm = T) > mean(corr_rgy_veg_rsds.list[[i]][x,y,], na.rm = T)){
          corr_rgy_veg_max.list[[i]][x,y,] <- corr_rgy_veg.list[[i]][x,y,]
          tas_or_rsds.array[x,y,i] <- 1 # temperature
          tas_cum.array[x,y] <- tas_cum.array[x,y] + 1
        }else{
          tas_or_rsds.array[x,y,i] <- 2 # shortwave incoming radiation
        }
      }
    }
  }
  dcorr_max.list[[i]] <- corr_wtr_veg.list[[i]] - corr_rgy_veg_max.list[[i]]
  print(paste0("source_id ", cmip6_data.df$source_id[i], " is done..."))
}


# options(warn=2) # find out why we get errors with calculating the month-of-year trend in maximum daily temperature
# from this we can infer the month of year trends
kendall_dcorr <-
  kendall_dcorr_max <- 
  kendall_corr_rgy_veg <-
  kendall_corr_wtr_veg <-
  kendall_dcorr_rsds <-
  kendall_corr_rgy_veg_rsds <-
  kendall_Tmaxbar_min <-
  kendall_av_tas <-
  kendall_av_EF <-
  kendall_av_rsds <-
  cor_ELI_Tmaxbar_min_source_id <-
  cor_ELI_max_Tmaxbar_min_source_id <- 
  cor_ELI_rsds_Tmaxbar_min_source_id <-
  pcor_ELI_Tmaxbar_min_source_id <-
  cor_EF_Tmaxbar_min_source_id <- 
  cor_ELI_EF_source_id <-
  cor_ELI_max_EF_source_id <-
  array(NaN,c(180,90,12,2))
mask_cmip6 <- array(NaN,c(180,90,12))
for(x in 1:180){
  for(y in 1:90){
    for(source_id in 1:length(cmip6_data.df$source_id)){
      if(!is.na(models_with_full_timeseries[x,y])){
        if(models_with_full_timeseries[x,y] > 8){
          if(sum(!is.na(dcorr.list[[source_id]][x,y,])) == 12 & sum(!is.na(Tmaxbar.list[[source_id]][x,y,])) == 12){
            mask_cmip6[x,y,source_id] <- 1
            kendall_dcorr[x,y,source_id,] <- c(unname(kendallTrendTest(dcorr.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(dcorr.list[[source_id]][x,y,])$p.value))
            kendall_dcorr_max[x,y,source_id,]  <- c(unname(kendallTrendTest(dcorr_max.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(dcorr.list[[source_id]][x,y,])$p.value))
            kendall_corr_rgy_veg[x,y,source_id,] <- c(unname(kendallTrendTest(corr_rgy_veg.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(corr_rgy_veg.list[[source_id]][x,y,])$p.value))
            kendall_corr_wtr_veg[x,y,source_id,] <- c(unname(kendallTrendTest(corr_wtr_veg.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(corr_wtr_veg.list[[source_id]][x,y,])$p.value))
            kendall_Tmaxbar_min[x,y,source_id,] <- c(unname(kendallTrendTest(Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,])$p.value))
            kendall_av_tas[x,y,source_id,] <- c(unname(kendallTrendTest(av_tas.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(av_tas.list[[source_id]][x,y,])$p.value))
            kendall_av_EF[x,y,source_id,] <- c(unname(kendallTrendTest(av_EF.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(av_EF.list[[source_id]][x,y,])$p.value))
            kendall_av_rsds[x,y,source_id,] <- c(unname(kendallTrendTest(av_rsds.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(av_rsds.list[[source_id]][x,y,])$p.value))
            cor_ELI_Tmaxbar_min_source_id[x,y,source_id,] <- c(cor(dcorr.list[[source_id]][x,y,],
                                                                   Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs"),
                                                               cor.test(dcorr.list[[source_id]][x,y,],
                                                                        Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs")$p.value)
            cor_ELI_max_Tmaxbar_min_source_id[x,y,source_id,] <- c(cor(dcorr_max.list[[source_id]][x,y,],
                                                                       Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs"),
                                                                   cor.test(dcorr_max.list[[source_id]][x,y,],
                                                                            Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs")$p.value)
            cor_ELI_rsds_Tmaxbar_min_source_id[x,y,source_id,] <- c(cor(dcorr_rsds.list[[source_id]][x,y,],
                                                                        Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs"),
                                                                    cor.test(dcorr_rsds.list[[source_id]][x,y,],
                                                                             Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs")$p.value)
            pcor_ELI_Tmaxbar_min_source_id[x,y,source_id,] <- c(pcor.test(dcorr.list[[source_id]][x,y,], Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], 
                                                                          av_rsds.list[[source_id]][x,y,], method = 'kendall')$estimate,
                                                                pcor.test(dcorr.list[[source_id]][x,y,], Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], 
                                                                          av_rsds.list[[source_id]][x,y,], method = 'kendall')$p.value)
            cor_EF_Tmaxbar_min_source_id[x,y,source_id,] <- c(cor(av_EF.list[[source_id]][x,y,],
                                                                  Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs"),
                                                              cor.test(av_EF.list[[source_id]][x,y,],
                                                                       Tmaxbar.list[[source_id]][x,y,] - av_tas.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs")$p.value)
            cor_ELI_EF_source_id[x,y,source_id,] <- c(cor(dcorr.list[[source_id]][x,y,],
                                                          av_EF.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs"),
                                                      cor.test(dcorr.list[[source_id]][x,y,],
                                                               av_EF.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs")$p.value)
            cor_ELI_max_EF_source_id[x,y,source_id,] <- c(cor(dcorr_max.list[[source_id]][x,y,],
                                                              av_EF.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs"),
                                                          cor.test(dcorr_max.list[[source_id]][x,y,],
                                                                   av_EF.list[[source_id]][x,y,], method = 'kendall', use = "pairwise.complete.obs")$p.value)
            kendall_dcorr_rsds[x,y,source_id,] <- c(unname(kendallTrendTest(dcorr_rsds.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(dcorr_rsds.list[[source_id]][x,y,])$p.value))
          }
          if(sum(!is.na(dcorr_rsds.list[[source_id]][x,y,])) == 12 & sum(!is.na(Tmaxbar.list[[source_id]][x,y,])) == 12){
            
            kendall_corr_rgy_veg_rsds[x,y,source_id,] <- c(unname(kendallTrendTest(corr_rgy_veg_rsds.list[[source_id]][x,y,])$estimate[2]), unname(kendallTrendTest(corr_rgy_veg_rsds.list[[source_id]][x,y,])$p.value))
          }
        }
      }
    }
  }
  print(x)
}


# mmmean_kendall
mmmean_prior_dcorr <- # prior_dcorr is the average dcorr over the first 30 years
  mmmean_prior_dcorr_max <- 
  mmmean_prior_EF <- 
  mmmean_kendall_corr_rgy_veg <-
  mmmean_kendall_corr_wtr_veg <-
  mmmean_prior_dcorr_rsds <- 
  mmmean_prior_corr_rgy_veg <-
  mmmean_prior_corr_rgy_veg_rsds <- 
  mmmean_prior_corr_wtr_veg <-
  mmmean_kendall_av_tas <-
  mmmean_kendall_av_EF <-
  mmmean_kendall_av_rsds <-
  array(NaN,c(180,90))
mmmean_kendall_dcorr <-
  mmmean_kendall_dcorr_max <-
  mmmean_kendall_dcorr_rsds <- 
  mmmean_kendall_Tmaxbar_min <-
  mmmean_cor_ELI_Tmaxbar_min <-
  mmmean_cor_ELI_max_Tmaxbar_min <-
  mmmean_cor_ELI_rsds_Tmaxbar_min <-
  mmmean_pcor_ELI_Tmaxbar_min <-
  mmmean_cor_EF_Tmaxbar_min <-
  mmmean_cor_ELI_EF <-
  mmmean_cor_ELI_max_EF <-
  array(NaN,c(180,90,3))
for(x in 1:180){
  for(y in 1:90){
    prior_dcorr <- prior_dcorr_max <- prior_EF <- prior_corr_rgy_veg <- prior_corr_rgy_veg_rsds <- prior_corr_wtr_veg <- prior_dcorr_rsds <- c()
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        for(source_id in 1:length(cmip6_data.df$source_id)){
          prior_dcorr[source_id] <- mean(dcorr.list[[source_id]][x,y,1:3])
          prior_dcorr_max[source_id] <- mean(dcorr_max.list[[source_id]][x,y,1:3])
          prior_EF[source_id] <- mean(av_EF.list[[source_id]][x,y,1:3])
          prior_dcorr_rsds[source_id] <- mean(dcorr_rsds.list[[source_id]][x,y,1:3])
          prior_corr_rgy_veg[source_id] <- mean(corr_rgy_veg.list[[source_id]][x,y,])
          prior_corr_rgy_veg_rsds[source_id] <- mean(corr_rgy_veg_rsds.list[[source_id]][x,y,])
          prior_corr_wtr_veg[source_id] <- mean(corr_wtr_veg.list[[source_id]][x,y,])
        }
        mmmean_prior_dcorr[x,y] <- mean(prior_dcorr, na.rm = T)
        mmmean_prior_dcorr_max[x,y] <- mean(prior_dcorr_max, na.rm = T)
        mmmean_prior_EF[x,y] <- mean(prior_EF, na.rm = T)
        mmmean_prior_dcorr_rsds[x,y] <- mean(prior_dcorr_rsds, na.rm = T)
        mmmean_prior_corr_rgy_veg[x,y] <- mean(prior_corr_rgy_veg, na.rm = T)
        mmmean_prior_corr_rgy_veg_rsds[x,y] <- mean(prior_corr_rgy_veg_rsds, na.rm = T)
        mmmean_prior_corr_wtr_veg[x,y] <- mean(prior_corr_wtr_veg, na.rm = T)
        mmmean_kendall_av_tas[x,y] <- mean(kendall_av_tas[x,y,,1], na.rm = T)
        mmmean_kendall_av_EF[x,y] <- mean(kendall_av_EF[x,y,,1], na.rm = T)
        mmmean_kendall_av_rsds[x,y] <- mean(kendall_av_rsds[x,y,,1], na.rm = T)
        mmmean_kendall_dcorr[x,y,1] <- mean(kendall_dcorr[x,y,,1], na.rm = T)
        mmmean_kendall_dcorr_max[x,y,1] <- mean(kendall_dcorr_max[x,y,,1], na.rm = T)
        mmmean_kendall_corr_rgy_veg[x,y] <- mean(kendall_corr_rgy_veg[x,y,,1], na.rm = T)
        mmmean_kendall_corr_wtr_veg[x,y] <- mean(kendall_corr_wtr_veg[x,y,,1], na.rm = T)
        mmmean_kendall_dcorr_rsds[x,y,1] <- mean(kendall_dcorr_rsds[x,y,,1], na.rm = T)
        mmmean_kendall_Tmaxbar_min[x,y,1] <- mean(kendall_Tmaxbar_min[x,y,,1], na.rm = T)
        mmmean_cor_ELI_Tmaxbar_min[x,y,1] <- mean(cor_ELI_Tmaxbar_min_source_id[x,y,,1], na.rm = T)
        mmmean_cor_ELI_max_Tmaxbar_min[x,y,1] <- mean(cor_ELI_max_Tmaxbar_min_source_id[x,y,,1], na.rm = T)
        mmmean_cor_ELI_rsds_Tmaxbar_min[x,y,1] <- mean(cor_ELI_rsds_Tmaxbar_min_source_id[x,y,,1], na.rm = T)
        mmmean_pcor_ELI_Tmaxbar_min[x,y,1] <- mean(pcor_ELI_Tmaxbar_min_source_id[x,y,,1], na.rm = T)
        mmmean_cor_EF_Tmaxbar_min[x,y,1] <- mean(cor_EF_Tmaxbar_min_source_id[x,y,,1], na.rm = T)
        mmmean_cor_ELI_EF[x,y,1] <- mean(cor_ELI_EF_source_id[x,y,,1], na.rm = T)
        mmmean_cor_ELI_max_EF[x,y,1] <- mean(cor_ELI_max_EF_source_id[x,y,,1], na.rm = T)
      }
    }
  }
}

# we need spatial images of 1) trend max_tasmax, 2) average highest max_tasmax, 3) trend max_tas, 4) trend max_tasmax - av_tas
mmmean_prior_dcorr_h3m_vec <- 
  mmmean_prior_dcorr_max_h3m_vec <- 
  mmmean_prior_EF_h3m_vec <- 
  mmmean_prior_corr_rgy_veg_h3m_vec <- 
  mmmean_prior_corr_wtr_veg_h3m_vec <-
  mmmean_kendall_dcorr_h3m_vec <- 
  mmmean_kendall_dcorr_max_h3m_vec <- 
  mmmean_kendall_corr_rgy_veg_h3m_vec <- 
  mmmean_kendall_corr_wtr_veg_h3m_vec <- 
  mmmean_kendall_av_tas_h3m_vec <-
  mmmean_kendall_av_EF_h3m_vec <-
  mmmean_kendall_av_rsds_h3m_vec <-
  mmmean_kendall_Tmaxbar_min_h3m_vec <-
  mmmean_cor_ELI_Tmaxbar_min_vec <-
  mmmean_cor_ELI_max_Tmaxbar_min_vec <-
  mmmean_cor_ELI_rsds_Tmaxbar_min_vec <-
  mmmean_pcor_ELI_Tmaxbar_min_vec <-
  mmmean_cor_EF_Tmaxbar_min_vec <-
  mmmean_cor_ELI_EF_vec <-
  mmmean_cor_ELI_max_EF_vec <-
  prior_dcorr_h3m_vec <- 
  kendall_dcorr_h3m_vec <- 
  prior_dcorr_max_h3m_vec <- 
  kendall_dcorr_max_h3m_vec <- 
  prior_corr_rgy_veg_h3m_vec <- 
  kendall_corr_rgy_veg_h3m_vec <- 
  prior_corr_wtr_veg_h3m_vec <- 
  kendall_corr_wtr_veg_h3m_vec <- 
  kendall_av_tas_h3m_vec <-
  kendall_av_EF_h3m_vec <-
  kendall_av_rsds_h3m_vec <-
  kendall_Tmaxbar_min_h3m_vec <-
  lon_vec <- 
  lat_vec <- 
  prior_dcorr <- prior_dcorr_max <- prior_corr_rgy_veg <- prior_corr_wtr_veg <- cor_ELI_Tmaxbar_min_source_id_vec <- cor_ELI_max_Tmaxbar_min_source_id_vec <- 
  sign_kendall_dcorr_vec <-
  sign_kendall_dcorr_max_vec <-
  sign_kendall_Tmaxbar_min_vec <-
  sign_kendall_av_rsds_vec <-
  sign_cor_ELI_Tmaxbar_min_source_id_vec <-
  sign_cor_ELI_max_Tmaxbar_min_source_id_vec <-
  sign_cor_ELI_rsds_Tmaxbar_min_source_id_vec <-
  sign_mmmean_kendall_dcorr_h3m_vec <- 
  sign_mmmean_kendall_dcorr_max_h3m_vec <- 
  sign_mmmean_kendall_Tmaxbar_min_h3m_vec <- 
  tas_cum_vec <-
  c()
for(x in 1:180){
  for(y in 1:90){
    if(models_with_full_timeseries[x,y] > 8){
      lon_vec <- c(lon_vec, lon[x])
      lat_vec <- c(lat_vec, lat[y])
      mmmean_prior_dcorr_h3m_vec <- c(mmmean_prior_dcorr_h3m_vec, mmmean_prior_dcorr[x,y])
      mmmean_prior_dcorr_max_h3m_vec <- c(mmmean_prior_dcorr_max_h3m_vec, mmmean_prior_dcorr_max[x,y])
      mmmean_prior_EF_h3m_vec <- c(mmmean_prior_EF_h3m_vec, mmmean_prior_EF[x,y])
      mmmean_prior_corr_rgy_veg_h3m_vec <- c(mmmean_prior_corr_rgy_veg_h3m_vec, mmmean_prior_corr_rgy_veg[x,y])
      mmmean_kendall_corr_rgy_veg_h3m_vec <- c(mmmean_kendall_corr_rgy_veg_h3m_vec,mmmean_kendall_corr_rgy_veg[x,y])
      mmmean_prior_corr_wtr_veg_h3m_vec <- c(mmmean_prior_corr_wtr_veg_h3m_vec, mmmean_prior_corr_wtr_veg[x,y])
      mmmean_kendall_corr_wtr_veg_h3m_vec <- c(mmmean_kendall_corr_wtr_veg_h3m_vec,mmmean_kendall_corr_wtr_veg[x,y])
      mmmean_kendall_av_tas_h3m_vec <- c(mmmean_kendall_av_tas_h3m_vec,mmmean_kendall_av_tas[x,y])
      mmmean_kendall_av_EF_h3m_vec <- c(mmmean_kendall_av_EF_h3m_vec,mmmean_kendall_av_EF[x,y])
      mmmean_kendall_av_rsds_h3m_vec <- c(mmmean_kendall_av_rsds_h3m_vec,mmmean_kendall_av_rsds[x,y])
      mmmean_kendall_dcorr_h3m_vec <- c(mmmean_kendall_dcorr_h3m_vec,mmmean_kendall_dcorr[x,y,1])
      mmmean_kendall_dcorr_max_h3m_vec <- c(mmmean_kendall_dcorr_max_h3m_vec,mmmean_kendall_dcorr_max[x,y,1])
      mmmean_kendall_Tmaxbar_min_h3m_vec <- c(mmmean_kendall_Tmaxbar_min_h3m_vec,mmmean_kendall_Tmaxbar_min[x,y,1])
      mmmean_cor_ELI_Tmaxbar_min_vec <- c(mmmean_cor_ELI_Tmaxbar_min_vec, mmmean_cor_ELI_Tmaxbar_min[x,y,1])
      mmmean_cor_ELI_max_Tmaxbar_min_vec <- c(mmmean_cor_ELI_max_Tmaxbar_min_vec, mmmean_cor_ELI_max_Tmaxbar_min[x,y,1])
      mmmean_cor_ELI_rsds_Tmaxbar_min_vec <- c(mmmean_cor_ELI_rsds_Tmaxbar_min_vec, mmmean_cor_ELI_rsds_Tmaxbar_min[x,y,1])
      mmmean_pcor_ELI_Tmaxbar_min_vec <- c(mmmean_pcor_ELI_Tmaxbar_min_vec, mmmean_pcor_ELI_Tmaxbar_min[x,y,1])
      mmmean_cor_EF_Tmaxbar_min_vec <- c(mmmean_cor_EF_Tmaxbar_min_vec, mmmean_cor_EF_Tmaxbar_min[x,y,1])
      mmmean_cor_ELI_EF_vec <- c(mmmean_cor_ELI_EF_vec, mmmean_cor_ELI_EF[x,y,1])
      mmmean_cor_ELI_max_EF_vec <- c(mmmean_cor_ELI_max_EF_vec, mmmean_cor_ELI_max_EF[x,y,1])
      for(source_id in 1:length(cmip6_data.df$source_id)){
        prior_dcorr[source_id] <- mean(dcorr.list[[source_id]][x,y,1:3])
        prior_dcorr_max[source_id] <- mean(dcorr_max.list[[source_id]][x,y,1:3])
        prior_corr_rgy_veg[source_id] <- mean(corr_rgy_veg.list[[source_id]][x,y,])
        prior_corr_wtr_veg[source_id] <- mean(corr_wtr_veg.list[[source_id]][x,y,])
      }
      cor_ELI_Tmaxbar_min_source_id_vec <- c(cor_ELI_Tmaxbar_min_source_id_vec, cor_ELI_Tmaxbar_min_source_id[x,y,,1])
      cor_ELI_max_Tmaxbar_min_source_id_vec <- c(cor_ELI_max_Tmaxbar_min_source_id_vec, cor_ELI_max_Tmaxbar_min_source_id[x,y,,1])
      sign_cor_ELI_Tmaxbar_min_source_id_vec <- c(sign_cor_ELI_Tmaxbar_min_source_id_vec, cor_ELI_Tmaxbar_min_source_id[x,y,,2])
      sign_cor_ELI_max_Tmaxbar_min_source_id_vec <- c(sign_cor_ELI_max_Tmaxbar_min_source_id_vec, cor_ELI_max_Tmaxbar_min_source_id[x,y,,2])
      sign_cor_ELI_rsds_Tmaxbar_min_source_id_vec <- c(sign_cor_ELI_rsds_Tmaxbar_min_source_id_vec, cor_ELI_rsds_Tmaxbar_min_source_id[x,y,,2])
      prior_dcorr_h3m_vec <- c(prior_dcorr_h3m_vec, prior_dcorr)
      kendall_dcorr_h3m_vec <- c(kendall_dcorr_h3m_vec, kendall_dcorr[x,y,,1])
      prior_dcorr_max_h3m_vec <- c(prior_dcorr_max_h3m_vec, prior_dcorr_max)
      kendall_dcorr_max_h3m_vec <- c(kendall_dcorr_max_h3m_vec, kendall_dcorr_max[x,y,,1])
      prior_corr_rgy_veg_h3m_vec <- c(prior_corr_rgy_veg_h3m_vec, prior_corr_rgy_veg)
      kendall_corr_rgy_veg_h3m_vec <- c(kendall_corr_rgy_veg_h3m_vec, kendall_corr_rgy_veg[x,y,,1])
      prior_corr_wtr_veg_h3m_vec <- c(prior_corr_wtr_veg_h3m_vec, prior_corr_wtr_veg)
      kendall_corr_wtr_veg_h3m_vec <- c(kendall_corr_wtr_veg_h3m_vec, kendall_corr_wtr_veg[x,y,,1])
      kendall_av_tas_h3m_vec <- c(kendall_av_tas_h3m_vec,kendall_av_tas[x,y,,1])
      kendall_av_EF_h3m_vec <- c(kendall_av_EF_h3m_vec,kendall_av_EF[x,y,,1])
      kendall_Tmaxbar_min_h3m_vec <- c(kendall_Tmaxbar_min_h3m_vec,kendall_Tmaxbar_min[x,y,,1])
      kendall_av_rsds_h3m_vec <- c(kendall_av_rsds_h3m_vec, kendall_av_rsds[x,y,,1])
      sign_kendall_dcorr_vec <- c(sign_kendall_dcorr_vec, kendall_dcorr[x,y,,2])
      sign_kendall_dcorr_max_vec <- c(sign_kendall_dcorr_max_vec, kendall_dcorr_max[x,y,,2])
      sign_kendall_Tmaxbar_min_vec <- c(sign_kendall_Tmaxbar_min_vec, kendall_Tmaxbar_min[x,y,,2])
      sign_kendall_av_rsds_vec <- c(sign_kendall_av_rsds_vec, kendall_av_rsds[x,y,,2])
      sign_mmmean_kendall_dcorr_h3m_vec <- c(sign_mmmean_kendall_dcorr_h3m_vec, mmmean_kendall_dcorr[x,y,3])
      sign_mmmean_kendall_dcorr_max_h3m_vec <- c(sign_mmmean_kendall_dcorr_max_h3m_vec, mmmean_kendall_dcorr_max[x,y,3])
      sign_mmmean_kendall_Tmaxbar_min_h3m_vec <- c(sign_mmmean_kendall_Tmaxbar_min_h3m_vec, mmmean_kendall_Tmaxbar_min[x,y,3])
      tas_cum_vec <- c(tas_cum_vec, tas_cum.array[x,y])
    }
  }
}

h3m_trends.df <- data.frame("mmmean_prior_dcorr_h3m" = mmmean_prior_dcorr_h3m_vec,
                            "mmmean_prior_dcorr_max_h3m" = mmmean_prior_dcorr_max_h3m_vec,
                            "mmmean_prior_EF_h3m" = mmmean_prior_EF_h3m_vec,
                            "mmmean_kendall_dcorr_h3m" = mmmean_kendall_dcorr_h3m_vec,
                            "mmmean_kendall_dcorr_max_h3m" = mmmean_kendall_dcorr_max_h3m_vec,
                            "mmmean_prior_corr_rgy_veg_h3m" = mmmean_prior_corr_rgy_veg_h3m_vec,
                            "mmmean_kendall_corr_rgy_veg_h3m" = mmmean_kendall_corr_rgy_veg_h3m_vec,
                            "mmmean_prior_corr_wtr_veg_h3m" = mmmean_prior_corr_wtr_veg_h3m_vec,
                            "mmmean_kendall_corr_wtr_veg_h3m" = mmmean_kendall_corr_wtr_veg_h3m_vec,
                            "mmmean_kendall_av_tas_h3m" = mmmean_kendall_av_tas_h3m_vec,
                            "mmmean_kendall_av_EF_h3m" = mmmean_kendall_av_EF_h3m_vec,
                            "mmmean_kendall_av_rsds_h3m" = mmmean_kendall_av_rsds_h3m_vec,
                            "mmmean_kendall_Tmaxbar_min_h3m" = mmmean_kendall_Tmaxbar_min_h3m_vec,
                            "mmmean_cor_ELI_Tmaxbar_min" = mmmean_cor_ELI_Tmaxbar_min_vec,
                            "mmmean_cor_ELI_max_Tmaxbar_min" = mmmean_cor_ELI_max_Tmaxbar_min_vec,
                            "mmmean_cor_ELI_rsds_Tmaxbar_min" = mmmean_cor_ELI_Tmaxbar_min_vec,
                            "mmmean_pcor_ELI_Tmaxbar_min" = mmmean_pcor_ELI_Tmaxbar_min_vec,
                            "mmmean_cor_EF_Tmaxbar_min" = mmmean_cor_EF_Tmaxbar_min_vec,
                            "mmmean_cor_ELI_EF" = mmmean_cor_ELI_EF_vec,
                            "mmmean_cor_ELI_max_EF" = mmmean_cor_ELI_max_EF_vec,
                            "tas_cum" = tas_cum_vec,
                            "lon" = lon_vec,
                            "lat" = lat_vec)
h3m_trends.df$tas_cum <- factor(h3m_trends.df$tas_cum, levels = c("0","1","2","3","4","5","6","7","8","9","10","11","12"))


h3m_trends_source_id.df <- data.frame("prior_dcorr_h3m" = prior_dcorr_h3m_vec,
                                      "kendall_dcorr_h3m" = kendall_dcorr_h3m_vec,
                                      "prior_dcorr_max_h3m" = prior_dcorr_max_h3m_vec,
                                      "kendall_dcorr_max_h3m" = kendall_dcorr_max_h3m_vec,
                                      "prior_corr_rgy_veg_h3m" = prior_corr_rgy_veg_h3m_vec,
                                      "kendall_corr_rgy_veg_h3m" = kendall_corr_rgy_veg_h3m_vec,
                                      "prior_corr_wtr_veg_h3m" = prior_corr_wtr_veg_h3m_vec,
                                      "kendall_corr_wtr_veg_h3m" = kendall_corr_wtr_veg_h3m_vec,
                                      "kendall_av_tas_h3m" = kendall_av_tas_h3m_vec,
                                      "kendall_av_EF_h3m" = kendall_av_EF_h3m_vec,
                                      "kendall_av_rsds_h3m" = kendall_av_rsds_h3m_vec,
                                      "kendall_Tmaxbar_min_h3m" = kendall_Tmaxbar_min_h3m_vec,
                                      "sign_kendall_dcorr" = sign_kendall_dcorr_vec,
                                      "sign_kendall_dcorr_max" = sign_kendall_dcorr_max_vec,
                                      "sign_kendall_Tmaxbar_min" = sign_kendall_Tmaxbar_min_vec,
                                      "sign_kendall_av_rsds" = sign_kendall_av_rsds_vec,
                                      "cor_ELI_Tmaxbar_min_source_id" = cor_ELI_Tmaxbar_min_source_id_vec,
                                      "cor_ELI_max_Tmaxbar_min_source_id" = cor_ELI_max_Tmaxbar_min_source_id_vec,
                                      "sign_cor_ELI_max_Tmaxbar_min_source_id" = sign_cor_ELI_max_Tmaxbar_min_source_id_vec,
                                      "lon" = rep(lon_vec, each = 12),
                                      "lat" = rep(lat_vec, each = 12),
                                      "source_id" = rep(cmip6_data.df$source_id,length(lon_vec)))

cols_EF <- rev(colorRamps::matlab.like(n = 6))
myvalues_EF <- c(seq(0,.5,.1),Inf)
dbar_EF <- plot_discrete_cbar(myvalues_EF,
                              colors = cols_EF,
                              spacing = "constant",
                              font_size = 6,
                              spacing_scaling = 2,
                              width = .2,
                              triangle_size = .175)
h3m_trends.df$cuts_mmmean_prior_EF_h3m <- cut(h3m_trends.df$mmmean_prior_EF_h3m, myvalues_EF, include.lowest = T)

cols_EF_h3m_trend <- brewer.pal(7,"RdYlBu")
myvalues_EF_h3m_trend <- c(-Inf,round(seq(-.015,.01,.005),3),Inf)
dbar_EF_h3m_trend <- plot_discrete_cbar(myvalues_EF_h3m_trend,
                                        colors = cols_EF_h3m_trend,
                                        legend_title = "(-/10yr)",
                                        spacing = "constant",
                                        font_size = 6,
                                        spacing_scaling = 2,
                                        width = .2,
                                        triangle_size = .175)
h3m_trends.df$cuts_mmmean_kendall_av_EF_h3m <- cut(h3m_trends.df$mmmean_kendall_av_EF_h3m, myvalues_EF_h3m_trend, include.lowest = T)

cols_rsds_h3m_trend <- rev(brewer.pal(6,"RdYlBu"))
myvalues_rsds_h3m_trend <- c(-Inf,round(seq(-2,2,1),3),Inf)
dbar_rsds_h3m_trend <- plot_discrete_cbar(myvalues_rsds_h3m_trend,
                                          colors = cols_rsds_h3m_trend,
                                          legend_title = "(W/m2/10yr)",
                                          spacing = "constant",
                                          font_size = 6,
                                          spacing_scaling = 2,
                                          width = .2,
                                          triangle_size = .175)
h3m_trends.df$cuts_mmmean_kendall_av_rsds_h3m <- cut(h3m_trends.df$mmmean_kendall_av_rsds_h3m, myvalues_rsds_h3m_trend, include.lowest = T)
h3m_trends_source_id.df$cuts_kendall_av_rsds_h3m <- cut(h3m_trends_source_id.df$kendall_av_rsds_h3m, myvalues_rsds_h3m_trend, include.lowest = T)

cols_Tmaxbar_min_h3m_trend <- rev(brewer.pal(9,"RdYlBu")[1:8])
myvalues_Tmaxbar_min_h3m_trend <- c(-Inf,round(seq(-.1,.2,.05),2),Inf)
dbar_Tmaxbar_min_h3m_trend <- plot_discrete_cbar(myvalues_Tmaxbar_min_h3m_trend,
                                                 colors = cols_Tmaxbar_min_h3m_trend,
                                                 legend_title = expression("(K/10yr)"),
                                                 spacing = "constant",
                                                 font_size = 6,
                                                 spacing_scaling = 2,
                                                 width = .2,
                                                 triangle_size = .175)
h3m_trends.df$cuts_mmmean_kendall_Tmaxbar_min_h3m <- cut(h3m_trends.df$mmmean_kendall_Tmaxbar_min_h3m, myvalues_Tmaxbar_min_h3m_trend, include.lowest = T)
h3m_trends_source_id.df$cuts_kendall_Tmaxbar_min_h3m <- cut(h3m_trends_source_id.df$kendall_Tmaxbar_min_h3m, myvalues_Tmaxbar_min_h3m_trend, include.lowest = T)

dcorrcol <- rev(brewer.pal(11,"BrBG"))
myvalues_prior_dcorr <- c(-Inf,-.8,-.6,-.4,-.2,-.0001,.0001,.2,.4,.6,.8,Inf)
h3m_trends.df$cuts_mmmean_prior_dcorr_h3m <- cut(h3m_trends.df$mmmean_prior_dcorr_h3m, myvalues_prior_dcorr, include.lowest = T)
h3m_trends_source_id.df$cuts_prior_dcorr_h3m <- cut(h3m_trends_source_id.df$prior_dcorr_h3m, myvalues_prior_dcorr, include.lowest = T)
h3m_trends.df$cuts_mmmean_prior_dcorr_max_h3m <- cut(h3m_trends.df$mmmean_prior_dcorr_max_h3m, myvalues_prior_dcorr, include.lowest = T)
h3m_trends_source_id.df$cuts_prior_dcorr_max_h3m <- cut(h3m_trends_source_id.df$prior_dcorr_max_h3m, myvalues_prior_dcorr, include.lowest = T)

dcorr_trendcol <- rev(brewer.pal(9,"RdYlBu"))[2:9]
myvalues_dcorr_trend <- c(-Inf,seq(-.03,.06,.015),Inf)
# plot the discrete colorbar for mean dcorr
dbar_dcorr_trend <- plot_discrete_cbar(breaks = c(myvalues_dcorr_trend),
                                       colors = c(dcorr_trendcol),
                                       legend_title = "(-/10yr)",
                                       spacing = "constant",
                                       font_size = 6,
                                       spacing_scaling = 2,
                                       width = .2,
                                       triangle_size = .175)
h3m_trends.df$cuts_mmmean_dcorr_trend_h3m <- cut(h3m_trends.df$mmmean_kendall_dcorr_h3m, myvalues_dcorr_trend, include.lowest = T)
h3m_trends_source_id.df$cuts_dcorr_trend_h3m <- cut(h3m_trends_source_id.df$kendall_dcorr_h3m, myvalues_dcorr_trend, include.lowest = T)
h3m_trends.df$cuts_mmmean_dcorr_max_trend_h3m <- cut(h3m_trends.df$mmmean_kendall_dcorr_max_h3m, myvalues_dcorr_trend, include.lowest = T)
h3m_trends_source_id.df$cuts_dcorr_max_trend_h3m <- cut(h3m_trends_source_id.df$kendall_dcorr_max_h3m, myvalues_dcorr_trend, include.lowest = T)

cols_cor_ELI_Tmaxbar_min <- c(brewer.pal(9,"PRGn")[1:4],brewer.pal(9,"PRGn")[6:9])
myvalues_cor_ELI_Tmaxbar_min <- round(c(-Inf,seq(-.45,.45,.15),Inf),2)
dbar_cor_ELI_Tmaxbar_min <- plot_discrete_cbar(myvalues_cor_ELI_Tmaxbar_min,
                                               colors = cols_cor_ELI_Tmaxbar_min,
                                               spacing = "constant",
                                               font_size = 6,
                                               spacing_scaling = 2,
                                               width = .2,
                                               triangle_size = .175)
h3m_trends.df$cuts_mmmean_cor_ELI_Tmaxbar_min <- cut(h3m_trends.df$mmmean_cor_ELI_Tmaxbar_min, myvalues_cor_ELI_Tmaxbar_min, include.lowest = T)
h3m_trends.df$cuts_mmmean_cor_ELI_max_Tmaxbar_min <- cut(h3m_trends.df$mmmean_cor_ELI_max_Tmaxbar_min, myvalues_cor_ELI_Tmaxbar_min, include.lowest = T)
h3m_trends.df$cuts_mmmean_cor_ELI_rsds_Tmaxbar_min <- cut(h3m_trends.df$mmmean_cor_ELI_rsds_Tmaxbar_min, myvalues_cor_ELI_Tmaxbar_min, include.lowest = T)
h3m_trends.df$cuts_mmmean_pcor_ELI_Tmaxbar_min <- cut(h3m_trends.df$mmmean_pcor_ELI_Tmaxbar_min, myvalues_cor_ELI_Tmaxbar_min, include.lowest = T)
h3m_trends_source_id.df$cuts_cor_ELI_Tmaxbar_min_source_id <- cut(h3m_trends_source_id.df$cor_ELI_Tmaxbar_min_source_id, myvalues_cor_ELI_Tmaxbar_min, include.lowest = T)
h3m_trends_source_id.df$cuts_cor_ELI_max_Tmaxbar_min_source_id <- cut(h3m_trends_source_id.df$cor_ELI_max_Tmaxbar_min_source_id, myvalues_cor_ELI_Tmaxbar_min, include.lowest = T)

cols_cor_ELI_EF <- rev(c(brewer.pal(9,"PRGn")[1:4],brewer.pal(9,"PRGn")[6:9]))
myvalues_cor_ELI_EF <- round(c(-Inf,seq(-.45,.45,.15),Inf),2)
dbar_cor_ELI_EF <- plot_discrete_cbar(myvalues_cor_ELI_EF,
                                      colors = cols_cor_ELI_EF,
                                      spacing = "constant",
                                      font_size = 6,
                                      spacing_scaling = 2,
                                      width = .2,
                                      triangle_size = .175)
h3m_trends.df$cuts_mmmean_cor_ELI_EF <- cut(h3m_trends.df$mmmean_cor_ELI_EF, myvalues_cor_ELI_EF, include.lowest = T)
h3m_trends.df$cuts_mmmean_cor_ELI_max_EF <- cut(h3m_trends.df$mmmean_cor_ELI_max_EF, myvalues_cor_ELI_EF, include.lowest = T)

cols_cor_EF_Tmaxbar_min <- rev(c(brewer.pal(9,"PRGn")[1:4],brewer.pal(9,"PRGn")[6:9]))
myvalues_cor_EF_Tmaxbar_min <- round(c(-Inf,seq(-.45,.45,.15),Inf),2)
dbar_cor_EF_Tmaxbar_min <- plot_discrete_cbar(myvalues_cor_EF_Tmaxbar_min,
                                              colors = cols_cor_EF_Tmaxbar_min,
                                              spacing = "constant",
                                              font_size = 6,
                                              spacing_scaling = 2,
                                              width = .2,
                                              triangle_size = .175)
h3m_trends.df$cuts_mmmean_cor_EF_Tmaxbar_min <- cut(h3m_trends.df$mmmean_cor_EF_Tmaxbar_min, myvalues_cor_EF_Tmaxbar_min, include.lowest = T)

# This is a data set from the maptools package
data(wrld_simpl)

# Create a data.frame object for ggplot. ggplot requires a data frame.
mymap <- fortify(wrld_simpl)

a <- ggplot(h3m_trends.df, aes(x=lon,y=lat,fill=cuts_mmmean_prior_dcorr_max_h3m)) +
  geom_tile() +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  scale_fill_manual(values = dcorrcol) + # watch how many unique col classes there are in dcorr.df$cuts
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("a) mean Ecosystem Limitation Index, 1980 - 2010")
a

# plot the discrete colorbar for mean dcorr
dbar <- plot_discrete_cbar(base::c(myvalues_prior_dcorr[1:5],0,myvalues_prior_dcorr[8:12]),
                           colors = base::c(dcorrcol[1:5],dcorrcol[7:11]),
                           spacing = "constant",
                           font_size = 6,
                           spacing_scaling = 2,
                           width = .2,
                           triangle_size = .175)

# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar <- dbar + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_smaller <- grid.arrange(empty, dbar, empty , ncol=3, widths = c(1,9,1))

plot_prior_dcorr <- grid.arrange(a,dbar_smaller, nrow = 2, heights = c(.8,.2))

a <- ggplot(h3m_trends_source_id.df[which(!is.na(h3m_trends_source_id.df$prior_dcorr_max_h3m)),], aes(x=lon,y=lat,fill=cuts_prior_dcorr_max_h3m)) +
  geom_tile() +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = dcorrcol) + # watch how many unique col classes there are in dcorr.df$cuts
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  facet_wrap(~source_id, ncol = 3) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24, hjust = 0.5),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("mean Ecosystem Limitation Index 1980 - 2010")
a

# plot the discrete colorbar for mean dcorr
dbar <- plot_discrete_cbar(breaks = c(myvalues_prior_dcorr[1:5],0,myvalues_prior_dcorr[8:12]),
                           colors = c(dcorrcol[1:5],dcorrcol[7:11]),
                           spacing = "constant",
                           font_size = 6,
                           spacing_scaling = 2,
                           width = .2,
                           triangle_size = .175)

# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar <- dbar + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_smaller <- grid.arrange(empty, dbar, empty , ncol=3, widths = c(1,6,1))

plot <- grid.arrange(a,dbar_smaller, nrow = 2, heights = c(.8,.2))
# ggsave("Figures/SFig9.png", plot = plot, width = 9*3/1.5, height = 9*3/1.5, units = "in")
ggsave("testdir/SFig9.png", plot = plot, width = 9*3/1.5, height = 9*3/1.5, units = "in")


a <- ggplot(h3m_trends.df, aes(x=lon,y=lat,fill=tas_cum)) +
  geom_tile() +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual("",
                    values = colorRampPalette(brewer.pal(9, "Oranges"))(13)) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle(expression("Sum of models for which cor(T"[a]*"',ET') > cor(SW"["in"]*"',ET')"))
a

# ggsave("Figures/SFig1.png", plot = a, width = 9, height = 6, units = "in")
ggsave("testdir/SFig1.png", plot = a, width = 9, height = 6, units = "in")

# agreement on sign of ELI
# this is the agreement between models
agreement.array <- 
  agreement_max.array <- 
  array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(kendall_dcorr[x,y,,1]) == 1)) > 8 | length(which(sign(kendall_dcorr[x,y,,1]) == -1)) > 8){
          agreement.array[x,y] <- 1
        }
        if(length(which(sign(kendall_dcorr_max[x,y,,1]) == 1)) > 8 | length(which(sign(kendall_dcorr_max[x,y,,1]) == -1)) > 8){
          agreement_max.array[x,y] <- 1
        }
      }
    }
  }
}

agreement.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                         c("agreement","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement.array[x,y])){
      agreement.df <- rbind(agreement.df,
                            data.frame("agreement" = agreement.array[x,y],
                                       "lon" = lon[x],
                                       "lat" = lat[y]))
    }
  }
}

agreement_max.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                             c("agreement_max","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement_max.array[x,y])){
      agreement_max.df <- rbind(agreement_max.df,
                                data.frame("agreement_max" = agreement_max.array[x,y],
                                           "lon" = lon[x],
                                           "lat" = lat[y]))
    }
  }
}
count_mmtrends.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_kendall_dcorr[,,1] > 0)])/total_land_area, sum(area.array[which(mmmean_kendall_dcorr[,,1] < 0)])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_kendall_dcorr[,,1] > 0 & !is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_kendall_dcorr[,,1] < 0 & !is.na(agreement.array))])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_kendall_dcorr[,,1] > 0 & is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_kendall_dcorr[,,1] < 0 & is.na(agreement.array))])/total_land_area)),
                                "trend" = rep(c("drying","wettening"),3),
                                "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends.df$trend_sign <- paste0(count_mmtrends.df$trend, "_",count_mmtrends.df$sign)
count_mmtrends.df$trend_f <- factor(count_mmtrends.df$trend, levels = c("wettening","drying"))

count_mmtrends_max.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_kendall_dcorr_max[,,1] > 0)])/total_land_area, sum(area.array[which(mmmean_kendall_dcorr_max[,,1] < 0)])/total_land_area),
                                               100*c(sum(area.array[which(mmmean_kendall_dcorr_max[,,1] > 0 & !is.na(agreement_max.array))])/total_land_area, sum(area.array[which(mmmean_kendall_dcorr_max[,,1] < 0 & !is.na(agreement_max.array))])/total_land_area),
                                               100*c(sum(area.array[which(mmmean_kendall_dcorr_max[,,1] > 0 & is.na(agreement_max.array))])/total_land_area, sum(area.array[which(mmmean_kendall_dcorr_max[,,1] < 0 & is.na(agreement_max.array))])/total_land_area)),
                                    "trend" = rep(c("drying","wettening"),3),
                                    "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends_max.df$trend_sign <- paste0(count_mmtrends_max.df$trend, "_",count_mmtrends_max.df$sign)
count_mmtrends_max.df$trend_f <- factor(count_mmtrends_max.df$trend, levels = c("wettening","drying"))


a <- ggplot(h3m_trends.df, aes(x=lon,y=lat,fill=cuts_mmmean_dcorr_max_trend_h3m)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = agreement_max.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = dcorr_trendcol) + # watch how many unique col classes there are in dcorr.df$cuts
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("b) Ecosystem Limitation Index trend, 1980 - 2100")
a
a <- ggplotGrob(a)

bar_mmtrends <- ggplot(count_mmtrends_max.df[which(count_mmtrends_max.df$sign != 'all'),], aes(x = trend_f, y = area, group = trend_f, fill = trend_sign)) +
  geom_bar(stat='identity',col = 'black', width = .5, position = "stack") +
  scale_x_discrete("") +
  geom_text(inherit.aes = F, data = count_mmtrends_max.df[which(count_mmtrends_max.df$sign == 'all'),], aes(x=trend_f,y=area-5,label=round(area,2)), size = 4) + 
  scale_y_continuous(expression(paste("land area-%")), expand = c(0,0)) +
  scale_fill_manual(values = c("drying_agree" = dcorr_trendcol[8],
                               "drying_disagree" = alpha(dcorr_trendcol[8], .25),
                               "wettening_disagree" = alpha(dcorr_trendcol[1], .25),
                               "wettening_agree" = dcorr_trendcol[1]),
                    drop = F) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=16),
        plot.title = element_text(size=10)
  )
bar_mmtrends
bar_mmtrends <- ggplotGrob(bar_mmtrends)

gt1 <- gtable(widths = unit(rep(2,32), c("null")), heights = unit(rep(2,32), "null"))
gt1 <- gtable_add_grob(gt1, a, t=1, b=32, l=1, r=32)
gt1 <- gtable_add_grob(gt1, bar_mmtrends, t = 29, l = 4, b = 20, r = 10)
grid.draw(gt1)

# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_dcorr_trend <- dbar_dcorr_trend + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_dcorr_trend_smaller <- grid.arrange(empty, dbar_dcorr_trend, empty , ncol=3, widths = c(1,6,1))

plot_kendall_dcorr <- grid.arrange(gt1,dbar_dcorr_trend_smaller, nrow = 2, heights = c(.8,.2))

sign_dcorr_max.df <- h3m_trends_source_id.df[which(h3m_trends_source_id.df$sign_kendall_dcorr_max < .05),]
f <- ggplot(h3m_trends_source_id.df[which(!is.na(h3m_trends_source_id.df$kendall_dcorr_max_h3m)),], aes(x=lon,y=lat,fill=cuts_dcorr_max_trend_h3m)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = sign_dcorr_max.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = dcorr_trendcol,
                    drop = F) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  facet_wrap(~source_id, ncol = 3) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24, hjust = 0.5),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("Ecosystem Limitation Index trend, 1980 - 2100")
f

# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_dcorr_trend <- dbar_dcorr_trend + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_dcorr_trend_smaller <- grid.arrange(empty, dbar_dcorr_trend, empty , ncol=3, widths = c(1,4,1))

plot <- grid.arrange(f,dbar_dcorr_trend_smaller, nrow = 2, heights = c(.8,.2))

# ggsave("Figures/SFig5.png", plot = plot, width = 9*3/1.5, height = 9*3/1.5, units = "in")
ggsave("testdir/SFig5.png", plot = plot, width = 9*3/1.5, height = 9*3/1.5, units = "in")

names.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                     c("lon","lat","label"))
reg <- c("SAM", "NAM", "CEU", "NAS")
jump_four <- seq(1,20,4)
for(i in 1:4){
  names.df <- rbind(names.df, data.frame("lon" = (min(hotspot_regs.df$x[jump_four[i]:(jump_four[i]+3)]) + max(hotspot_regs.df$x[jump_four[i]:(jump_four[i]+3)]))/2,
                                         "lat" = max(hotspot_regs.df$y[jump_four[i]:(jump_four[i]+3)]),
                                         "label" = reg[i]))
}

agreement.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(kendall_Tmaxbar_min[x,y,,1]) == 1)) > 8 | length(which(sign(kendall_Tmaxbar_min[x,y,,1]) == -1)) > 8){
          agreement.array[x,y] <- 1
        }
      }
    }
  }
}

agreement.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                         c("agreement","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement.array[x,y])){
      agreement.df <- rbind(agreement.df,
                            data.frame("agreement" = agreement.array[x,y],
                                       "lon" = lon[x],
                                       "lat" = lat[y]))
    }
  }
}

count_mmtrends.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_kendall_Tmaxbar_min[,,1] > 0)])/total_land_area, sum(area.array[which(mmmean_kendall_Tmaxbar_min[,,1] < 0)])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_kendall_Tmaxbar_min[,,1] > 0 & !is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_kendall_Tmaxbar_min[,,1] < 0 & !is.na(agreement.array))])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_kendall_Tmaxbar_min[,,1] > 0 & is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_kendall_Tmaxbar_min[,,1] < 0 & is.na(agreement.array))])/total_land_area)),
                                "trend" = rep(c("warming","cooling"),3),
                                "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends.df$trend_sign <- paste0(count_mmtrends.df$trend, "_",count_mmtrends.df$sign)

f <- ggplot(h3m_trends.df, aes(x=lon,y=lat,fill=cuts_mmmean_kendall_Tmaxbar_min_h3m)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = agreement.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_text(inherit.aes = F, data = names.df, aes(x=lon,y=lat+3.5,label=label), size = 5, col = 1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = cols_Tmaxbar_min_h3m_trend,
                    drop = F) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle(expression("a) Temperature excess trend, 1980 - 2100"))
f
f <- ggplotGrob(f)

bar_mmtrends <- ggplot(count_mmtrends.df[which(count_mmtrends.df$sign != 'all'),], aes(x = trend, y = area, group = trend, fill = trend_sign)) +
  geom_bar(stat='identity',col = 'black', width = .5, position = "stack") +
  scale_x_discrete("") +
  geom_text(inherit.aes = F, data = count_mmtrends.df[which(count_mmtrends.df$sign == 'all'),], aes(x=trend,y=area-5,label=round(area,2)), size = 4) + 
  scale_y_continuous(expression(paste("land area-%")), expand = c(0,0)) +
  scale_fill_manual(values = c("warming_agree" = cols_Tmaxbar_min_h3m_trend[8],
                               "warming_disagree" = alpha(cols_Tmaxbar_min_h3m_trend[8], .25),
                               "cooling_disagree" = alpha(cols_Tmaxbar_min_h3m_trend[1], .25),
                               "cooling_agree" = cols_Tmaxbar_min_h3m_trend[1]),
                    drop = F) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=16),
        plot.title = element_text(size=10)
  )
bar_mmtrends
bar_mmtrends <- ggplotGrob(bar_mmtrends)

gt1 <- gtable(widths = unit(rep(2,32), c("null")), heights = unit(rep(2,32), "null"))
gt1 <- gtable_add_grob(gt1, f, t=1, b=32, l=1, r=32)
gt1 <- gtable_add_grob(gt1, bar_mmtrends, t = 29, l = 4, b = 20, r = 10)
grid.draw(gt1)


# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_Tmaxbar_min_h3m_trend <- dbar_Tmaxbar_min_h3m_trend + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_Tmaxbar_min_h3m_trend_smaller <- grid.arrange(empty, dbar_Tmaxbar_min_h3m_trend, empty , ncol=3, widths = c(1,6,1))

plot_kendall_Tmaxbar_min <- grid.arrange(gt1,dbar_Tmaxbar_min_h3m_trend_smaller, nrow = 2, heights = c(.8,.2))


sign_Tmaxbar_min.df <- h3m_trends_source_id.df[which(h3m_trends_source_id.df$sign_kendall_Tmaxbar_min < .05),]
f <- ggplot(h3m_trends_source_id.df[which(!is.na(h3m_trends_source_id.df$kendall_Tmaxbar_min_h3m)),], aes(x=lon,y=lat,fill=cuts_kendall_Tmaxbar_min_h3m)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = sign_Tmaxbar_min.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = cols_Tmaxbar_min_h3m_trend,
                    drop = F) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  facet_wrap(~source_id, ncol = 3) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24, hjust = 0.5),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("Temperature excess trend, 1980 - 2100")
f

# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_Tmaxbar_min_h3m_trend <- dbar_Tmaxbar_min_h3m_trend + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_Tmaxbar_min_h3m_trend_smaller <- grid.arrange(empty, dbar_Tmaxbar_min_h3m_trend, empty , ncol=3, widths = c(1,4,1))

plot <- grid.arrange(f,dbar_Tmaxbar_min_h3m_trend_smaller, nrow = 2, heights = c(.8,.2))

# ggsave("Figures/SFig3.png", plot = plot, width = 9*3/1.5, height = 9*3/1.5, units = "in")
ggsave("testdir/SFig3.png", plot = plot, width = 9*3/1.5, height = 9*3/1.5, units = "in")









agreement.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(kendall_av_rsds[x,y,,1]) == 1)) > 8 | length(which(sign(kendall_av_rsds[x,y,,1]) == -1)) > 8){
          agreement.array[x,y] <- 1
        }
      }
    }
  }
}
agreement.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                         c("agreement","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement.array[x,y])){
      agreement.df <- rbind(agreement.df,
                            data.frame("agreement" = agreement.array[x,y],
                                       "lon" = lon[x],
                                       "lat" = lat[y]))
    }
  }
}

count_mmtrends.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_kendall_av_rsds > 0)])/total_land_area, sum(area.array[which(mmmean_kendall_av_rsds < 0)])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_kendall_av_rsds > 0 & !is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_kendall_av_rsds < 0 & !is.na(agreement.array))])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_kendall_av_rsds > 0 & is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_kendall_av_rsds < 0 & is.na(agreement.array))])/total_land_area)),
                                "trend" = rep(c("more","less"),3),
                                "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends.df$trend_sign <- paste0(count_mmtrends.df$trend, "_",count_mmtrends.df$sign)

f <- ggplot(h3m_trends.df, aes(x=lon,y=lat,fill=cuts_mmmean_kendall_av_rsds_h3m)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = agreement.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = cols_rsds_h3m_trend,
                    drop = F) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle(expression("Incoming shortwave radiation trend, 1980 - 2100"))
f
f <- ggplotGrob(f)

bar_mmtrends <- ggplot(count_mmtrends.df[which(count_mmtrends.df$sign != 'all'),], aes(x = trend, y = area, group = trend, fill = trend_sign)) +
  geom_bar(stat='identity',col = 'black', width = .5, position = "stack") +
  scale_x_discrete("") +
  geom_text(inherit.aes = F, data = count_mmtrends.df[which(count_mmtrends.df$sign == 'all'),], aes(x=trend,y=area-5,label=round(area,2)), size = 4) + 
  scale_y_continuous(expression(paste("land area-%")), expand = c(0,0)) +
  scale_fill_manual(values = c("more_agree" = cols_rsds_h3m_trend[6],
                               "more_disagree" = alpha(cols_rsds_h3m_trend[6], .25),
                               "less_disagree" = alpha(cols_rsds_h3m_trend[1], .25),
                               "less_agree" = cols_rsds_h3m_trend[1]),
                    drop = F) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=16),
        plot.title = element_text(size=10)
  )
bar_mmtrends
bar_mmtrends <- ggplotGrob(bar_mmtrends)

gt1 <- gtable(widths = unit(rep(2,32), c("null")), heights = unit(rep(2,32), "null"))
gt1 <- gtable_add_grob(gt1, f, t=1, b=32, l=1, r=32)
gt1 <- gtable_add_grob(gt1, bar_mmtrends, t = 29, l = 4, b = 20, r = 10)
grid.draw(gt1)


# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_rsds_h3m_trend <- dbar_rsds_h3m_trend + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_rsds_h3m_trend_smaller <- grid.arrange(empty, dbar_rsds_h3m_trend, empty , ncol=3, widths = c(1,6,1))

plot_kendall_rsds <- grid.arrange(gt1,dbar_rsds_h3m_trend_smaller, nrow = 2, heights = c(.8,.2))
# ggsave("Figures/SFig4.png", plot = plot_kendall_rsds, width = 9*1.3, height = 7*1.3, units = "in")
ggsave("testdir/SFig4.png", plot = plot_kendall_rsds, width = 9*1.3, height = 7*1.3, units = "in")

sign_kendall_av_rsds.df <- h3m_trends_source_id.df[which(h3m_trends_source_id.df$sign_kendall_av_rsds < .05),]
f <- ggplot(h3m_trends_source_id.df[which(!is.na(h3m_trends_source_id.df$kendall_av_rsds_h3m)),], aes(x=lon,y=lat,fill=cuts_kendall_av_rsds_h3m)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = sign_kendall_av_rsds.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = cols_rsds_h3m_trend,
                    drop = F) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  facet_wrap(~source_id, ncol = 3) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24, hjust = 0.5),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("Incoming shortwave radiation trend, 1980 - 2100")
f

# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_dcorr_trend <- dbar_dcorr_trend + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_dcorr_trend_smaller <- grid.arrange(empty, dbar_dcorr_trend, empty , ncol=3, widths = c(1,4,1))

plot <- grid.arrange(f,dbar_dcorr_trend_smaller, nrow = 2, heights = c(.8,.2))



















agreement.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(cor_ELI_Tmaxbar_min_source_id[x,y,,1]) == 1)) > 8 | length(which(sign(cor_ELI_Tmaxbar_min_source_id[x,y,,1]) == -1)) > 8){
          agreement.array[x,y] <- 1
        }
      }
    }
  }
}

agreement.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                         c("agreement","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement.array[x,y])){
      agreement.df <- rbind(agreement.df,
                            data.frame("agreement" = agreement.array[x,y],
                                       "lon" = lon[x],
                                       "lat" = lat[y]))
    }
  }
}

count_mmtrends.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] >= 0)])/total_land_area, sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] < 0)])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] >= 0 & !is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] < 0 & !is.na(agreement.array))])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] >= 0 & is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] < 0 & is.na(agreement.array))])/total_land_area)),
                                "trend" = rep(c("positive","negative"),3),
                                "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends.df$trend_sign <- paste0(count_mmtrends.df$trend, "_",count_mmtrends.df$sign)

agreement.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(cor_ELI_Tmaxbar_min_source_id[x,y,,1]) == 1)) > 8 | length(which(sign(cor_ELI_Tmaxbar_min_source_id[x,y,,1]) == -1)) > 8){
          agreement.array[x,y] <- 1
        }
      }
    }
  }
}

agreement.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                         c("agreement","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement.array[x,y])){
      agreement.df <- rbind(agreement.df,
                            data.frame("agreement" = agreement.array[x,y],
                                       "lon" = lon[x],
                                       "lat" = lat[y]))
    }
  }
}

count_mmtrends.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] >= 0)])/total_land_area, sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] < 0)])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] >= 0 & !is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] < 0 & !is.na(agreement.array))])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] >= 0 & is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_Tmaxbar_min[,,1] < 0 & is.na(agreement.array))])/total_land_area)),
                                "trend" = rep(c("positive","negative"),3),
                                "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends.df$trend_sign <- paste0(count_mmtrends.df$trend, "_",count_mmtrends.df$sign)

agreement_max.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                             c("agreement_max","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement_max.array[x,y])){
      agreement_max.df <- rbind(agreement_max.df,
                                data.frame("agreement_max" = agreement_max.array[x,y],
                                           "lon" = lon[x],
                                           "lat" = lat[y]))
    }
  }
}

count_mmtrends_max.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_cor_ELI_max_Tmaxbar_min[,,1] >= 0)])/total_land_area, sum(area.array[which(mmmean_cor_ELI_max_Tmaxbar_min[,,1] < 0)])/total_land_area),
                                               100*c(sum(area.array[which(mmmean_cor_ELI_max_Tmaxbar_min[,,1] >= 0 & !is.na(agreement_max.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_max_Tmaxbar_min[,,1] < 0 & !is.na(agreement_max.array))])/total_land_area),
                                               100*c(sum(area.array[which(mmmean_cor_ELI_max_Tmaxbar_min[,,1] >= 0 & is.na(agreement_max.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_max_Tmaxbar_min[,,1] < 0 & is.na(agreement_max.array))])/total_land_area)),
                                    "trend" = rep(c("positive","negative"),3),
                                    "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends_max.df$trend_sign <- paste0(count_mmtrends_max.df$trend, "_",count_mmtrends_max.df$sign)

f <- ggplot(h3m_trends.df, aes(x=lon,y=lat,fill=cuts_mmmean_cor_ELI_max_Tmaxbar_min)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = agreement_max.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = cols_cor_ELI_Tmaxbar_min,
                    drop = F) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle(expression("c) cor(Temperature excess,ELI), 1980 - 2100"))
f
f <- ggplotGrob(f)

bar_mmtrends <- ggplot(count_mmtrends_max.df[which(count_mmtrends_max.df$sign != 'all'),], aes(x = trend, y = area, group = trend, fill = trend_sign)) +
  geom_bar(stat='identity',col = 'black', width = .5, position = "stack") +
  scale_x_discrete("") +
  geom_text(inherit.aes = F, data = count_mmtrends_max.df[which(count_mmtrends_max.df$sign == 'all'),], aes(x=trend,y=area-5,label=round(area,2)), size = 4) + 
  scale_y_continuous(expression(paste("land area-%")), expand = c(0,0)) +
  scale_fill_manual(values = c("positive_agree" = cols_cor_ELI_Tmaxbar_min[8],
                               "positive_disagree" = alpha(cols_cor_ELI_Tmaxbar_min[8], .25),
                               "negative_disagree" = alpha(cols_cor_ELI_Tmaxbar_min[1], .25),
                               "negative_agree" = cols_cor_ELI_Tmaxbar_min[1]),
                    drop = F) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=16),
        plot.title = element_text(size=10)
  )
bar_mmtrends
bar_mmtrends <- ggplotGrob(bar_mmtrends)

gt1 <- gtable(widths = unit(rep(2,32), c("null")), heights = unit(rep(2,32), "null"))
gt1 <- gtable_add_grob(gt1, f, t=1, b=32, l=1, r=32)
gt1 <- gtable_add_grob(gt1, bar_mmtrends, t = 29, l = 4, b = 20, r = 10)
grid.draw(gt1)


# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_cor_ELI_Tmaxbar_min <- dbar_cor_ELI_Tmaxbar_min + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_cor_ELI_Tmaxbar_min_smaller <- grid.arrange(empty, dbar_cor_ELI_Tmaxbar_min, empty , ncol=3, widths = c(1,6,1))

plot_cor_ELI_Tmaxbar_min <- grid.arrange(gt1,dbar_cor_ELI_Tmaxbar_min_smaller, nrow = 2, heights = c(.8,.2))

plot_spat <- grid.arrange(plot_kendall_Tmaxbar_min, plot_kendall_dcorr, plot_cor_ELI_Tmaxbar_min)
# ggsave("Figures/Fig1.png", plot = plot_spat, width = 9*1.3, height = 21*1.3, units = "in")
ggsave("testdir/Fig1.png", plot = plot_spat, width = 9*1.3, height = 21*1.3, units = "in")





sign_cor_ELI_max_Tmaxbar_min.df <- h3m_trends_source_id.df[which(h3m_trends_source_id.df$sign_cor_ELI_max_Tmaxbar_min_source_id < .05),]
f <- ggplot(h3m_trends_source_id.df[which(!is.na(h3m_trends_source_id.df$kendall_Tmaxbar_min_h3m)),], aes(x=lon,y=lat,fill=cuts_cor_ELI_max_Tmaxbar_min_source_id)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = sign_cor_ELI_max_Tmaxbar_min.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = cols_cor_ELI_Tmaxbar_min,
                    drop = F) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  facet_wrap(~source_id, ncol = 3) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24, hjust = 0.5),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle(expression("cor(Temperature excess,ELI), 1980 - 2100"))
f

# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_cor_ELI_Tmaxbar_min <- dbar_cor_ELI_Tmaxbar_min + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_cor_ELI_Tmaxbar_min_smaller <- grid.arrange(empty, dbar_cor_ELI_Tmaxbar_min, empty , ncol=3, widths = c(1,4,1))

plot <- grid.arrange(f,dbar_cor_ELI_Tmaxbar_min_smaller, nrow = 2, heights = c(.8,.2))

# ggsave("Figures/SFig6.png", plot = plot, width = 9*3/1.5, height = 9*3/1.5, units = "in")
ggsave("testdir/SFig6.png", plot = plot, width = 9*3/1.5, height = 9*3/1.5, units = "in")


# agreement on sign of ELI
# this is the agreement between models
agreement.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(kendall_av_EF[x,y,,1]) == 1)) > 8 | length(which(sign(kendall_av_EF[x,y,,1]) == -1)) > 8){
          agreement.array[x,y] <- 1
        }
      }
    }
  }
}

agreement.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                         c("agreement","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement.array[x,y])){
      agreement.df <- rbind(agreement.df,
                            data.frame("agreement" = agreement.array[x,y],
                                       "lon" = lon[x],
                                       "lat" = lat[y]))
    }
  }
}

count_mmtrends.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_kendall_av_EF > 0)])/total_land_area, sum(area.array[which(mmmean_kendall_av_EF < 0)])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_kendall_av_EF > 0 & !is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_kendall_av_EF < 0 & !is.na(agreement.array))])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_kendall_av_EF > 0 & is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_kendall_av_EF < 0 & is.na(agreement.array))])/total_land_area)),
                                "trend" = rep(c("wettening","drying"),3),
                                "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends.df$trend_sign <- paste0(count_mmtrends.df$trend, "_",count_mmtrends.df$sign)
count_mmtrends.df$trend_f <- factor(count_mmtrends.df$trend, levels = c("drying","wettening"))


a <- ggplot(h3m_trends.df, aes(x=lon,y=lat,fill=cuts_mmmean_kendall_av_EF_h3m)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = agreement.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = cols_EF_h3m_trend) + # watch how many unique col classes there are in dcorr.df$cuts
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("a) EF trend, 1980 - 2100")
a
a <- ggplotGrob(a)

bar_mmtrends <- ggplot(count_mmtrends.df[which(count_mmtrends.df$sign != 'all'),], aes(x = trend_f, y = area, group = trend_f, fill = trend_sign)) +
  geom_bar(stat='identity',col = 'black', width = .5, position = "stack") +
  scale_x_discrete("") +
  geom_text(inherit.aes = F, data = count_mmtrends.df[which(count_mmtrends.df$sign == 'all'),], aes(x=trend_f,y=area-5,label=round(area,2)), size = 4) + 
  scale_y_continuous(expression(paste("land area-%")), expand = c(0,0)) +
  scale_fill_manual(values = c("drying_agree" = dcorr_trendcol[8],
                               "drying_disagree" = alpha(dcorr_trendcol[8], .25),
                               "wettening_disagree" = alpha(dcorr_trendcol[1], .25),
                               "wettening_agree" = dcorr_trendcol[1]),
                    drop = F) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=16),
        plot.title = element_text(size=10)
  )
bar_mmtrends
bar_mmtrends <- ggplotGrob(bar_mmtrends)

gt1 <- gtable(widths = unit(rep(2,32), c("null")), heights = unit(rep(2,32), "null"))
gt1 <- gtable_add_grob(gt1, a, t=1, b=32, l=1, r=32)
gt1 <- gtable_add_grob(gt1, bar_mmtrends, t = 29, l = 4, b = 20, r = 10)
grid.draw(gt1)

# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_EF_h3m_trend <- dbar_EF_h3m_trend + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_EF_h3m_trend_smaller <- grid.arrange(empty, dbar_EF_h3m_trend, empty , ncol=3, widths = c(1,6,1))

plot_kendall_av_EF <- grid.arrange(gt1,dbar_EF_h3m_trend_smaller, nrow = 2, heights = c(.8,.2))

agreement.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(cor_ELI_EF_source_id[x,y,,1]) == 1)) > 8 | length(which(sign(cor_ELI_EF_source_id[x,y,,1]) == -1)) > 8){
          agreement.array[x,y] <- 1
        }
      }
    }
  }
}

agreement.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                         c("agreement","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement.array[x,y])){
      agreement.df <- rbind(agreement.df,
                            data.frame("agreement" = agreement.array[x,y],
                                       "lon" = lon[x],
                                       "lat" = lat[y]))
    }
  }
}

count_mmtrends.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_cor_ELI_EF[,,1] >= 0)])/total_land_area, sum(area.array[which(mmmean_cor_ELI_EF[,,1] < 0)])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_cor_ELI_EF[,,1] >= 0 & !is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_EF[,,1] < 0 & !is.na(agreement.array))])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_cor_ELI_EF[,,1] >= 0 & is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_EF[,,1] < 0 & is.na(agreement.array))])/total_land_area)),
                                "trend" = rep(c("positive","negative"),3),
                                "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends.df$trend_sign <- paste0(count_mmtrends.df$trend, "_",count_mmtrends.df$sign)

agreement_max.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(cor_ELI_max_EF_source_id[x,y,,1]) == 1)) > 8 | length(which(sign(cor_ELI_max_EF_source_id[x,y,,1]) == -1)) > 8){
          agreement_max.array[x,y] <- 1
        }
      }
    }
  }
}

agreement_max.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                             c("agreement_max","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement_max.array[x,y])){
      agreement_max.df <- rbind(agreement_max.df,
                                data.frame("agreement_max" = agreement_max.array[x,y],
                                           "lon" = lon[x],
                                           "lat" = lat[y]))
    }
  }
}

count_mmtrends_max.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_cor_ELI_max_EF[,,1] >= 0)])/total_land_area, sum(area.array[which(mmmean_cor_ELI_max_EF[,,1] < 0)])/total_land_area),
                                               100*c(sum(area.array[which(mmmean_cor_ELI_max_EF[,,1] >= 0 & !is.na(agreement_max.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_max_EF[,,1] < 0 & !is.na(agreement_max.array))])/total_land_area),
                                               100*c(sum(area.array[which(mmmean_cor_ELI_max_EF[,,1] >= 0 & is.na(agreement_max.array))])/total_land_area, sum(area.array[which(mmmean_cor_ELI_max_EF[,,1] < 0 & is.na(agreement_max.array))])/total_land_area)),
                                    "trend" = rep(c("positive","negative"),3),
                                    "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends_max.df$trend_sign <- paste0(count_mmtrends_max.df$trend, "_",count_mmtrends_max.df$sign)


f <- ggplot(h3m_trends.df, aes(x=lon,y=lat,fill=cuts_mmmean_cor_ELI_max_EF)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = agreement_max.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = cols_cor_ELI_EF,
                    drop = F) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle(expression("c) cor(EF,ELI), 1980 - 2100"))
f
f <- ggplotGrob(f)

bar_mmtrends <- ggplot(count_mmtrends_max.df[which(count_mmtrends_max.df$sign != 'all'),], aes(x = trend, y = area, group = trend, fill = trend_sign)) +
  geom_bar(stat='identity',col = 'black', width = .5, position = "stack") +
  scale_x_discrete("") +
  geom_text(inherit.aes = F, data = count_mmtrends_max.df[which(count_mmtrends_max.df$sign == 'all'),], aes(x=trend,y=area-5,label=round(area,2)), size = 4) +
  scale_y_continuous(expression(paste("land area-%")), expand = c(0,0)) +
  scale_fill_manual(values = c("positive_agree" = cols_cor_ELI_EF[8],
                               "positive_disagree" = alpha(cols_cor_ELI_EF[8], .25),
                               "negative_disagree" = alpha(cols_cor_ELI_EF[1], .25),
                               "negative_agree" = cols_cor_ELI_EF[1]),
                    drop = F) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=16),
        plot.title = element_text(size=10)
  )
bar_mmtrends
bar_mmtrends <- ggplotGrob(bar_mmtrends)

gt1 <- gtable(widths = unit(rep(2,32), c("null")), heights = unit(rep(2,32), "null"))
gt1 <- gtable_add_grob(gt1, f, t=1, b=32, l=1, r=32)
gt1 <- gtable_add_grob(gt1, bar_mmtrends, t = 29, l = 4, b = 20, r = 10)
grid.draw(gt1)


# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_cor_ELI_EF <- dbar_cor_ELI_EF + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_cor_ELI_EF_smaller <- grid.arrange(empty, dbar_cor_ELI_EF, empty , ncol=3, widths = c(1,6,1))

plot_cor_ELI_EF <- grid.arrange(gt1,dbar_cor_ELI_EF_smaller, nrow = 2, heights = c(.8,.2))

agreement.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(cor_EF_Tmaxbar_min_source_id[x,y,,1]) == 1)) > 8 | length(which(sign(cor_EF_Tmaxbar_min_source_id[x,y,,1]) == -1)) > 8){
          agreement.array[x,y] <- 1
        }
      }
    }
  }
}

agreement.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                         c("agreement","lon","lat"))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(agreement.array[x,y])){
      agreement.df <- rbind(agreement.df,
                            data.frame("agreement" = agreement.array[x,y],
                                       "lon" = lon[x],
                                       "lat" = lat[y]))
    }
  }
}

count_mmtrends.df <- data.frame("area" = c(100*c(sum(area.array[which(mmmean_cor_EF_Tmaxbar_min[,,1] >= 0)])/total_land_area, sum(area.array[which(mmmean_cor_EF_Tmaxbar_min[,,1] < 0)])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_cor_EF_Tmaxbar_min[,,1] >= 0 & !is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_cor_EF_Tmaxbar_min[,,1] < 0 & !is.na(agreement.array))])/total_land_area),
                                           100*c(sum(area.array[which(mmmean_cor_EF_Tmaxbar_min[,,1] >= 0 & is.na(agreement.array))])/total_land_area, sum(area.array[which(mmmean_cor_EF_Tmaxbar_min[,,1] < 0 & is.na(agreement.array))])/total_land_area)),
                                "trend" = rep(c("positive","negative"),3),
                                "sign" = rep(c("all","agree","disagree"),each=2))
count_mmtrends.df$trend_sign <- paste0(count_mmtrends.df$trend, "_",count_mmtrends.df$sign)

f <- ggplot(h3m_trends.df, aes(x=lon,y=lat,fill=cuts_mmmean_cor_EF_Tmaxbar_min)) +
  geom_tile() +
  geom_point(inherit.aes = F, data = agreement.df, aes(x=lon,y=lat), size = .05) +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  # geom_text(inherit.aes = F, data = names.df, aes(x=lon,y=lat+3.5,label=label), size = 5, col = 1) +
  geom_segment(data = hotspot_regs.df, inherit.aes = F, aes(x = x, y = y, xend = xend, yend = yend), lty = 'dashed', size = 1.25) +
  scale_fill_manual(values = cols_cor_EF_Tmaxbar_min,
                    drop = F) +
  scale_x_continuous("",
                     limits=c(-180,180),
                     expand=c(0,0)) +
  scale_y_continuous("",
                     limits=c(-60,70),
                     expand=c(0,0)) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        # axis.title = element_text(size=20),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
        # ) + ggtitle(expression("b) cor("*bar("T"[max])*" divergence,EF), 1980 - 2100"))
  ) + ggtitle(expression("b) cor(Temperature excess,EF), 1980 - 2100"))
f
f <- ggplotGrob(f)

bar_mmtrends <- ggplot(count_mmtrends.df[which(count_mmtrends.df$sign != 'all'),], aes(x = trend, y = area, group = trend, fill = trend_sign)) +
  geom_bar(stat='identity',col = 'black', width = .5, position = "stack") +
  scale_x_discrete("") +
  geom_text(inherit.aes = F, data = count_mmtrends.df[which(count_mmtrends.df$sign == 'all'),], aes(x=trend,y=area-5,label=round(area,2)), size = 4) + 
  scale_y_continuous(expression(paste("land area-%")), expand = c(0,0)) +
  scale_fill_manual(values = c("positive_agree" = cols_cor_EF_Tmaxbar_min[8],
                               "positive_disagree" = alpha(cols_cor_EF_Tmaxbar_min[8], .25),
                               "negative_disagree" = alpha(cols_cor_EF_Tmaxbar_min[1], .25),
                               "negative_agree" = cols_cor_EF_Tmaxbar_min[1]),
                    drop = F) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.text = element_text(size=16),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title = element_text(size=16),
        plot.title = element_text(size=10)
  )
bar_mmtrends
bar_mmtrends <- ggplotGrob(bar_mmtrends)

gt1 <- gtable(widths = unit(rep(2,32), c("null")), heights = unit(rep(2,32), "null"))
gt1 <- gtable_add_grob(gt1, f, t=1, b=32, l=1, r=32)
gt1 <- gtable_add_grob(gt1, bar_mmtrends, t = 29, l = 4, b = 20, r = 10)
grid.draw(gt1)


# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_cor_EF_Tmaxbar_min <- dbar_cor_EF_Tmaxbar_min + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_cor_EF_Tmaxbar_min_smaller <- grid.arrange(empty, dbar_cor_EF_Tmaxbar_min, empty , ncol=3, widths = c(1,6,1))

plot_cor_EF_Tmaxbar_min <- grid.arrange(gt1,dbar_cor_EF_Tmaxbar_min_smaller, nrow = 2, heights = c(.8,.2))

plot_spat <- grid.arrange(plot_kendall_av_EF, plot_cor_EF_Tmaxbar_min, plot_cor_ELI_EF, ncol = 1)
# ggsave("Figures/Fig2.png", plot = plot_spat, width = 9*1.3, height = 21*1.3, units = "in")
ggsave("testdir/Fig2.png", plot = plot_spat, width = 9*1.3, height = 21*1.3, units = "in")







# make a mask since the observations
mask_all <- array(0,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    for(source_id in 1:length(cmip6_data.df$source_id)){
      if(sum(!is.na(dcorr.list[[source_id]][x,y,])) == 12){ # check if all the models have a value there
        mask_all[x,y] <- mask_all[x,y] + 1
      }
    }
  }
}

mask_full <- array(NaN,c(180,90)); mask_full[which(mask_all == 12)] <- 1

save(mask_full, file = paste0(path_RData, "mask_full_acccmip6.RData"))

center_year <- seq(1985,2095,10)
global_tseries.df <- setNames(data.frame(matrix(ncol = 10, nrow = 0)),
                              c("ELI", "ELI_max", "corr_rgy_veg", "corr_wtr_veg", "Tmeanbar", "Tmaxbar", "Tmaxbar_min","center_year","year","source_id"))
for(source_id in 1:length(cmip6_data.df$source_id)){
  for(i in 1:12){
    global_tseries.df <- rbind(global_tseries.df,
                               data.frame("ELI" = weighted.mean(x = dcorr.list[[source_id]][,,i] * mask_full, w = area.array, na.rm = T),
                                          "ELI_max" = weighted.mean(x = dcorr_max.list[[source_id]][,,i] * mask_full, w = area.array, na.rm = T),
                                          "corr_rgy_veg" = weighted.mean(x = corr_rgy_veg.list[[source_id]][,,i] * mask_full, w = area.array, na.rm = T),
                                          "corr_wtr_veg" = weighted.mean(x = corr_wtr_veg.list[[source_id]][,,i] * mask_full, w = area.array, na.rm = T),
                                          "Tmeanbar" = weighted.mean(x = av_tas.list[[source_id]][,,i] * mask_full, w = area.array, na.rm = T),
                                          "Tmaxbar" = weighted.mean(x = Tmaxbar.list[[source_id]][,,i] * mask_full, w = area.array, na.rm = T),
                                          "Tmaxbar_min" = weighted.mean(x = (Tmaxbar.list[[source_id]][,,i] - av_tas.list[[source_id]][,,i]) * mask_full, w = area.array, na.rm = T),
                                          "center_year" = center_year[i],
                                          "year" = center_year[i]-5,
                                          "source_id" = cmip6_data.df$source_id[source_id]
                               ))
  }
}
global_tseries.df$reg <- "Global"
global_tseries_lim.df <- global_tseries.df

agreement.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(models_with_full_timeseries[x,y])){
      if(models_with_full_timeseries[x,y] > 8){
        if(length(which(sign(kendall_Tmaxbar_min[x,y,,1]) == 1)) > 8 | length(which(sign(kendall_Tmaxbar_min[x,y,,1]) == -1)) > 8){
          agreement.array[x,y] <- 1
        }
      }
    }
  }
}

global_tseries_ELI.df <- setNames(data.frame(matrix(ncol = 11, nrow = 0)),
                                  c("ELI_all", "ELI_agree", "ELI_pos", "ELI_neg", 
                                    "ELI_max_all", "ELI_max_agree", "ELI_max_pos", "ELI_max_neg", 
                                    "center_year","year","source_id"))
for(source_id in 1:length(cmip6_data.df$source_id)){
  for(i in 1:12){
    global_tseries_ELI.df <- rbind(global_tseries_ELI.df,
                                   data.frame("ELI_all" = weighted.mean(x = dcorr.list[[source_id]][,,i] * mask_full, w = area.array, na.rm = T),
                                              "ELI_agree" = weighted.mean(x = dcorr.list[[source_id]][,,i] * mask_full * agreement.array, w = area.array, na.rm = T),
                                              "ELI_pos" = weighted.mean(x = (dcorr.list[[source_id]][,,i] * mask_full * agreement.array)[which(mmmean_kendall_Tmaxbar_min > 0)], w = area.array[which(mmmean_kendall_Tmaxbar_min > 0)], na.rm = T),
                                              "ELI_neg" = weighted.mean(x = (dcorr.list[[source_id]][,,i] * mask_full * agreement.array)[which(mmmean_kendall_Tmaxbar_min < 0)], w = area.array[which(mmmean_kendall_Tmaxbar_min < 0)], na.rm = T),
                                              "ELI_max_all" = weighted.mean(x = dcorr_max.list[[source_id]][,,i] * mask_full, w = area.array, na.rm = T),
                                              "ELI_max_agree" = weighted.mean(x = dcorr_max.list[[source_id]][,,i] * mask_full * agreement.array, w = area.array, na.rm = T),
                                              "ELI_max_pos" = weighted.mean(x = (dcorr_max.list[[source_id]][,,i] * mask_full * agreement.array)[which(mmmean_kendall_Tmaxbar_min > 0)], w = area.array[which(mmmean_kendall_Tmaxbar_min > 0)], na.rm = T),
                                              "ELI_max_neg" = weighted.mean(x = (dcorr_max.list[[source_id]][,,i] * mask_full * agreement.array)[which(mmmean_kendall_Tmaxbar_min < 0)], w = area.array[which(mmmean_kendall_Tmaxbar_min < 0)], na.rm = T),
                                              "center_year" = center_year[i],
                                              "year" = center_year[i]-5,
                                              "source_id" = cmip6_data.df$source_id[source_id]
                                   ))
  }
}


prior_dcorr_mask <- array(NaN,c(180,90,12,3))
for(x in 1:180){
  for(y in 1:90){
    for(source_id in 1:length(cmip6_data.df$source_id)){
      inter_mean_ELI <- mean(dcorr_max.list[[source_id]][x,y,1:3])
      if(!is.na(inter_mean_ELI)){
        if(inter_mean_ELI < -.2){
          prior_dcorr_mask[x,y,source_id,1] <- 1
        }else if(inter_mean_ELI < 0.2){
          prior_dcorr_mask[x,y,source_id,2] <- 1
        }else{
          prior_dcorr_mask[x,y,source_id,3] <- 1
        }
      }
    }
  }
}


reg <- c("SAM", "NAM", "CEU", "NAS")
reg_lim <- c("rgy_lim","trans","wtr_lim")
count_j <- 1
for(j in seq(1,16,4)){ # loop over all regions
  lonmin <- min(hotspot_regs.df$x[j:(j+3)])
  lonmax <- max(hotspot_regs.df$x[j:(j+3)])
  latmin <- min(hotspot_regs.df$y[j:(j+3)])
  latmax <- max(hotspot_regs.df$y[j:(j+3)])
  reg_mask <- array(NaN,c(180,90)); reg_mask[which(lon >= lonmin & lon <= lonmax), which(lat >= latmin & lat <= latmax)] <- 1
  for(source_id in 1:length(cmip6_data.df$source_id)){
    for(t in 1:12){ # loop over all decades                                                                                                                                                                                          which(lat >= latmin & lat <= latmax),t])[2])){
      global_tseries.df <- rbind(global_tseries.df,
                                 data.frame("ELI" = weighted.mean(x = (dcorr.list[[source_id]][,,t] * reg_mask * mask_full),
                                                                  w = area.array,
                                                                  na.rm = T),
                                            "ELI_max" = weighted.mean(x = (dcorr_max.list[[source_id]][,,t] * reg_mask * mask_full),
                                                                      w = area.array,
                                                                      na.rm = T),
                                            "corr_rgy_veg" = weighted.mean(x = (corr_rgy_veg.list[[source_id]][,,t] * reg_mask * mask_full),
                                                                           w = area.array,
                                                                           na.rm = T),
                                            "corr_wtr_veg" = weighted.mean(x = (corr_wtr_veg.list[[source_id]][,,t] * reg_mask * mask_full),
                                                                           w = area.array,
                                                                           na.rm = T),
                                            "Tmeanbar" = weighted.mean(x = (av_tas.list[[source_id]][,,t] * reg_mask * mask_full),
                                                                       w = area.array,
                                                                       na.rm = T),
                                            "Tmaxbar" = weighted.mean(x = (Tmaxbar.list[[source_id]][,,t] * reg_mask * mask_full),
                                                                      w = area.array,
                                                                      na.rm = T),
                                            "Tmaxbar_min" = weighted.mean(x = (Tmaxbar.list[[source_id]][,,t] -
                                                                                 av_tas.list[[source_id]][,,t]) * reg_mask * mask_full,
                                                                          w = area.array,
                                                                          na.rm = T),
                                            "center_year" = center_year[t],
                                            "year" = center_year[t]-5,
                                            "source_id" = cmip6_data.df$source_id[source_id],
                                            "reg" = reg[count_j]))
    }
  }
  count_j <- count_j + 1
}

count_j <- 1
for(j in 1:3){ # loop over all reg_lim
  for(source_id in 1:length(cmip6_data.df$source_id)){
    for(t in 1:12){ # loop over all decades
      global_tseries_lim.df <- rbind(global_tseries_lim.df,
                                     data.frame("ELI" = weighted.mean(x = (dcorr.list[[source_id]][,,t] * prior_dcorr_mask[,,source_id,j] * mask_full),
                                                                      w = area.array,
                                                                      na.rm = T),
                                                "ELI_max" = weighted.mean(x = (dcorr_max.list[[source_id]][,,t] * prior_dcorr_mask[,,source_id,j] * mask_full),
                                                                          w = area.array,
                                                                          na.rm = T),
                                                "corr_rgy_veg" = weighted.mean(x = (corr_rgy_veg.list[[source_id]][,,t] * prior_dcorr_mask[,,source_id,j] * mask_full),
                                                                               w = area.array,
                                                                               na.rm = T),
                                                "corr_wtr_veg" = weighted.mean(x = (corr_wtr_veg.list[[source_id]][,,t] * prior_dcorr_mask[,,source_id,j] * mask_full),
                                                                               w = area.array,
                                                                               na.rm = T),
                                                "Tmeanbar" = weighted.mean(x = (av_tas.list[[source_id]][,,t] * prior_dcorr_mask[,,source_id,j] * mask_full),
                                                                           w = area.array,
                                                                           na.rm = T),
                                                "Tmaxbar" = weighted.mean(x = (Tmaxbar.list[[source_id]][,,t] * prior_dcorr_mask[,,source_id,j] * mask_full),
                                                                          w = area.array,
                                                                          na.rm = T),
                                                "Tmaxbar_min" = weighted.mean(x = (Tmaxbar.list[[source_id]][,,t] -
                                                                                     av_tas.list[[source_id]][,,t]) * prior_dcorr_mask[,,source_id,j] * mask_full,
                                                                              w = area.array,
                                                                              na.rm = T),
                                                "center_year" = center_year[t],
                                                "year" = center_year[t]-5,
                                                "source_id" = cmip6_data.df$source_id[source_id],
                                                "reg" = reg_lim[count_j]))
      
    }
  }
  count_j <- count_j + 1
}

for(i in 1:length(global_tseries.df$ELI)){
  global_tseries.df$ELI_max_rel_1980[i] <- (global_tseries.df$ELI_max[i] - global_tseries.df$ELI_max[which(global_tseries.df$reg == global_tseries.df$reg[i] & global_tseries.df$source_id == global_tseries.df$source_id[i])][1])
  global_tseries.df$corr_rgy_veg_rel_1980[i] <- (global_tseries.df$corr_rgy_veg[i] - global_tseries.df$corr_rgy_veg[which(global_tseries.df$reg == global_tseries.df$reg[i] & global_tseries.df$source_id == global_tseries.df$source_id[i])][1])
  global_tseries.df$corr_wtr_veg_rel_1980[i] <- (global_tseries.df$corr_wtr_veg[i] - global_tseries.df$corr_wtr_veg[which(global_tseries.df$reg == global_tseries.df$reg[i] & global_tseries.df$source_id == global_tseries.df$source_id[i])][1])
  global_tseries.df$Tmeanbar_rel_1980[i] <- (global_tseries.df$Tmeanbar[i] - global_tseries.df$Tmeanbar[which(global_tseries.df$reg == global_tseries.df$reg[i] & global_tseries.df$source_id == global_tseries.df$source_id[i])][1])
  global_tseries.df$Tmaxbar_rel_1980[i] <- (global_tseries.df$Tmaxbar[i] - global_tseries.df$Tmaxbar[which(global_tseries.df$reg == global_tseries.df$reg[i] & global_tseries.df$source_id == global_tseries.df$source_id[i])][1])
  global_tseries.df$Tmaxbar_min_rel_1980[i] <- (global_tseries.df$Tmaxbar_min[i] - global_tseries.df$Tmaxbar_min[which(global_tseries.df$reg == global_tseries.df$reg[i] & global_tseries.df$source_id == global_tseries.df$source_id[i])][1])
}

for(i in 1:length(global_tseries_lim.df$ELI)){
  global_tseries_lim.df$ELI_rel_1980[i] <- (global_tseries_lim.df$ELI[i] - global_tseries_lim.df$ELI[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i] & global_tseries_lim.df$source_id == global_tseries_lim.df$source_id[i])][1])
  global_tseries_lim.df$ELI_max_rel_1980[i] <- (global_tseries_lim.df$ELI_max[i] - global_tseries_lim.df$ELI_max[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i] & global_tseries_lim.df$source_id == global_tseries_lim.df$source_id[i])][1])
  global_tseries_lim.df$corr_rgy_veg_rel_1980[i] <- (global_tseries_lim.df$corr_rgy_veg[i] - global_tseries_lim.df$corr_rgy_veg[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i] & global_tseries_lim.df$source_id == global_tseries_lim.df$source_id[i])][1])
  global_tseries_lim.df$corr_wtr_veg_rel_1980[i] <- (global_tseries_lim.df$corr_wtr_veg[i] - global_tseries_lim.df$corr_wtr_veg[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i] & global_tseries_lim.df$source_id == global_tseries_lim.df$source_id[i])][1])
  global_tseries_lim.df$Tmeanbar_rel_1980[i] <- (global_tseries_lim.df$Tmeanbar[i] - global_tseries_lim.df$Tmeanbar[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i] & global_tseries_lim.df$source_id == global_tseries_lim.df$source_id[i])][1])
  global_tseries_lim.df$Tmaxbar_rel_1980[i] <- (global_tseries_lim.df$Tmaxbar[i] - global_tseries_lim.df$Tmaxbar[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i] & global_tseries_lim.df$source_id == global_tseries_lim.df$source_id[i])][1])
  global_tseries_lim.df$Tmaxbar_min_rel_1980[i] <- (global_tseries_lim.df$Tmaxbar_min[i] - global_tseries_lim.df$Tmaxbar_min[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i] & global_tseries_lim.df$source_id == global_tseries_lim.df$source_id[i])][1])
}

col_reg <- c("Global" = "black",
             "CEU" = brewer.pal(12,"Paired")[2],
             "NAM" = brewer.pal(12,"Paired")[4],
             "NAS" = brewer.pal(12,"Paired")[6],
             "SAM" = brewer.pal(12,"Paired")[8])

col_reg_lim <- c("Global" = "black",
                 "rgy_lim" = brewer.pal(9,"BrBG")[8],
                 "trans" = 'slategrey',
                 "wtr_lim" = brewer.pal(9,"BrBG")[2])

cols_source_id <- c("ACCESS-ESM1-5" = colorRampPalette(brewer.pal(9, "Set1"))(12)[1],
                    "BCC-CSM2-MR" = colorRampPalette(brewer.pal(9, "Set1"))(12)[2],
                    "CMCC-ESM2" = colorRampPalette(brewer.pal(9, "Set1"))(12)[3],
                    "CNRM-CM6-1" = colorRampPalette(brewer.pal(9, "Set1"))(12)[4],
                    "CNRM-ESM2-1" = colorRampPalette(brewer.pal(9, "Set1"))(12)[5],
                    "EC-Earth3-CC" = colorRampPalette(brewer.pal(9, "Set1"))(12)[6],
                    "GFDL-ESM4" = colorRampPalette(brewer.pal(9, "Set1"))(12)[7],
                    "HadGEM3-GC31-LL" = colorRampPalette(brewer.pal(9, "Set1"))(12)[8],
                    "MPI-ESM1-2-HR" = colorRampPalette(brewer.pal(9, "Set1"))(12)[9],
                    "MPI-ESM1-2-LR" = colorRampPalette(brewer.pal(9, "Set1"))(12)[10],
                    "MRI-ESM2-0" = colorRampPalette(brewer.pal(9, "Set1"))(12)[11],
                    "UKESM1-0-LL" = colorRampPalette(brewer.pal(9, "Set1"))(12)[12])

global_tseries.df$reg_f <- factor(global_tseries.df$reg, levels = c("Global","CEU","NAM","NAS","SAM"))


sign.df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
                    c("sign_Tmaxbar_min","sign_ELI_max","reg","source_id"))
for(r in unique(global_tseries.df$reg)){
  for(s in unique(global_tseries.df$source_id)){
    sign.df <- rbind(sign.df,
                     data.frame("sign_Tmaxbar_min" = unname(kendallTrendTest(global_tseries.df$Tmaxbar_min[which(global_tseries.df$reg == r & global_tseries.df$source_id == s)])$p.value),
                                "sign_ELI_min" = unname(kendallTrendTest(global_tseries.df$ELI_max[which(global_tseries.df$reg == r & global_tseries.df$source_id == s)])$p.value),
                                "reg" = r,
                                "source_id" = s))
  }
}
insign_Tmaxbar_min <- sign.df[which(sign.df$sign_Tmaxbar_min > .05),]
insign_ELI_max <- sign.df[which(sign.df$sign_ELI_min > .05),]
sign_Tmaxbar_min <- sign.df[which(sign.df$sign_Tmaxbar_min < .05),]
sign_ELI_max <- sign.df[which(sign.df$sign_ELI_min < .05),]

global_tseries_sign_Tmaxbar_min.df <- setNames(data.frame(matrix(ncol = ncol(global_tseries.df), nrow = 0)),
                                               names(global_tseries.df))
for(i in 1:length(sign_Tmaxbar_min$reg)){
  global_tseries_sign_Tmaxbar_min.df <- rbind(global_tseries_sign_Tmaxbar_min.df,
                                              global_tseries.df[which(global_tseries.df$reg == sign_Tmaxbar_min$reg[i] & global_tseries.df$source_id == sign_Tmaxbar_min$source_id[i]),])
}

global_tseries_sign_ELI_max.df <- setNames(data.frame(matrix(ncol = ncol(global_tseries.df), nrow = 0)),
                                           names(global_tseries.df))
for(i in 1:length(sign_ELI_max$reg)){
  global_tseries_sign_ELI_max.df <- rbind(global_tseries_sign_ELI_max.df,
                                          global_tseries.df[which(global_tseries.df$reg == sign_ELI_max$reg[i] & global_tseries.df$source_id == sign_ELI_max$source_id[i]),])
}

global_tseries_insign_Tmaxbar_min.df <- setNames(data.frame(matrix(ncol = ncol(global_tseries.df), nrow = 0)),
                                                 names(global_tseries.df))
for(i in 1:length(insign_Tmaxbar_min$reg)){
  global_tseries_insign_Tmaxbar_min.df <- rbind(global_tseries_insign_Tmaxbar_min.df,
                                                global_tseries.df[which(global_tseries.df$reg == insign_Tmaxbar_min$reg[i] & global_tseries.df$source_id == insign_Tmaxbar_min$source_id[i]),])
}

global_tseries_insign_ELI_max.df <- setNames(data.frame(matrix(ncol = ncol(global_tseries.df), nrow = 0)),
                                             names(global_tseries.df))
for(i in 1:length(insign_ELI_max$reg)){
  global_tseries_insign_ELI_max.df <- rbind(global_tseries_insign_ELI_max.df,
                                            global_tseries.df[which(global_tseries.df$reg == insign_ELI_max$reg[i] & global_tseries.df$source_id == insign_ELI_max$source_id[i]),])
}



# abs 1980
qa4 <- ggplot(global_tseries_sign_Tmaxbar_min.df, aes(x=year,y=Tmaxbar_min_rel_1980,col=source_id)) +
  geom_line(linetype = 'solid', size = 1.25) +
  geom_line(data = global_tseries_insign_Tmaxbar_min.df, linetype = 'dotted', size = 1.25) +
  geom_hline(yintercept=0) +
  scale_x_continuous("",
                     breaks = seq(1980,2090,20),
                     labels = seq(1980,2090,20), expand = c(0,0)) +
  scale_y_continuous(expression("change since 1980 (K)"), expand = c(0,0)) +
  scale_color_manual("",
                     values = cols_source_id) +
  facet_wrap(~reg_f) +
  theme(legend.position = c(.825,.25),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        panel.spacing = unit(3, "lines"),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24, hjust = 0.5),
        plot.tag.position = c(.55,0.03)
  ) + 
  ggtitle(expression("a) Temperature excess"))
qa4

# abs 1980
qb4 <- ggplot(global_tseries_sign_ELI_max.df, aes(x=year,y=ELI_max_rel_1980,col=source_id)) +
  geom_line(linetype = 'solid', size = 1.25) +
  geom_line(data = global_tseries_insign_ELI_max.df, linetype = 'dotted', size = 1.25) +
  geom_hline(yintercept=0) +
  scale_x_continuous("",
                     breaks = seq(1980,2090,20),
                     labels = seq(1980,2090,20), expand = c(0,0)) +
  scale_y_continuous(expression("change since 1980 (-)"), expand = c(0,0)) +
  scale_color_manual("",
                     values = cols_source_id) +
  facet_wrap(~reg_f) +
  theme(legend.position = "none",
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        panel.spacing = unit(3, "lines"),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24, hjust = 0.5),
        plot.tag.position = c(.55,0.03)
  ) + 
  ggtitle(expression("b) Ecosystem Limitation Index"))
qb4

plots <- plot_grid(qa4, qb4, nrow=2, align = 'v')
# ggsave("Figures/SFig7.png", plot = plots, width = 9*2/1.5, height = 9*2/1.5, units = "in")
ggsave("testdir/SFig7.png", plot = plots, width = 9*2/1.5, height = 9*2/1.5, units = "in")


global_tseries_ELI.df$ELI_all_rel_1980 <-
  global_tseries_ELI.df$ELI_agree_rel_1980 <-
  global_tseries_ELI.df$ELI_pos_rel_1980 <-
  global_tseries_ELI.df$ELI_neg_rel_1980 <-
  NaN
for(i in 1:length(global_tseries_ELI.df$ELI_all)){
  global_tseries_ELI.df$ELI_all_rel_1980[i] <- (global_tseries_ELI.df$ELI_all[i] - global_tseries_ELI.df$ELI_all[which(global_tseries_ELI.df$source_id == global_tseries_ELI.df$source_id[i])][1])
  global_tseries_ELI.df$ELI_agree_rel_1980[i] <- (global_tseries_ELI.df$ELI_agree[i] - global_tseries_ELI.df$ELI_agree[which(global_tseries_ELI.df$source_id == global_tseries_ELI.df$source_id[i])][1])
  global_tseries_ELI.df$ELI_pos_rel_1980[i] <- (global_tseries_ELI.df$ELI_pos[i] - global_tseries_ELI.df$ELI_pos[which(global_tseries_ELI.df$source_id == global_tseries_ELI.df$source_id[i])][1])
  global_tseries_ELI.df$ELI_neg_rel_1980[i] <- (global_tseries_ELI.df$ELI_neg[i] - global_tseries_ELI.df$ELI_neg[which(global_tseries_ELI.df$source_id == global_tseries_ELI.df$source_id[i])][1])
}

global_tseries_ELI.df$ELI_max_all_rel_1980 <-
  global_tseries_ELI.df$ELI_max_agree_rel_1980 <-
  global_tseries_ELI.df$ELI_max_pos_rel_1980 <-
  global_tseries_ELI.df$ELI_max_neg_rel_1980 <-
  NaN
for(i in 1:length(global_tseries_ELI.df$ELI_max_all)){
  global_tseries_ELI.df$ELI_max_all_rel_1980[i] <- (global_tseries_ELI.df$ELI_max_all[i] - global_tseries_ELI.df$ELI_max_all[which(global_tseries_ELI.df$source_id == global_tseries_ELI.df$source_id[i])][1])
  global_tseries_ELI.df$ELI_max_agree_rel_1980[i] <- (global_tseries_ELI.df$ELI_max_agree[i] - global_tseries_ELI.df$ELI_max_agree[which(global_tseries_ELI.df$source_id == global_tseries_ELI.df$source_id[i])][1])
  global_tseries_ELI.df$ELI_max_pos_rel_1980[i] <- (global_tseries_ELI.df$ELI_max_pos[i] - global_tseries_ELI.df$ELI_max_pos[which(global_tseries_ELI.df$source_id == global_tseries_ELI.df$source_id[i])][1])
  global_tseries_ELI.df$ELI_max_neg_rel_1980[i] <- (global_tseries_ELI.df$ELI_max_neg[i] - global_tseries_ELI.df$ELI_max_neg[which(global_tseries_ELI.df$source_id == global_tseries_ELI.df$source_id[i])][1])
}

# mmmean and mmsd per i) region and ii) center_year
global_mmmtseries.df <- setNames(data.frame(matrix(ncol = 24, nrow = 0)),
                                 c("mmmean_ELI", "mmminsd_ELI", "mmplssd_ELI",
                                   "mmmean_ELI_max", "mmminsd_ELI_max", "mmplssd_ELI_max",
                                   "mmmean_corr_rgy_veg", "mmminsd_corr_rgy_veg", "mmplssd_corr_rgy_veg",
                                   "mmmean_corr_wtr_veg", "mmminsd_corr_wtr_veg", "mmplssd_corr_wtr_veg",
                                   "mmmean_Tmeanbar", "mmminsd_Tmeanbar", "mmplssd_Tmeanbar",
                                   "mmmean_Tmaxbar", "mmminsd_Tmaxbar", "mmplssd_Tmaxbar",
                                   "mmmean_Tmaxbar_min", "mmminsd_Tmaxbar_min", "mmplssd_Tmaxbar_min",
                                   "center_year","year","reg"))
for(i_reg in unique(global_tseries.df$reg)){
  for(i in unique(global_tseries.df$year)){
    global_mmmtseries.df <- rbind(global_mmmtseries.df,
                                  data.frame("mmmean_ELI" = mean(global_tseries.df$ELI_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmminsd_ELI" = mean(global_tseries.df$ELI_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) - sd(global_tseries.df$ELI_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmplssd_ELI" = mean(global_tseries.df$ELI_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) + sd(global_tseries.df$ELI_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmmean_ELI_max" = mean(global_tseries.df$ELI_max_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmminsd_ELI_max" = mean(global_tseries.df$ELI_max_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) - sd(global_tseries.df$ELI_max_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmplssd_ELI_max" = mean(global_tseries.df$ELI_max_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) + sd(global_tseries.df$ELI_max_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmmean_corr_rgy_veg" = mean(global_tseries.df$corr_rgy_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmminsd_corr_rgy_veg" = mean(global_tseries.df$corr_rgy_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) - sd(global_tseries.df$corr_rgy_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmplssd_corr_rgy_veg" = mean(global_tseries.df$corr_rgy_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) + sd(global_tseries.df$corr_rgy_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmmean_corr_wtr_veg" = mean(global_tseries.df$corr_wtr_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmminsd_corr_wtr_veg" = mean(global_tseries.df$corr_wtr_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) - sd(global_tseries.df$corr_wtr_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmplssd_corr_wtr_veg" = mean(global_tseries.df$corr_wtr_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) + sd(global_tseries.df$corr_wtr_veg_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmmean_Tmeanbar" = mean(global_tseries.df$Tmeanbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmminsd_Tmeanbar" = mean(global_tseries.df$Tmeanbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) - sd(global_tseries.df$Tmeanbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmplssd_Tmeanbar" = mean(global_tseries.df$Tmeanbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) + sd(global_tseries.df$Tmeanbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmmean_Tmaxbar" = mean(global_tseries.df$Tmaxbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmminsd_Tmaxbar" = mean(global_tseries.df$Tmaxbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) - sd(global_tseries.df$Tmaxbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmplssd_Tmaxbar" = mean(global_tseries.df$Tmaxbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) + sd(global_tseries.df$Tmaxbar_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmmean_Tmaxbar_min" = mean(global_tseries.df$Tmaxbar_min_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmminsd_Tmaxbar_min" = mean(global_tseries.df$Tmaxbar_min_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) - sd(global_tseries.df$Tmaxbar_min_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "mmplssd_Tmaxbar_min" = mean(global_tseries.df$Tmaxbar_min_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T) + sd(global_tseries.df$Tmaxbar_min_rel_1980[which(global_tseries.df$reg == i_reg & global_tseries.df$year == i)], na.rm = T),
                                             "center_year" = i+5,
                                             "year" = i,
                                             "reg" = i_reg))
  }
}

# mmmean and mmsd per i) region and ii) center_year
global_mmmtseries_lim.df <- setNames(data.frame(matrix(ncol = 24, nrow = 0)),
                                     c("mmmean_ELI", "mmminsd_ELI", "mmplssd_ELI",
                                       "mmmean_ELI_max", "mmminsd_ELI_max", "mmplssd_ELI_max",
                                       "mmmean_corr_rgy_veg", "mmminsd_corr_rgy_veg", "mmplssd_corr_rgy_veg",
                                       "mmmean_corr_wtr_veg", "mmminsd_corr_wtr_veg", "mmplssd_corr_wtr_veg",
                                       "mmmean_Tmeanbar", "mmminsd_Tmeanbar", "mmplssd_Tmeanbar",
                                       "mmmean_Tmaxbar", "mmminsd_Tmaxbar", "mmplssd_Tmaxbar",
                                       "mmmean_Tmaxbar_min", "mmminsd_Tmaxbar_min", "mmplssd_Tmaxbar_min",
                                       "center_year","year","reg"))
for(i_reg in unique(global_tseries_lim.df$reg)){
  for(i in unique(global_tseries_lim.df$year)){
    global_mmmtseries_lim.df <- rbind(global_mmmtseries_lim.df,
                                      data.frame("mmmean_ELI" = mean(global_tseries_lim.df$ELI_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmminsd_ELI" = mean(global_tseries_lim.df$ELI_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) - sd(global_tseries_lim.df$ELI_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmplssd_ELI" = mean(global_tseries_lim.df$ELI_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) + sd(global_tseries_lim.df$ELI_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmmean_ELI_max" = mean(global_tseries_lim.df$ELI_max_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmminsd_ELI_max" = mean(global_tseries_lim.df$ELI_max_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) - sd(global_tseries_lim.df$ELI_max_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmplssd_ELI_max" = mean(global_tseries_lim.df$ELI_max_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) + sd(global_tseries_lim.df$ELI_max_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmmean_corr_rgy_veg" = mean(global_tseries_lim.df$corr_rgy_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmminsd_corr_rgy_veg" = mean(global_tseries_lim.df$corr_rgy_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) - sd(global_tseries_lim.df$corr_rgy_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmplssd_corr_rgy_veg" = mean(global_tseries_lim.df$corr_rgy_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) + sd(global_tseries_lim.df$corr_rgy_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmmean_corr_wtr_veg" = mean(global_tseries_lim.df$corr_wtr_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmminsd_corr_wtr_veg" = mean(global_tseries_lim.df$corr_wtr_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) - sd(global_tseries_lim.df$corr_wtr_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmplssd_corr_wtr_veg" = mean(global_tseries_lim.df$corr_wtr_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) + sd(global_tseries_lim.df$corr_wtr_veg_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmmean_Tmeanbar" = mean(global_tseries_lim.df$Tmeanbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmminsd_Tmeanbar" = mean(global_tseries_lim.df$Tmeanbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) - sd(global_tseries_lim.df$Tmeanbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmplssd_Tmeanbar" = mean(global_tseries_lim.df$Tmeanbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) + sd(global_tseries_lim.df$Tmeanbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmmean_Tmaxbar" = mean(global_tseries_lim.df$Tmaxbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmminsd_Tmaxbar" = mean(global_tseries_lim.df$Tmaxbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) - sd(global_tseries_lim.df$Tmaxbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmplssd_Tmaxbar" = mean(global_tseries_lim.df$Tmaxbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) + sd(global_tseries_lim.df$Tmaxbar_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmmean_Tmaxbar_min" = mean(global_tseries_lim.df$Tmaxbar_min_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmminsd_Tmaxbar_min" = mean(global_tseries_lim.df$Tmaxbar_min_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) - sd(global_tseries_lim.df$Tmaxbar_min_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "mmplssd_Tmaxbar_min" = mean(global_tseries_lim.df$Tmaxbar_min_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T) + sd(global_tseries_lim.df$Tmaxbar_min_rel_1980[which(global_tseries_lim.df$reg == i_reg & global_tseries_lim.df$year == i)], na.rm = T),
                                                 "center_year" = i+5,
                                                 "year" = i,
                                                 "reg" = i_reg))
  }
}

# mmmean and mmsd per i) region and ii) center_year
global_mmmtseries_ELI.df <- setNames(data.frame(matrix(ncol = 27, nrow = 0)),
                                     c("mmmean_ELI_all", "mmminsd_ELI_all", "mmplssd_ELI_all",
                                       "mmmean_ELI_agree", "mmminsd_ELI_agree", "mmplssd_ELI_agree",
                                       "mmmean_ELI_pos", "mmminsd_ELI_pos", "mmplssd_ELI_pos",
                                       "mmmean_ELI_neg", "mmminsd_ELI_neg", "mmplssd_ELI_neg",
                                       "mmmean_ELI_max_all", "mmminsd_ELI_max_all", "mmplssd_ELI_max_all",
                                       "mmmean_ELI_max_agree", "mmminsd_ELI_max_agree", "mmplssd_ELI_max_agree",
                                       "mmmean_ELI_max_pos", "mmminsd_ELI_max_pos", "mmplssd_ELI_max_pos",
                                       "mmmean_ELI_max_neg", "mmminsd_ELI_max_neg", "mmplssd_ELI_max_neg",
                                       "center_year","year","reg"))

for(i in unique(global_tseries_ELI.df$year)){
  global_mmmtseries_ELI.df <- rbind(global_mmmtseries_ELI.df,
                                    data.frame("mmmean_ELI_all" = mean(global_tseries_ELI.df$ELI_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmminsd_ELI_all" = mean(global_tseries_ELI.df$ELI_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) - sd(global_tseries_ELI.df$ELI_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmplssd_ELI_all" = mean(global_tseries_ELI.df$ELI_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) + sd(global_tseries_ELI.df$ELI_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmmean_ELI_agree" = mean(global_tseries_ELI.df$ELI_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmminsd_ELI_agree" = mean(global_tseries_ELI.df$ELI_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) - sd(global_tseries_ELI.df$ELI_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmplssd_ELI_agree" = mean(global_tseries_ELI.df$ELI_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) + sd(global_tseries_ELI.df$ELI_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmmean_ELI_pos" = mean(global_tseries_ELI.df$ELI_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmminsd_ELI_pos" = mean(global_tseries_ELI.df$ELI_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) - sd(global_tseries_ELI.df$ELI_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmplssd_ELI_pos" = mean(global_tseries_ELI.df$ELI_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) + sd(global_tseries_ELI.df$ELI_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmmean_ELI_neg" = mean(global_tseries_ELI.df$ELI_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmminsd_ELI_neg" = mean(global_tseries_ELI.df$ELI_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) - sd(global_tseries_ELI.df$ELI_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmplssd_ELI_neg" = mean(global_tseries_ELI.df$ELI_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) + sd(global_tseries_ELI.df$ELI_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               
                                               "mmmean_ELI_max_all" = mean(global_tseries_ELI.df$ELI_max_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmminsd_ELI_max_all" = mean(global_tseries_ELI.df$ELI_max_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) - sd(global_tseries_ELI.df$ELI_max_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmplssd_ELI_max_all" = mean(global_tseries_ELI.df$ELI_max_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) + sd(global_tseries_ELI.df$ELI_max_all_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmmean_ELI_max_agree" = mean(global_tseries_ELI.df$ELI_max_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmminsd_ELI_max_agree" = mean(global_tseries_ELI.df$ELI_max_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) - sd(global_tseries_ELI.df$ELI_max_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmplssd_ELI_max_agree" = mean(global_tseries_ELI.df$ELI_max_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) + sd(global_tseries_ELI.df$ELI_max_agree_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmmean_ELI_max_pos" = mean(global_tseries_ELI.df$ELI_max_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmminsd_ELI_max_pos" = mean(global_tseries_ELI.df$ELI_max_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) - sd(global_tseries_ELI.df$ELI_max_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmplssd_ELI_max_pos" = mean(global_tseries_ELI.df$ELI_max_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) + sd(global_tseries_ELI.df$ELI_max_pos_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmmean_ELI_max_neg" = mean(global_tseries_ELI.df$ELI_max_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmminsd_ELI_max_neg" = mean(global_tseries_ELI.df$ELI_max_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) - sd(global_tseries_ELI.df$ELI_max_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "mmplssd_ELI_max_neg" = mean(global_tseries_ELI.df$ELI_max_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T) + sd(global_tseries_ELI.df$ELI_max_neg_rel_1980[which(global_tseries_ELI.df$year == i)], na.rm = T),
                                               "center_year" = i+5,
                                               "year" = i))
}

global_mmmtseries_ELI.df_stacked <- data.frame("mmmean_ELI" = c(global_mmmtseries_ELI.df$mmmean_ELI_all,
                                                                global_mmmtseries_ELI.df$mmmean_ELI_agree,
                                                                global_mmmtseries_ELI.df$mmmean_ELI_pos,
                                                                global_mmmtseries_ELI.df$mmmean_ELI_neg),
                                               "mmminsd_ELI" = c(global_mmmtseries_ELI.df$mmminsd_ELI_all,
                                                                 global_mmmtseries_ELI.df$mmminsd_ELI_agree,
                                                                 global_mmmtseries_ELI.df$mmminsd_ELI_pos,
                                                                 global_mmmtseries_ELI.df$mmminsd_ELI_neg),
                                               "mmplssd_ELI" = c(global_mmmtseries_ELI.df$mmplssd_ELI_all,
                                                                 global_mmmtseries_ELI.df$mmplssd_ELI_agree,
                                                                 global_mmmtseries_ELI.df$mmplssd_ELI_pos,
                                                                 global_mmmtseries_ELI.df$mmplssd_ELI_neg),
                                               "mmmean_ELI_max" = c(global_mmmtseries_ELI.df$mmmean_ELI_max_all,
                                                                    global_mmmtseries_ELI.df$mmmean_ELI_max_agree,
                                                                    global_mmmtseries_ELI.df$mmmean_ELI_max_pos,
                                                                    global_mmmtseries_ELI.df$mmmean_ELI_max_neg),
                                               "mmminsd_ELI_max" = c(global_mmmtseries_ELI.df$mmminsd_ELI_max_all,
                                                                     global_mmmtseries_ELI.df$mmminsd_ELI_max_agree,
                                                                     global_mmmtseries_ELI.df$mmminsd_ELI_max_pos,
                                                                     global_mmmtseries_ELI.df$mmminsd_ELI_max_neg),
                                               "mmplssd_ELI_max" = c(global_mmmtseries_ELI.df$mmplssd_ELI_max_all,
                                                                     global_mmmtseries_ELI.df$mmplssd_ELI_max_agree,
                                                                     global_mmmtseries_ELI.df$mmplssd_ELI_max_pos,
                                                                     global_mmmtseries_ELI.df$mmplssd_ELI_max_neg),
                                               "center_year" = rep(global_mmmtseries_ELI.df$center_year, 4),
                                               "year" = rep(global_mmmtseries_ELI.df$year, 4),
                                               "type" = rep(c("all","agree","pos","neg"),each=12))


col_ELI <- c("all" = "black",
             "agree" = "black",
             "pos" = brewer.pal(9, "Reds")[7],
             "neg" = brewer.pal(9, "Blues")[7])
lty_ELI <- c("all" = "solid",
             "agree" = "dashed",
             "pos" = "dashed",
             "neg" = "dashed")


# abs 1980
m4 <- ggplot(global_mmmtseries_ELI.df_stacked, aes(x=year,y=mmmean_ELI_max,col=type,linetype=type)) +
  geom_line(size = 1.25) +
  geom_hline(yintercept=0) +
  scale_x_continuous("",
                     breaks = seq(1980,2090,20),
                     labels = seq(1980,2090,20), expand = c(0,0)) +
  scale_y_continuous(expression(paste("change since 1980")), expand = c(0,0)) +
  scale_color_manual("",
                     values = col_ELI,
                     labels = c("all",
                                expression("only where at least 8/12 models agree on temperature excess trends"),
                                expression("only where at least 8/12 models agree on positive temperature excess trends"),
                                expression("only where at least 8/12 models agree on negative temperature excess trends"))) +
  scale_fill_manual("",
                    values = c("all" = "black",
                               "agree" = "black",
                               "pos" = "black",
                               "neg" = "black"),
                    labels = c("all",
                               expression("only where at least 8/12 models agree on temperature excess trends"),
                               expression("only where at least 8/12 models agree on positive temperature excess trends"),
                               expression("only where at least 8/12 models agree on negative temperature excess trends"))) +
  scale_linetype_manual("",
                        values = lty_ELI,
                        labels = c("all",
                                   expression("only where at least 8/12 models agree on temperature excess trends"),
                                   expression("only where at least 8/12 models agree on positive temperature excess trends"),
                                   expression("only where at least 8/12 models agree on negative temperature excess trends"))) +
  theme(legend.position = c(.3,.875),
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.width = unit(2, "cm"),
        legend.background = element_rect(fill = 'transparent'),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24),
        axis.line.y.right = element_line(color = "grey80"),
        axis.ticks.y.right = element_line(color = "grey80"),
        axis.text.y.right = element_text(color = "grey80"),
        axis.title.y.right = element_text(color = "grey80")
        
  ) + ggtitle("Ecosystem Limitation Index")
m4

# ggsave("Figures/SFig8.png", plot = m4, width = 10*1.25, height = 5*1.25, units = "in")
ggsave("testdir/SFig8.png", plot = m4, width = 10*1.25, height = 5*1.25, units = "in")



sign.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                    c("sign_mmmean_Tmaxbar_min","sign_mmmean_ELI_max","reg"))
for(r in unique(global_mmmtseries.df$reg)){
  sign.df <- rbind(sign.df,
                   data.frame("sign_mmmean_Tmaxbar_min" = unname(kendallTrendTest(global_mmmtseries.df$mmmean_Tmaxbar_min[which(global_mmmtseries.df$reg == r)])$p.value),
                              "sign_mmmean_ELI_min" = unname(kendallTrendTest(global_mmmtseries.df$mmmean_ELI_max[which(global_mmmtseries.df$reg == r)])$p.value),
                              "reg" = r))
}

# abs 1980
qa4 <- ggplot(global_mmmtseries.df, aes(x=year,y=mmmean_Tmaxbar_min,col=reg)) +
  geom_line(linetype = 'solid', size = 1.25) +
  geom_hline(yintercept=0) +
  scale_x_continuous("",
                     breaks = seq(1980,2090,20),
                     labels = seq(1980,2090,20), expand = c(0,0)) +
  scale_y_continuous(expression("change since 1980 (K)"), expand = c(0,0)) +
  scale_color_manual("region",
                     values = col_reg) +
  scale_fill_manual("region",
                    values = c("Global" = "black",
                               # "NEA" = "black",
                               "CEU" = "black",
                               "NAM" = "black",
                               "NAS" = "black",
                               "SAM" = "black")) +
  guides(fill = guide_legend(order = 1, title.position = "top", title.hjust = .5, label.position = "right", nrow=2,
                             label.theme = element_text(angle = 0))) +
  guides(col = guide_legend(order = 1, title.position = "top", title.hjust = .5, label.position = "right", nrow=2,
                            label.theme = element_text(angle = 0))) +
  theme(legend.position = c(.275,.875),
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24),
        axis.line.y.right = element_line(color = "grey80"),
        axis.ticks.y.right = element_line(color = "grey80"),
        axis.text.y.right = element_text(color = "grey80"),
        axis.title.y.right = element_text(color = "grey80")
        
  ) +
  ggtitle(expression("a) Temperature excess"))
qa4

# abs 1980
m4 <- ggplot(global_mmmtseries.df, aes(x=year,y=mmmean_ELI_max,col=reg)) +
  geom_line(linetype = 'solid', size = 1.25) +
  geom_hline(yintercept=0) +
  scale_x_continuous("",
                     breaks = seq(1980,2090,20),
                     labels = seq(1980,2090,20), expand = c(0,0)) +
  scale_y_continuous(expression(paste("change since 1980")), expand = c(0,0)) +
  scale_color_manual("region",
                     values = col_reg) +
  scale_fill_manual("region",
                    values = c("Global" = "black",
                               # "NEA" = "black",
                               "CEU" = "black",
                               "NAM" = "black",
                               "NAS" = "black",
                               "SAM" = "black")) +
  guides(fill = guide_legend(order = 1, title.position = "top", title.hjust = .5, label.position = "right", nrow=2,
                             label.theme = element_text(angle = 0))) +
  guides(col = guide_legend(order = 1, title.position = "top", title.hjust = .5, label.position = "right", nrow=2,
                            label.theme = element_text(angle = 0))) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24),
        axis.line.y.right = element_line(color = "grey80"),
        axis.ticks.y.right = element_line(color = "grey80"),
        axis.text.y.right = element_text(color = "grey80"),
        axis.title.y.right = element_text(color = "grey80")
        
  ) + ggtitle("b) Ecosystem Limitation Index")
m4

plots <- plot_grid(qa4, m4, ncol=2, align = 'h')
# ggsave("Figures/Fig3.png", plot = plots, width = 10*1.25, height = 5*1.25, units = "in")
ggsave("testdir/Fig3.png", plot = plots, width = 10*1.25, height = 5*1.25, units = "in")


load(paste0(path_RData, "global_tseries_ERA5.RData"))

global_tseries_CMIP6.df <- setNames(data.frame(matrix(ncol = 6, nrow = 0)),
                                    c("ELI_max","Tmaxbar_min", "center_year", "year", "source_id", "reg"))

# For the minimum
for(i in c('Global',reg)){
  print(i)
  for(t in seq(1985,2015,10)){
    global_tseries_CMIP6.df <- rbind(global_tseries_CMIP6.df,
                                     data.frame("ELI_max" = min(global_tseries.df$ELI_max_rel_1980[which(global_tseries.df$center_year == t & global_tseries.df$reg == i)]),
                                                "Tmaxbar_min" = min(global_tseries.df$Tmaxbar_min_rel_1980[which(global_tseries.df$center_year == t & global_tseries.df$reg == i)]),
                                                "center_year" = t,
                                                "year" = t - 5,
                                                "source_id" = "min",
                                                "reg" = i))
  }
}
# For the minimum
for(i in c('Global',reg)){
  print(i)
  for(t in seq(1985,2015,10)){
    global_tseries_CMIP6.df <- rbind(global_tseries_CMIP6.df,
                                     data.frame("ELI_max" = max(global_tseries.df$ELI_max_rel_1980[which(global_tseries.df$center_year == t & global_tseries.df$reg == i)]),
                                                "Tmaxbar_min" = max(global_tseries.df$Tmaxbar_min_rel_1980[which(global_tseries.df$center_year == t & global_tseries.df$reg == i)]),
                                                "center_year" = t,
                                                "year" = t - 5,
                                                "source_id" = "max",
                                                "reg" = i))
  }
}
# Add the mmmean
for(i in c('Global',reg)){
  print(i)
  for(t in seq(1985,2015,10)){
    global_tseries_CMIP6.df <- rbind(global_tseries_CMIP6.df,
                                     data.frame("ELI_max" = global_mmmtseries.df$mmmean_ELI_max[which(global_mmmtseries.df$center_year == t & global_mmmtseries.df$reg == i)],
                                                "Tmaxbar_min" = global_mmmtseries.df$mmmean_Tmaxbar_min[which(global_mmmtseries.df$center_year == t & global_mmmtseries.df$reg == i)],
                                                "center_year" = t,
                                                "year" = t - 5,
                                                "source_id" = "mean",
                                                "reg" = i))
  }
}
# Add the ERA5 values
for(i in c('Global',reg)){
  print(i)
  for(t in seq(1985,2015,10)){
    global_tseries_CMIP6.df <- rbind(global_tseries_CMIP6.df,
                                     data.frame("ELI_max" = global_tseries_ERA5.df$ELI_max_rel_1980[which(global_tseries_ERA5.df$center_year == t & global_tseries_ERA5.df$reg == i)],
                                                "Tmaxbar_min" = global_tseries_ERA5.df$Tmaxbar_min_rel_1980[which(global_tseries_ERA5.df$center_year == t & global_tseries_ERA5.df$reg == i)],
                                                "center_year" = t,
                                                "year" = t - 5,
                                                "source_id" = "ERA5",
                                                "reg" = i))
  }
}

global_tseries_CMIP6.df$reg_f <- factor(global_tseries_CMIP6.df$reg, levels = c("Global","CEU","NAM","NAS","SAM"))

cols_CMIP6_ERA5 = c("mean" = "black",
                    "ERA5" = "black")
lty_CMIP6_ERA5 <- c("mean" = "solid",
                    "ERA5" = "dashed")

global_tseries_CMIP6_minmax.df <- global_tseries_CMIP6.df[c("center_year", "year", "source_id", "reg", "reg_f")]
global_tseries_CMIP6_minmax.df$ELI_max_min <- global_tseries_CMIP6.df$ELI_max[which(global_tseries_CMIP6.df$source_id == "min")]
global_tseries_CMIP6_minmax.df$ELI_max_max <- global_tseries_CMIP6.df$ELI_max[which(global_tseries_CMIP6.df$source_id == "max")]
global_tseries_CMIP6_minmax.df$Tmaxbar_min_min <- global_tseries_CMIP6.df$Tmaxbar_min[which(global_tseries_CMIP6.df$source_id == "min")]
global_tseries_CMIP6_minmax.df$Tmaxbar_min_max <- global_tseries_CMIP6.df$Tmaxbar_min[which(global_tseries_CMIP6.df$source_id == "max")]

# abs 1980
qa4 <- ggplot(global_tseries_CMIP6.df[which(global_tseries_CMIP6.df$source_id == "mean" | global_tseries_CMIP6.df$source_id == "ERA5"),], aes(x=year,y=Tmaxbar_min,col=source_id,linetype=source_id)) +
  geom_ribbon(inherit.aes = F, data = global_tseries_CMIP6_minmax.df, aes(x=year,ymin=Tmaxbar_min_min,ymax=Tmaxbar_min_max, group = reg_f), col = "grey", fill = "grey", alpha = 0.5) +
  geom_line(size = 1.25) +
  geom_hline(yintercept=0) +
  scale_x_continuous("",
                     breaks = seq(1980,2000,10),
                     labels = seq(1980,2000,10), expand = c(0,0)) +
  scale_y_continuous(expression("change since 1980 (K)"), expand = c(0,0)) +
  scale_color_manual("",
                     labels = c("CMIP6 multi-model mean",
                                "ERA5-Land"),
                     values = cols_CMIP6_ERA5) +
  scale_linetype_manual("",
                        labels = c("CMIP6 multi-model mean",
                                   "ERA5-Land"),
                        values = lty_CMIP6_ERA5) +
  facet_wrap(~reg_f) +
  theme(legend.position = c(.825,.25),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.width = unit(3, "lines"),
        legend.background = element_rect(colour = NA, fill = NA),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        panel.spacing = unit(3, "lines"),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24, hjust = 0.5),
        plot.tag.position = c(.55,0.03)
  ) +
  ggtitle(expression("a) Temperature excess"))
qa4

# abs 1980
qb4 <- ggplot(global_tseries_CMIP6.df[which(global_tseries_CMIP6.df$source_id == "mean" | global_tseries_CMIP6.df$source_id == "ERA5"),], aes(x=year,y=ELI_max,col=source_id,linetype=source_id)) +
  geom_ribbon(inherit.aes = F, data = global_tseries_CMIP6_minmax.df, aes(x=year,ymin=ELI_max_min,ymax=ELI_max_max, group = reg_f), col = "grey", fill = "grey", alpha = 0.5) +
  geom_line(size = 1.25) +
  geom_hline(yintercept=0) +
  scale_x_continuous("",
                     breaks = seq(1980,2000,10),
                     labels = seq(1980,2000,10), expand = c(0,0)) +
  scale_y_continuous(expression("change since 1980 (-)"), expand = c(0,0)) +
  scale_color_manual("",
                     labels = c("CMIP6 multi-model mean",
                                "ERA5-Land"),
                     values = cols_CMIP6_ERA5) +
  scale_linetype_manual("",
                        labels = c("CMIP6 multi-model mean",
                                   "ERA5-Land"),
                        values = lty_CMIP6_ERA5) +
  facet_wrap(~reg_f) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        panel.spacing = unit(3, "lines"),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24, hjust = 0.5),
        plot.tag.position = c(.55,0.03)
  ) +
  ggtitle(expression("b) Ecosystem Limitation Index"))
qb4

plots <- plot_grid(qa4, qb4, nrow=2, align = 'v')
# ggsave("Figures/Fig4.png", plot = plots, width = 9*2/1.5, height = 9*2/1.5, units = "in")
ggsave("testdir/Fig4.png", plot = plots, width = 9*2/1.5, height = 9*2/1.5, units = "in")

# abs 1980
qa4_init_ELI <- ggplot(global_mmmtseries_lim.df, aes(x=year,y=mmmean_Tmaxbar_min,col=reg)) +
  geom_line(linetype = 'solid', size = 1.25) +
  geom_hline(yintercept=0) +
  scale_x_continuous("",
                     breaks = seq(1980,2090,20),
                     labels = seq(1980,2090,20), expand = c(0,0)) +
  scale_y_continuous(expression("change since 1980 (K)"), expand = c(0,0)) +
  scale_color_manual("mean ELI (1980 - 2010)",
                     values = col_reg_lim,
                     labels = c("Global","energy limited","transitional","water limited")) +
  scale_fill_manual("mean ELI (1980 - 2010)",
                    values = c("Global" = "black",
                               "rgy_lim" = "black",
                               "trans" = "black",
                               "wtr_lim" = "black"),
                    labels = c("Global","energy limited","transitional","water limited")) +
  theme(legend.position = c(.3,.75),
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24),
        axis.line.y.right = element_line(color = "grey80"),
        axis.ticks.y.right = element_line(color = "grey80"),
        axis.text.y.right = element_text(color = "grey80"),
        axis.title.y.right = element_text(color = "grey80")
        
  ) +
  ggtitle(expression("b) Temperature excess"))
qa4_init_ELI

# abs 1980
m4_init_ELI <- ggplot(global_mmmtseries_lim.df, aes(x=year,y=mmmean_ELI_max,col=reg)) +
  geom_line(linetype = 'solid', size = 1.25) +
  geom_hline(yintercept=0) +
  scale_x_continuous("",
                     breaks = seq(1980,2090,20),
                     labels = seq(1980,2090,20), expand = c(0,0)) +
  scale_y_continuous(expression(paste("change since 1980")), expand = c(0,0)) +
  scale_color_manual("region",
                     values = col_reg_lim,
                     labels = c("Global","energy limited","transitional","water limited")) +
  scale_fill_manual("region",
                    values = c("Global" = "black",
                               "rgy_lim" = "black",
                               "trans" = "black",
                               "wtr_lim" = "black"),
                    labels = c("Global","energy limited","transitional","water limited")) +
  guides(fill = guide_legend(order = 1, title.position = "top", title.hjust = .5, label.position = "right", nrow=2,
                             label.theme = element_text(angle = 0))) +
  guides(col = guide_legend(order = 1, title.position = "top", title.hjust = .5, label.position = "right", nrow=2,
                            label.theme = element_text(angle = 0))) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24),
        axis.line.y.right = element_line(color = "grey80"),
        axis.ticks.y.right = element_line(color = "grey80"),
        axis.text.y.right = element_text(color = "grey80"),
        axis.title.y.right = element_text(color = "grey80")
        
  ) + ggtitle("c) Ecosystem Limitation Index")
m4_init_ELI

plots <- plot_grid(qa4, m4, ncol=2, align = 'h')

# abs 1980
sa4_lim_reg <- ggplot(global_mmmtseries_lim.df, aes(x=mmmean_ELI_max,y=mmmean_Tmaxbar_min,col=reg,fill=reg)) +
  geom_point(size = 1.25) +
  geom_hline(yintercept=0) +
  geom_vline(xintercept=0) +
  stat_smooth(method = "lm", formula=y~x) +
  scale_x_continuous("ELI change since 1980 (-)", expand = c(0,0)) +
  scale_y_continuous(expression("Temperature excess change since 1980 (K)"), expand = c(0,0)) +
  scale_color_manual("region",
                     values = col_reg_lim,
                     labels = c("Global","energy limited","transitional","water limited")) +
  scale_fill_manual("region",
                    values = col_reg_lim,
                    labels = c("Global","energy limited","transitional","water limited")) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        legend.key = element_rect(colour = NA, fill = NA),
        legend.background = element_rect(colour = NA, fill = NA),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24),
        axis.line.y.right = element_line(color = "grey80"),
        axis.ticks.y.right = element_line(color = "grey80"),
        axis.text.y.right = element_text(color = "grey80"),
        axis.title.y.right = element_text(color = "grey80")
        
  ) +
  ggtitle("d) Linear regression")
sa4_lim_reg

plots3 <- plot_grid(qa4_init_ELI, m4_init_ELI, ncol = 2, align='h')
plots4 <- grid.arrange(grid.arrange(empty, plot_prior_dcorr, empty, ncol = 3, widths = c(1.25,7.5,1.25)), plots3, grid.arrange(empty, sa4_lim_reg, empty, ncol = 3, widths = c(1.25,7.5,1.25)), nrow = 3)
# ggsave("Figures/Fig5.png", plot = plots4, width = 10*1.25, height = 15*1.25, units = "in")
ggsave("testdir/Fig5.png", plot = plots4, width = 10*1.25, height = 15*1.25, units = "in")





#######################################################################################################################
##############################################                                     ####################################
##############################################   cor: source_id separate & reg av  ####################################
##############################################                                     ####################################
#######################################################################################################################
global_tseries_boot.df <- setNames(data.frame(matrix(ncol = 8, nrow = 0)),
                                   c("ELI", "ELI_max", "Tmeanbar", "Tmaxbar_min", "center_year", "year", "source_id", "bootstrap_cycle"))
for(source_id in 1:length(cmip6_data.df$source_id)){
  # the original values
  for(t in 1:12){
    global_tseries_boot.df <- rbind(global_tseries_boot.df,
                                    data.frame("ELI" = weighted.mean(x = dcorr.list[[source_id]][,,t] * mask_full, w = area.array, na.rm = T),
                                               "ELI_max" = weighted.mean(x = dcorr_max.list[[source_id]][,,t] * mask_full, w = area.array, na.rm = T),
                                               "Tmeanbar" = weighted.mean(x = av_tas.list[[source_id]][,,t] * mask_full, w = area.array, na.rm = T),
                                               "Tmaxbar_min" = weighted.mean(x = (Tmaxbar.list[[source_id]][,,t] - av_tas.list[[source_id]][,,t]) * mask_full, w = area.array, na.rm = T),
                                               "center_year" = center_year[t],
                                               "year" = center_year[t]-5,
                                               "source_id" = cmip6_data.df$source_id[source_id],
                                               "bootstrap_cycle" = 0
                                    ))
  }
}

# bootstrap over the Globally generated pairs of lons and lats
for(bootstrap_cycle in 1:100){ # how many times bootstrapping is needed
  resample_index <- sample(x = seq(1:length(which(!is.na(mask_full)))), size = length(which(!is.na(mask_full))), replace = T)
  for(source_id in 1:length(cmip6_data.df$source_id)){
    for(t in 1:12){
      dcorr_vec <- c(dcorr.list[[source_id]][,,t][which(!is.na(mask_full))])[resample_index]
      dcorr_max_vec <- c(dcorr_max.list[[source_id]][,,t][which(!is.na(mask_full))])[resample_index]
      av_tas_vec <- c(av_tas.list[[source_id]][,,t][which(!is.na(mask_full))])[resample_index]
      Tmaxbar_vec <- c(Tmaxbar.list[[source_id]][,,t][which(!is.na(mask_full))])[resample_index]
      area_vec <- c(area.array[which(!is.na(mask_full))])[resample_index]
      global_tseries_boot.df <- rbind(global_tseries_boot.df,
                                      data.frame("ELI" = weighted.mean(x = dcorr_vec, w = area_vec, na.rm = T),
                                                 "ELI_max" = weighted.mean(x = dcorr_max_vec, w = area_vec, na.rm = T),
                                                 "Tmeanbar" = weighted.mean(x = av_tas_vec, w = area_vec, na.rm = T),
                                                 "Tmaxbar_min" = weighted.mean(x = (Tmaxbar_vec - av_tas_vec), w = area_vec, na.rm = T),
                                                 "center_year" = center_year[t],
                                                 "year" = center_year[t]-5,
                                                 "source_id" = cmip6_data.df$source_id[source_id],
                                                 "bootstrap_cycle" = bootstrap_cycle
                                      ))
    }
  }
  print(paste0("bootstrap_cycle ",bootstrap_cycle, " is done..."))
}
global_tseries_boot.df$reg <- "Global"
global_tseries_lim_boot.df <- global_tseries_boot.df


reg <- c("SAM", "NAM", "CEU", "NAS")
count_j <- 1
for(j in seq(1,16,4)){ # loop over all regions
  lonmin <- min(hotspot_regs.df$x[j:(j+3)])
  lonmax <- max(hotspot_regs.df$x[j:(j+3)])
  latmin <- min(hotspot_regs.df$y[j:(j+3)])
  latmax <- max(hotspot_regs.df$y[j:(j+3)])
  reg_mask <- array(NaN,c(180,90)); reg_mask[which(lon >= lonmin & lon <= lonmax), which(lat >= latmin & lat <= latmax)] <- 1
  for(source_id in 1:length(cmip6_data.df$source_id)){
    for(t in 1:12){ # loop over all decades
      dcorr_vec <- c(dcorr.list[[source_id]][,,t][which(!is.na(reg_mask))])
      dcorr_max_vec <- c(dcorr_max.list[[source_id]][,,t][which(!is.na(reg_mask))])
      av_tas_vec <- c(av_tas.list[[source_id]][,,t][which(!is.na(reg_mask))])
      Tmaxbar_vec <- c(Tmaxbar.list[[source_id]][,,t][which(!is.na(reg_mask))])
      area_vec <- c(area.array[which(!is.na(reg_mask))])
      global_tseries_boot.df <- rbind(global_tseries_boot.df,
                                      data.frame("ELI" = weighted.mean(x = dcorr_vec, w = area_vec, na.rm = T),
                                                 "ELI_max" = weighted.mean(x = dcorr_max_vec, w = area_vec, na.rm = T),
                                                 "Tmeanbar" = weighted.mean(x = av_tas_vec, w = area_vec, na.rm = T),
                                                 "Tmaxbar_min" = weighted.mean(x = (Tmaxbar_vec - av_tas_vec), w = area_vec, na.rm = T),
                                                 "center_year" = center_year[t],
                                                 "year" = center_year[t]-5,
                                                 "source_id" = cmip6_data.df$source_id[source_id],
                                                 "reg" = reg[count_j],
                                                 "bootstrap_cycle" = 0))
    }
  }
  count_j <- count_j + 1
}

count_j <- 1
for(bootstrap_cycle in 1:100){ # how many times bootstrapping is needed
  for(j in seq(1,16,4)){ # loop over all regions
    lonmin <- min(hotspot_regs.df$x[j:(j+3)])
    lonmax <- max(hotspot_regs.df$x[j:(j+3)])
    latmin <- min(hotspot_regs.df$y[j:(j+3)])
    latmax <- max(hotspot_regs.df$y[j:(j+3)])
    reg_mask <- array(NaN,c(180,90)); reg_mask[which(lon >= lonmin & lon <= lonmax), which(lat >= latmin & lat <= latmax)] <- 1
    resample_index <- sample(x = seq(1:length(which(!is.na(reg_mask)))), size = length(which(!is.na(reg_mask))), replace = T)
    for(source_id in 1:length(cmip6_data.df$source_id)){
      for(t in 1:12){ # loop over all decades
        dcorr_vec <- c(dcorr.list[[source_id]][,,t][which(!is.na(reg_mask))])[resample_index]
        dcorr_max_vec <- c(dcorr_max.list[[source_id]][,,t][which(!is.na(reg_mask))])[resample_index]
        av_tas_vec <- c(av_tas.list[[source_id]][,,t][which(!is.na(reg_mask))])[resample_index]
        Tmaxbar_vec <- c(Tmaxbar.list[[source_id]][,,t][which(!is.na(reg_mask))])[resample_index]
        area_vec <- c(area.array[which(!is.na(reg_mask))])[resample_index]
        global_tseries_boot.df <- rbind(global_tseries_boot.df,
                                        data.frame("ELI" = weighted.mean(x = dcorr_vec, w = area_vec, na.rm = T),
                                                   "ELI_max" = weighted.mean(x = dcorr_max_vec, w = area_vec, na.rm = T),
                                                   "Tmeanbar" = weighted.mean(x = av_tas_vec, w = area_vec, na.rm = T),
                                                   "Tmaxbar_min" = weighted.mean(x = (Tmaxbar_vec - av_tas_vec), w = area_vec, na.rm = T),
                                                   "center_year" = center_year[t],
                                                   "year" = center_year[t]-5,
                                                   "source_id" = cmip6_data.df$source_id[source_id],
                                                   "reg" = reg[count_j],
                                                   "bootstrap_cycle" = bootstrap_cycle))
      }
    }
    count_j <- count_j + 1
  }
  print(paste0("bootstrap_cycle ",bootstrap_cycle, " is done..."))
  count_j <- 1
}

cor.df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)),
                   c("reg", "source_id",
                     "cor_Tmaxbar_min_vs_ELI", "cor_Tmaxbar_min_vs_ELI_max", "bootstrap_cycle"))
# So let's try attribution with a global time series & per region & per source_id
for(bootstrap_cycle in 0:100){
  for(i in unique(global_tseries_boot.df$reg)){
    for(source_id in unique(global_tseries_boot.df$source_id)){
      boot_reg_source_id_index <- which(global_tseries_boot.df$reg == i & global_tseries_boot.df$source_id == source_id & global_tseries_boot.df$bootstrap_cycle == bootstrap_cycle)
      cor.df <- rbind(cor.df,
                      data.frame("reg" = i,
                                 "source_id" = source_id,
                                 "cor_Tmaxbar_min_vs_ELI" = cor(global_tseries_boot.df$Tmaxbar_min[boot_reg_source_id_index],
                                                                # global_tseries_boot.df$ELI[boot_reg_source_id_index], use = "pairwise.complete.obs"),
                                                                global_tseries_boot.df$ELI[boot_reg_source_id_index], method = 'kendall', use = "pairwise.complete.obs"),
                                 "cor_Tmaxbar_min_vs_ELI_max" = cor(global_tseries_boot.df$Tmaxbar_min[boot_reg_source_id_index],
                                                                    # global_tseries_boot.df$ELI[boot_reg_source_id_index], use = "pairwise.complete.obs"),
                                                                    global_tseries_boot.df$ELI_max[boot_reg_source_id_index], method = 'kendall', use = "pairwise.complete.obs"),
                                 "bootstrap_cycle" = bootstrap_cycle))
    }
  }
  print(paste0("bootstrap_cycle ",bootstrap_cycle, " is done..."))
}
cor.df$reg_f <- factor(cor.df$reg, levels = c("Global","CEU","NAM","NAS","SAM"))



# lim
global_tseries_lim_boot.df$reg <- "Global"
reg_lim <- c("rgy_lim","trans","wtr_lim")
count_j <- 1
for(j in 1:3){ # loop over all regions
  for(source_id in 1:length(cmip6_data.df$source_id)){
    for(t in 1:12){ # loop over all decades
      dcorr_vec <- c(dcorr.list[[source_id]][,,t][which(!is.na(prior_dcorr_mask[,,source_id,j]))])
      dcorr_max_vec <- c(dcorr_max.list[[source_id]][,,t][which(!is.na(prior_dcorr_mask[,,source_id,j]))])
      av_tas_vec <- c(av_tas.list[[source_id]][,,t][which(!is.na(prior_dcorr_mask[,,source_id,j]))])
      Tmaxbar_vec <- c(Tmaxbar.list[[source_id]][,,t][which(!is.na(prior_dcorr_mask[,,source_id,j]))])
      area_vec <- c(area.array[which(!is.na(prior_dcorr_mask[,,source_id,j]))])
      global_tseries_lim_boot.df <- rbind(global_tseries_lim_boot.df,
                                          data.frame("ELI" = weighted.mean(x = dcorr_vec, w = area_vec, na.rm = T),
                                                     "ELI_max" = weighted.mean(x = dcorr_max_vec, w = area_vec, na.rm = T),
                                                     "Tmeanbar" = weighted.mean(x = av_tas_vec, w = area_vec, na.rm = T),
                                                     "Tmaxbar_min" = weighted.mean(x = (Tmaxbar_vec - av_tas_vec), w = area_vec, na.rm = T),
                                                     "center_year" = center_year[t],
                                                     "year" = center_year[t]-5,
                                                     "source_id" = cmip6_data.df$source_id[source_id],
                                                     "reg" = reg_lim[count_j],
                                                     "bootstrap_cycle" = 0))
    }
  }
  count_j <- count_j + 1
}

count_j <- 1
for(bootstrap_cycle in 1:100){ # how many times bootstrapping is needed
  for(j in 1:3){ # loop over all regions
    for(source_id in 1:length(cmip6_data.df$source_id)){
      resample_index <- sample(x = seq(1:length(which(!is.na(prior_dcorr_mask[,,source_id,j])))), size = length(which(!is.na(prior_dcorr_mask[,,source_id,j]))), replace = T)
      for(t in 1:12){ # loop over all decades
        dcorr_vec <- c(dcorr.list[[source_id]][,,t][which(!is.na(prior_dcorr_mask[,,source_id,j]))])[resample_index]
        dcorr_max_vec <- c(dcorr_max.list[[source_id]][,,t][which(!is.na(prior_dcorr_mask[,,source_id,j]))])[resample_index]
        av_tas_vec <- c(av_tas.list[[source_id]][,,t][which(!is.na(prior_dcorr_mask[,,source_id,j]))])[resample_index]
        Tmaxbar_vec <- c(Tmaxbar.list[[source_id]][,,t][which(!is.na(prior_dcorr_mask[,,source_id,j]))])[resample_index]
        area_vec <- c(area.array[which(!is.na(prior_dcorr_mask[,,source_id,j]))])[resample_index]
        global_tseries_lim_boot.df <- rbind(global_tseries_lim_boot.df,
                                            data.frame("ELI" = weighted.mean(x = dcorr_vec, w = area_vec, na.rm = T),
                                                       "ELI_max" = weighted.mean(x = dcorr_max_vec, w = area_vec, na.rm = T),
                                                       "Tmeanbar" = weighted.mean(x = av_tas_vec, w = area_vec, na.rm = T),
                                                       "Tmaxbar_min" = weighted.mean(x = (Tmaxbar_vec - av_tas_vec), w = area_vec, na.rm = T),
                                                       "center_year" = center_year[t],
                                                       "year" = center_year[t]-5,
                                                       "source_id" = cmip6_data.df$source_id[source_id],
                                                       "reg" = reg_lim[count_j],
                                                       "bootstrap_cycle" = bootstrap_cycle))
      }
    }
    count_j <- count_j + 1
  }
  print(paste0("bootstrap_cycle ",bootstrap_cycle, " is done..."))
  count_j <- 1
}

cor_lim.df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)),
                       c("reg", "source_id",
                         "cor_Tmaxbar_min_vs_ELI", "cor_Tmaxbar_min_vs_ELI_max", "bootstrap_cycle"))
# So let's try attribution with a global time series & per region & per source_id
for(bootstrap_cycle in 0:100){
  for(i in unique(global_tseries_lim_boot.df$reg)){
    for(source_id in unique(global_tseries_lim_boot.df$source_id)){
      boot_reg_source_id_index <- which(global_tseries_lim_boot.df$reg == i & global_tseries_lim_boot.df$source_id == source_id & global_tseries_lim_boot.df$bootstrap_cycle == bootstrap_cycle)
      cor_lim.df <- rbind(cor_lim.df,
                          data.frame("reg" = i,
                                     "source_id" = source_id,
                                     "cor_Tmaxbar_min_vs_ELI" = cor(global_tseries_lim_boot.df$Tmaxbar_min[boot_reg_source_id_index],
                                                                    global_tseries_lim_boot.df$ELI[boot_reg_source_id_index], method = 'kendall', use = "pairwise.complete.obs"),
                                     "cor_Tmaxbar_min_vs_ELI_max" = cor(global_tseries_lim_boot.df$Tmaxbar_min[boot_reg_source_id_index],
                                                                        global_tseries_lim_boot.df$ELI_max[boot_reg_source_id_index], method = 'kendall', use = "pairwise.complete.obs"),
                                     "bootstrap_cycle" = bootstrap_cycle))
    }
  }
  print(paste0("bootstrap_cycle ",bootstrap_cycle, " is done..."))
}

cols_source_id <- c("ACCESS-ESM1-5" = colorRampPalette(brewer.pal(9, "Set1"))(12)[1],
                    "BCC-CSM2-MR" = colorRampPalette(brewer.pal(9, "Set1"))(12)[2],
                    "CMCC-ESM2" = colorRampPalette(brewer.pal(9, "Set1"))(12)[3],
                    "CNRM-CM6-1" = colorRampPalette(brewer.pal(9, "Set1"))(12)[4],
                    "CNRM-ESM2-1" = colorRampPalette(brewer.pal(9, "Set1"))(12)[5],
                    "EC-Earth3-CC" = colorRampPalette(brewer.pal(9, "Set1"))(12)[6],
                    "GFDL-ESM4" = colorRampPalette(brewer.pal(9, "Set1"))(12)[7],
                    "HadGEM3-GC31-LL" = colorRampPalette(brewer.pal(9, "Set1"))(12)[8],
                    "MPI-ESM1-2-HR" = colorRampPalette(brewer.pal(9, "Set1"))(12)[9],
                    "MPI-ESM1-2-LR" = colorRampPalette(brewer.pal(9, "Set1"))(12)[10],
                    "MRI-ESM2-0" = colorRampPalette(brewer.pal(9, "Set1"))(12)[11],
                    "UKESM1-0-LL" = colorRampPalette(brewer.pal(9, "Set1"))(12)[12])


boxplot_lim_reg <- ggplot(cor_lim.df, aes(x = reg, y = cor_Tmaxbar_min_vs_ELI_max, fill = source_id)) +
  geom_hline(yintercept=0, linetype = 'dashed', size=1.25) +
  geom_boxplot() + 
  scale_x_discrete("Region", labels = c("Global","Energy limited","Transitional","Water limited")) +
  scale_y_continuous(expression("cor(Temperature excess,ELI)"),expand=c(0,0)) +
  scale_fill_manual("",values = cols_source_id,
                    drop = F) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("c) bootstrapped x 100")
boxplot_lim_reg

bar_lim_reg <- ggplot(cor_lim.df[which(cor_lim.df$bootstrap_cycle == 0),], aes(x = reg, y = cor_Tmaxbar_min_vs_ELI_max, fill = source_id)) +
  geom_bar(stat='identity',col = 'black', width = .75, position = "dodge") +
  geom_hline(yintercept=0, linetype = 'dashed', size=1.25) +
  scale_x_discrete("Region", labels = c("Global","Energy limited","Transitional","Water limited")) +
  scale_y_continuous(expression("cor(Temperature excess,ELI)"),expand=c(0,0)) +
  scale_fill_manual("",values = cols_source_id,
                    # labels = unique(combined_cor_lim.df$source_id_combined)[c(1,6,2,7,3,8,4,9,5,10)],
                    drop = F) +
  theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.text.x = element_text(angle=45,vjust=0.5),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("b)")
bar_lim_reg
plots2 <- grid.arrange(bar_lim_reg, boxplot_lim_reg, nrow = 2)

#######################################################################################################################
##############################################                                     ####################################
##############################################          cor: mmm & reg av          ####################################
##############################################                                     ####################################
#######################################################################################################################


for(i in 1:length(global_tseries_boot.df$ELI)){
  global_tseries_boot.df$ELI_rel_1980[i] <- (global_tseries_boot.df$ELI[i] - global_tseries_boot.df$ELI[which(global_tseries_boot.df$reg == global_tseries_boot.df$reg[i] & global_tseries_boot.df$source_id == global_tseries_boot.df$source_id[i] & global_tseries_boot.df$bootstrap_cycle == global_tseries_boot.df$bootstrap_cycle[i])][1])
  global_tseries_boot.df$ELI_max_rel_1980[i] <- (global_tseries_boot.df$ELI_max[i] - global_tseries_boot.df$ELI_max[which(global_tseries_boot.df$reg == global_tseries_boot.df$reg[i] & global_tseries_boot.df$source_id == global_tseries_boot.df$source_id[i] & global_tseries_boot.df$bootstrap_cycle == global_tseries_boot.df$bootstrap_cycle[i])][1])
  global_tseries_boot.df$Tmeanbar_rel_1980[i] <- (global_tseries_boot.df$Tmeanbar[i] - global_tseries_boot.df$Tmeanbar[which(global_tseries_boot.df$reg == global_tseries_boot.df$reg[i] & global_tseries_boot.df$source_id == global_tseries_boot.df$source_id[i])][1])
  global_tseries_boot.df$Tmaxbar_min_rel_1980[i] <- (global_tseries_boot.df$Tmaxbar_min[i] - global_tseries_boot.df$Tmaxbar_min[which(global_tseries_boot.df$reg == global_tseries_boot.df$reg[i] & global_tseries_boot.df$source_id == global_tseries_boot.df$source_id[i])][1])
}

for(i in 1:length(global_tseries_lim_boot.df$ELI)){
  global_tseries_lim_boot.df$ELI_rel_1980[i] <- (global_tseries_lim_boot.df$ELI[i] - global_tseries_lim_boot.df$ELI[which(global_tseries_lim_boot.df$reg == global_tseries_lim_boot.df$reg[i] & global_tseries_lim_boot.df$source_id == global_tseries_lim_boot.df$source_id[i] & global_tseries_lim_boot.df$bootstrap_cycle == global_tseries_lim_boot.df$bootstrap_cycle[i])][1])
  global_tseries_lim_boot.df$ELI_max_rel_1980[i] <- (global_tseries_lim_boot.df$ELI_max[i] - global_tseries_lim_boot.df$ELI_max[which(global_tseries_lim_boot.df$reg == global_tseries_lim_boot.df$reg[i] & global_tseries_lim_boot.df$source_id == global_tseries_lim_boot.df$source_id[i] & global_tseries_lim_boot.df$bootstrap_cycle == global_tseries_lim_boot.df$bootstrap_cycle[i])][1])
  global_tseries_lim_boot.df$Tmeanbar_rel_1980[i] <- (global_tseries_lim_boot.df$Tmeanbar[i] - global_tseries_lim_boot.df$Tmeanbar[which(global_tseries_lim_boot.df$reg == global_tseries_lim_boot.df$reg[i] & global_tseries_lim_boot.df$source_id == global_tseries_lim_boot.df$source_id[i])][1])
  global_tseries_lim_boot.df$Tmaxbar_min_rel_1980[i] <- (global_tseries_lim_boot.df$Tmaxbar_min[i] - global_tseries_lim_boot.df$Tmaxbar_min[which(global_tseries_lim_boot.df$reg == global_tseries_lim_boot.df$reg[i] & global_tseries_lim_boot.df$source_id == global_tseries_lim_boot.df$source_id[i])][1])
}

# mmmean and mmsd per i) region and ii) center_year
global_mmmtseries_boot.df <- setNames(data.frame(matrix(ncol = 16, nrow = 0)),
                                      c("mmmean_ELI", "mmminsd_ELI", "mmplssd_ELI",
                                        "mmmean_ELI_max", "mmminsd_ELI_max", "mmplssd_ELI_max",
                                        "mmmean_Tmeanbar", "mmminsd_Tmeanbar", "mmplssd_Tmeanbar",
                                        "mmmean_Tmaxbar_min", "mmminsd_Tmaxbar_min", "mmplssd_Tmaxbar_min",
                                        "center_year","year","reg","bootstrap_cycle"))
for(bootstrap_cycle in 0:100){
  for(i_reg in unique(global_tseries_boot.df$reg)){
    for(i in unique(global_tseries_boot.df$year)){
      boot_reg_source_id_index <- which(global_tseries_boot.df$reg == i_reg & global_tseries_boot.df$year == i & global_tseries_boot.df$bootstrap_cycle == bootstrap_cycle)
      global_mmmtseries_boot.df <- rbind(global_mmmtseries_boot.df,
                                         data.frame("mmmean_ELI" = mean(global_tseries_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmminsd_ELI" = mean(global_tseries_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T) - sd(global_tseries_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmplssd_ELI" = mean(global_tseries_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T) + sd(global_tseries_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmmean_ELI_max" = mean(global_tseries_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmminsd_ELI_max" = mean(global_tseries_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T) - sd(global_tseries_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmplssd_ELI_max" = mean(global_tseries_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T) + sd(global_tseries_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmmean_Tmeanbar" = mean(global_tseries_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmminsd_Tmeanbar" = mean(global_tseries_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T) - sd(global_tseries_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmplssd_Tmeanbar" = mean(global_tseries_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T) + sd(global_tseries_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmmean_Tmaxbar_min" = mean(global_tseries_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmminsd_Tmaxbar_min" = mean(global_tseries_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T) - sd(global_tseries_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "mmplssd_Tmaxbar_min" = mean(global_tseries_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T) + sd(global_tseries_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                    "center_year" = i+5,
                                                    "year" = i,
                                                    "reg" = i_reg,
                                                    "bootstrap_cycle" = bootstrap_cycle))
    }
  }
  print(paste0("bootstrap_cycle ",bootstrap_cycle, " is done..."))
}

# mmmean and mmsd per i) region and ii) center_year
global_mmmtseries_lim_boot.df <- setNames(data.frame(matrix(ncol = 16, nrow = 0)),
                                          c("mmmean_ELI", "mmminsd_ELI", "mmplssd_ELI",
                                            "mmmean_ELI_max", "mmminsd_ELI_max", "mmplssd_ELI_max",
                                            "mmmean_Tmeanbar", "mmminsd_Tmeanbar", "mmplssd_Tmeanbar",
                                            "mmmean_Tmaxbar_min", "mmminsd_Tmaxbar_min", "mmplssd_Tmaxbar_min",
                                            "center_year","year","reg","bootstrap_cycle"))
for(bootstrap_cycle in 0:100){
  for(i_reg in unique(global_tseries_lim_boot.df$reg)){
    for(i in unique(global_tseries_lim_boot.df$year)){
      boot_reg_source_id_index <- which(global_tseries_lim_boot.df$reg == i_reg & global_tseries_lim_boot.df$year == i & global_tseries_lim_boot.df$bootstrap_cycle == bootstrap_cycle)
      global_mmmtseries_lim_boot.df <- rbind(global_mmmtseries_lim_boot.df,
                                             data.frame("mmmean_ELI" = mean(global_tseries_lim_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmminsd_ELI" = mean(global_tseries_lim_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T) - sd(global_tseries_lim_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmplssd_ELI" = mean(global_tseries_lim_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T) + sd(global_tseries_lim_boot.df$ELI_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmmean_ELI_max" = mean(global_tseries_lim_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmminsd_ELI_max" = mean(global_tseries_lim_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T) - sd(global_tseries_lim_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmplssd_ELI_max" = mean(global_tseries_lim_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T) + sd(global_tseries_lim_boot.df$ELI_max_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmmean_Tmeanbar" = mean(global_tseries_lim_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmminsd_Tmeanbar" = mean(global_tseries_lim_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T) - sd(global_tseries_lim_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmplssd_Tmeanbar" = mean(global_tseries_lim_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T) + sd(global_tseries_lim_boot.df$Tmeanbar_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmmean_Tmaxbar_min" = mean(global_tseries_lim_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmminsd_Tmaxbar_min" = mean(global_tseries_lim_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T) - sd(global_tseries_lim_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "mmplssd_Tmaxbar_min" = mean(global_tseries_lim_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T) + sd(global_tseries_lim_boot.df$Tmaxbar_min_rel_1980[boot_reg_source_id_index], na.rm = T),
                                                        "center_year" = i+5,
                                                        "year" = i,
                                                        "reg" = i_reg,
                                                        "bootstrap_cycle" = bootstrap_cycle))
    }
  }
  print(paste0("bootstrap_cycle ",bootstrap_cycle, " is done..."))
}

mmmcor_boot.df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
                           c("reg",
                             "cor_Tmaxbar_min_vs_ELI",
                             "cor_Tmaxbar_min_vs_ELI_max",
                             "bootstrap_cycle"))
for(bootstrap_cycle in 0:100){
  for(i in unique(global_tseries.df$reg)){
    boot_reg_index <- which(global_mmmtseries_boot.df$reg == i & global_mmmtseries_boot.df$bootstrap_cycle == bootstrap_cycle)
    mmmcor_boot.df <- rbind(mmmcor_boot.df,
                            data.frame("reg" = i,
                                       "cor_Tmaxbar_min_vs_ELI" = cor(global_mmmtseries_boot.df$mmmean_Tmaxbar_min[boot_reg_index],
                                                                      global_mmmtseries_boot.df$mmmean_ELI[boot_reg_index], method = 'kendall', use = "pairwise.complete.obs"),
                                       "cor_Tmaxbar_min_vs_ELI_max" = cor(global_mmmtseries_boot.df$mmmean_Tmaxbar_min[boot_reg_index],
                                                                          global_mmmtseries_boot.df$mmmean_ELI_max[boot_reg_index], method = 'kendall', use = "pairwise.complete.obs"),
                                       "bootstrap_cycle" = bootstrap_cycle))
    
  }
  print(paste0("bootstrap_cycle ",bootstrap_cycle, " is done..."))
}
mmmcor_boot.df$reg_f <- factor(mmmcor_boot.df$reg, levels = c("Global","CEU","NAM","NAS","SAM"))

mmmcor_lim_boot.df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
                               c("reg",
                                 "cor_Tmaxbar_min_vs_ELI",
                                 "cor_Tmaxbar_min_vs_ELI_max",
                                 "bootstrap_cycle"))
for(bootstrap_cycle in 0:100){
  for(i in unique(global_tseries_lim.df$reg)){
    boot_reg_index <- which(global_mmmtseries_lim_boot.df$reg == i & global_mmmtseries_lim_boot.df$bootstrap_cycle == bootstrap_cycle)
    mmmcor_lim_boot.df <- rbind(mmmcor_lim_boot.df,
                                data.frame("reg" = i,
                                           "cor_Tmaxbar_min_vs_ELI" = cor(global_mmmtseries_lim_boot.df$mmmean_Tmaxbar_min[boot_reg_index],
                                                                          global_mmmtseries_lim_boot.df$mmmean_ELI[boot_reg_index], method = 'kendall', use = "pairwise.complete.obs"),
                                           "cor_Tmaxbar_min_vs_ELI_max" = cor(global_mmmtseries_lim_boot.df$mmmean_Tmaxbar_min[boot_reg_index],
                                                                              global_mmmtseries_lim_boot.df$mmmean_ELI_max[boot_reg_index], method = 'kendall', use = "pairwise.complete.obs"),
                                           "bootstrap_cycle" = bootstrap_cycle))
    
  }
  print(paste0("bootstrap_cycle ",bootstrap_cycle, " is done..."))
}
mmmcor_lim_boot.df$reg_f <- factor(mmmcor_lim_boot.df$reg, levels = c("Global","rgy_lim","trans","wtr_lim"))


# boxplot_reg <- ggplot(mmmcor_boot.df, aes(x = reg_f, y = cor_Tmaxbar_min_vs_ELI, fill = reg_f)) +
#   geom_boxplot() +
#   geom_point(data = mmmcor_boot.df[which(mmmcor_boot.df$bootstrap_cycle == 0),], shape = 4, size  = 5) +
#   scale_x_discrete("Region") +
#   scale_y_continuous(expression("cor(Temperature excess,ELI)"),expand=c(0,0)) +
#   scale_fill_manual("",values = col_reg,
#                     drop = F) +
#   theme(legend.position = "none",
#         legend.text = element_text(size=12),
#         legend.title = element_text(size=20),
#         strip.background = element_blank(),
#         axis.line.x = element_line(colour = "black"),
#         axis.line.y = element_line(colour = "black"),
#         strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_rect(color = 'black', fill = NA, size = 1),
#         panel.background = element_blank(),
#         axis.text = element_text(size=18),
#         axis.title = element_text(size=20),
#         plot.title = element_text(size=24),
#         plot.tag.position = c(.55,0.03)
#   ) + ggtitle("b)")
# boxplot_reg
# 
# plots3 <- grid.arrange(sa4_reg, boxplot_reg, ncol = 2)

mmmboxplot_lim_reg <- ggplot(mmmcor_lim_boot.df, aes(x = reg_f, y = cor_Tmaxbar_min_vs_ELI, fill = reg_f)) +
  geom_boxplot() +
  geom_point(data = mmmcor_lim_boot.df[which(mmmcor_lim_boot.df$bootstrap_cycle == 0),], shape = 4, size  = 5) +
  scale_x_discrete("",
                   labels = c("rgy_lim" = "energy limited",
                              "trans" = "transitional",
                              "wtr_lim" = "water limited")) +
  scale_y_continuous(expression("cor(Temperature excess, ELI)"),expand=c(0,0)) +
  scale_fill_manual("",values = col_reg_lim,
                    # labels = unique(combined_mmmcor_lim_boot.df$linreg)[c(1,2)],
                    drop = F) +
  theme(legend.position = "none",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        strip.background = element_blank(),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        strip.text.x = element_text(size=18, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.text.x = element_text(angle=45,vjust=0.5),
        axis.title = element_text(size=20),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("a)")
mmmboxplot_lim_reg

plots3 <- plot_grid(mmmboxplot_lim_reg, bar_lim_reg, ncol = 2, rel_widths = c(1,2), align='h')
plots4 <- grid.arrange(plots3, boxplot_lim_reg, nrow = 2)
# ggsave("Figures/SFig10.png", plot = plots4, width = 10*1.25, height = 8*1.25, units = "in")
ggsave("testdir/SFig10.png", plot = plots4, width = 10*1.25, height = 8*1.25, units = "in")



#######################################################################################################################
##############################################                                     ####################################
##############################################          bars: per source_id        ####################################
##############################################                                     ####################################
#######################################################################################################################
# bounds for bar charts
dcorr_trend_bounds <- seq(-.05,.15,.05)
prior_dcorr_bounds <- c(-1,-.2,.2,1)
min_points <- 10

prior_ELI <- prior_ELI_rsds <- prior_ELI_max <- array(NaN,c(180,90,12))
for(x in 1:180){
  for(y in 1:90){
    for(source_id in 1:length(cmip6_data.df$source_id)){
      prior_ELI[x,y,source_id] <- mean(dcorr.list[[source_id]][x,y,1:3])
      prior_ELI_max[x,y,source_id] <- mean(dcorr_max.list[[source_id]][x,y,1:3])
      prior_ELI_rsds[x,y,source_id] <- mean(dcorr_rsds.list[[source_id]][x,y,1:3])
    }
  }
}

# make a models_with_full_timeseries_mask with c(180,90,9)
models_with_full_timeseries_mask <- array(NaN,c(180,90))
models_with_full_timeseries_mask[which(models_with_full_timeseries > 5)] <- 1
models_with_full_timeseries_mask_rep <- area.array_rep <- array(NaN,c(180,90,12))
for(i in 1:length(cmip6_data.df$source_id)){
  models_with_full_timeseries_mask_rep[,,i] <- models_with_full_timeseries_mask
  area.array_rep[,,i] <- area.array
}


ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id.df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
                                                             c("Tmaxbar_min","kendall_ELI","reg","source_id"))
ELI_rsds_Tmaxbar_min_bars_10yr_1980_2100_source_id.df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
                                                                  c("Tmaxbar_min","kendall_ELI","reg","source_id"))
ELI_max_Tmaxbar_min_bars_10yr_1980_2100_source_id.df <- setNames(data.frame(matrix(ncol = 4, nrow = 0)),
                                                                 c("Tmaxbar_min","kendall_ELI","reg","source_id"))
# source_id
for(source_id in 1:length(cmip6_data.df$source_id)){
  Tmaxbar_min_lim_per_ELI_trend <- calc_boxes_wts(x_array = kendall_dcorr[,,source_id,1]*models_with_full_timeseries_mask, y_array = prior_ELI[,,source_id]*models_with_full_timeseries_mask, col_array = kendall_Tmaxbar_min[,,source_id,1]*models_with_full_timeseries_mask,
                                                  func = "weighted mean", x_bounds = dcorr_trend_bounds, y_bounds = prior_dcorr_bounds, min_points = min_points, min_models = 0, wts_array = area.array)
  Tmaxbar_min_lim_per_ELI_rsds_trend <- calc_boxes_wts(x_array = kendall_dcorr_rsds[,,source_id,1]*models_with_full_timeseries_mask, y_array = prior_ELI_rsds[,,source_id]*models_with_full_timeseries_mask, col_array = kendall_Tmaxbar_min[,,source_id,1]*models_with_full_timeseries_mask,
                                                       func = "weighted mean", x_bounds = dcorr_trend_bounds, y_bounds = prior_dcorr_bounds, min_points = min_points, min_models = 0, wts_array = area.array)
  Tmaxbar_min_lim_per_ELI_max_trend <- calc_boxes_wts(x_array = kendall_dcorr_max[,,source_id,1]*models_with_full_timeseries_mask, y_array = prior_ELI_max[,,source_id]*models_with_full_timeseries_mask, col_array = kendall_Tmaxbar_min[,,source_id,1]*models_with_full_timeseries_mask,
                                                      func = "weighted mean", x_bounds = dcorr_trend_bounds, y_bounds = prior_dcorr_bounds, min_points = min_points, min_models = 0, wts_array = area.array)
  
  Tmaxbar_min_glob_per_ELI_trend <- 
    Tmaxbar_min_glob_per_ELI_rsds_trend <- 
    Tmaxbar_min_glob_per_ELI_max_trend <- 
    c()
  for(i in 1:4){
    Tmaxbar_min_glob_per_ELI_trend[i] <- weighted.mean(x = (kendall_Tmaxbar_min[,,source_id,1]*models_with_full_timeseries_mask)[which(kendall_dcorr[,,source_id,1]*models_with_full_timeseries_mask > dcorr_trend_bounds[i] &
                                                                                                                                         kendall_dcorr[,,source_id,1]*models_with_full_timeseries_mask < dcorr_trend_bounds[i+1])],
                                                       w = (area.array)[which(kendall_dcorr[,,source_id,1]*models_with_full_timeseries_mask > dcorr_trend_bounds[i] &
                                                                                kendall_dcorr[,,source_id,1]*models_with_full_timeseries_mask < dcorr_trend_bounds[i+1])],
                                                       na.rm = T)
    
    Tmaxbar_min_glob_per_ELI_rsds_trend[i] <- weighted.mean(x = (kendall_Tmaxbar_min[,,source_id,1]*models_with_full_timeseries_mask)[which(kendall_dcorr_rsds[,,source_id,1]*models_with_full_timeseries_mask > dcorr_trend_bounds[i] &
                                                                                                                                              kendall_dcorr_rsds[,,source_id,1]*models_with_full_timeseries_mask < dcorr_trend_bounds[i+1])],
                                                            w = (area.array)[which(kendall_dcorr_rsds[,,source_id,1]*models_with_full_timeseries_mask > dcorr_trend_bounds[i] &
                                                                                     kendall_dcorr_rsds[,,source_id,1]*models_with_full_timeseries_mask < dcorr_trend_bounds[i+1])],
                                                            na.rm = T)
    
    Tmaxbar_min_glob_per_ELI_max_trend[i] <- weighted.mean(x = (kendall_Tmaxbar_min[,,source_id,1]*models_with_full_timeseries_mask)[which(kendall_dcorr_max[,,source_id,1]*models_with_full_timeseries_mask > dcorr_trend_bounds[i] &
                                                                                                                                             kendall_dcorr_max[,,source_id,1]*models_with_full_timeseries_mask < dcorr_trend_bounds[i+1])],
                                                           w = (area.array)[which(kendall_dcorr_max[,,source_id,1]*models_with_full_timeseries_mask > dcorr_trend_bounds[i] &
                                                                                    kendall_dcorr_max[,,source_id,1]*models_with_full_timeseries_mask < dcorr_trend_bounds[i+1])],
                                                           na.rm = T)
  }
  
  ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id.df <- rbind(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id.df, 
                                                            data.frame("Tmaxbar_min" = Tmaxbar_min_glob_per_ELI_trend,
                                                                       "kendall_ELI" = c(1:4),
                                                                       "reg" = rep("Global",4),
                                                                       "source_id" = cmip6_data.df$source_id[source_id]))
  ELI_rsds_Tmaxbar_min_bars_10yr_1980_2100_source_id.df <- rbind(ELI_rsds_Tmaxbar_min_bars_10yr_1980_2100_source_id.df, 
                                                                 data.frame("Tmaxbar_min" = Tmaxbar_min_glob_per_ELI_rsds_trend,
                                                                            "kendall_ELI" = c(1:4),
                                                                            "reg" = rep("Global",4),
                                                                            "source_id" = cmip6_data.df$source_id[source_id]))
  ELI_max_Tmaxbar_min_bars_10yr_1980_2100_source_id.df <- rbind(ELI_max_Tmaxbar_min_bars_10yr_1980_2100_source_id.df, 
                                                                data.frame("Tmaxbar_min" = Tmaxbar_min_glob_per_ELI_max_trend,
                                                                           "kendall_ELI" = c(1:4),
                                                                           "reg" = rep("Global",4),
                                                                           "source_id" = cmip6_data.df$source_id[source_id]))
  for(x in 1:4){
    for(y in 1:3){
      ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id.df <- rbind(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id.df,
                                                                data.frame("Tmaxbar_min" = Tmaxbar_min_lim_per_ELI_trend$val_per_box[x,y],
                                                                           "kendall_ELI" = x,
                                                                           "reg" = reg_lim[y],
                                                                           "source_id" = cmip6_data.df$source_id[source_id]))
      ELI_rsds_Tmaxbar_min_bars_10yr_1980_2100_source_id.df <- rbind(ELI_rsds_Tmaxbar_min_bars_10yr_1980_2100_source_id.df,
                                                                     data.frame("Tmaxbar_min" = Tmaxbar_min_lim_per_ELI_rsds_trend$val_per_box[x,y],
                                                                                "kendall_ELI" = x,
                                                                                "reg" = reg_lim[y],
                                                                                "source_id" = cmip6_data.df$source_id[source_id]))
      ELI_max_Tmaxbar_min_bars_10yr_1980_2100_source_id.df <- rbind(ELI_max_Tmaxbar_min_bars_10yr_1980_2100_source_id.df,
                                                                    data.frame("Tmaxbar_min" = Tmaxbar_min_lim_per_ELI_max_trend$val_per_box[x,y],
                                                                               "kendall_ELI" = x,
                                                                               "reg" = reg_lim[y],
                                                                               "source_id" = cmip6_data.df$source_id[source_id]))
    }
  }
}

count <- 1
for(source_id in 1:length(cmip6_data.df$source_id)){
  ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id.df$area[(count):(count+15)] <- c(100*(sum(area.array_rep[which(kendall_dcorr[,,source_id,1] > -.05 & kendall_dcorr[,,source_id,1] < 0)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(kendall_dcorr[,,source_id,1] > 0 & kendall_dcorr[,,source_id,1] < .05)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(kendall_dcorr[,,source_id,1] > .05 & kendall_dcorr[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(kendall_dcorr[,,source_id,1] > 0.1 & kendall_dcorr[,,source_id,1] < 0.15)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > -2 & prior_ELI[,,source_id] < -.2 & kendall_dcorr[,,source_id,1] > -.05 & kendall_dcorr[,,source_id,1] < 0)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > -.2 & prior_ELI[,,source_id] < .2 & kendall_dcorr[,,source_id,1] > -.05 & kendall_dcorr[,,source_id,1] < 0)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > .2 & prior_ELI[,,source_id] < 2 & kendall_dcorr[,,source_id,1] > -.05 & kendall_dcorr[,,source_id,1] < 0)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > -2 & prior_ELI[,,source_id] < -.2 & kendall_dcorr[,,source_id,1] > 0 & kendall_dcorr[,,source_id,1] < .05)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > -.2 & prior_ELI[,,source_id] < .2 & kendall_dcorr[,,source_id,1] > 0 & kendall_dcorr[,,source_id,1] < .05)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > .2 & prior_ELI[,,source_id] < 2 & kendall_dcorr[,,source_id,1] > 0 & kendall_dcorr[,,source_id,1] < .05)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > -2 & prior_ELI[,,source_id] < -.2 & kendall_dcorr[,,source_id,1] > .05 & kendall_dcorr[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > -.2 & prior_ELI[,,source_id] < .2 & kendall_dcorr[,,source_id,1] > .05 & kendall_dcorr[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > .2 & prior_ELI[,,source_id] < 2 & kendall_dcorr[,,source_id,1] > .05 & kendall_dcorr[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > -2 & prior_ELI[,,source_id] < -.2 & kendall_dcorr[,,source_id,1] > 0.1 & kendall_dcorr[,,source_id,1] < 0.15)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > -.2 & prior_ELI[,,source_id] < .2 & kendall_dcorr[,,source_id,1] > 0.1 & kendall_dcorr[,,source_id,1] < 0.15)])/(total_land_area)),
                                                                                 100*(sum(area.array_rep[which(prior_ELI[,,source_id] > .2 & prior_ELI[,,source_id] < 2 & kendall_dcorr[,,source_id,1] > 0.1 & kendall_dcorr[,,source_id,1] < 0.15)])/(total_land_area)))
  count <- count + 16
}

count <- 1
for(source_id in 1:length(cmip6_data.df$source_id)){
  ELI_rsds_Tmaxbar_min_bars_10yr_1980_2100_source_id.df$area[(count):(count+15)] <- c(100*(sum(area.array_rep[which(kendall_dcorr_rsds[,,source_id,1] > -.05 & kendall_dcorr_rsds[,,source_id,1] < 0)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(kendall_dcorr_rsds[,,source_id,1] > 0 & kendall_dcorr_rsds[,,source_id,1] < .05)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(kendall_dcorr_rsds[,,source_id,1] > .05 & kendall_dcorr_rsds[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(kendall_dcorr_rsds[,,source_id,1] > 0.1 & kendall_dcorr_rsds[,,source_id,1] < 0.15)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > -2 & prior_ELI_rsds[,,source_id] < -.2 & kendall_dcorr_rsds[,,source_id,1] > -.05 & kendall_dcorr_rsds[,,source_id,1] < 0)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > -.2 & prior_ELI_rsds[,,source_id] < .2 & kendall_dcorr_rsds[,,source_id,1] > -.05 & kendall_dcorr_rsds[,,source_id,1] < 0)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > .2 & prior_ELI_rsds[,,source_id] < 2 & kendall_dcorr_rsds[,,source_id,1] > -.05 & kendall_dcorr_rsds[,,source_id,1] < 0)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > -2 & prior_ELI_rsds[,,source_id] < -.2 & kendall_dcorr_rsds[,,source_id,1] > 0 & kendall_dcorr_rsds[,,source_id,1] < .05)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > -.2 & prior_ELI_rsds[,,source_id] < .2 & kendall_dcorr_rsds[,,source_id,1] > 0 & kendall_dcorr_rsds[,,source_id,1] < .05)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > .2 & prior_ELI_rsds[,,source_id] < 2 & kendall_dcorr_rsds[,,source_id,1] > 0 & kendall_dcorr_rsds[,,source_id,1] < .05)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > -2 & prior_ELI_rsds[,,source_id] < -.2 & kendall_dcorr_rsds[,,source_id,1] > .05 & kendall_dcorr_rsds[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > -.2 & prior_ELI_rsds[,,source_id] < .2 & kendall_dcorr_rsds[,,source_id,1] > .05 & kendall_dcorr_rsds[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > .2 & prior_ELI_rsds[,,source_id] < 2 & kendall_dcorr_rsds[,,source_id,1] > .05 & kendall_dcorr_rsds[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > -2 & prior_ELI_rsds[,,source_id] < -.2 & kendall_dcorr_rsds[,,source_id,1] > 0.1 & kendall_dcorr_rsds[,,source_id,1] < 0.15)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > -.2 & prior_ELI_rsds[,,source_id] < .2 & kendall_dcorr_rsds[,,source_id,1] > 0.1 & kendall_dcorr_rsds[,,source_id,1] < 0.15)])/(total_land_area)),
                                                                                      100*(sum(area.array_rep[which(prior_ELI_rsds[,,source_id] > .2 & prior_ELI_rsds[,,source_id] < 2 & kendall_dcorr_rsds[,,source_id,1] > 0.1 & kendall_dcorr_rsds[,,source_id,1] < 0.15)])/(total_land_area)))
  count <- count + 16
}

count <- 1
for(source_id in 1:length(cmip6_data.df$source_id)){
  ELI_max_Tmaxbar_min_bars_10yr_1980_2100_source_id.df$area[(count):(count+15)] <- c(100*(sum(area.array_rep[which(kendall_dcorr_max[,,source_id,1] > -.05 & kendall_dcorr_max[,,source_id,1] < 0)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(kendall_dcorr_max[,,source_id,1] > 0 & kendall_dcorr_max[,,source_id,1] < .05)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(kendall_dcorr_max[,,source_id,1] > .05 & kendall_dcorr_max[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(kendall_dcorr_max[,,source_id,1] > 0.1 & kendall_dcorr_max[,,source_id,1] < 0.15)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > -2 & prior_ELI_max[,,source_id] < -.2 & kendall_dcorr_max[,,source_id,1] > -.05 & kendall_dcorr_max[,,source_id,1] < 0)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > -.2 & prior_ELI_max[,,source_id] < .2 & kendall_dcorr_max[,,source_id,1] > -.05 & kendall_dcorr_max[,,source_id,1] < 0)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > .2 & prior_ELI_max[,,source_id] < 2 & kendall_dcorr_max[,,source_id,1] > -.05 & kendall_dcorr_max[,,source_id,1] < 0)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > -2 & prior_ELI_max[,,source_id] < -.2 & kendall_dcorr_max[,,source_id,1] > 0 & kendall_dcorr_max[,,source_id,1] < .05)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > -.2 & prior_ELI_max[,,source_id] < .2 & kendall_dcorr_max[,,source_id,1] > 0 & kendall_dcorr_max[,,source_id,1] < .05)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > .2 & prior_ELI_max[,,source_id] < 2 & kendall_dcorr_max[,,source_id,1] > 0 & kendall_dcorr_max[,,source_id,1] < .05)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > -2 & prior_ELI_max[,,source_id] < -.2 & kendall_dcorr_max[,,source_id,1] > .05 & kendall_dcorr_max[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > -.2 & prior_ELI_max[,,source_id] < .2 & kendall_dcorr_max[,,source_id,1] > .05 & kendall_dcorr_max[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > .2 & prior_ELI_max[,,source_id] < 2 & kendall_dcorr_max[,,source_id,1] > .05 & kendall_dcorr_max[,,source_id,1] < 0.1)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > -2 & prior_ELI_max[,,source_id] < -.2 & kendall_dcorr_max[,,source_id,1] > 0.1 & kendall_dcorr_max[,,source_id,1] < 0.15)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > -.2 & prior_ELI_max[,,source_id] < .2 & kendall_dcorr_max[,,source_id,1] > 0.1 & kendall_dcorr_max[,,source_id,1] < 0.15)])/(total_land_area)),
                                                                                     100*(sum(area.array_rep[which(prior_ELI_max[,,source_id] > .2 & prior_ELI_max[,,source_id] < 2 & kendall_dcorr_max[,,source_id,1] > 0.1 & kendall_dcorr_max[,,source_id,1] < 0.15)])/(total_land_area)))
  count <- count + 16
}


save(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id.df, col_reg_lim,
     ELI_rsds_Tmaxbar_min_bars_10yr_1980_2100_source_id.df, 
     ELI_max_Tmaxbar_min_bars_10yr_1980_2100_source_id.df, 
     file = paste0(path_RData, "Fig6a_source_id.RData"))


#######################################################################################################################
##############################################                                     ####################################
##############################################  bars: per source_id but together   ####################################
##############################################                                     ####################################
#######################################################################################################################
# I'll try to make bar charts that show Tmaxbar_min per class of ELI trend. The black bars should show the average across all regimes and then 
# maybe to add are rgy-lim (darkgreen), trans (grey) and brown (wtr-lim) bars according to prior ELI

# bounds for bar charts
dcorr_trend_bounds <- seq(-.05,.15,.05)
prior_dcorr_bounds <- c(-2,-.2,.2,2)
min_points <- 10

prior_ELI <- prior_ELI_rsds <- prior_ELI_max <- array(NaN,c(180,90,12))
for(x in 1:180){
  for(y in 1:90){
    for(source_id in 1:length(cmip6_data.df$source_id)){
      prior_ELI[x,y,source_id] <- mean(dcorr.list[[source_id]][x,y,1:3])
      prior_ELI_max[x,y,source_id] <- mean(dcorr_max.list[[source_id]][x,y,1:3])
      prior_ELI_rsds[x,y,source_id] <- mean(dcorr_rsds.list[[source_id]][x,y,1:3])
    }
  }
}

# make a models_with_full_timeseries_mask with c(180,90,9) 
models_with_full_timeseries_mask <- array(NaN,c(180,90))
models_with_full_timeseries_mask[which(models_with_full_timeseries > 5)] <- 1
models_with_full_timeseries_mask_rep <- area.array_rep <- array(NaN,c(180,90,12))
for(i in 1:length(cmip6_data.df$source_id)){
  models_with_full_timeseries_mask_rep[,,i] <- models_with_full_timeseries_mask
  area.array_rep[,,i] <- area.array
}

# source_id
Tmaxbar_min_lim_per_ELI_trend <- calc_boxes_wts(x_array = kendall_dcorr[,,,1]*models_with_full_timeseries_mask_rep, y_array = prior_ELI*models_with_full_timeseries_mask_rep, col_array = kendall_Tmaxbar_min[,,,1]*models_with_full_timeseries_mask_rep,
                                                func = "weighted mean", x_bounds = dcorr_trend_bounds, y_bounds = prior_dcorr_bounds, min_points = min_points, min_models = 0, wts_array = area.array_rep)
Tmaxbar_min_lim_per_ELI_trend_rsds <- calc_boxes_wts(x_array = kendall_dcorr_rsds[,,,1]*models_with_full_timeseries_mask_rep, y_array = prior_ELI_rsds*models_with_full_timeseries_mask_rep, col_array = kendall_Tmaxbar_min[,,,1]*models_with_full_timeseries_mask_rep,
                                                     func = "weighted mean", x_bounds = dcorr_trend_bounds, y_bounds = prior_dcorr_bounds, min_points = min_points, min_models = 0, wts_array = area.array_rep)
Tmaxbar_min_lim_per_ELI_trend_max <- calc_boxes_wts(x_array = kendall_dcorr_max[,,,1]*models_with_full_timeseries_mask_rep, y_array = prior_ELI_max*models_with_full_timeseries_mask_rep, col_array = kendall_Tmaxbar_min[,,,1]*models_with_full_timeseries_mask_rep,
                                                    func = "weighted mean", x_bounds = dcorr_trend_bounds, y_bounds = prior_dcorr_bounds, min_points = min_points, min_models = 0, wts_array = area.array_rep)

Tmaxbar_min_glob_per_ELI_trend <- Tmaxbar_min_glob_per_ELI_trend_rsds <- Tmaxbar_min_glob_per_ELI_trend_max <- c()
for(i in 1:4){
  Tmaxbar_min_glob_per_ELI_trend[i] <- weighted.mean(x = (kendall_Tmaxbar_min[,,,1]*models_with_full_timeseries_mask_rep)[which(kendall_dcorr[,,,1]*models_with_full_timeseries_mask_rep > dcorr_trend_bounds[i] &
                                                                                                                                  kendall_dcorr[,,,1]*models_with_full_timeseries_mask_rep < dcorr_trend_bounds[i+1])],
                                                     w = (area.array_rep)[which(kendall_dcorr[,,,1]*models_with_full_timeseries_mask_rep > dcorr_trend_bounds[i] &
                                                                                  kendall_dcorr[,,,1]*models_with_full_timeseries_mask_rep < dcorr_trend_bounds[i+1])],
                                                     na.rm = T)
  Tmaxbar_min_glob_per_ELI_trend_rsds[i] <- weighted.mean(x = (kendall_Tmaxbar_min[,,,1]*models_with_full_timeseries_mask_rep)[which(kendall_dcorr_rsds[,,,1]*models_with_full_timeseries_mask_rep > dcorr_trend_bounds[i] &
                                                                                                                                       kendall_dcorr_rsds[,,,1]*models_with_full_timeseries_mask_rep < dcorr_trend_bounds[i+1])],
                                                          w = (area.array_rep)[which(kendall_dcorr_rsds[,,,1]*models_with_full_timeseries_mask_rep > dcorr_trend_bounds[i] &
                                                                                       kendall_dcorr_rsds[,,,1]*models_with_full_timeseries_mask_rep < dcorr_trend_bounds[i+1])],
                                                          na.rm = T)
  Tmaxbar_min_glob_per_ELI_trend_max[i] <- weighted.mean(x = (kendall_Tmaxbar_min[,,,1]*models_with_full_timeseries_mask_rep)[which(kendall_dcorr_max[,,,1]*models_with_full_timeseries_mask_rep > dcorr_trend_bounds[i] &
                                                                                                                                      kendall_dcorr_max[,,,1]*models_with_full_timeseries_mask_rep < dcorr_trend_bounds[i+1])],
                                                         w = (area.array_rep)[which(kendall_dcorr_max[,,,1]*models_with_full_timeseries_mask_rep > dcorr_trend_bounds[i] &
                                                                                      kendall_dcorr_max[,,,1]*models_with_full_timeseries_mask_rep < dcorr_trend_bounds[i+1])],
                                                         na.rm = T)
}

ELI_Tmaxbar_min_bars_10yr_1980_2100.df <- data.frame("Tmaxbar_min" = c(Tmaxbar_min_glob_per_ELI_trend,Tmaxbar_min_glob_per_ELI_trend_rsds, Tmaxbar_min_glob_per_ELI_trend_max),
                                                     "kendall_ELI" = rep(c(1:4),3),
                                                     "reg" = rep("Global",12),
                                                     "dcorr_with" = rep(c("tas","rsds","max"),each=4))
reg_lim <- c("rgy_lim","trans","wtr_lim")
for(x in 1:4){
  for(y in 1:3){
    ELI_Tmaxbar_min_bars_10yr_1980_2100.df <- rbind(ELI_Tmaxbar_min_bars_10yr_1980_2100.df,
                                                    data.frame("Tmaxbar_min" = c(Tmaxbar_min_lim_per_ELI_trend$val_per_box[x,y],Tmaxbar_min_lim_per_ELI_trend_rsds$val_per_box[x,y],Tmaxbar_min_lim_per_ELI_trend_max$val_per_box[x,y]),
                                                               "kendall_ELI" = rep(x,3),
                                                               "reg" = rep(reg_lim[y],3),
                                                               "dcorr_with" = c("tas","rsds","max")))
  }
}

col_reg_lim <- c("Global" = "black",
                 "rgy_lim" = brewer.pal(9,"BrBG")[8],
                 "trans" = 'slategrey',
                 "wtr_lim" = brewer.pal(9,"BrBG")[2])

ELI_Tmaxbar_min_bars_10yr_1980_2100.df$area <- c(100*(sum(area.array_rep[which(kendall_dcorr[,,,1] > -.05 & kendall_dcorr[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr[,,,1] > 0 & kendall_dcorr[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr[,,,1] > .05 & kendall_dcorr[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr[,,,1] > 0.1 & kendall_dcorr[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr_rsds[,,,1] > -.05 & kendall_dcorr_rsds[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr_rsds[,,,1] > 0 & kendall_dcorr_rsds[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr_rsds[,,,1] > .05 & kendall_dcorr_rsds[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr_rsds[,,,1] > 0.1 & kendall_dcorr_rsds[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr_max[,,,1] > -.05 & kendall_dcorr_max[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr_max[,,,1] > 0 & kendall_dcorr_max[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr_max[,,,1] > .05 & kendall_dcorr_max[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(kendall_dcorr_max[,,,1] > 0.1 & kendall_dcorr_max[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > -2 & prior_ELI < -.2 & kendall_dcorr[,,,1] > -.05 & kendall_dcorr[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > -2 & prior_ELI_rsds < -.2 & kendall_dcorr_rsds[,,,1] > -.05 & kendall_dcorr_rsds[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > -2 & prior_ELI_max < -.2 & kendall_dcorr_max[,,,1] > -.05 & kendall_dcorr_max[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > -.2 & prior_ELI < .2 & kendall_dcorr[,,,1] > -.05 & kendall_dcorr[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > -.2 & prior_ELI_rsds < .2 & kendall_dcorr_rsds[,,,1] > -.05 & kendall_dcorr_rsds[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > -.2 & prior_ELI_max < .2 & kendall_dcorr_max[,,,1] > -.05 & kendall_dcorr_max[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > .2 & prior_ELI < 2 & kendall_dcorr[,,,1] > -.05 & kendall_dcorr[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > .2 & prior_ELI_rsds < 2 & kendall_dcorr_rsds[,,,1] > -.05 & kendall_dcorr_rsds[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > .2 & prior_ELI_max < 2 & kendall_dcorr_max[,,,1] > -.05 & kendall_dcorr_max[,,,1] < 0)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > -2 & prior_ELI < -.2 & kendall_dcorr[,,,1] > 0 & kendall_dcorr[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > -2 & prior_ELI_rsds < -.2 & kendall_dcorr_rsds[,,,1] > 0 & kendall_dcorr_rsds[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > -2 & prior_ELI_max < -.2 & kendall_dcorr_max[,,,1] > 0 & kendall_dcorr_max[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > -.2 & prior_ELI < .2 & kendall_dcorr[,,,1] > 0 & kendall_dcorr[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > -.2 & prior_ELI_rsds < .2 & kendall_dcorr_rsds[,,,1] > 0 & kendall_dcorr_rsds[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > -.2 & prior_ELI_max < .2 & kendall_dcorr_max[,,,1] > 0 & kendall_dcorr_max[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > .2 & prior_ELI < 2 & kendall_dcorr[,,,1] > 0 & kendall_dcorr[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > .2 & prior_ELI_rsds < 2 & kendall_dcorr_rsds[,,,1] > 0 & kendall_dcorr_rsds[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > .2 & prior_ELI_max < 2 & kendall_dcorr_max[,,,1] > 0 & kendall_dcorr_max[,,,1] < .05)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > -2 & prior_ELI < -.2 & kendall_dcorr[,,,1] > .05 & kendall_dcorr[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > -2 & prior_ELI_rsds < -.2 & kendall_dcorr_rsds[,,,1] > .05 & kendall_dcorr_rsds[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > -2 & prior_ELI_max < -.2 & kendall_dcorr_max[,,,1] > .05 & kendall_dcorr_max[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > -.2 & prior_ELI < .2 & kendall_dcorr[,,,1] > .05 & kendall_dcorr[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > -.2 & prior_ELI_rsds < .2 & kendall_dcorr_rsds[,,,1] > .05 & kendall_dcorr_rsds[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > -.2 & prior_ELI_max < .2 & kendall_dcorr_max[,,,1] > .05 & kendall_dcorr_max[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > .2 & prior_ELI < 2 & kendall_dcorr[,,,1] > .05 & kendall_dcorr[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > .2 & prior_ELI_rsds < 2 & kendall_dcorr_rsds[,,,1] > .05 & kendall_dcorr_rsds[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > .2 & prior_ELI_max < 2 & kendall_dcorr_max[,,,1] > .05 & kendall_dcorr_max[,,,1] < 0.1)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > -2 & prior_ELI < -.2 & kendall_dcorr[,,,1] > 0.1 & kendall_dcorr[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > -2 & prior_ELI_rsds < -.2 & kendall_dcorr_rsds[,,,1] > 0.1 & kendall_dcorr_rsds[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > -2 & prior_ELI_max < -.2 & kendall_dcorr_max[,,,1] > 0.1 & kendall_dcorr_max[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > -.2 & prior_ELI < .2 & kendall_dcorr[,,,1] > 0.1 & kendall_dcorr[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > -.2 & prior_ELI_rsds < .2 & kendall_dcorr_rsds[,,,1] > 0.1 & kendall_dcorr_rsds[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > -.2 & prior_ELI_max < .2 & kendall_dcorr_max[,,,1] > 0.1 & kendall_dcorr_max[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI > .2 & prior_ELI < 2 & kendall_dcorr[,,,1] > 0.1 & kendall_dcorr[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_rsds > .2 & prior_ELI_rsds < 2 & kendall_dcorr_rsds[,,,1] > 0.1 & kendall_dcorr_rsds[,,,1] < 0.15)])/(12*total_land_area)),
                                                 100*(sum(area.array_rep[which(prior_ELI_max > .2 & prior_ELI_max < 2 & kendall_dcorr_max[,,,1] > 0.1 & kendall_dcorr_max[,,,1] < 0.15)])/(12*total_land_area)))



save(ELI_Tmaxbar_min_bars_10yr_1980_2100.df, col_reg_lim, 
     file = paste0(path_RData, "Fig6a.RData"))

