# 2022-01-21
# by: Jasper Denissen
# Calculating the warm land area from the ensemble of CMIP6 projections
pdf(NULL)

##############################################################################################################################
##############################################################################################################################
######################## !!! Don't forget to reset the working directory to a directory of your choosing!!! ##################
##############################################################################################################################
##############################################################################################################################
setwd('/Net/Groups/BGI/work_3/HydroBioClim/archive/Denissen_etal_2024_ESD/')

# from /RData
load("RData/202401_dcorr_cmip6_10yr_h3m_combined_mask_acccmip6_rsds.RData")
path_RData <- "RData/"
# # from /testdir
# load("testdir/202401_dcorr_cmip6_10yr_h3m_combined_mask_acccmip6_rsds.RData")
# path_RData <- "testdir/"

# Get the proper packages
source('Scripts/to_be_loaded_packages.R')
lon <- seq(-179,179,2)
lat <- seq(-89,89,2)
count <- 1
for(i in 1:12){
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- dcorr.list[[i]][91:180,,]; test[91:180,,] <- dcorr.list[[i]][1:90,,]; dcorr.list[[i]] <- test; dcorr.list[[i]] <- -1*dcorr.list[[i]]
  test <- array(NaN,c(180,90,12)); test[1:90,,] <- dcorr_netrad.list[[i]][91:180,,]; test[91:180,,] <- dcorr_netrad.list[[i]][1:90,,]; dcorr_netrad.list[[i]] <- test; dcorr_netrad.list[[i]] <- -1*dcorr_netrad.list[[i]]
  
  count <- count + 1
}


# select the right models and put in array
dcorr_all.array <- 
  dcorr_netrad.array <-
  array(NaN,c(180,90,12*12))
count_all <- 1
for(i in 1:12){
  dcorr_all.array[,,count_all:(count_all+11)] <- dcorr.list[[i]]
  dcorr_netrad.array[,,count_all:(count_all+11)] <- dcorr_netrad.list[[i]]
  count_all <- count_all + 12
  print(paste(i, " is done...",sep=''))
}

models_with_full_timeseries <- 
  models_with_full_timeseries_mask<- array(0,c(180,90))
count_i <- 1
for(i in seq(1,(12*12),12)){
  for(x in 1:180){
    for(y in 1:90){
      if(sum(!is.na(dcorr_all.array[x,y,(i:(i+11))])) == 12){ # only calculate if the model has full time series
        models_with_full_timeseries[x,y] <- models_with_full_timeseries[x,y] + 1
        if(models_with_full_timeseries[x,y] > 8){
          models_with_full_timeseries_mask[x,y] <- 1
        }
      }
    }
  }
  print(i)
  count_i <- count_i + 1
}

models_with_full_timeseries_netrad <- 
  models_with_full_timeseries_netrad_mask<- array(0,c(180,90))
count_i <- 1
for(i in seq(1,(12*12),12)){
  for(x in 1:180){
    for(y in 1:90){
      if(sum(!is.na(dcorr_netrad.array[x,y,(i:(i+11))])) == 12){ # only calculate if the model has full time series
        models_with_full_timeseries_netrad[x,y] <- models_with_full_timeseries_netrad[x,y] + 1
        if(models_with_full_timeseries_netrad[x,y] > 8){
          models_with_full_timeseries_netrad_mask[x,y] <- 1
        }
      }
    }
  }
  print(i)
  count_i <- count_i + 1
}

# maybe now per surface area instead of % of grid cells?
r <- raster()  # by default 1 by 1 degree
res(r) <- 2 # so change the resolution
a <- raster::area(r) # calculate the area of a 2x2 degree grid since N - S, as area varies only by latitude, not longitude
area <- a[,1]
area.array <- array(NaN,c(180,90))
for(x in 1:180){
  area.array[x,] <- area
}
total_land_area <- sum(models_with_full_timeseries_mask[,which(lat > -60 & lat < 70)]*area.array[,which(lat > -60 & lat < 70)], na.rm = T)
total_land_area_netrad <- sum(models_with_full_timeseries_netrad_mask[,which(lat > -60 & lat < 70)]*area.array[,which(lat > -60 & lat < 70)], na.rm = T)

models_with_full_timeseries[,which(lat < -60 | lat > 70)] <- 0
models_with_full_timeseries_netrad[,which(lat < -60 | lat > 70)] <- 0

save(models_with_full_timeseries, total_land_area,
     models_with_full_timeseries_netrad, total_land_area_netrad, 
     file = paste0(path_RData, "total_land_area_CMIP6_climex_acccmip6.RData"))
