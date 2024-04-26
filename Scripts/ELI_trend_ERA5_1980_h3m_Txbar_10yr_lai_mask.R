# # By: Jasper Denissen
# 2021-02-04
# Script to process ERA5(-Land) data
# 1980-2019, 2.0x2.0 grid cell resolution WITH THE RASTER PACKAGE!
pdf(NULL)

#####################################################################################################################
#####################################################################################################################
######################## !!! Don't forget to reset the working directory to a directory of your choosing!!! #########
#####################################################################################################################
#####################################################################################################################
setwd('/Net/Groups/BGI/work_3/HydroBioClim/archive/Denissen_etal_2024_ESD/')

# packages
source('Scripts/to_be_loaded_packages.R')

##########################################################################################
##########################################################################################
########################################## 2d00 ##########################################
##########################################################################################
##########################################################################################


# RData
load("RData/202207_dcorr_ERA5_10yr_h3m_daily_combined_mask_terra_tas_ssrd.RData")
load('RData/total_land_area_CMIP6_climex_acccmip6.RData')
# Load mask based on CMIP6 model data availability for fair comparison...
load("RData/mask_full_acccmip6.RData")
path_RData <- "RData/"
# # testdir
# load("testdir/202207_dcorr_ERA5_10yr_h3m_daily_combined_mask_terra_tas_ssrd.RData")
# load('testdir/total_land_area_CMIP6_climex_acccmip6.RData')
# # Load mask based on CMIP6 model data availability for fair comparison...
# load("testdir/mask_full_acccmip6.RData")
# path_RData <- "testdir/"

# source functions
source('/Net/Groups/BGI/people/jdenis/scripts/functions/plot_discrete_cbar.R')

# maybe now per surface area instead of % of grid cells?
r <- raster()  # by default 1 by 1 degree
res(r) <- 2 # so change the resolution
a <- raster::area(r) # calculate the area of a 2x2 degree grid from N - S, as area varies only by latitude, not longitude
area <- a[,1]
area.array <- array(NaN,c(180,90))
for(x in 1:180){
  area.array[x,] <- area
}

lon <- seq(-179,179,2)
lat <- seq(-89,89,2)

corr_rgy_veg_max.array <- corr_rgy_veg_ssrd.array
dcorr_max.array <- array(NaN,c(180,90,12))
# now the strongest average energy proxy
t2m_or_ssrd.array <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    for(i in 1:4){}
    if(sum(!is.na(corr_rgy_veg.array[x,y,i])) > 0){
      if(mean(corr_rgy_veg.array[x,y,i], na.rm = T) > mean(corr_rgy_veg_ssrd.array[x,y,i], na.rm = T)){
        corr_rgy_veg_max.array[x,y,i] <- corr_rgy_veg.array[x,y,i]
        t2m_or_ssrd.array[x,y] <- 1
      }else{
        t2m_or_ssrd.array[x,y] <- 2
        corr_rgy_veg_max.array[x,y,i] <- corr_rgy_veg_ssrd.array[x,y,i]
      }
    }
  }
}
dcorr_max.array <- corr_wtr_veg.array - corr_rgy_veg_max.array
dcorr_max_1m.array <- corr_wtr_veg_swvl123.array - corr_rgy_veg_max.array

# prior_ELI <- 
prior_ELI_max <- prior_ELI_max_1m <- array(NaN,c(180,90))
for(x in 1:180){
  for(y in 1:90){
    prior_ELI_max[x,y] <- mean(dcorr_max.array[x,y,4:7])
    prior_ELI_max_1m[x,y] <- mean(dcorr_max_1m.array[x,y,4:7])
  }
}

prior_ELI_max_vec <- 
  prior_ELI_max_1m_vec <- 
  lon_vec <- lat_vec <-
  c()
for(x in 1:180){
  for(y in 1:90){
    lon_vec <- c(lon_vec, lon[x])
    lat_vec <- c(lat_vec, lat[y])
    prior_ELI_max_vec <- c(prior_ELI_max_vec, prior_ELI_max[x,y])
    prior_ELI_max_1m_vec <- c(prior_ELI_max_1m_vec, prior_ELI_max_1m[x,y])
  }
}

prior_ELI.df <- data.frame("ELI_max" = prior_ELI_max_vec,
                           "ELI_max_1m" = prior_ELI_max_1m_vec,
                           # "ELI" = prior_ELI_vec,
                           "lon" = lon_vec,
                           "lat" = lat_vec)
dcorrcol <- rev(brewer.pal(11,"BrBG"))
myvalues_prior_ELI <- c(-Inf,-.8,-.6,-.4,-.2,-.0001,.0001,.2,.4,.6,.8,Inf)
prior_ELI.df$cuts_prior_ELI_max <- cut(prior_ELI.df$ELI_max, myvalues_prior_ELI, include.lowest = T)
prior_ELI.df$cuts_prior_ELI_max_1m <- cut(prior_ELI.df$ELI_max_1m, myvalues_prior_ELI, include.lowest = T)

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

center_year <- seq(1955,2015,10)
global_tseries_ELI.df <- setNames(data.frame(matrix(ncol = 6, nrow = 0)),
                                  c("ELI_max","ELI_max_1m","T2maxbar_min","center_year","year","reg"))
for(i in 4:7){
  global_tseries_ELI.df <- rbind(global_tseries_ELI.df,
                                 data.frame("ELI_max" = weighted.mean(x = dcorr_max.array[,,i] * mask_full, w = area.array, na.rm = T),
                                            "ELI_max_1m" = weighted.mean(x = dcorr_max_1m.array[,,i] * mask_full, w = area.array, na.rm = T),
                                            "T2maxbar_min" = weighted.mean(x = (T2maxbar.array[,,i] - av_t2m.array[,,i]) * mask_full, w = area.array, na.rm = T),
                                            "center_year" = center_year[i],
                                            "year" = center_year[i]-5,
                                            "reg" = "Global"
                                 ))
}
global_tseries_lim.df <- global_tseries_ELI.df

prior_ELI_mask <- array(NaN,c(180,90,3))
for(x in 1:180){
  for(y in 1:90){
    if(!is.na(prior_ELI_max_1m[x,y])){
      if(prior_ELI_max_1m[x,y] < -.2){
        prior_ELI_mask[x,y,1] <- 1
      }else if(prior_ELI_max_1m[x,y] < 0.2){
        prior_ELI_mask[x,y,2] <- 1
      }else{
        prior_ELI_mask[x,y,3] <- 1
      }
    }
  }
}

reg_lim <- c("rgy_lim","trans","wtr_lim")
for(j in 1:3){ # loop over all reg_lim
  for(t in 4:7){ # loop over all decades
    global_tseries_lim.df <- rbind(global_tseries_lim.df,
                                   data.frame("ELI_max" = weighted.mean(x = (dcorr_max.array[,,t] * prior_ELI_mask[,,j]),
                                                                        w = area.array,
                                                                        na.rm = T),
                                              "ELI_max_1m" = weighted.mean(x = (dcorr_max_1m.array[,,t] * prior_ELI_mask[,,j]),
                                                                           w = area.array,
                                                                           na.rm = T),
                                              "T2maxbar_min" = weighted.mean(x = (T2maxbar.array[,,t] -
                                                                                    av_t2m.array[,,t]) * prior_ELI_mask[,,j],
                                                                             w = area.array,
                                                                             na.rm = T),
                                              "center_year" = center_year[t],
                                              "year" = center_year[t]-5,
                                              "reg" = reg_lim[j]))
    
  }
}

for(i in 1:length(global_tseries_lim.df$ELI_max)){
  global_tseries_lim.df$ELI_max_rel_1980[i] <- (global_tseries_lim.df$ELI_max[i] - global_tseries_lim.df$ELI_max[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i])][1])
  global_tseries_lim.df$ELI_max_1m_rel_1980[i] <- (global_tseries_lim.df$ELI_max_1m[i] - global_tseries_lim.df$ELI_max_1m[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i])][1])
  global_tseries_lim.df$T2maxbar_min_rel_1980[i] <- (global_tseries_lim.df$T2maxbar_min[i] - global_tseries_lim.df$T2maxbar_min[which(global_tseries_lim.df$reg == global_tseries_lim.df$reg[i])][1])
}

col_reg_lim <- c("Global" = "black",
                 "rgy_lim" = brewer.pal(9,"BrBG")[8],
                 "trans" = 'slategrey',
                 "wtr_lim" = brewer.pal(9,"BrBG")[2])

# Define hotspots
hotspot_regs.df <- data.frame("x" = c(lon[53]-1, lon[69]+1, lon[53]-1, lon[53]-1,
                                      lon[28]-1, lon[45]+1, lon[28]-1, lon[28]-1,
                                      lon[91]-1, lon[114]+1, lon[91]-1, lon[91]-1,
                                      lon[118]-1, lon[143]+1, lon[118]-1, lon[118]-1),
                              "y" = c(lat[34]-1, lat[34]-1, lat[34]-1, lat[46]+1,
                                      lat[67]-1, lat[67]-1, lat[67]-1, lat[73]+1,
                                      lat[64]-1, lat[64]-1, lat[64]-1, lat[72]+1,
                                      lat[72]-1, lat[72]-1, lat[72]-1, lat[77]+1),
                              "xend" = c(lon[53]-1, lon[69]+1, lon[69]+1, lon[69]+1,
                                         lon[28]-1, lon[45]+1, lon[45]+1, lon[45]+1,
                                         lon[91]-1, lon[114]+1, lon[114]+1, lon[114]+1,
                                         lon[118]-1, lon[143]+1, lon[143]+1, lon[143]+1),
                              "yend" = c(lat[46]+1, lat[46]+1, lat[34]-1, lat[46]+1,
                                         lat[73]+1, lat[73]+1, lat[67]-1, lat[73]+1,
                                         lat[72]+1, lat[72]+1, lat[64]-1, lat[72]+1,
                                         lat[77]+1, lat[77]+1, lat[72]-1, lat[77]+1))

global_tseries_ERA5.df <- global_tseries_ELI.df
colnames(global_tseries_ERA5.df) <- c("ELI","ELI_max","Tmaxbar_min", "center_year", "year", "reg")
global_tseries_ERA5.df['source_id'] = "ERA5"
global_tseries_ERA5.df <- global_tseries_ERA5.df[c("ELI_max","Tmaxbar_min", "center_year", "year", "source_id", "reg")]
# Saved into the same dataframe structure as the individual CMIP6 models, within which dcorr_max_1m.array from ERA5 will be usd as ELI_max
reg <- c("SAM", "NAM", "CEU", "NAS")
count_j <- 1
for(j in seq(1,16,4)){ # loop over all regions
  lonmin <- min(hotspot_regs.df$x[j:(j+3)])
  lonmax <- max(hotspot_regs.df$x[j:(j+3)])
  latmin <- min(hotspot_regs.df$y[j:(j+3)])
  latmax <- max(hotspot_regs.df$y[j:(j+3)])
  reg_mask <- array(NaN,c(180,90)); reg_mask[which(lon >= lonmin & lon <= lonmax), which(lat >= latmin & lat <= latmax)] <- 1
  for(t in 4:7){ # loop over all decades
    global_tseries_ERA5.df <- rbind(global_tseries_ERA5.df,
                               data.frame("ELI_max" = weighted.mean(x = dcorr_max_1m.array[,,t] * reg_mask * mask_full, 
                                                                    w = area.array, 
                                                                    na.rm = T),
                                          "Tmaxbar_min" = weighted.mean(x = (T2maxbar.array[,,t] - av_t2m.array[,,t]) * reg_mask * mask_full, 
                                                                        w = area.array, 
                                                                        na.rm = T),
                                          "center_year" = center_year[t],
                                          "year" = center_year[t]-5,
                                          "source_id" = "ERA5",
                                          "reg" = reg[count_j]))
  }
  count_j <- count_j + 1
}

for(i in 1:length(global_tseries_ERA5.df$ELI_max)){
  global_tseries_ERA5.df$ELI_max_rel_1980[i] <- (global_tseries_ERA5.df$ELI_max[i] - global_tseries_ERA5.df$ELI_max[which(global_tseries_ERA5.df$reg == global_tseries_ERA5.df$reg[i])][1])
  global_tseries_ERA5.df$Tmaxbar_min_rel_1980[i] <- (global_tseries_ERA5.df$Tmaxbar_min[i] - global_tseries_ERA5.df$Tmaxbar_min[which(global_tseries_ERA5.df$reg == global_tseries_ERA5.df$reg[i])][1])
}

save(global_tseries_ERA5.df,
     file = paste0(path_RData, "global_tseries_ERA5.RData"))
