# By: Jasper Denissen
# 2024/01/06
# script to 
# 1) regrid data downloaded with acccmip6 to 2.0x2.0 degree grid resolution using bilinear interpolation and 
# 2) only keep the values between 1980-2100 and
# 3) compute weighted average of top meter soil moisture from total water content per soil layer and
# 4) write as .nc
pdf(NULL)

#####################################################################################################################
#####################################################################################################################
######################## !!! Don't forget to reset the working directory to a directory of your choosing!!! #########
#####################################################################################################################
#####################################################################################################################
setwd('/Net/Groups/BGI/work_3/HydroBioClim/archive/Denissen_etal_2024_ESD/')

# packages
source('Scripts/to_be_loaded_packages.R')

path_cmip6 <- "Data/cmip6_202312_climex_ESD_acccmip6/"
path_cmip6_2.0x2.0 <- "Data/cmip6_202312_climex_ESD_acccmip6_2.0x2.0/"
# # testdir
# path_cmip6_2.0x2.0 <- "testdir/cmip6_202312_climex_ESD_acccmip6_2.0x2.0/"
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

# mrsol SHOULD BE PART OF THE VARIABLE LIST, AS IT IS USED TO MAKE A LAND MASK!
list_vars <- list_all[grep(paste0("*",unique(source_id[1]),"*"),list_all)]
vars <- c()
for(i in 1:length(list_vars)){
  vars <- c(vars, strsplit(list_vars[i],"_")[[1]][1])
}
uvars <- unique(vars)
uvars[which(uvars == 'tas')] <- "tas_" # otherwise "tas" returns both "tas" and "tasmax" in the string search.
uvars <- uvars[c(3,1,2,4:9)]
uvars_mrso <- uvars
uvars_mrso[which(uvars == 'mrsol')] <- 'mrso'

# Rotate to center on Europe & Africa
target_rast <- rotate(rast('/Data/hfls_historical_Amon_ACCESS-ESM1-5_r10i1p1f1_1980_2015_2.0x2.0.nc'), left = T)

for(source in cmip6_data.df$source_id){
  list_source <- list_all[grep(paste0("*",source,"*"),list_all)]
  if(cmip6_data.df$SM[which(cmip6_data.df$source_id == source)] == 'mrsol'){
    for(var in uvars){
      list_source_var <- list_source[grep(paste(var,"*"),list_source)]
      stack_var <- rotate(rast(paste0(path_cmip6, list_source_var)), left = T)
      time_var <- time(stack_var)
      if(var != 'mrsol'){
        stack_var_t <- stack_var[[which(time_var == '1981-01-16'):which(time_var == '2100-12-16')]]
      }else{
        ndepths <- length(which(time_var == '1981-01-16'))
        depths <- c()
        pre_depths_char <- names(stack_var)[1:ndepths]
        for(i in 1:ndepths){
          pre1_depths_char <- strsplit(pre_depths_char[i],"_")[[1]][2]
          pre2_depths_char <- strsplit(pre1_depths_char,"=")[[1]][2]
          depths <- c(depths, as.numeric(pre2_depths_char))
        }
        stack_var_t <- stack_var[[which(time_var == '1981-01-16')[1]:which(time_var == '2100-12-16')[ndepths]]]
        dweights <- depths[1]
        for(j in 2:which(depths > 1)[1]){
          if(j == which(depths > 1)[1]){
            dweights <- c(dweights,1 - sum(dweights))
          }else{
            dweights <- c(dweights,depths[j] - depths[j-1])
          }
        }
        wmeans <- stack_var_t[[seq(1,length(time_var),ndepths)]]*dweights[1]
        for(k in 2:length(dweights)){
          wmeans <- wmeans + stack_var_t[[seq(k,length(time_var),ndepths)]]*dweights[k]
          print(k)
        }
        stack_var_t <- wmeans
      }
      if(var == 'mrsol'){
        # use 1-m average soil moisture as a mask
        # first, two models (EC-Earth and GFDL) do not have NA's in the ocean, but zero. So replace 0 with NA
        copy <- stack_var_t
        copy[copy == 0] <- NA
        # "aggregate NAs"
        # make a raster object with 1 if NA and 0 if not NA
        x_NA <- is.na(copy)
        # resample this NA-mask to the target resolution and/or origin using nearest neighbor
        x_NA2 <- resample(x_NA, target_rast, method = "near")
      }else{
        stack_var_t <- mask(stack_var_t, x_NA, maskvalue = 1)
      }
      stack_var_t_rs <- resample(stack_var_t, target_rast, method = 'bilinear')
      # and use as mask
      stack_var_t_rs <- mask(stack_var_t_rs, x_NA2, maskvalue = 1)
      test <- writeCDF(rotate(stack_var_t_rs, left = F), paste0(path_cmip6_2.0x2.0,strsplit(var,"_")[[1]][1],"_",source,"_",cmip6_data.df$member_id[which(cmip6_data.df$source_id == source)],"_1980_2100_2.0x2.0.nc"), overwrite=T, varname = var)
      print(paste0(var, " for ", source, " is done..."))
      
      

    }
  }else{
    # routine for mrso
    # might skip altogether, as with mrsol we already get 12 models. 
    print(paste0("skipping ", source, " as it doesn't have mrsol..."))
  }
}

