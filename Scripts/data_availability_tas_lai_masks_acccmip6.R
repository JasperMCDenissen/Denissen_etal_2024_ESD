# 2022-01-21
# by: Jasper Denissen
# Calculating the warm land area from the ensemble of CMIP6 projections
pdf(NULL)

##############################################################################################################################
##############################################################################################################################
######################## !!! Don't forget to reset the working directory to a directory of your choosing!!! ##################
##############################################################################################################################
###############################################################################################################
setwd('/Net/Groups/BGI/work_3/HydroBioClim/archive/Denissen_etal_2024_ESD/')

# RData
load("RData/hottest_3_moy_acccmip6.RData")
load("RData/202401_dcorr_cmip6_10yr_h3m_combined_mask_acccmip6_rsds.RData")
# # testdir
# load("testdir/hottest_3_moy_acccmip6.RData")
# load("testdir/202401_dcorr_cmip6_10yr_h3m_combined_mask_acccmip6_rsds.RData")

# Get the proper packages
source('Scripts/to_be_loaded_packages.R')
# plot discrete color bar
source('Scripts/plot_discrete_cbar.R')

lon <- seq(-179,179,2)
lat <- seq(-89,89,2)


hottest_3_moy <- aperm(as.array(hottest_3_moy), c(2,1,3,4))[,90:1,,]
# so hottest_3_moy has a 1 if that months belongs to the hottest three months over 120 years.
hottest_3_moy_tas_mask <- 
  hottest_3_moy_lai_mask <- 
  hottest_3_moy_combined_mask <- 
  array(NaN,c(180,90,12,1440))
i_rep <- rep(1:12,120)
for(source_id in 1:length(cmip6_data.df$source_id)){
  for(i in 1:1440){
    hottest_3_moy_tas_mask[,,source_id,i] <- hottest_3_moy[,,source_id,i_rep[i]]*tas_mask.list[[source_id]][,,i]
    hottest_3_moy_lai_mask[,,source_id,i] <- hottest_3_moy[,,source_id,i_rep[i]]*lai_mask.list[[source_id]][,,i]
    hottest_3_moy_combined_mask[,,source_id,i] <- hottest_3_moy[,,source_id,i_rep[i]]*combined_mask.list[[source_id]][,,i]
  }
}

sum_hottest_3_moy_tas_mask <- 
  sum_hottest_3_moy_lai_mask <- 
  sum_hottest_3_moy_combined_mask <- 
  array(NaN,c(180,90,12))
for(x in 1:180){
  for(y in 1:90){
    for(source_id in 1:length(cmip6_data.df$source_id)){
      sum_hottest_3_moy_tas_mask[x,y,source_id] <- sum(hottest_3_moy_tas_mask[x,y,source_id,], na.rm = T)
      sum_hottest_3_moy_lai_mask[x,y,source_id] <- sum(hottest_3_moy_lai_mask[x,y,source_id,], na.rm = T)
      sum_hottest_3_moy_combined_mask[x,y,source_id] <- sum(hottest_3_moy_combined_mask[x,y,source_id,], na.rm = T)
    }
  }
}
sum_hottest_3_moy_tas_mask[which(sum_hottest_3_moy_tas_mask == 0)] <- NaN
sum_hottest_3_moy_lai_mask[which(sum_hottest_3_moy_lai_mask == 0)] <- NaN
sum_hottest_3_moy_combined_mask[which(sum_hottest_3_moy_combined_mask == 0)] <- NaN

copy_sum_hottest_3_moy_tas_mask <- sum_hottest_3_moy_tas_mask
copy_sum_hottest_3_moy_lai_mask <- sum_hottest_3_moy_lai_mask
copy_sum_hottest_3_moy_combined_mask <- sum_hottest_3_moy_combined_mask

sum_hottest_3_moy_tas_mask[1:90,,] <- copy_sum_hottest_3_moy_tas_mask[91:180,,]; sum_hottest_3_moy_tas_mask[91:180,,] <- copy_sum_hottest_3_moy_tas_mask[1:90,,]
sum_hottest_3_moy_lai_mask[1:90,,] <- copy_sum_hottest_3_moy_lai_mask[91:180,,]; sum_hottest_3_moy_lai_mask[91:180,,] <- copy_sum_hottest_3_moy_lai_mask[1:90,,]
sum_hottest_3_moy_combined_mask[1:90,,] <- copy_sum_hottest_3_moy_combined_mask[91:180,,]; sum_hottest_3_moy_combined_mask[91:180,,] <- copy_sum_hottest_3_moy_combined_mask[1:90,,]

list_of_masks <- list(sum_hottest_3_moy_tas_mask, sum_hottest_3_moy_lai_mask, sum_hottest_3_moy_combined_mask)
masks.df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)),
                     c("data_count","lon","lat","mask_type","source_id"))
mask_types <- c("tas","lai","combined")
for(i_mask in 1:3){
  for(source_id in 1:length(cmip6_data.df$source_id)){
    for(x in 1:180){
      for(y in 1:90){
        if(!is.na(list_of_masks[[i_mask]][x,y,source_id])){
          masks.df <- rbind(masks.df,
                            data.frame("data_count" = list_of_masks[[i_mask]][x,y,source_id],
                                       "lon" = lon[x]+180,
                                       "lat" = lat[y],
                                       "mask_type" = mask_types[i_mask],
                                       "source_id" = cmip6_data.df$source_id[source_id]))
        }
      }
    }
    print(paste0("source_id ",cmip6_data.df$source_id[source_id]," is done..."))
  }
}
masks.df$lon <- masks.df$lon - 180
masks.df$mask_type_f = factor(masks.df$mask_type, levels=c('tas','lai','combined'))

cols_data_count <- colorRamps::matlab.like(n = 9)
myvalues_data_count <- c(1,seq(40,360,40))
dbar_data_count <- plot_discrete_cbar(myvalues_data_count,
                                      colors = cols_data_count,
                                      spacing = "constant",
                                      font_size = 6,
                                      spacing_scaling = 2,
                                      width = .2,
                                      triangle_size = .175)
masks.df$cuts_data_count <- cut(masks.df$data_count, myvalues_data_count, include.lowest = T)

# This is a data set from the maptools package
data(wrld_simpl)

# Create a data.frame object for ggplot. ggplot requires a data frame.
mymap <- fortify(wrld_simpl)

# New facet label names for mask_type variable
mask_type.labs <- c("Ta > 10°C", "LAI > 0.2", "T > 10°C & LAI > 0.2")

names(mask_type.labs) <- c("tas", "lai", "combined")

a <- ggplot(masks.df, aes(x=lon,y=lat,fill=cuts_data_count)) +
  geom_tile() +
  geom_map(data = mymap, map = mymap, aes(x = long, y = lat, map_id = id),
           inherit.aes = F, fill = NA, color = "black", size = .1) +
  facet_grid(source_id ~ mask_type_f,
             labeller = labeller(mask_type_f = mask_type.labs)) + 
  scale_fill_manual(values = cols_data_count) + 
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
        strip.text = element_text(size=14, face="bold", margin = margin(0, 0, .5, 0, "cm")),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA, size = 1),
        panel.background = element_blank(),
        axis.text = element_text(size=18),
        axis.title = element_blank(),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("Values retained after masking")
a

# reduce top and bottom margins
empty <- ggplot() + theme_void()
dbar_data_count <- dbar_data_count + theme(plot.margin = unit(c(-35, 10, -30, 10), "pt"))
dbar_data_count_smaller <- grid.arrange(empty, dbar_data_count, empty , ncol=3, widths = c(1,9,1))
plot <- grid.arrange(a,dbar_data_count_smaller, nrow = 2, heights = c(.975,.025))

# ggsave("Figures/SFig2.png", plot = plot, width = 10*1.3, height = 20*1.3, units = "in")
ggsave("testdir/SFig2.png", plot = plot, width = 10*1.3, height = 20*1.3, units = "in")
