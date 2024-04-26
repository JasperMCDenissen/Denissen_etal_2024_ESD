# By: Jasper Denissen
# 2024-03-12
# Script to Fig6
pdf(NULL)

############################################################################################################################
############################################################################################################################
######################## !!! Don't forget to reset the working directory to a directory of your choosing!!! ################
############################################################################################################################
###############################################################################################################
setwd('/Net/Groups/BGI/work_3/HydroBioClim/archive/Denissen_etal_2024_ESD/')

# Get the proper packages
source('Scripts/to_be_loaded_packages.R')

# RData
load("RData/Fig6a.RData")
load("RData/Fig6a_source_id.RData")
# # testdir
# load("testdir/Fig6a.RData")
# load("testdir/Fig6a_source_id.RData")

ELI_Tmaxbar_min_bars_10yr_1980_2100_copymax.df <- ELI_Tmaxbar_min_bars_10yr_1980_2100.df[which(ELI_Tmaxbar_min_bars_10yr_1980_2100.df$dcorr_with == 'max'),]
ELI_Tmaxbar_min_bars_10yr_1980_2100_copymax.df <- ELI_Tmaxbar_min_bars_10yr_1980_2100_copymax.df[c(1,2,3,5)]
ELI_Tmaxbar_min_bars_10yr_1980_2100_copymax.df$source_id <- "Multi-model mean"

ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df <- rbind(ELI_max_Tmaxbar_min_bars_10yr_1980_2100_source_id.df, ELI_Tmaxbar_min_bars_10yr_1980_2100_copymax.df)
ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$source_id_f <- factor(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$source_id, levels = c("Multi-model mean", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3-CC", "GFDL-ESM4", "HadGEM3-GC31-LL", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", "UKESM1-0-LL"))

ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$Tmaxbar_min_for_area <- ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$Tmaxbar_min
ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$Tmaxbar_min_for_area[which(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$Tmaxbar_min_for_area < 0)] <- .005

ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$area_round <- NaN
for(i in 1:length(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$area_round)){
  if(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$area[i] < 1){
    ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$area_round[i] <- round(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$area[i], 1)
  }else{
    ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$area_round[i] <- round(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$area[i], 0)
  }
}

cols_source_id <- c("Multi-model mean" = "black",
                    "ACCESS-ESM1-5" = colorRampPalette(brewer.pal(9, "Set1"))(12)[1],
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


a4 <- ggplot(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df[which(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$reg == 'Global'),], aes(x = kendall_ELI, y = Tmaxbar_min, col = source_id_f, fill = source_id_f, label = round(area,1))) +
  geom_bar(stat='identity',col = 'black', width = .75, position = "dodge") +
  geom_hline(yintercept=0) +
  geom_text(size = 3.5, position = position_dodge(width=.75), vjust = -1,
            inherit.aes = F, data = ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df[which(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$reg == 'Global'),], aes(x = kendall_ELI, y = Tmaxbar_min_for_area, col = source_id_f, label = area_round)) +
  scale_x_continuous("ELI trend (-/10yr)",
                     breaks = seq(1,4,1),
                     labels = c("-0.05 < ELI < 0", "0 < ELI < 0.05", "0.05 < ELI < 0.1", "0.1 < ELI < 0.15")) + 
  scale_y_continuous(expression("Temperature excess trend (K/10yr)"),
                     limits = c(min(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$Tmaxbar_min[which(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$reg == 'Global')], na.rm = T),
                                max(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$Tmaxbar_min[which(ELI_Tmaxbar_min_bars_10yr_1980_2100_source_id_combined.df$reg == 'Global')], na.rm = T) + .01),
                     expand=c(0,0)) +
  scale_colour_manual("",
                      values = cols_source_id) +
  scale_fill_manual("",
                    values = cols_source_id) +
  theme(legend.position = c(.125,.7),
        # theme(legend.position = "right",
        legend.text = element_text(size=12),
        legend.title = element_text(size=20),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = NA, colour = NA),
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
        # axis.ticks.x = element_line(hjust = 1),
        plot.title = element_text(size=24),
        plot.tag.position = c(.55,0.03)
  ) + ggtitle("CMIP6 1980 - 2100 (per decade)")
a4

# ggsave("Figures/Fig6.png", plot = a4, width = 3.33*3.5, height = 3.33*2, units = "in")
ggsave("testdir/Fig6.png", plot = a4, width = 3.33*3.5, height = 3.33*2, units = "in")

