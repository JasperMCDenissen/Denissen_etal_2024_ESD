# Denissen_etal_2024_ESD

2024-04-22
by: Jasper Denissen

This is a text message describing the use of all the files included in this directory for the paper "Intensified future heat extremes linked with increasing ecosystem water limitation", Denissen et al., 2024. 

This main directory should have 4 subfolders:
- Data: This can be downloaded from https://zenodo.org/doi/10.5281/zenodo.11072826
- RData: This can be downloaded from https://zenodo.org/doi/10.5281/zenodo.11072826
- Figures: A folder with the original paper figures.
- testdir: An empty folder meant for output of all scripts. Should be made by the user. 

Any .R script can be executed by opening R from the terminal by simply typing 'R' in the command prompt. In R, execute 'source('name_script.R')'. Take care to check where any output (.nc files, .RData files, figures) is stored. To prevent overwriting existing files, change any output path to '/testdir/'.

[NOTE] It is not necessary to download and aggregate the CMIP6 runs and pre-process the data, the paper figures can be made directly in [3].
- [0] states the R environment originally used.
- [1] download and aggregate the CMIP6 and ERA5(-Land) data
- [2] Data processing for plot figures
- [3] Making the paper figures

[0] R environment
The following version of R has been used to compile the scripts
R version 4.3.1 (2023-06-16) -- "Beagle Scouts"
Copyright (C) 2023 The R Foundation for Statistical Computing
Platform: x86_64-redhat-linux-gnu (64-bit)

Install all necessary packages from within the R terminal: source('Scripts/to_be_installed_packages.R')

Make sure to re-set the working directory to a directory of your choosing, preferably where you downloaded the data & scripts, at the start of every script!

[1] download and aggregate the CMIP6 runs
The python script to download CMIP6 runs accompanied by a tutorial is available from https://github.com/TaufiqHassan/acccmip6
For details on variable_id, institution_id, source_id and member_id, look up Table 1 in the Materials and Methods section of the manuscript. 
If any original CMIP6 data is downloaded, it should be added in /Data/cmip6_202312_climex_ESD_acccmip6! 
Aggregating the CMIP6 data is not necessary, as regridded data has been stored in Data/cmip6_202312_climex_ESD_acccmip6_2.0x2.0
Further, daily temperature data and monthly averaged variables from ERA5 and ERA5Land have been saved in Data/ERA5_t2m_2d00_daily and Data/era5land.

Data/
- ERA5_t2m_2d00_daily                                                   # folder containing daily temperature (t2m) from ERA5 at 2 degree grid resolution
- cmip6_202312_climex_ESD_acccmip6_2.0x2.0                              # folder containing regridded data from 12 CMIP6 models at 2 degree grid resolution
- era5land                                                              # folder containing regridded variables from era5 and era5land at 2 degree grid resolution
- hfls_historical_Amon_ACCESS-ESM1-5_r10i1p1f1_1980_2015_2.0x2.0.nc     # netCDF file that containst he desired spatial grid resolution

[2] Data pre-processing
It is not necessary to compute the pre-processed data, as all the .RData files and functions necessary to make the paper figures can be found here:

RData/
- 202207_dcorr_ERA5_10yr_h3m_daily_combined_mask_terra_tas_ssrd.RData   # .RData file containing spatially (2.0x2.0 degree) and temporally (decennial) aggregated ERA5(-Land) variables and the ELI from 198
- 202401_dcorr_cmip6_10yr_h3m_combined_mask_acccmip6_rsds.RData         # .RData file containing spatially (2.0x2.0 degree) and temporally (decennial) aggregated CMIP6 variables and the ELI from 1981-2100
- 202401_max_tasmax_moy_terra.RData                                     # .RData file containing a list with the maximum daily temperature per month of year per decade (more accurate description in Materi
- Fig6a.RData                                                           # .RData file containing input R data.frames with ELI per CMIP6 model for Fig6
- Fig6a_source_id.RData                                                 # .RData file containing input R data.frames with temperature excess per CMIP6 model for Fig6
- global_tseries_ERA5.RData                                             # .RData file containing input R data.frames with ERA5-Land ELI and temperature excess for Fig4
- hottest_3_moy_acccmip6.RData                                          # .RData file containing a mask (list) that gives 1 if that month of year per decade is one of the three hottest months of year for
- mask_full_acccmip6.RData                                              # .RData file containing a mask(array) that gives 1 when all 12 CMIP6 models have full time series in that respective grid cell
- total_land_area_CMIP6_climex_acccmip6.RData                           # .RData file containing an array with the number of CMIP6 models that have full time series in that respective grid cell and the to

Scripts/
- calc_boxes.R                                                          # Function to average a target variable across bins along axes of two other variables (side panels in Figure 3)
- parallel_raster_functions.R                                           # Collection of functions to apply to raster objects in a parallel computing environment
- plot_discrete_cbar.R                                                  # Function that compiles a discrete color bar

If there is interest to pre-process the data, it can be done with the following scripts in this sequence:

Scripts/
- check_regrid_acccmip6_models.R                                        # Script to read and regrid CMIP6 model data to 2.0 x 2.0 degree grid resolution and compute weighted average of soil moisture over
- read_cmip6_2d00_1980_h3m_Txbar_lai_mask_terra_acccmip6.R              # Script to process all CMIP6 model data to temporally (decennial) aggregate variables and compute the ELI from 1981-2100. This scri
- total_land_area_acccmip6.R                                            # Script to output total_land_area_CMIP6_climex_acccmip6.RData
- read_ERA5_2d00_1950_h3m_Txbar_10yr_lai_mask_terra_ssrd.R              # Script to output RData/202207_dcorr_ERA5_10yr_h3m_daily_combined_mask_terra_tas_ssrd.RData
- ELI_trend_seascycle_climex_mean_sourceid_Txbar_h3m_rsds_acccmip6.R    # Script to output RData/202207_dcorr_ERA5_10yr_h3m_daily_combined_mask_terra_tas_ssrd.RData, RData/Fig6a.RData, RData/Fig6a_source_
- ELI_trend_ERA5_1980_h3m_Txbar_10yr_lai_mask.R                         # Script to output RData/global_tseries_ERA5.RData

[3] Making the paper figures
All the original paper figures can be found in Figures/. 

With the following scripts one can remake all the original paper and supplementary material figures (Fig 1-4 and SFig 1-20):

Scripts/
- ELI_trend_seascycle_climex_mean_sourceid_Txbar_h3m_rsds_acccmip6.R   # Script to output Fig1-5 and SFig2-10
- Fig6_acccmip6.R                                                      # Script to output Fig6
- data_availability_tas_lai_masks_acccmip6.R                           # Script to output SFig1
