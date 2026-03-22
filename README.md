# README

## Requirements
* seawater
* gsw
* Python 3.8.1
* Numpy
* Matplotlib
* OpenCV
* pyTMD
* Xarray
* Pandas
* cmcrameri

## Code description
ADCP/calibration.py: reads .mat file, applies Purdie 2012 calibration and saves data in netCDF format
ADCP/Wind_rose_plot.ipynb: plots Figure 10
ADCP/TMD_currents.py: used to generate the timeseries of modelled tidal current speeds

CT_2021/read.py: parses CNV files and stores netCDF file
CT_2021/CT_arm_processing.ipynb: runs the parser, calculates supercooling and stores netCDF file
CTD_2021/plot_CTD.py: visualisation routines
CTD_2021/error_analysis.py: routine to add error bars to CTD plots'
CTD_2021/SBE_cast2_3.ipynb: plots all good SBE CTD data
Icefin/analyze_frazil_2Nov.ipynb: processes 4k videos into frames and plots depth of Icefin during videos
Icefin/Depth_of_Field.ipynb: calculates depth of field for SubC 4
Icefin/frazil_size_3Nov.ipynb: extracts apparent diameters and apparent concentrations from Icefin videos
Icefin/analyze_videos.py: routines for frazil ice apparent diameter and concentrations
Icefin/read.py: parser for the pressure data from CSV files
ADCP_and_Icefin.ipynb: comparison plots of apparent frazil ice concentration (Icefin) and frazil ice proxy (ADCP)
ADCP_timeseries_analysis.ipynb: timeseries plots Figure 11
all_supercooling_2021.ipynb: supercooling timeseries using both CT and CTD data
FIC_SC.ipynb: combines apparent frazil ice concentrations from Icefin with CTD data
frazil_size_and_concentration_uncertainty.ipynb: estimates for uncertainty apparent frazil ice concentration
supercooling_uncertainty.ipynb: error estimates for CTD data

