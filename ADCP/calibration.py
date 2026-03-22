import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
import xarray as xr
import os
import sys
import pandas as pd
import datetime
import os
from TMD_currents import*
os.chdir("../..")
savedir = os.path.join(os.getcwd(),"DATA")
path = os.getcwd()
path = os.path.join(path,"DATA/SIOS21/adcp_data_and_analysis_lars_smedsrud/")

# matfile = loadmat(path+"Nortek_ADCP_currents_Oct20_to_Nov4_avgd.mat")
matfile = loadmat(path+"/Nortek_ADCP_currents_Oct20_to_Nov4_avgd.mat")

# depth_profile is the average depth of the timeseries, for usage with imshow
# depth_bins is the true depth of the bins
ds = xr.Dataset(coords=dict(time = matfile["time"][:,0], \
                            depth_profile = matfile["bin_centres"][1:,0].mean()+np.arange(matfile["amplitude_average"].shape[1])))

ds["depth_bins"] = (["depth_profile", "time"], matfile["bin_centres"].transpose())
ds["amplitude"] = (["depth_profile", "time"], matfile["amplitude_average"].transpose())
ds["instrument_depth"] = ("time", matfile["depth"][:,0])
ds.coords["u"] = ("time", matfile["u"].mean(axis=1)*100)
# ds.coords["u"] = ("time", matfile["u"][:,0]*100)
ds.coords["v"] = ("time", matfile["v"].mean(axis=1)*100)
ds.coords["w"] = ("time", matfile["w"].mean(axis=1)*100)
# ds.coords["v"] = ("time", matfile["v"][:,0]*100)
ds.coords["temperature"] = ("time", matfile["temperature"][:,0])

# plot the tidal height, temperature, velocity and echogram as subplots
# add the packages from AZFP
# convert matlab datenum to numpy datetime64[ns]
time = matfile["time"][:,0]
python_datetime = np.zeros_like(time).astype("datetime64[ns]")
for i in range(len(time)):
    python_datetime[i] = pd.Timestamp(datetime.datetime.fromordinal(int(matfile["time"][i,0])) + datetime.timedelta(days=matfile["time"][i,0]%1) - datetime.timedelta(days = 366))
# convert NZDT to UTC
ds.coords["time"] = python_datetime + np.timedelta64(13, "h")

# evaluate tidal currents
# camp 2021 coordinates
lon = 166.23333333333332
lat = -77.86666666666666
dattime = ds.time.data
# model the tidal current by using the CATS2008 wrapper for single point tidal timeseries in McMurdo Sound
if not os.path.exists(path+"tide.npz"):
    tide = tidalCurrentAZFP(dattime, lon, lat)
    tide["z"] = tidalHeightAZFP(dattime, lon, lat)
# save the tidal current
    np.savez(path+"tide.npz", tide=tide, allow_pickle=True)
# or load the tidal current to save time
else:
    npz_file = np.load(path+"tide.npz", allow_pickle=True)
    tide = npz_file['tide'].item()   
    
# add the tides to the data set
ds.coords["tide_z"] = ("time", tide["z"])
ds.coords["tide_u"] = ("time", tide["u"].reshape(tide["u"].shape[2]))
ds.coords["tide_v"] = ("time", tide["v"].reshape(tide["v"].shape[2]))

# save the data for further analysis in MMS2021_CATS2008_ADCP_comparison or other notebooks
if os.path.exists(path+"Nortek_ADCP_currents_Oct20_to_Nov4_avgd.nc"):
    os.remove(path+"Nortek_ADCP_currents_Oct20_to_Nov4_avgd.nc")
    
ds.to_netcdf(path+"Nortek_ADCP_currents_Oct20_to_Nov4_avgd.nc", engine="netcdf4")

## CALIBRATION
#ADCP 1000 kHz calibration (Purdie 2012)
#$$BS_v = EL + 10\log_{10} (R^2) + 2 \alpha_w R + 20 R \int \alpha_p dr$$
#$$= 0.43 ADCP_{counts}+20\log_{10}(R) + 2 \alpha_w R + 20 R \int \alpha_p dr$$

#where $\alpha_w$=0.56 dB m$^{-1}$ (http://resource.npl.co.uk/acoustics/techguides/seaabsorption/ using salinity 34.5 and temperature -1.9)

#ignore $\alpha_p$ for first order approximation, because the water was clear

#R is depth-depth[0], because the ADCP is pointing down

if os.path.exists(path+"Nortek_ADCP_currents_Oct20_to_Nov4_avgd_cal.nc"):
    os.remove(path+"Nortek_ADCP_currents_Oct20_to_Nov4_avgd_cal.nc")

alpha_w = 0.067 # using salinity 34.64 psu; temperature -1.94 degrees C; depth 0.01 km; pH 8.03 (Matson et al., 2011) and the tool mentioned above
R = ds.depth_bins-ds.depth_bins[0,:]
BSv = ds.amplitude + 10*np.log10(R**2) + 2 * alpha_w * R
ds["BSv"] = BSv

ds.to_netcdf(path+"Nortek_ADCP_currents_Oct20_to_Nov4_avgd_cal.nc", engine="netcdf4")
    
