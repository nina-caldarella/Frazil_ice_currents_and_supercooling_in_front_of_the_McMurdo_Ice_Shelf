import os
import pyproj
import datetime
import numpy as np
import matplotlib
matplotlib.rcParams['axes.linewidth'] = 2.0
matplotlib.rcParams["animation.html"] = "jshtml"
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import ipywidgets as widgets
from IPython.display import HTML

# import tide programs
import pyTMD.io
import pyTMD.time
import pyTMD.predict
import pyTMD.tools
import pyTMD.utilities

# imports for tidal heights
# from __future__ import print_function
import IPython.display


def tidalCurrentAZFP(dattime, lon, lat):
    """_summary_

    Args:
        dattime (array): time in datetime.datetime format or np.datetime64
        lon (_type_): _description_
        lat (_type_): _description_

    Returns:
        dict: tide
    """
    # available model list
    model_list = sorted(pyTMD.io.model.global_current() + pyTMD.io.model.antarctic_current())
    # display widgets for setting directory and model
    TMDwidgets = pyTMD.tools.widgets()
    TMDwidgets.model.options = model_list
    TMDwidgets.model.value = 'CATS2008'
    widgets.VBox([
        TMDwidgets.directory,
        TMDwidgets.model,
        TMDwidgets.atlas,
        TMDwidgets.compress,
        TMDwidgets.datepick
    ])

    # get model parameters
    model = pyTMD.io.model('/home/nina/CODE',
        format=TMDwidgets.atlas.value,
        compressed=TMDwidgets.compress.value
    ).current(TMDwidgets.model.value)

    # # convert from calendar date to days relative to Jan 1, 1992 (48622 MJD)
    # use the AZFP time instead
    # convert to np.datetime64 to datetime.datetime
    if type(dattime[0])==np.datetime64:
        unix_epoch = np.datetime64(0, 's')
        one_second = np.timedelta64(1, 's')
        seconds_since_epoch = (dattime - unix_epoch) / one_second
        seconds_since_epoch = seconds_since_epoch.astype("int")
        dt = []
        for i in range(len(seconds_since_epoch)):
            dt.append(datetime.datetime.utcfromtimestamp(seconds_since_epoch[i]))
    else:
        dt = dattime
    tide_time = np.zeros(len(dattime))
    for i in range(len(dattime)):
        tide_time[i] = timescale.convert_calendar_dates(dt[i].year, dt[i].month, dt[i].day, dt[i].hour, dt[i].minute)
    # lon = 166.6472
    # lat = -77.8946
    # save tide currents
    tide_samples = np.arange(0, len(dattime), 1)
    tide = {}
    # iterate over u and v currents
    for TYPE in model.type:
        # read tidal constants and interpolate to grid points
        if model.format in ('OTIS','ATLAS','ESR'):
            amp,ph,D,c = pyTMD.io.OTIS.extract_constants(lon, lat, model.grid_file,
                model.model_file['u'], model.projection, type=TYPE,
                method='spline', grid=model.format)
            DELTAT = np.zeros_like(tide_time)
        elif (model.format == 'netcdf'):
            amp,ph,D,c = pyTMD.io.ATLAS.extract_constants(lon, lat, model.grid_file,
                model.model_file[TYPE], type=TYPE, method='spline',
                scale=model.scale, compressed=model.compressed)
            DELTAT = np.zeros_like(tide_time)
        elif (model.format == 'GOT'):
            amp,ph,c = pyTMD.io.GOT.extract_constants(lon, lat, model.model_file[TYPE],
                method='spline', scale=model.scale,
                compressed=model.compressed)
            # interpolate delta times from calendar dates to tide time
            DELTAT = pyTMD.time.interpolate_delta_time(delta_file, tide_time)
        elif (model.format == 'FES'):
            amp,ph = pyTMD.io.FES.extract_constants(lon, lat, model.model_file[TYPE],
                type=TYPE, version=model.version, method='spline',
                scale=model.scale, compressed=model.compressed)
            c = model.constituents
            # interpolate delta times from calendar dates to tide time
            DELTAT = pyTMD.time.interpolate_delta_time(delta_file, tide_time)
                # calculate complex phase in radians for Euler's
        cph = -1j*ph*np.pi/180.0
        # calculate constituent oscillation
        hc = amp*np.exp(cph)

        # allocate for tide current map calculated every hour
        tide[TYPE] = np.ma.zeros((1,1,len(tide_samples)))
        k=0
        for i in tide_samples:
            # predict tidal elevations at time and infer minor corrections
            TIDE = pyTMD.predict.map(tide_time[i], hc, c, deltat=DELTAT[i],
                corrections=model.format)
            MINOR = pyTMD.predict.infer_minor(tide_time[i], hc, c,
                deltat=DELTAT[i], corrections=model.format)
            # add major and minor components and reform grid
            tide[TYPE][:,:,k] = np.reshape((TIDE+MINOR),(1,1))
            k+=1
    return tide

def tidalHeightAZFP(dattime, lon, lat):
    """_summary_

    Args:
        dattime (array): time in datetime.datetime format or np.datetime64
        lon (_type_): _description_
        lat (_type_): _description_

    Returns:
        dict: tide
    """
    # available model list
    model_list = sorted(pyTMD.io.model.ocean_elevation())
    # display widgets for setting directory and model
    TMDwidgets = pyTMD.tools.widgets()
    TMDwidgets.model.options = model_list
    TMDwidgets.model.value = 'CATS2008'
    widgets.VBox([
        TMDwidgets.directory,
        TMDwidgets.model,
        TMDwidgets.atlas,
        TMDwidgets.compress,
        TMDwidgets.datepick])
    # default coordinates to use
    LAT,LON = (lat,lon)
    # get model parameters
    model = pyTMD.io.model('/home/nina/CODE',
        format=TMDwidgets.atlas.value,
        compressed=TMDwidgets.compress.value
    ).elevation(TMDwidgets.model.value)

    if type(dattime[0])==np.datetime64:
        unix_epoch = np.datetime64(0, 's')
        one_second = np.timedelta64(1, 's')
        seconds_since_epoch = (dattime - unix_epoch) / one_second
        seconds_since_epoch = seconds_since_epoch.astype("int")
        dt = []
        for i in range(len(seconds_since_epoch)):
            dt.append(datetime.datetime.utcfromtimestamp(seconds_since_epoch[i]))
    else:
        dt = dattime
    tide_time = np.zeros(len(dattime))
    for i in range(len(dattime)):
        tide_time[i] = pyTMD.time.convert_calendar_dates(dt[i].year, dt[i].month, dt[i].day, dt[i].hour, dt[i].minute)

    # delta time (TT - UT1) file
    delta_file = pyTMD.utilities.get_data_path(['data','merged_deltat.data'])

    # read tidal constants and interpolate to leaflet points
    if model.format in ('OTIS','ATLAS','TMD3'):
        constituents = pyTMD.io.OTIS.read_constants(
            model.grid_file, model.model_file,
            model.projection, type=model.type,
            grid=model.format)
        c = constituents.fields
        DELTAT = np.zeros_like(tide_time)
    elif (model.format == 'netcdf'):
        constituents = pyTMD.io.ATLAS.read_constants(
            model.grid_file, model.model_file,
            type=model.type, compressed=model.compressed)
        c = constituents.fields
        DELTAT = np.zeros_like(tide_time)
    elif (model.format == 'GOT'):
        constituents = pyTMD.io.GOT.read_constants(
            model.model_file, compressed=model.compressed)
        c = constituents.fields
        # interpolate delta times from calendar dates to tide time
        DELTAT = pyTMD.time.interpolate_delta_time(delta_file, tide_time)
    elif (model.format == 'FES'):
        constituents = pyTMD.io.FES.read_constants(model.model_file,
            type=model.type, version=model.version,
            compressed=model.compressed)
        c = model.constituents
        # interpolate delta times from calendar dates to tide time
        DELTAT = pyTMD.time.interpolate_delta_time(delta_file, tide_time)

    # update the tide prediction and plot
    def update_tide_prediction(*args):
        if model.format in ('OTIS','ATLAS','TMD3'):
            amp,ph,D = pyTMD.io.OTIS.interpolate_constants(
                np.atleast_1d(LON), np.atleast_1d(LAT),
                constituents, model.projection, type=model.type,
                method='spline', extrapolate=True)
        elif (model.format == 'netcdf'):
            amp,ph,D = pyTMD.io.ATLAS.interpolate_constants(
                np.atleast_1d(LON), np.atleast_1d(LAT),
                constituents, type=model.type, scale=model.scale,
                method='spline', extrapolate=True)
        elif (model.format == 'GOT'):
            amp,ph = pyTMD.io.GOT.interpolate_constants(
                np.atleast_1d(LON), np.atleast_1d(LAT),
                constituents, scale=model.scale,
                method='spline', extrapolate=True)
        elif (model.format == 'FES'):
            amp,ph = pyTMD.io.FES.interpolate_constants(
                np.atleast_1d(LON), np.atleast_1d(LAT),
                constituents, scale=model.scale,
                method='spline', extrapolate=True)
        # calculate complex phase in radians for Euler's
        cph = -1j*ph*np.pi/180.0
        # calculate constituent oscillation
        hc = amp*np.exp(cph)
        # predict tidal elevations at time 1 and infer minor corrections
        TIDE = pyTMD.predict.time_series(tide_time, hc, c,
            deltat=DELTAT, corrections=model.format)
        MINOR = pyTMD.predict.infer_minor(tide_time, hc, c,
            deltat=DELTAT, corrections=model.format)
        TIDE.data[:] += MINOR.data[:]
        # convert to centimeters
        TIDE.data[:] *= 100.0

        # differentiate to calculate high and low tides
        diff = np.zeros_like(tide_time, dtype=np.float64)
        # forward differentiation for starting point
        diff[0] = TIDE.data[1] - TIDE.data[0]
        # backward differentiation for end point
        diff[-1] = TIDE.data[-1] - TIDE.data[-2]
        # centered differentiation for all others
        diff[1:-1] = (TIDE.data[2:] - TIDE.data[0:-2])/2.0
        return TIDE
    # run tide prediction at initial location
    TIDE = update_tide_prediction()
    return TIDE.data

if __name__ == "__main__":
    from socket import gethostname
    from scipy.io import loadmat
    import xarray as xr
    # define path
    hn = gethostname()
    if hn=='funkaholic-koala':
        path = "/media/nina/K892A_Nina/DATA/adcp_data_and_analysis_lars_smedsrud/"
    elif hn=='flake20sm1':
        path = "/home/nina/DATA/adcp_data_and_analysis_lars_smedsrud/"
    matfile = loadmat(path+"/Nortek_ADCP_currents_Oct20_to_Nov4_avgd.mat")
    # convert some of the matlab variables to Xarray
    ds = xr.Dataset(coords=dict(time = matfile["time"][:,0], bins = np.arange(20)))
    ds["SvArrays"] = (["bins", "time"], matfile["amplitude_average"].transpose())
    ds["depth"] = ("time", matfile["depth"][:,0])
    ds["u"] = ("time", matfile["u"][:,0]*100)
    ds["v"] = ("time", matfile["v"][:,0]*100)
    ds["temperature"] = ("time", matfile["temperature"][:,0])
    # plot the tidal height, temperature, velocity and echogram as subplots
    # add the packages from AZFP
    import sys
    sys.path.append("/home/nina/CODE/ninas_PhD")
    from analysis.plotfrazil import *
    from azfp_data_processing.TMD_currents import *
    # convert matlab datenum to numpy datetime64[ns]
    import datetime
    time = matfile["time"][:,0]
    python_datetime = np.zeros_like(time).astype("datetime64[ns]")
    for i in range(len(time)):
        python_datetime[i] = pd.Timestamp(datetime.datetime.fromordinal(int(matfile["time"][i,0])) + datetime.timedelta(days=matfile["time"][i,0]%1) - datetime.timedelta(days = 366))

    ds.coords["time"] = python_datetime
    # evaluate tidal currents
    # camp 2021 coordinates
    lon = 166.23333333333332
    lat = -77.86666666666666
    dattime = ds.time.data
    TIDE = tidalHeightAZFP(dattime, lon, lat)
    print(" ")
