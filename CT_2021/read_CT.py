import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import re

def lines_that_contain(string, fp):
    k=0
    for line in fp:
        if string in line:
            break
        else:
            k+=1
    return k

def read_single_cnv(path,file):
    filename = path+file
    # read individual lines of CNV file
    # Using readlines()
    file1 = open(filename, 'r')
    Lines = file1.readlines()

    # read line 19-26 to load the variable names
    # find line starting with "# name 0 ="
    var_name_line_0 = lines_that_contain("# name 0 =",Lines)
    var_name_lines = np.arange(var_name_line_0,var_name_line_0+9)+1
    var_names = list()
    k=0
    for i in var_name_lines-1:
        var_names.append(Lines[i])
        # select last part of the line to name variable
        cut_string_here = var_names[k].index(": ")+2
        var_names[k] = var_names[k][cut_string_here:-1]
        k+=1
    # overwrite variable names with manual entries to allow saving NetCDF
    var_names[0] = "Conductivity_S_m"
    var_names[1] = "Depth_m"
    var_names[2] = "Salinity_PSU"
    var_names[3] = "Scan_count"
    var_names[4] = "Temperature_DegC"
    var_names[5] = "Time_Elapsed_s"
    var_names[6] = "Time_System_s"
    var_names[7] = "Pressure_db"
    var_names[8] = "oo"

    # fill in the columns with data from the CNV file
    # select line 126 until last line
    # find line containing *END* and start reading on next line
    data_line_0 = lines_that_contain("*END*", Lines)+1
    selected = Lines[data_line_0:]
    # split with one or more spaces as delimiter
    data = list()
    for i in range(len(selected)):
        data.append(re.split(r" {1,}", selected[i].strip()))
    # convert to array with float64
    data = np.array(data).astype("float64")
    return data, var_names

def read_write_all_cnv(path, netcdf=True):
    # select only CNV and sort file names
    files = [ f for f in os.listdir(path) if f.endswith('.cnv') ]
    files.sort(key=lambda f: int(re.sub('\D', '', f)))
    # read all files and store in dataframe
    if not files:
        print("The file list is empty.")
    k=0
    for file in files:
        if k==0: # first data set
            [data, var_names] = read_single_cnv(path,file)
            df = pd.DataFrame(columns = var_names, data = data)
            k+=1
        else: # next data sets
            [data, var_names] = read_single_cnv(path,file)
            df = pd.concat([df, pd.DataFrame(columns = var_names, data = data)], ignore_index=True)
    if netcdf:
        import xarray as xr
        # add time coordinate
        time = pd.to_datetime(df.Time_System_s, unit="s").values
        ds = df.to_xarray()
        ds = ds.assign_coords(time = ("index", time))
        ds.to_netcdf(path+"CT_arm.nc", mode="w", engine='netcdf4')
        return ds
    else:
        df.to_csv(path+"CT_arm.csv")
        return df
