import xarray as xr
import numpy as np
import pandas as pd
import os
import glob
import subprocess


def read_pressure(path, file):
    """
    Read Icefin's pressure log and save in Xarray format
    """
    if "PRESSURE_STAT" not in file:
        raise ValueError(f"The specified file is incorrect.")
    print("Reading Icefin's pressure log...")

    if not os.path.exists(path + file[:-4] + "_test.csv"):
        csv = pd.read_csv(path + file, skip_blank_lines=True, skiprows=4, header=0)
        ds = csv.to_xarray()
        # workaround for making datetime the main dimension for sub-setting data
        # add the datetime to the dataframe
        csv["datetime"] = pd.to_datetime(ds.UNIX_timestamp, unit="s")
        # overwrite the csv file with the datetime column
        csv.to_csv(path + file[:-4] + "_test.csv")

    # parse the csv file again
    from datetime import datetime

    dateparse = lambda x: datetime.strptime(x, "%Y-%m-%d %H:%M:%S.%f")
    csv = pd.read_csv(
        path + file[:-4] + "_test.csv", parse_dates=True, date_parser=dateparse
    )
    # now store it as an Xarray dataset and save
    ds = csv.to_xarray()
    ds = ds.set_coords("datetime").swap_dims({"index": "datetime"})
    # update data type to np.datetime64 so sub-setting is possible
    ds["datetime"] = pd.to_datetime(ds.datetime)
    if not os.path.exists(path + file[:-4] + ".nc"):
        ds.to_netcdf(path + file[:-4] + ".nc", engine="netcdf4")
    return ds


def read_position(path, file):
    """
    Read Icefin's position log and save in Xarray format
    """
    if "NORBIT_POSITION" not in file:
        raise ValueError(f"The specified file is incorrect.")
    print("Reading Icefin's position log...")
    from datetime import datetime

    dateparse = lambda x: datetime.strptime(x, "%Y-%m-%d %H:%M:%S")
    csv = pd.read_csv(
        path + file,
        skip_blank_lines=True,
        skiprows=4,
        header=0,
        parse_dates=True,
        date_parser=dateparse,
    )
    # now store it as an Xarray dataset and save
    ds = csv.to_xarray()
    ds = ds.set_coords("datetime").swap_dims({"index": "datetime"})
    # update data type to np.datetime64 so sub-setting is possible
    ds["datetime"] = pd.to_datetime(ds.datetime)
    if not os.path.exists(path + file[:-4] + ".nc"):
        ds.to_netcdf(path + file[:-4] + ".nc", engine="netcdf4")
    return ds


def get_length(filename):
    result = subprocess.run(
        [
            "ffprobe",
            "-v",
            "error",
            "-show_entries",
            "format=duration",
            "-of",
            "default=noprint_wrappers=1:nokey=1",
            filename,
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return float(result.stdout)


def time_video_UTC(path):
    """
    Extract the video timing from the file names for all the boxes
    """
    print("Reading Icefin's video file metadata...")
    df = pd.DataFrame(
        np.nan, index=range(100), columns=[]
    )  # arbitrary size of 100 is used, because the number of videos in each box in a single deployment is expected to be well below 100
    duration = pd.DataFrame()
    boxes = ["box1", "box2", "box3", "box4", "box5", "box6"]
    for box in boxes:
        path_box = path + box
        files = glob.glob(path_box + "/MARS_ICE*.mov", recursive=True)
        unix_timestamp = np.ones(len(files))
        video_length_s = np.ones(len(files))
        basenames = [os.path.basename(filepath) for filepath in files]
        for i in range(len(files)):
            file = files[i]
            if len(file) >= 29:  # Ensures file has enough characters
                unix_timestamp[i] = float(files[i][-29:-10].replace("_", "."))
                video_length_s[i] = float(get_length(files[i]))
        # sort the table in chronological order
        args = np.argsort(unix_timestamp)
        sorted_names = [basenames[i] for i in args]
        df[box + " start (UTC)"] = pd.Series(unix_timestamp[args])
        df[box + " start (UTC)"] = pd.to_datetime(df[box + " start (UTC)"], unit="s")
        df[box + " duration (minutes)"] = (
            pd.Series(video_length_s[args], dtype="float") / 60
        )
        df[box + " file name"] = pd.Series(sorted_names, dtype="object")
    df = df.dropna(how="all")
    df.to_csv(path + "video_timing_UTC.csv")
    return df


if __name__ == "__main__":

    from socket import gethostname

    hn = gethostname()
    if hn == "funkaholic-koala":
        path = "/media/nina/K892A_Nina/"
    elif hn == "flake20sm1":
        path = "/home/nina/"

    #     # process the pressure and position data and video timing from 1 November 2021
    #     read_pressure(
    #         path + "MARS_ICE02_032/MARS_ICE02_032_CSV/",
    #         "MARS_ICE02_032_telemetry.gssbin_OPENINS_PRESSURE_STAT.csv",
    #     )
    #     read_position(
    #         path + "MARS_ICE02_032/MARS_ICE02_032_CSV/",
    #         "MARS_ICE02_032_telemetry.gssbin_NORBIT_POSITION.csv",
    #     )
    #     time_video_UTC(path + "MARS_ICE02_032/MARS_ICE02_032_VIDEO/")
    # process the pressure and position data and video timing from 2 November 2021
    read_pressure(
        path + "MARS_ICE02_033/MARS_ICE02_033_CSV/",
        "MARS_ICE02_033_telemetry.gssbin_OPENINS_PRESSURE_STAT.csv",
    )
    read_position(
        path + "MARS_ICE02_033/MARS_ICE02_033_CSV/",
        "MARS_ICE02_033_telemetry.gssbin_NORBIT_POSITION.csv",
    )
    time_video_UTC(path + "MARS_ICE02_033/MARS_ICE02_033_VIDEO/")
    # process the pressure and position data and video timing from 3 November 2021
    read_pressure(
        path + "MARS_ICE02_034/MARS_ICE02_034_CSV/",
        "MARS_ICE02_034_1_telemetry.gssbin_OPENINS_PRESSURE_STAT.csv",
    )
    read_position(
        path + "MARS_ICE02_034/MARS_ICE02_034_CSV/",
        "MARS_ICE02_034_1_telemetry.gssbin_NORBIT_POSITION.csv",
    )
    time_video_UTC(path + "MARS_ICE02_034/MARS_ICE02_034_VIDEO/")
