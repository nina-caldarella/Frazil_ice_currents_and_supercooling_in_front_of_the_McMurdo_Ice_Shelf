import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import os

def plot_up(ds, uT_C, res, fig, ax, mid_time, end_time, save_fig=None, color = "black"):
    up_start = mid_time
    up_end = end_time
    time_up = up_start

    up_cast = ds.sel(time=slice(up_start, up_end))
    u_rho = res.sel(time=slice(up_start, up_end)).u_rho.data

    ax.plot(
        up_cast.Temperature,
        up_cast.depth, 
        label="up on " + time_up.strftime("%d-%m") + " at " + time_up.strftime("%H:%M"),
        color = color,
        linewidth=3
    )
    ax.fill_betweenx(up_cast.depth, up_cast.Temperature - uT_C, up_cast.Temperature + uT_C, alpha=0.18, color=color)
    ax.set_xlim([-1.96, -1.90])
    ax.set_ylim([0, 550])
    ax.set_xlabel("In-situ temperature t90 ($^o$C) \n (a)")
    ax.set_ylabel("Depth (m)")
    ax.invert_yaxis()
    ax.tick_params(axis="x", labelrotation=45)
    # ax.legend(loc="upper right")

#    ax2.plot(
#        up_cast.rho - 1000,
#        up_cast.depth,
#        label="up on " + time_up.strftime("%d-%m") + " at " + time_up.strftime("%H:%M"),
#        color = color
#    )
#    ax2.fill_betweenx(up_cast.depth, up_cast.rho - 1000 - u_rho, up_cast.rho - 1000 + u_rho, alpha=0.18, color=color)
#    ax2.set_xlim([27.8, 27.9])
#    ax2.set_ylim([0, 550])
#    ax2.set_xlabel(r"Potential density (kg m$^{-3}$)")
#    # ax2.set_ylabel("Pressure (dbar)")
#    ax2.invert_yaxis()
#    ax2.tick_params(axis="x", labelrotation=45)
#    ax2.legend()
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.tight_layout()
    plt.rcParams["font.size"] = 22

    return fig, ax


def plot_down_up(ds, start_time, mid_time, end_time, save_fig=None):
    down_start = start_time
    down_end = mid_time
    time_down = down_start
    up_start = mid_time
    up_end = end_time
    time_up = up_start

    down_cast = ds.sel(time=slice(down_start, down_end))
    up_cast = ds.sel(time=slice(up_start, up_end))

    fig, ((ax, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 15))
    ax.plot(
        down_cast.Temperature,
        down_cast.depth,
        label="down at " + time_down.strftime("%m-%d %H:%M"),
        color="black",
    )
    ax.plot(
        up_cast.Temperature,
        up_cast.depth,
        label="up at " + time_up.strftime("%m-%d %H:%M"),
        color="blue",
    )
    ax.set_xlim([-1.96, -1.90])
    ax.set_ylim([0, 550])
    ax.set_xlabel("Temperature t90 (deg C) \n (a)")
    ax.set_ylabel("Depth (m)")
    ax.invert_yaxis()
    ax.tick_params(axis="x", labelrotation=45)
    # ax.legend(loc="upper right")

    ax2.plot(
        down_cast.SP,
        down_cast.depth,
        color="black",
        label="down at " + time_down.strftime("%H:%M"),
    )
    ax2.plot(
        up_cast.SP,
        up_cast.depth,
        color="blue",
        label="up at " + time_up.strftime("%H:%M"),
    )
    ax2.set_xlim([34.61, 34.75])
    ax2.set_ylim([0, 550])
    ax2.set_xlabel("Practical salinity (psu) \n (b)")
    # ax2.set_ylabel("Pressure (dbar)")
    ax2.invert_yaxis()
    ax2.tick_params(axis="x", labelrotation=45)

    ax3.plot(
        down_cast.rho,
        down_cast.depth,
        color="black",
        label="down at " + time_down.strftime("%m-%d %H:%M"),
    )
    ax3.plot(
        up_cast.rho,
        up_cast.depth,
        color="blue",
        label="up at " + time_up.strftime("%m-%d %H:%M"),
    )
    ax3.set_xlim([1027.85, 1028.0])
    ax3.set_ylim([0, 550])
    ax3.set_xlabel(r"Potential density (kg m$^{-3}$)"
                    "\n"
                    "(c)")
    # ax2.set_xlabel(r"Practical salinity (psu)")
    ax3.set_ylabel("Depth (m)")
    ax3.invert_yaxis()
    ax3.tick_params(axis="x", labelrotation=45)
    ax3.legend(loc="best")
    # conductivity plot
    ax4.plot(
        down_cast.con,
        down_cast.depth,
        color="black",
        label="down at " + time_down.strftime("%m-%d %H:%M"),
    )
    ax4.plot(
        up_cast.con,
        up_cast.depth,
        color="blue",
        label="up at " + time_up.strftime("%m-%d %H:%M"),
    )
    ax4.set_xlim([27.1, 27.5])
    ax4.set_ylim([0, 550])
    ax4.set_xlabel(r"Conductivity (S m$^{-1}$)"
                    "\n"
                    "(d)")
    ax4.invert_yaxis()
    ax4.tick_params(axis="x", labelrotation=45)
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.tight_layout()
    plt.rcParams["font.size"] = 22
#    ax.set_title("cast " + str(start_time) + " to " + str(end_time))

    if save_fig != None:
        fig.savefig(save_fig, bbox_inches="tight")
