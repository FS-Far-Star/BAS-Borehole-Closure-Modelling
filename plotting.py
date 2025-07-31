import pandas as pd
import matplotlib.pyplot as plt
from global_settings import *

def plot_graphs():
    # --- Load data
    depth           = pd.read_csv("data/depth.csv")["depth"].values
    bore_diameter   = pd.read_csv("data/bore_diameter.csv").values
    delta_bore      = pd.read_csv("data/delta_bore.csv").values
    temp            = pd.read_csv("data/temperature.csv").values
    p_t             = pd.read_csv("data/p_t.csv").values
    p_fluid         = pd.read_csv("data/p_fluid.csv").values
    p_ice           = pd.read_csv("data/p_ice.csv").values
    rho             = pd.read_csv("data/rho.csv").values

    timeseries      = pd.read_csv("data/timeseries.csv")
    t               = timeseries["time"]
    fluid_volume    = timeseries["fluid_volume"]
    fluid_height    = timeseries["fluid_height"]
    bore_depth      = timeseries["bore_depth"]

    # --- Plot borehole profiles every X timesteps
    if plot_bore_profile_over_time:
        plt.figure()
        plt.ylabel('Depth [m]')
        plt.gca().invert_yaxis()

        for i in range(0, bore_diameter.shape[1], plot_spacing):
            nonzero_indices = np.nonzero(bore_diameter[:, i])[0]
            if len(nonzero_indices) > 0:
                last_index = nonzero_indices[-1]
                plt.plot(bore_diameter[:, i][:last_index+1]/2, depth[:last_index+1], '.')
                plt.plot(-bore_diameter[:, i][:last_index+1]/2, depth[:last_index+1], '.')
        nonzero_indices = np.nonzero(bore_diameter[:, -1])[0]
        if len(nonzero_indices) > 0:
                last_index = nonzero_indices[-1]
                plt.plot(bore_diameter[:, -1][:last_index+1]/2, depth[:last_index+1], color='red', linewidth=2)
                plt.plot(-bore_diameter[:, -1][:last_index+1]/2, depth[:last_index+1], color='red', linewidth=2)

        plt.title("Borehole Profiles Over Time")
        plt.grid(True)
        plt.tight_layout()

    # --- Plot borehole minimum non-zero diameter
    if plot_mim_bore_diameter_over_time:
        plt.figure()
        # Mask zero values with NaN (so they don't affect the min)
        masked_diameters = np.where(bore_diameter > 0, bore_diameter, np.nan)

        # Take the minimum ignoring NaNs (i.e., ignoring zeros)
        min_nonzero_diameter = np.nanmin(masked_diameters, axis=0)

        plt.plot(t, min_nonzero_diameter, label='Minimum Non-Zero Diameter')
        plt.ylabel('Borehole Minimum Diameter [mm]')
        plt.xlabel('Time [s]')
        plt.title("Borehole Diameter Over Time")
        plt.grid(True)
        plt.tight_layout()

    # --- Plot borehole closure rate every 10 timesteps
    if plot_borehole_closure_rate_over_time:
        plt.figure()
        plt.xlabel('Diameter Change[mm]')
        plt.ylabel('Depth [m]')
        plt.gca().invert_yaxis()

        for i in range(0, delta_bore.shape[1], plot_spacing):
            nonzero_indices = np.nonzero(delta_bore[:, i])[0]
            if len(nonzero_indices) > 0:
                last_index = nonzero_indices[-1]
                plt.plot(delta_bore[:, i][:last_index+1], depth[:last_index+1], '.')

        plt.title("Borehole Closure Rate Over Time")
        plt.grid(True)
        plt.tight_layout()

    # --- Multi-axis plot for timeseries data
    if plot_bore_status_over_time:
        fig1, ax1 = plt.subplots()

        # Fluid Volume
        color1 = 'tab:blue'
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Fluid Volume [m³]', color=color1)
        ax1.plot(t, fluid_volume, color=color1, label='Fluid Volume')
        ax1.tick_params(axis='y', labelcolor=color1)

        # Fluid Height
        ax2 = ax1.twinx()
        color2 = 'tab:red'
        ax2.set_ylabel('Fluid Height [m]', color=color2)
        ax2.plot(t, fluid_height, color=color2, linestyle='--', label='Fluid Height')
        ax2.tick_params(axis='y', labelcolor=color2)

        # Bore Depth
        ax3 = ax1.twinx()
        color3 = 'tab:green'
        ax3.spines["right"].set_position(("axes", 1.15))
        ax3.set_frame_on(True)
        ax3.patch.set_visible(False)
        ax3.set_ylabel('Bore Depth [m]', color=color3)
        ax3.plot(t, bore_depth, color=color3, linestyle='-.', label='Bore Depth')
        ax3.tick_params(axis='y', labelcolor=color3)

        # Combined legend
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines3, labels3 = ax3.get_legend_handles_labels()
        ax1.legend(lines1 + lines2 + lines3, labels1 + labels2 + labels3, loc='lower right')

        plt.title('Borehole Status Over Time')
        plt.grid(True)
        plt.tight_layout()

    # --- Plot temperature and density profile
    if plot_temperature_density_profile:
        fig2, ax1 = plt.subplots()
        ax1.plot(temp, depth, color='tab:red', label='Temperature')
        ax1.set_xlabel('Temperature [°C]')
        ax1.set_ylabel('Depth [m]')
        ax1.invert_yaxis()
        ax1.tick_params(axis='x', labelcolor='tab:red')
        ax1.set_title('Temperature and Density Profile')

        # Create second y-axis (same depth, different x-scale)
        ax2 = ax1.twiny()
        ax2.plot(rho, depth, color='tab:blue', label='Density')
        ax2.set_xlabel('Density [kg/m³]')
        ax2.tick_params(axis='x', labelcolor='tab:blue')
        plt.tight_layout()

    # --- Plot pressure profiles
    if plot_pressure_profile_over_time:
        plt.figure()
        for i in range(1, p_t.shape[1], plot_spacing):  # pressure is not calculated for 0th time step
            plt.plot(p_t[:int(bore_depth[i]//dh), i]/1e6, depth[:int(bore_depth[i]//dh)])
            plt.plot(p_fluid[:int(bore_depth[i]//dh), i]/1e6, depth[:int(bore_depth[i]//dh)])
        plt.plot(p_ice/1e6, depth,'--', label='p_ice', color='black')

        plt.legend()
        plt.title("Pressure Profiles Over Time")
        plt.grid(True)
        plt.tight_layout()
        plt.xlabel('Pressure [MPa]')
        plt.ylabel('Depth [m]')
        plt.gca().invert_yaxis()

        # show all plots
        plt.show()
    return None
# plot_graphs()