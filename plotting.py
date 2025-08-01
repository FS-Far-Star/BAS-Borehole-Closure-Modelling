import pandas as pd
import matplotlib.pyplot as plt
from global_settings import *
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np

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

    # --- Shared color settings
    norm = mcolors.Normalize(vmin=0, vmax=bore_diameter.shape[1])
    colormap = cm.get_cmap('viridis')
    sm = cm.ScalarMappable(cmap=colormap, norm=norm)
    sm.set_array([])

    # --- Plot borehole profiles
    if plot_bore_profile_over_time:
        fig, ax = plt.subplots()
        ax.set_ylabel('Depth [m]')
        ax.set_xlabel('Borehole Walls [mm]')
        ax.invert_yaxis()

        for i in range(0, bore_diameter.shape[1], plot_spacing):
            last_index = int(bore_depth.iloc[i] // dh) - 1
            color = colormap(norm(i))
            ax.plot(bore_diameter[:last_index+1, i] / 2, depth[:last_index+1], '.', color=color)
            ax.plot(-bore_diameter[:last_index+1, i] / 2, depth[:last_index+1], '.', color=color)

        last_index = int(bore_depth.iloc[-1] // dh) - 1
        ax.plot(bore_diameter[:last_index+1, -1] / 2, depth[:last_index+1], color='red', linewidth=2, label='Final')
        ax.plot(-bore_diameter[:last_index+1, -1] / 2, depth[:last_index+1], color='red', linewidth=2)

        ax.set_title("Borehole Profiles Over Time")
        ax.grid(True)
        fig.tight_layout()
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Timestep Index")

    # --- Plot minimum bore diameter
    if plot_mim_bore_diameter_over_time:
        fig, ax = plt.subplots()
        masked_diameters = np.where(bore_diameter > 0, bore_diameter, np.nan)
        min_nonzero_diameter = np.nanmin(masked_diameters, axis=0)
        ax.plot(t, min_nonzero_diameter)
        ax.set_ylabel('Borehole Minimum Diameter [mm]')
        ax.set_xlabel('Time [s]')
        ax.set_title("Borehole Diameter Over Time")
        ax.grid(True)
        fig.tight_layout()

    # --- Plot closure rate
    if plot_borehole_closure_rate_over_time:
        fig, ax = plt.subplots()
        ax.set_xlabel('Diameter Change [mm]')
        ax.set_ylabel('Depth [m]')
        ax.invert_yaxis()

        for i in range(0, delta_bore.shape[1], plot_spacing):
            last_index = int(bore_depth.iloc[i] // dh) - 1
            color = colormap(norm(i))
            ax.plot(delta_bore[:last_index+1, i], depth[:last_index+1], '.', color=color)

        ax.set_title("Borehole Closure Rate Over Time")
        ax.grid(True)
        fig.tight_layout()
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Timestep Index")

    # --- Plot bore status
    if plot_bore_status_over_time:
        fig, ax1 = plt.subplots()

        color1 = 'tab:blue'
        ax1.set_xlabel('Time')
        ax1.set_ylabel('Fluid Volume [m³]', color=color1)
        ax1.plot(t, fluid_volume, color=color1, label='Fluid Volume')
        ax1.tick_params(axis='y', labelcolor=color1)

        ax2 = ax1.twinx()
        color2 = 'tab:red'
        ax2.set_ylabel('Fluid Height [m]', color=color2)
        ax2.plot(t, fluid_height, color=color2, linestyle='--', label='Fluid Height')
        ax2.tick_params(axis='y', labelcolor=color2)
        ax2.invert_yaxis()

        ax3 = ax1.twinx()
        color3 = 'tab:green'
        ax3.spines["right"].set_position(("axes", 1.15))
        ax3.set_frame_on(True)
        ax3.patch.set_visible(False)
        ax3.set_ylabel('Bore Depth [m]', color=color3)
        ax3.plot(t, bore_depth, color=color3, linestyle='-.', label='Bore Depth')
        ax3.tick_params(axis='y', labelcolor=color3)
        ax3.invert_yaxis()

        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        lines3, labels3 = ax3.get_legend_handles_labels()
        ax1.legend(lines1 + lines2 + lines3, labels1 + labels2 + labels3, loc='lower right')

        ax1.set_title('Borehole Status Over Time')
        ax1.grid(True)
        fig.tight_layout()

    # --- Plot temperature and density profile
    if plot_temperature_density_profile:
        fig, ax1 = plt.subplots()
        ax1.plot(temp, depth, color='tab:red', label='Temperature')
        ax1.set_xlabel('Temperature [°C]')
        ax1.set_ylabel('Depth [m]')
        ax1.invert_yaxis()
        ax1.tick_params(axis='x', labelcolor='tab:red')
        ax1.set_title('Temperature and Density Profile')

        ax2 = ax1.twiny()
        ax2.plot(rho, depth, color='tab:blue', label='Density')
        ax2.set_xlabel('Density [kg/m³]')
        ax2.tick_params(axis='x', labelcolor='tab:blue')
        fig.tight_layout()

    # --- Plot pressure profiles
    if plot_pressure_profile_over_time:
        fig, ax = plt.subplots()
        for i in range(1, p_t.shape[1], plot_spacing):
            last_index = int(bore_depth.iloc[i] // dh) - 1
            color = colormap(norm(i))
            ax.plot(p_t[:last_index+1, i]/1e6, depth[:last_index+1], '-', color=color, alpha=0.8, label=None)
            ax.plot(p_fluid[:last_index+1, i]/1e6, depth[:last_index+1], '--', color=color, alpha=0.5, label=None)
        ax.plot(p_ice/1e6, depth, '--', label='p_ice', color='black')
        ax.set_title("Pressure Profiles Over Time")
        ax.set_xlabel('Pressure [MPa]')
        ax.set_ylabel('Depth [m]')
        ax.invert_yaxis()
        ax.grid(True)
        fig.tight_layout()
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Timestep Index")

    # --- Merged plot: 3 key plots side by side with shared colorbar
    if merged_plot:
        fig, axs = plt.subplots(1, 3, figsize=(18, 6), sharey=True, constrained_layout=True)
        fig.suptitle("Borehole Evolution Over Time (Shared Color Scale)", fontsize=14)

        # 1. Borehole Profile
        axs[0].set_title("Borehole Profile")
        axs[0].set_xlabel("Radius [mm]")
        axs[0].set_ylabel("Depth [m]")
        axs[0].invert_yaxis()
        for i in range(0, bore_diameter.shape[1], plot_spacing):
            last_index = int(bore_depth.iloc[i] // dh) - 1
            color = colormap(norm(i))
            axs[0].plot(bore_diameter[:last_index+1, i]/2, depth[:last_index+1], '.', color=color)
            axs[0].plot(-bore_diameter[:last_index+1, i]/2, depth[:last_index+1], '.', color=color)
        last_index = int(bore_depth.iloc[-1] // dh) - 1
        axs[0].plot(bore_diameter[:last_index+1, -1] / 2, depth[:last_index+1], color='red', linewidth=2)
        axs[0].plot(-bore_diameter[:last_index+1, -1] / 2, depth[:last_index+1], color='red', linewidth=2)

        # 2. Closure Rate
        axs[1].set_title("Closure Rate")
        axs[1].set_xlabel("ΔDiameter [mm]")
        for i in range(0, delta_bore.shape[1], plot_spacing):
            last_index = int(bore_depth.iloc[i] // dh) - 1
            color = colormap(norm(i))
            axs[1].plot(delta_bore[:last_index+1, i], depth[:last_index+1], '.', color=color)
        axs[1].plot(delta_bore[:last_index+1, -1], depth[:last_index+1], color='red', linewidth=2)

        # 3. Pressure
        axs[2].set_title("Pressure")
        axs[2].set_xlabel("Pressure [MPa]")
        for i in range(1, p_t.shape[1], plot_spacing):
            last_index = int(bore_depth.iloc[i] // dh) - 1
            color = colormap(norm(i))
            axs[2].plot(p_t[:last_index+1, i]/1e6, depth[:last_index+1], '-', color=color, alpha=0.8)
            axs[2].plot(p_fluid[:last_index+1, i]/1e6, depth[:last_index+1], '--', color=color, alpha=0.4)
        axs[2].plot(p_ice/1e6, depth, '--', color='black', label='p_ice')
        axs[2].legend()

        # Shared colorbar
        sm = cm.ScalarMappable(cmap=colormap, norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axs.ravel().tolist(), orientation='vertical', pad=0.02)
        cbar.set_label("Timestep Index")
    
    # --- Show all plots
    plt.show()
    return None
