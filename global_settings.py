import numpy as np
import scipy
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import math
import warnings
import os
import pandas as pd

# Ensure data directory exists
os.makedirs("data", exist_ok=True)

# === Ice settings ===
H = 800                            # ice thickness (m)
acc = 0.37/0.910                   # surface accumulation (m/yr)
Ts = -22.3                         # Surface temp (°C)
QG = 50 * pow(10,-3)               # geothermal flux (W/m^2)
melt = 0.0                         # meltrate (m/yr)

# === Simulation parameters ===
max_time = 60*60*24*365             # seconds
dt = 24*60*60                       # time step in seconds, this setting only works if NOT loading for spreadsheet 
dh = 1                              # depth step in meters

# Load existing data?
loading = True                     # Set to True if loading existing data
current_bore_diameter = np.ones(int(H//dh))*100                                 # example: 100 mm bore diameter all along
# current_bore_diameter = pd.read_csv("data/bore_diameter.csv").values[:,-1]    # load from csv, dh step size must match, must be 1D
current_fluid_height = 0                                                        # enter current fluid height, fluid_volume will be calculated
filename = "Detailed Scenario Planning.xlsx"            # <-- file name here
sheet_name = "04-07-25-Scen3-Good,24-16hr"              # <-- sheet name here

# === Flow model parameters ===
p_exp = 3                          # LLiboutry
rhoi = 909                         # ice density (kg/m^3)
rhos = 404                         # surface density (kg/m^3)
Lrho = 33                          # half depth decay (m) - assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho)
g = 9.81                           # gravitational acceleration (m/s²)
ke = 1                             # enhancement coefficient for borehole closure
fluid_density = 793                # density of drilling fluid (kg/m³)

# === Drilling parameters ===
refill_limit = 90
refill_level = 80
leak_threshold = 40
drill_type = 'shallow'
drilling_speed = 0.0002                         # meters per second. This is not active if loading data.

# === Plotting settings ===
plot_spacing = 20                               # plot bore profile/closure rate every X timesteps

merged_plot = True                              # Merges the three following plots into one figure
plot_bore_profile_over_time = False   
plot_borehole_closure_rate_over_time = False  
plot_pressure_profile_over_time = False

plot_mim_bore_diameter_over_time = True
plot_bore_status_over_time = True
plot_temperature_density_profile = False
