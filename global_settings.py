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

# === Simulation parameters ===
max_time = 60*60*24*365             # seconds
max_depth = 2000                    # meters
dt = 24*60*60                       # time step in seconds
dh = 1                              # depth step in meters

# Load existing data?
loading = False                     # Set to True if loading existing data
current_bore_depth = 1500                                                 # enter current bore depth
current_bore_diameter = np.ones(int(max_depth//dh))*100                         # example: 100 mm bore diameter all along
# current_bore_diameter = pd.read_csv("data/bore_diameter.csv").values[:,-1]    # load from csv, dh step size must match, must be 1D
current_fluid_height = 80                                               # enter current fluid height, fluid_volume will be calculated

# === Ice settings ===
H = 800                            # ice thickness (m)
acc = 0.37/0.910                   # surface accumulation (m/yr)
Ts = -22.3                         # Surface temp (°C)
QG = 50 * pow(10,-3)               # geothermal flux (W/m^2)
melt = 0.0                         # meltrate (m/yr)

# === Flow model parameters ===
p_exp = 3                          # LLiboutry
rhoi = 909                         # ice density (kg/m^3)
rhos = 404                         # surface density (kg/m^3)
Lrho = 33                          # half depth decay (m) - assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho)
g = 9.81                           # gravitational acceleration (m/s²)
ke = 1                             # enhancement coefficient for borehole closure
fluid_density = 793                # density of drilling fluid (kg/m³)

# === Drilling parameters ===
drilling_speed = 0.0001            # meters per second

# === Plotting settings ===
plot_spacing = 20                               # plot bore profile/closure rate every X timesteps
plot_bore_profile_over_time = True                
plot_mim_bore_diameter_over_time = True
plot_borehole_closure_rate_over_time = True        
plot_bore_status_over_time = True
plot_temperature_density_profile = True
plot_pressure_profile_over_time = True
