import numpy as np
import scipy
import matplotlib.pyplot as plt
import math
import warnings


# Simulation parameters
max_time = 60*60*24*365  # seconds
max_depth = 2000  # meters
dt = 24*60*60  # time step in seconds
dh = 1  # depth step in meters

t = np.arange(0 , max_time + dt , dt)
h = np.arange(0 , max_depth + dh , dh)

fluid_volume = np.zeros_like(t)       # fluid volume over time
bore_depth = np.zeros_like(t)         # borehole depth over time

bore_diameter = np.zeros((len(h), len(t)))  # bore diameter as a function of depth and time
temp = np.zeros((len(h), len(t)))           # temperature as a function of depth and time

# print(t.shape, bore_diameter.shape)

# Test values
H = 800   # ice thickness (m)
acc = 0.37/0.910    #surface accumulation (m/yr)
Ts = -22.3  # Surface temp (°C)
QG = 50 * pow(10,-3)    # geothermal flux (W/m^2)
melt = 0.0  # meltrate (m/yr)

# Flow model parameters
p_exp = 3 #LLiboutry
# Density (an 'ok' fit Carlos did to "Palmer density.xlsx")
rhoi = 909  # ice density (kg/m^3)
rhos = 404  # surface density (kg/m^3)
Lrho = 33   # half depth decay (m) - assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho)
g = 9.81  # gravitational acceleration (m/s²)
ke = 1 # enhancement coefficient for borehole closure

# drill parameters
drilling_speed = 0.001  # meters per second
flow_rate = 1e-3  # cubic meters per second

def drill_diameter(drill_type):
    if drill_type == 'shallow':
        return 0.1  # shallow bore diameter in meters
    elif drill_type == 'deep':
        return 0.2  # deep bore diameter in meters
    else:
        raise ValueError("Unknown drill type: {}".format(drill_type))