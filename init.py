import numpy as np
import scipy
import matplotlib.pyplot as plt


# Simulation parameters
max_time = 60*60*24*365  # seconds
max_depth = 2000  # meters
dt = 100  # time step in seconds
dh = 1  # depth step in meters

t = np.arange(0 , max_time + dt , dt)
h = np.arange(0 , max_depth + dh , dh)

fluid_volume = np.zeros_like(t)       # fluid volume over time
bore_depth = np.zeros_like(t)         # borehole depth over time

bore_diameter = np.zeros((len(h), len(t)))  # bore diameter as a function of depth and time

print(t.shape, bore_diameter.shape)