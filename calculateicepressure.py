import numpy as np
import scipy
import math
#import matplotlib.pyplot as plt

def calculateicepressure(z, rho, g):
    """Takes rho ordered with z=0 at bottom of ice, and +ve z up, and outputs a density for z=0 at surface and +ve z down"""
    rhonew = np.copy(rho)
    shift = math.floor(np.size(z)/2)
    rhonew = np.roll(np.flip(rhonew), shift)
    rhonew[z<=0] = 0

    p = scipy.integrate.cumulative_trapezoid(rhonew * g, z, initial=0)

    """fig, ax1 = plt.subplots()
    ax1.set_xlim(xmin=0,xmax=np.max(z))
    ax1.set_xlabel('z (m)')
    ax1.set_ylabel('rho', color="tab:red")
    ax1.plot(z,rho, color = ("tab:red", 0.3))
    ax1.plot(z,rhonew, color = "tab:red")
    ax1.tick_params(axis='y', labelcolor="tab:red")
    ax2 = ax1.twinx()
    ax2.set_ylabel('p', color="tab:blue")
    ax2.plot(z,p, color = "tab:blue")
    ax2.tick_params(axis='y', labelcolor="tab:blue")"""

    return p