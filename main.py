# James Veale 2024, adapted from Carlos Martin Garcia's Matlab
# Borehole Temperature, Pressure, and Closure

import numpy as np
import math
from calculateTSS import calculateTSS
from calculateexponentialdensity import calculateexponentialdensity
from createlliboutryshapefunction import createlliboutryshapefunction
from calculateicepressure import calculateicepressure
from calculateboreholeclosure import calculateboreholeclosure
from calculatefluidpressure import calculatefluidpressure
import matplotlib.pyplot as plt

# Temperatures
#Carlos' values
"""# glaciological
H = 700   # ice thickness (m)
acc = 0.37/0.910    #surface accumulation (m/yr)
Ts = -22.3  # Surface temp (°C)
QG = 50 * pow(10,-3)    # geothermal flux (W/m^2)
melt = 0.0  # meltrate (m/yr)

# Cold scenario
HCold = 600
accCold = 0.69 / 0.910
TsCold = -23.3
QGCold = 40 * pow(10,-3)

# Warm Scenario
HWarm = 800
accWarm = 0.25/0.910
TsWarm = -21.3
QGWarm = 60 * pow(10,-3)"""

# Test values
H = 800   # ice thickness (m)
acc = 0.37/0.910    #surface accumulation (m/yr)
Ts = -22.3  # Surface temp (°C)
QG = 50 * pow(10,-3)    # geothermal flux (W/m^2)
melt = 0.0  # meltrate (m/yr)

# Flow model parameters
p_ice = 3 #LLiboutry
# Density (an 'ok' fit Carlos did to "Palmer density.xlsx")
rhoi = 909  # ice density (kg/m^3)
rhos = 404  # surface density (kg/m^3)
Lrho = 33   # half depth decay (m) - assumes rho=rhoi+(rhos-rhoi)*exp(-(H-z)/Lrho)

# Geometry
zinterval = 1   # resolution in m
nz = int(2*(H/zinterval) + 1)
print(nz)
#nz = 101
#nz = 1401
z = np.linspace(-H, H, nz)
#zCold = np.linspace(-HCold, HCold, nz)
#zWarm = np.linspace(-HWarm, HWarm, nz)

## Flow Model

# Exponential density
#rho = rhoi + (rhos-rhoi) * np.power(math.e,-(H-z)/Lrho)
rho = calculateexponentialdensity(z,rhoi,rhos,H,Lrho)

#Lliboutry shape function
#eta = 1 - (p+2) / (p+1) * (1- z/H) + (1- z/H)**(p+2)/(p+1)
#eta = eta*rhoi/rho
#eta[z<=0] = 0
eta = createlliboutryshapefunction(z, p_ice, H, rhoi, rho)

# SS Temperature
T = calculateTSS(z, acc, 0.0, eta, rho, rhoi, 0.0, Ts, QG)
#TCold = calculateTSS(zCold, acc, 0.0, eta, rho, rhoi, 0.0, TsCold, QGCold)
#TWarm = calculateTSS(zWarm, acc, 0.0, eta, rho, rhoi, 0.0, TsWarm, QGWarm)


zT = np.column_stack((z[z>=0],np.flip(T[z>=0])))
print("zT",zT)

g = 9.83 # gravitational acceleration at poles (m/s^2)
fluidlevel = 70 # close-off expected at ~65m at this site
densityexxsold60 = 793  # kg/m^3 # TODO: vary at depths with temperature

p_ice = calculateicepressure(z,rho,g)
p_fluid = calculatefluidpressure(z[z>=0], fluidlevel, densityexxsold60, g)
p_total = p_ice[z>=0] - p_fluid # +ve pressure closes borehole
print('p_t',p_total)

# Borehole
D0_s = 143    # initial diameter shallow drilling (mm)
D0_d = 132    # initial diameter deep drilling (mm)
ke = 1  # enhancement coefficient
t = 365 # time over which borehole left (days)
Dcore_d = 98    # core diameter deep drill (mm)
Dbarrel_d = 115 # outerbarrel deep drill (mm)

D_s_nofluid = calculateboreholeclosure(z[z>=0], D0_s, ke, np.flip(T[z>=0]), p_ice[z>=0], t)
print(np.column_stack((z[z>=0],D0_s+D_s_nofluid)))

D_d_nofluid = calculateboreholeclosure(z[z>=0], D0_d, ke, np.flip(T[z>=0]), p_ice[z>=0], t)

D_s = calculateboreholeclosure(z[z>=0], D0_s, ke, np.flip(T[z>=0]), p_total, t)
print(np.column_stack((z[z>=0],D0_s+D_s)))
D_d = calculateboreholeclosure(z[z>=0], D0_d, ke, np.flip(T[z>=0]), p_total, t)

Dnew_s_nofluid = D0_s + D_s_nofluid
Dnew_d_nofluid = D0_d + D_d_nofluid
Dnew_s = D0_s + D_s
Dnew_d = D0_d + D_d

depthswitch = 90 # m

Dnew_resultant_nofluid = np.copy(Dnew_s_nofluid)
Dnew_resultant_nofluid[z[z>=0]>depthswitch] = Dnew_d_nofluid[z[z>=0]>depthswitch]

Dnew_resultant = np.copy(Dnew_s)
Dnew_resultant[z[z>=0]>depthswitch] = Dnew_d[z[z>=0]>depthswitch]
#Dnew_s[z[z>=0]>depthswitch] = Dnew_d[z[z>=0]>depthswitch]

fig, ax = plt.subplots(2,2)
ax[0,0].plot(z[z>=0],np.roll(np.flip(rho), math.floor(nz/2))[z>=0])
ax[0,0].set_title("rho[z>=0]")
ax[0,0].set_ylabel("Density [kg.m^-3]")
ax[0,1].plot(z[z>=0],np.flip(T[z>=0]))
ax[0,1].set_title("np.flip(T[z>=0])")
ax[0,1].set_ylabel("Temperature [°C]]")
ax[1,0].plot(z[z>=0],p_ice[z>=0])
ax[1,0].plot(z[z>=0], p_fluid)
ax[1,0].plot(z[z>=0], p_total)
ax[1,0].set_title("p[z>=0]")
ax[1,0].set_ylabel("Pressure [Pa]")
ax[1,1].plot(z[z>=0],Dnew_s_nofluid)
ax[1,1].plot(z[z>=0],Dnew_d_nofluid)
ax[1,1].plot(z[z>=0],Dnew_s)
ax[1,1].plot(z[z>=0],Dnew_d)
#ax[1,1].plot(z[z>=0],Dnew_s)
ax[1,1].set_title("D0 + D")
ax[1,1].set_ylabel("New Borehole Diameter after " + str(t) + " days, mm")
ax[1,1].plot(z[z>=0], np.full(np.size(z[z>=0]),D0_d))
ax[1,1].plot(z[z>=0], np.full(np.size(z[z>=0]),D0_s))

fig2, ax2 = plt.subplots(1,1)
ax2.set_title("Borehole Diameters")
ax2.set_ylabel("New Borehole Diameter after " + str(t) + " days, mm")
ax2.set_xlabel("Depth from surface, m")
ax2.plot(z[z>=0], np.full(np.size(z[z>=0]),D0_s), linestyle="dashed", label="Shallow drill, initial")
ax2.plot(z[z>=0], np.full(np.size(z[z>=0]),D0_d), linestyle="dashed", label="Deep drill, initial")
ax2.plot(z[z>=0], np.full(np.size(z[z>=0]),Dbarrel_d), linestyle="dashed", label="Deep drill, barrel diameter")
ax2.plot(z[z>=0], np.full(np.size(z[z>=0]),Dcore_d), linestyle="dashed", label="Deep drill, core diameter")
ax2.plot(z[z>=0],Dnew_s_nofluid, linestyle="dashed", label="Shallow drill, no fluid")
ax2.plot(z[z>=0],Dnew_d_nofluid, linestyle="dashed", label="Deep drill, no fluid")
ax2.plot(z[z>=0],Dnew_s, linestyle="dashed", label="Shallow drill, fluid to " + str(fluidlevel) + "m")
ax2.plot(z[z>=0],Dnew_d, linestyle="dashed", label="Deep drill, fluid to " + str(fluidlevel) + "m")
ax2.plot(z[z>=0],Dnew_resultant_nofluid, label="borehole, no fluid, switched to deep drill at " + str(depthswitch) + "m")
ax2.plot(z[z>=0],Dnew_resultant, label="borehole, fluid to " + str(fluidlevel) + "m, switched to deep drill at " + str(depthswitch) + "m")
ax2.legend()
ax2.set_yticks(np.arange(0, 150, 5))
ax2.set_xticks(np.arange(0, 800, 25))
plt.grid()
plt.show()