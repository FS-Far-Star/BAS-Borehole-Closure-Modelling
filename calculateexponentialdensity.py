import numpy as np
import math

def calculateexponentialdensity(z, rhoi, rhos, H, Lrho):
    return rhoi + (rhos-rhoi)* np.power(math.e, -(H-z)/Lrho)
