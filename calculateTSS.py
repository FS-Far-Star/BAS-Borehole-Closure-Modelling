# James Veale 2024, from Carlos Martin Garcia's Matlab

import numpy as np
import math
import scipy
import warnings

def calculateTSS(z, acc, melt, eta, rho, rhoi, Tm, Ts, QG):
    """
    tss(z, w, Ts, kc, kd, QG)
    
    Solves:
        w dTdz - kd d2T/dz2 = 0
        T(z(end)) = Ts
        dt(z(1))/dz = -QG/kc

    by trapezoidal integration
    T(z) = Ts - int{u dz}_z^H
    where u(z) = -QG/kc*exp(int{w/kd dz})
    """

    numiter = 100
    Tss0 = np.zeros(np.size(z))

    for i in range(numiter):
        print("iter=",i)
        kc = 2 * (9.828*np.power(math.e, -0.0057*(273.15+Tss0))) * rho/rhoi / (3-rho/rhoi)
        Cp = (152.5 + 7.122*(273.15 + Tss0)) / 31557600

        kc[z<0] = (9.828 * np.power(math.e, -0.0057 * (273.15)))
        Cp[z<0] = (152.5 + 7.122 * (273.15)) / 31557600

        kd = kc / (rho * Cp)    # Diffusivity (m^2/yr)

        w = -(acc-melt)*eta - melt

        #I1 = int{w/kd dz}_0^z
        I1 = scipy.integrate.cumulative_trapezoid(w/kd, z, initial=0)
        u = -QG / kc * np.power(math.e, I1)

        #I2 = int{u dz}_z^H
        I2 = np.zeros(np.size(z))
        I2 = scipy.integrate.cumulative_trapezoid(np.flip(u), np.flip(z), initial=0)

        #T = Ts - int{u dz}_z^H
        Tss = Ts + np.flip(I2)

        if sum(abs(Tss0-Tss)) < pow(10, -6):
            print("Converged at i = ", i)
            #break
            return Tss
        else:
            Tss0 = Tss
        
    if i == numiter:
        warnings.warn("CalculateSS did not converge!", UserWarning)
