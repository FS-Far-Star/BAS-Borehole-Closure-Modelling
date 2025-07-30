from init import *

import numpy as np
import scipy.integrate
import warnings

def TSS(h, acc, melt, eta, rho, rhoi, Ts, QG):
    """
    Computes steady-state temperature profile in glacier ice.

    Solves:
        w dT/dh - kd d²T/dh² = 0
        T(h=0) = Ts (surface)
        dT/dh at base = -QG/kc

    Parameters:
        h    : np.ndarray, depth array [m], increasing positively downward
        acc  : float, surface accumulation rate [m/yr]
        melt : float, basal melt rate [m/yr]
        eta  : np.ndarray, Lliboutry shape function (same length as h)
        rho  : np.ndarray, density profile at h [kg/m³]
        rhoi : float, reference ice density [kg/m³]
        Ts   : float, surface temperature [°C]
        QG   : float, geothermal heat flux [W/m²]

    Returns:
        Tss  : np.ndarray, steady-state temperature profile [°C]
    """
    numiter = 100
    Tss0 = np.zeros_like(h)

    for i in range(numiter):
        kc = 2 * (9.828 * np.exp(-0.0057 * (273.15 + Tss0))) * rho / rhoi / (3 - rho / rhoi)
        Cp = (152.5 + 7.122 * (273.15 + Tss0)) / 31557600  # convert to per second

        kd = kc / (rho * Cp)  # thermal diffusivity [m²/yr]

        w = -(acc - melt) * eta - melt  # vertical advection velocity [m/yr]

        # I1 = ∫ (w / kd) dh from surface to h
        I1 = scipy.integrate.cumulative_trapezoid(w / kd, h, initial=0)
        u = -QG / kc * np.exp(I1)

        # I2 = ∫ u dh from h to surface (since h increases downward)
        I2 = scipy.integrate.cumulative_trapezoid(u, h, initial=0)

        # T(h) = Ts + ∫ u dh from h to surface
        Tss = Ts + I2

        if np.sum(np.abs(Tss0 - Tss)) < 1e-6:
            return Tss
        else:
            Tss0 = Tss

    warnings.warn("TSS did not converge after {} iterations".format(numiter), UserWarning)
    return Tss0


def delta_d_bore(D0, ke, T, p, t):
    """Calculates change in borehole diameter in time t (days) for range of depths, z (m) given
        initialdiameter, D0 (mm)
        enhancement coefficient, ke
        ice temperatures at z, T (C)
        overburden pressure at depth at z, p (Pa)

        Using Talalay, deltaD = D0(1-exp[(6*10^-21)*ke*(e^(0.12T))*(p^3)*t])

        N.B. When D0, T, p are given as arrays, each value must be at the same depths as z as the desired output
        """
    delta = np.zeros_like(D0)

    for i in range(len(D0)):
        if D0[i] > 0:
            eq = ke * math.exp(0.12 * T[i]) * t * (p[i] ** 3) * 6e-21
            d = D0[i] * (1 - math.exp(eq))

            # Clamp to physical limit
            delta[i] = max(d, -D0[i])
        else:
            delta[i] = 0  # no borehole → no change

    return delta

def fluid_pressure(h, fluid_volume, bore_depth, bore_diameter, drill_diameter, drilling, fluid_density=993, g=9.81):
    """
    Calculate hydrostatic fluid pressure profile in a borehole.
    
    Parameters:
        h (np.ndarray): 1D array of vertical depth grid points [m], increasing positively downward
        fluid_volume (float): Total fluid volume in borehole [m³]
        bore_depth (float): Current depth of borehole [m]
        bore_diameter (np.ndarray): Borehole diameter at each depth in h [m]
        drill_diameter (float): Diameter at current drill tip [m]
        fluid_density (float): Density of borehole fluid [kg/m³]
        g (float): Gravitational acceleration [m/s²]
    
    Returns:
        np.ndarray: Fluid pressure profile [Pa], same shape as h
    """
    
    p_fluid = np.zeros_like(h)
    dh = abs(h[1] - h[0])  # assumed uniform spacing

    # Find bottom index i such that h[i] < bore_depth <= h[i+1]
    bore_index = np.searchsorted(h, bore_depth, side='right') - 1

    if bore_index < 0:
        return p_fluid  # borehole hasn't reached first grid point yet

    # Volume of partial cell between h[i] and bore_depth
    dz_partial = bore_depth - h[bore_index]
    D_partial = 0.5 * (bore_diameter[bore_index] + drill_diameter)
    V_partial = (np.pi / 4) * D_partial**2 * dz_partial

    # Remaining fluid volume to distribute upward
    remaining_volume = fluid_volume - V_partial

    # handling start of drilling
    if remaining_volume < 0:
        remaining_volume = 0

    # Fill upward until volume runs out
    fluid_top_index = 0
    dz_fill = 0

    for i in range(bore_index - 1, -1, -1):
        dv = (np.pi / 4) * bore_diameter[i]**2 * dh
        if remaining_volume >= dv:
            remaining_volume -= dv
        else:
            fluid_top_index = i
            if fluid_top_index > 1:
                dz_fill = remaining_volume / ((np.pi / 4) * ((bore_diameter[i]+bore_diameter[i-1])/2)**2)
            else:
                dz_fill = 0
                print("Warning: Fluid volume exceeds borehole capacity")
            break

    # Compute true fluid surface height
    fluid_height = h[fluid_top_index] - dz_fill

    # Refill condition: fluid dropped too deep during drilling
    if fluid_height > 90 and drilling:
        refill_level = 80
        refill_index = np.searchsorted(h, refill_level, side='left')  # h[refill_index] >= 80
        print('Refill')

        # Top partial cell (from 80 m to next grid point)
        if h[refill_index] > refill_level and refill_index > 0:
            dz_top_partial = h[refill_index] - refill_level
            V_top_partial = (np.pi / 4) * ((bore_diameter[refill_index] + bore_diameter[refill_index - 1])/2)**2 * dz_top_partial
        else:
            V_top_partial = 0

        # Fully filled cells from refill_index down to i
        V_full = 0
        for i in range(refill_index, bore_index):
            dv = (np.pi / 4) * bore_diameter[i]**2 * dh
            V_full += dv

        # Total fluid volume including partial at bottom
        fluid_volume = V_partial + V_full + V_top_partial

        # Reset fluid height to refill level
        fluid_height = refill_level
    else:
        pass    # no refill

    # Compute hydrostatic pressure from fluid_height down to bore_depth
    for i in range(fluid_top_index, bore_index + 1):
        if h[i] >= fluid_height:
            p_fluid[i] = fluid_density * g * (h[i] - fluid_height)  # will be negative unless corrected

    return p_fluid,fluid_volume

def ice_pressure(h, rho, g=9.81):
    """
    Calculates the overburden pressure in ice from surface to depth.

    Parameters:
        h   : array-like, depth [m], increasing positively downward
        rho : array-like, ice density at depth [kg/m³]
        g   : float, gravitational acceleration [m/s²]

    Returns:
        p   : array-like, overburden pressure at each depth [Pa]
    """
    p = scipy.integrate.cumulative_trapezoid(rho * g, h, initial=0)
    return p


def shape_function(h, p_exp, H, rhoi, rho):
    """
    Computes the Lliboutry shape function η(h) over depth h,
    assuming h = 0 at surface and h = H at bed.

    Parameters:
        h     : array-like, depth [m], increasing positively downward
        p_exp : Glen's law exponent (usually 3)
        H     : ice thickness [m]
        rhoi  : reference ice density [kg/m³]
        rho   : array-like, depth-dependent ice density [kg/m³]

    Returns:
        eta   : shape function values at each depth
    """
    s = 1 - h / H  # normalized vertical coordinate, 1 at surface, 0 at bed

    eta = 1 - (p_exp + 2) / (p_exp + 1) * (1 - s) + s**(p_exp + 2) / (p_exp + 1)
    eta *= (rhoi / rho)

    # Zero out below bed if needed (e.g., if h exceeds H due to padding)
    eta = np.where(h > H, 0, eta)

    return eta

def density(h, rhoi, rhos, H, Lrho):
    """
    Returns vertical density profile using an exponential fit.

    Parameters:
        h     : array-like, depth [m], increasing positively downward
        rhoi  : float, ice density [kg/m³] at depth
        rhos  : float, surface density [kg/m³]
        H     : float, total ice thickness [m]
        Lrho  : float, e-folding depth scale [m]

    Returns:
        rho   : array-like, density profile [kg/m³]
    """
    return rhoi + (rhos - rhoi) * np.exp(-(H - h) / Lrho)   

def check_status(time):
    # Placeholder for drilling strategy
    drilling, drill_type = True, 'shallow'
    return [drilling, drill_type]