import numpy as np
import math

def calculateboreholeclosure(z, D0, ke, T, p, t):
    """Calculates change in borehole diameter in time t (days) for range of depths, z (m) given
        initialdiameter, D0 (mm)
        enhancement coefficient, ke
        ice temperatures at z, T (C)
        overburden pressure at depth at z, p (Pa)

        Using Talalay, deltaD = D0(1-exp[(6*10^-21)*ke*(e^(0.12T))*(p^3)*t])

        N.B. When D0, T, p are given as arrays, each value must be at the same depths as z as the desired output
        """
    
    eq1 = ke * np.power(math.e, 0.12*T) * t * np.power(p, 3) * 6 * pow(10, -21)
    D = D0 * (1 - np.power(math.e, eq1))

    D[D<-D0] = -D0

    return D
    
