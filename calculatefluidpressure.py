def calculatefluidpressure(z, fluiddepth, rhofluid, g):
    """Returns an array of fluid pressures at depths below ice surface given in array z, given:
        fluiddepth, metres below z the top of the fluid is,
        rhofluid, density of the drilling fluid,
        g, local gravitational constant"""
    
    # TODO: rhofluid vs temperature!
    
    #rho.g.h
    p = rhofluid * g * (z - fluiddepth)
    
    #fix p above fluid level to 0
    p[z<=fluiddepth] = 0

    return p
