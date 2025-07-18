def createlliboutryshapefunction(z, p, H, rhoi, rho):
    eta = 1 - (p+2) / (p+1) * (1- z/H) + (1- z/H)**(p+2)/(p+1)
    eta = eta*rhoi/rho
    eta[z<=0] = 0
    return eta