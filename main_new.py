from init import *
from functions import *

# Ice properties
rho = density(h,rhoi,rhos,H,Lrho)
eta = shape_function(h, p_exp, H, rhoi, rho)
p_ice = ice_pressure(h, rho, g)
print('p_ice', p_ice)
temp = TSS(h, acc, melt, eta, rho, rhoi, Ts, QG)  #ignoring drilling fluid contribution?
print('temp', temp)

for i in range(1,3):
    print('Timestep:', i)
    # update bore status with drilling status
    bore_diameter[:,i] = bore_diameter[:,i-1]  # copy previous bore diameter
    drilling, drill_type = check_status(t[i])
    if drilling == True:
        last_depth = bore_depth[i-1]//dh
        bore_depth[i] = bore_depth[i-1] + drilling_speed * dt
    else:
        bore_depth[i] = bore_depth[i-1]

    print('bore depth',bore_depth[i])
    h_index = bore_depth[i]//dh

    # new drilling in this time step
    if h_index > last_depth:
        # print(last_depth, h_index)
        bore_diameter[last_depth:h_index,i] = drill_diameter(drill_type)
        print(bore_diameter[:,i])

    # update pressure distribution
    p_fluid, fluid_volume[i] = fluid_pressure(h, fluid_volume[i-1], bore_depth[i], bore_diameter[:,i], drill_diameter(drill_type), drilling, fluid_density=993, g=9.81)
    print('fluid_volume', fluid_volume[i])
    p_t = p_ice - p_fluid   # positive value means bore is closing

    print("p_t", p_t)
    print('before',bore_diameter[:,i])
    print('change', delta_d_bore(bore_diameter[:,i], ke, temp, p_t, dt))
    bore_diameter[:,i] = bore_diameter[:,i] + delta_d_bore(bore_diameter[:,i], ke, temp, p_t, dt)
    print('after',bore_diameter[:,i])
