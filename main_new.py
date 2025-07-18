from init import *
from functions import *
drilling = False
adding_fluid = False

for i in range(1,len(t)+1):
    # update bore status with drilling status
    drilling, adding_fluid = check_status(i)
    if drilling == True:
        pre_drill = bore_depth[i-1]//dh
        bore_depth[i] = bore_depth[i-1] + drilling_speed * dt
    else:
        bore_depth[i] = bore_depth[i-1]

    h_index = bore_depth[i]//dh

    if h_index > pre_drill:
        bore_diameter[pre_drill:h_index,i] = shallow_drill_diameter
    
    if adding_fluid == True:
        fluid_volume[i] = fluid_volume[i-1] + flow_rate * dt
    else:
        fluid_volume[i] = fluid_volume[i-1]

    p_ice = ice_pressure()
    p_fluid = fluid_pressure()
    p_t = p_ice - p_fluid   # positive value means bore is closing
    bore_diameter[0:h_index,i] = bore_diameter[0:h_index,i-1] + delta_d_bore(fluid_volume[i], bore_depth[i], bore_diameter[0:h_index,i-1])
    


