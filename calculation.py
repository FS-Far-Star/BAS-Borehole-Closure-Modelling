from global_settings import *
from functions import *

def run_calculation():
    "Main calculation function"
    global dt,drilling_speed
    if loading == True:
        dt = 24*60*60                                           # time step set to 1 day to match spreadsheet

    # Initialize 
    t = np.arange(0 , max_time + dt , dt)                       # time 
    h = np.arange(0 , H + dh , dh)                              # depth 
    bore_diameter  = np.zeros((len(h), len(t)))                 # bore diameter as a function of depth and time
    delta_bore  = np.zeros((len(h), len(t)))                    # bore diameter change every time step
    bore_depth = np.ones(len(t)) * current_bore_depth           # borehole depth over time

    if loading == False:
        # --- initialize arrays ---
        fluid_height   = np.zeros(len(t))                       # height of fluid in borehole over time
        fluid_volume   = np.zeros(len(t))                       # fluid volume over time
    elif loading == True:
        # --- Load data, skip first 2 rows ---
        df = pd.read_excel(filename, sheet_name=sheet_name, skiprows=2)
        df = df[df.iloc[:, 2].notna()]
        last_valid_index = df[df.iloc[:, 2] == 'Camp Closed'].last_valid_index()
        df = df.loc[:last_valid_index]
        working_status = df.iloc[:, 2].to_numpy()               # Can also use column E instead if 'Drilling' is to be changed to other stuff
        EOD_depth = np.nan_to_num(df.iloc[:, 8].to_numpy(), nan=0.0)    # nan to 0.0, to handle the empty cells
        # --- initialize arrays ---
        bore_depth[0:len(EOD_depth)] = EOD_depth                # borehole depth over time
        bore_diameter[0:len(current_bore_diameter),0] = current_bore_diameter  # load data
        fluid_height   = np.ones(len(t)) * current_fluid_height # height of fluid in borehole over time
        current_fluid_volume = calc_volume(bore_depth[0], bore_diameter[:,0], fluid_height[0], dh, h)
        fluid_volume   = np.ones(len(t)) * current_fluid_volume # fluid volume over time

    p_fluid    = np.zeros((len(h), len(t)))                     # fluid pressure 
    p_t        = np.zeros((len(h), len(t)))                     # total pressure

    # Ice properties
    rho  = density(h, rhoi, rhos, Lrho)
    eta  = shape_function(h, p_exp, H, rhoi, rho)
    p_ice = ice_pressure(h, rho, g)                             # hydrostatic pressure profile of ice is constant over time
    temp  = TSS(h, acc, melt, eta, rho, rhoi, Ts, QG)           # ignoring drilling fluid contribution and transient, ice temperature is constant over time

    # Iterating
    for i in range(1, len(t)):
        print('Timestep:', i)
        # update bore status with drilling status
        bore_diameter[:, i] = bore_diameter[:, i-1]             # copy previous bore diameter
        if loading == True:
            if i < len(working_status) and working_status[i] == 'Drilling':
                drilling = True
            else:
                drilling = False
        elif loading == False:
            drilling = True     # unless loading data, assume drilling is always True

        last_depth      = int(bore_depth[i-1] // dh)

        if loading == False:  
            drilling = check_status(t[i],max_time)
            if drilling == True:
                bore_depth[i]   = bore_depth[i-1] + drilling_speed * dt
                if bore_depth[i] > H:
                    bore_depth[i] = H                       # max depth reached
            else:
                bore_depth[i]   = bore_depth[i-1]
        elif loading == True:
            if drilling == True:
                pass
            elif drilling == False:
                bore_depth[i] = bore_depth[i-1]            # keep the same depth if not drilling

        print(f"bore depth {bore_depth[i]:.2f}m")
        print('')

        h_index = int(bore_depth[i] // dh)                 

        # set diameter of new drilling in this time step
        if h_index > last_depth:
            bore_diameter[last_depth:h_index, i] = drill_diameter(drill_type)

        # update pressure distribution
        p_fluid[:, i], fluid_volume[i], fluid_height[i] = fluid_pressure(h, fluid_volume[i-1], bore_depth[i], bore_diameter[:, i], drilling, fluid_density, g)
        p_t[:, i] = p_ice - p_fluid[:, i]                # Sign notation: positive pressure will close the borehole

        # apply borehole closure
        delta_bore[:, i]     = delta_d_bore(bore_diameter[:, i], ke, temp, p_t[:, i], dt)
        bore_diameter[:, i] += delta_bore[:, i]
        
    # Save 2D arrays: each column is a timestep
    pd.DataFrame(bore_diameter).to_csv("data/bore_diameter.csv", index=False)
    pd.DataFrame(delta_bore).to_csv("data/delta_bore.csv", index=False)
    pd.DataFrame(temp).to_csv("data/temperature.csv", index=False)
    pd.DataFrame(p_t).to_csv("data/p_t.csv", index=False)
    pd.DataFrame(p_fluid).to_csv("data/p_fluid.csv", index=False)
    pd.DataFrame(p_ice).to_csv("data/p_ice.csv", index=False)
    pd.DataFrame(rho).to_csv("data/rho.csv", index=False)

    # Save 1D time series
    pd.DataFrame({
        "time": t,
        "fluid_volume": fluid_volume,
        "fluid_height": fluid_height,
        "bore_depth": bore_depth
    }).to_csv("data/timeseries.csv", index=False)

    # Save spatial axis (e.g., depth)
    pd.DataFrame({"depth": h}).to_csv("data/depth.csv", index=False)
    return None

if __name__ == "__main__":
    run_calculation()
