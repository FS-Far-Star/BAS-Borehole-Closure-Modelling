# Borehole Closure Modeling Tool

This tool simulates borehole closure over time due to pressure and temperature in glacial ice, with drilling progression and fluid refill mechanics. It is based on time-stepped modeling and incorporates empirical drill performance, ice physics, and fluid pressure balance.



## Project Structure
```
BAS-Borehole-Closure-Modelling/
├── run.py                          # Top-level script to run the entire simulation
├── plotting.py                     # All visualization logic
├── global_settings.py              # Simulation setting
├── functions.py                    # Core computation functions (physics, drilling, fluid)
├── Detailed Scenario Planning.xlsx # Planning spreadsheet that drives the simulation
├── data/                           # Saved output from each simulation run
│ ├── bore_diameter.csv             # Bore diameter
│ ├── delta_bore.csv                # Bore diameter change per timestep
│ ├── temperature.csv               # Extrapolated ice temperature
│ ├── p_t.csv                       # Total pressure as a function of time and depth
│ ├── p_fluid.csv                   # Fluid pressure as a function of time and depth
│ ├── p_ice.csv                     # Total pressure as a function of depth
│ ├── rho.csv                       # Extrapolated ice density
│ ├── timeseries.csv                # Other time-dependent variables
│ └── depth.csv                     # Depth series used in modelling
└── README.md                       # This file
```

## User instructions
Unless you are attempting to develop the code further, calculation.py, plotting.py, and functions.py should NOT be altered. 
1. Spreadsheet driven simulation (default setting)
    - In ```global_settings.py```
        - specify spreadsheet file name and sheet name. 
        - Set loading = True
        - Modify ice settings to match drill site 
        - Modify simulation parameters like time and depth step sizes
        - Modify plot settings to decide which plots you want in output
    - Execute ```run.py```

2. Run bore closure simulation
    - In ```global_settings.py```
        - Set ```loading = False```
        - Modify ice settings to match drill site 
        - Modify simulation parameters including ```max_time```, ```dt```, ```dh```
        - Specify ```current_bore_diameter``` (import from csv or specify otherwise)
        - Specify ```current_bore_depth``` and ```current_fluid_height```
        - Set function output of ```check_status(time)``` to ```False```, this disables drilling
        - Modify plot settings to decide which plots you want in output
    - Execute ```run.py```

3. Run drilling and bore closure simulation
    - In ```global_settings.py```
        - Set ```loading = False```
        - Modify ice settings to match drill site 
        - Modify simulation parameters including ```max_time```, ```dt```, ```dh```
        - Set ```current bore diameter``` to ```np.zeros(int(H//dh))```
        - Set ```current_bore_depth``` and ```current_fluid_height``` to 0
        - Set function output of ```check_status(time)``` to ```drilling``` and specify duration of drilling within the function
        - Modify plot settings to decide which plots you want in output
    - Execute ```run.py```
4. Plot the data of the last simulation
    - In ```run.py```, comment out ```run_calculation()```, and execute ```run.py```

The default unit of time in this simulation is seconds. It's recommended to use ```dt``` of at least 1 day (86400 seconds). The default unit of length is meters, except for diameters where millimeters are used. ```h``` is postive downward (i.e. depth from surface). 