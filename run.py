# Alex Jiang 04/08/2025
# Transient Borehole Geometry Modeling
# Based on James Veale's drilling planning spreadsheet and previous work

from calculation import run_calculation
from plotting import plot_graphs
import os

# Clear the console
os.system('cls' if os.name == 'nt' else 'clear')

run_calculation()
plot_graphs()