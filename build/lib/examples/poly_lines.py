"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to use the high level polyLinesToVTK function.

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
from VTKwrite.interface import polyLinesToVTK
import numpy as np

FILE_PATH = "poly_lines"
def clean():
    try:
        os.remove(FILE_PATH + ".vtu")
    except:
        pass
        
def run():
    print("Running poly_lines...")
    # Positions of points that define the lines
    npoints = 7
    x = np.zeros(npoints)
    y = np.zeros(npoints)
    z = np.zeros(npoints)

    # First line
    x[0], y[0], z[0] = 0.0, 0.0, 0.0
    x[1], y[1], z[1] = 1.0, 1.0, 0.0
    x[2], y[2], z[2] = 2.0, 0.0, 0.0
    x[3], y[3], z[3] = 3.0, -1.0, 0.0

    # Second line
    x[4], y[4], z[4] = 0.0, 0.0, 3.0
    x[5], y[5], z[5] = 1.0, 1.0, 3.0
    x[6], y[6], z[6] = 2.0, 0.0, 3.0

    # Connectivity of the lines
    pointsPerLine = np.zeros(2)
    pointsPerLine[0] = 4
    pointsPerLine[1] = 3
    ncells = 2

    # Some variables
    pressure = np.random.rand(npoints)
    temp = np.random.rand(npoints)
    vel = np.zeros(ncells*3)
    vel[0:3]  = [1.0, -1.0, 2.0]
    vel[3:6]  = [3.0, 4.0, -5.0]

    # initialize the data structure
    all_cell_data = {"scalars" : None, "vectors" : None}
    # vectors
    cellData_vc = {"vel" : vel}
    all_cell_data["vectors"] = cellData_vc
    
    # initialize the data structure
    all_point_data = {"scalars" : None, "vectors" : None}
    # scalars
    pointData_sc = {"temp" : temp}
    pointData_sc["pressure"] = pressure
    all_point_data["scalars"] = pointData_sc

    polyLinesToVTK(FILE_PATH, x, y, z, pointsPerLine = pointsPerLine, all_cell_data = all_cell_data, all_point_data = all_point_data)

if __name__ == "__main__":
    run()
