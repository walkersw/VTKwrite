"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to use the high level linesToVTK function.

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
from VTKwrite.interface import linesToVTK
import numpy as np

FILE_PATH = "./lines"
def clean():
    try:
        os.remove(FILE_PATH + ".vtu")
    except:
        pass
        
def run():
    print("Running lines...")
    # Positions of points that define the lines
    npoints = 4
    x = np.zeros(npoints)
    y = np.zeros(npoints)
    z = np.zeros(npoints)
    temp = np.random.rand(npoints)
    flux = np.random.rand(npoints*3)
    
    ncells = 2
    pressure = np.random.rand(ncells)
    vel = np.zeros(ncells*3)
    vel[0] = 1.0
    vel[1] = 5.0
    vel[2] = -3.0
    vel[3] = 7.0
    vel[4] = -6.0
    vel[5] = 1.0

    x[0], y[0], z[0] = 0.0, 0.0, 0.0
    x[1], y[1], z[1] = 1.0, 1.0, 1.0
    x[2], y[2], z[2] = 0.0, 0.0, 0.0
    x[3], y[3], z[3] = -1.0, 1.0, 1.0
    
    # data defined on cells
    all_cell_data = []
    all_cell_data.append(["pressure", "scalars", pressure])
    all_cell_data.append(["velocity", "vectors", vel])

    # data defined at points
    all_point_data = []
    all_point_data.append(["temp", "scalars", temp])
    all_point_data.append(["flux", "vectors", flux])
    
    comments = [ "comment 1", "comment 2" ]
    linesToVTK(FILE_PATH, x, y, z, all_cell_data = all_cell_data, all_point_data = all_point_data, comments = comments)

if __name__ == "__main__":
    run()

