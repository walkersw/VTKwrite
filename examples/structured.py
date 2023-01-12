"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to export a structured grid.

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
from VTKwrite.interface import structuredToVTK
import numpy as np
import random as rnd

FILE_PATH = "./structured"
def clean():
    try:
        os.remove(FILE_PATH + ".vts")
    except:
        pass

def run():
    print("Running structured...")

    # Dimensions
    nx, ny, nz = 6, 6, 2
    lx, ly, lz = 1.0, 1.0, 1.0
    dx, dy, dz = lx/nx, ly/ny, lz/nz

    ncells = nx * ny * nz
    npoints = (nx + 1) * (ny + 1) * (nz + 1)

    # Coordinates
    X = np.arange(0, lx + 0.1*dx, dx, dtype='float64')
    Y = np.arange(0, ly + 0.1*dy, dy, dtype='float64')
    Z = np.arange(0, lz + 0.1*dz, dz, dtype='float64')

    x = np.zeros((nx + 1, ny + 1, nz + 1))
    y = np.zeros((nx + 1, ny + 1, nz + 1))
    z = np.zeros((nx + 1, ny + 1, nz + 1))

    # We add some random fluctuation to make the grid
    # more interesting
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                x[i,j,k] = X[i] + (0.5 - rnd.random()) * 0.1 * dx 
                y[i,j,k] = Y[j] + (0.5 - rnd.random()) * 0.1 * dy
                z[i,j,k] = Z[k] + (0.5 - rnd.random()) * 0.1 * dz

    # Variables
    pressure = np.random.rand(ncells).reshape( (nx, ny, nz))
    temp = np.random.rand(npoints).reshape( (nx + 1, ny + 1, nz + 1))

    # data defined on cells
    all_cell_data = []
    all_cell_data.append(["pressure", "scalars", pressure])

    # data defined at points
    all_point_data = []
    all_point_data.append(["temp", "scalars", temp])

    comments = [ "comment 1", "comment 2" ]
    structuredToVTK(FILE_PATH, x, y, z, all_cell_data = all_cell_data, all_point_data = all_point_data, comments = comments)

if __name__ == "__main__":
    run()
