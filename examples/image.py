"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to use the high level imageToVTK function.

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
from VTKwrite.interface import imageToVTK
import numpy as np

FILE_PATH = "./image"
def clean():
    try:
        os.remove(FILE_PATH + ".vti")
    except:
        pass
        
def run():
    print("Running image...")

    # Dimensions
    nx, ny, nz = 6, 6, 2
    ncells = nx * ny * nz
    npoints = (nx + 1) * (ny + 1) * (nz + 1)

    # Variables
    pressure = np.random.rand(ncells).reshape( (nx, ny, nz), order = 'C')
    temp = np.random.rand(npoints).reshape( (nx + 1, ny + 1, nz + 1))

    # data defined on cells
    all_cell_data = []
    all_cell_data.append(["pressure", "scalars", pressure])

    # data defined at points
    all_point_data = []
    all_point_data.append(["temp", "scalars", temp])

    comments = [ "comment 1", "comment 2" ]
    imageToVTK(FILE_PATH, all_cell_data = all_cell_data, all_point_data = all_point_data, comments = comments )

if __name__ == "__main__":
    run()

