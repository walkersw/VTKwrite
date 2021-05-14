"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to use the high level unstructuredGridToVTK function.
This example shows how to export an unstructured grid given its vertices and
topology through connectivity and offset lists.
Check the VTK file format for details of the unstructured grid.

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
from VTKwrite.interface import unstructuredGridToVTK
from VTKwrite.vtkbin import VtkTriangle, VtkQuad
import numpy as np

FILE_PATH = "./unstructured"
def clean():
    try:
        os.remove(FILE_PATH + ".vtu")
    except:
        pass
        
def run():
    print("running unstructured...")

    # Define vertices
    x = np.zeros(6)
    y = np.zeros(6)
    z = np.zeros(6)

    x[0], y[0], z[0] = 0.0, 0.0, 0.0
    x[1], y[1], z[1] = 1.0, 0.0, 0.0
    x[2], y[2], z[2] = 2.0, 0.0, 0.0
    x[3], y[3], z[3] = 0.0, 1.0, 0.0
    x[4], y[4], z[4] = 1.0, 1.0, 0.0
    x[5], y[5], z[5] = 2.0, 1.0, 0.0

    # Define connectivity or vertices that belongs to each element
    conn = np.zeros(10)

    conn[0], conn[1], conn[2] = 0, 1, 3              # first triangle
    conn[3], conn[4], conn[5] = 1, 4, 3              # second triangle
    conn[6], conn[7], conn[8], conn[9] = 1, 2, 5, 4  # rectangle

    # Define offset of last vertex of each element
    offset = np.zeros(3)
    offset[0] = 3
    offset[1] = 6
    offset[2] = 10

    # Define cell types

    ctype = np.zeros(3)
    ctype[0], ctype[1] = VtkTriangle.tid, VtkTriangle.tid
    ctype[2] = VtkQuad.tid
    
    # data defined on cells
    all_cell_data = []
    all_cell_data.append(["pressure0", "scalars", np.random.rand(3)])
    all_cell_data.append(["pressure1", "scalars", np.random.rand(3)])
    all_cell_data.append(["vel0", "vectors", np.random.rand(3*3)])
    all_cell_data.append(["vel1", "vectors", np.random.rand(3*3)])

    # data defined at points
    all_point_data = []
    all_point_data.append(["potential0", "scalars", np.random.rand(6)])
    all_point_data.append(["potential1", "scalars", np.random.rand(6)])
    all_point_data.append(["flux0", "vectors", np.random.rand(6*3)])
    all_point_data.append(["flux1", "vectors", np.random.rand(6*3)])
    
    comments = [ "comment 1", "comment 2" ]
    unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity = conn, offsets = offset, cell_types = ctype, all_cell_data = all_cell_data, all_point_data = all_point_data, comments = comments)

if __name__ == "__main__":
    run()
