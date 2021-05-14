"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to use the high level unstructuredGridToVTK function.
This example shows how to export *only* an unstructured grid given its
vertices and topology through connectivity and offset lists.
It also shows how to store a static grid and a "time-sequence" of data
all in one .vtu file.  In addition, the time values are stored separately
using a "fake" grid.
Check the VTK file format for details of the unstructured grid.

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
from VTKwrite.interface import unstructuredGridToVTK
from VTKwrite.interface import timeseries_unstructuredGrid
from VTKwrite.vtkbin import VtkTriangle, VtkQuad
from VTKwrite.vtkbin import VtkGroup
import numpy as np

FILE_PATH_GRID = "./unstructured_timedep_grid_only"
FILE_PATH_ALL  = "./unstructured_timedep_all"
FILE_PATH_TIME = "./unstructured_timedep_time_values"
def clean():
    try:
        os.remove(FILE_PATH_GRID + ".vtu")
    except:
        pass
    try:
        os.remove(FILE_PATH_ALL + ".vtu")
    except:
        pass
    try:
        os.remove(FILE_PATH_TIME + ".vtu")
    except:
        pass

def run():
    print("running unstructured_timedep...")

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
    
    comments = [ "comment 1", "comment 2" ]
    # write grid only
    unstructuredGridToVTK(FILE_PATH_GRID, x, y, z, connectivity = conn, offsets = offset, cell_types = ctype, all_cell_data = None, all_point_data = None, comments = comments)

    # create time-dependent data
    p0_ref = np.zeros(3)
    p0_ref[0] = 0
    p0_ref[1] = 1
    p0_ref[2] = 2
    p1_ref = np.zeros(3)
    p1_ref[0] = -1
    p1_ref[1] = 4
    p1_ref[2] = -3
    v0_ref = np.zeros(3*3)
    v0_ref[0:3] = [1, 1, 1]
    v0_ref[3:6] = [0, -2, 3]
    v0_ref[6:9] = [-2, 0, 1]
    v1_ref = np.zeros(3*3)
    v1_ref[0:3] = [-1, 4, -1]
    v1_ref[3:6] = [1, 1, 6]
    v1_ref[6:9] = [2, -8, 3]
    
    pot0_ref = x + y + z
    pot1_ref = np.sin(2*x + y)
    f0_ref = np.zeros((3,6))
    f0_ref[0,:] = x
    f0_ref[1,:] = -z
    f0_ref[2,:] = 2*y
    f0_ref = f0_ref.flatten('F')
    f1_ref = np.zeros((3,6))
    f1_ref[0,:] = np.cos(x)
    f1_ref[1,:] = -y
    f1_ref[2,:] = 3*(x + z)
    f1_ref = f1_ref.flatten('F')
    
    p0 = p0_ref
    p1 = p1_ref
    v0 = v0_ref
    v1 = v1_ref
    pot0 = pot0_ref
    pot1 = pot1_ref
    f0 = f0_ref
    f1 = f1_ref
    
    # data defined on cells
    all_cell_data = []
    all_cell_data.append(["pressure0", "scalars", p0])
    all_cell_data.append(["pressure1", "scalars", p1])
    all_cell_data.append(["vel0", "vectors", v0])
    all_cell_data.append(["vel1", "vectors", v1])

    # data defined at points
    all_point_data = []
    all_point_data.append(["potential0", "scalars", pot0])
    all_point_data.append(["potential1", "scalars", pot1])
    all_point_data.append(["flux0", "vectors", f0])
    all_point_data.append(["flux1", "vectors", f1])
    
    # vector of time values
    tv_vec = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    
    # now begin writing the grid that will include all time-dependent data
    ts_ugrid = timeseries_unstructuredGrid(FILE_PATH_ALL, tv_vec)
    ts_ugrid.init_unstructuredGridToVTK(x, y, z, connectivity = conn, offsets = offset, cell_types = ctype, all_cell_data = all_cell_data, all_point_data = all_point_data, comments = comments)
    
    # BEGIN: create time-dependent data
    p0_tv = [p0] * 5
    p1_tv = [p1] * 5
    v0_tv = [v0] * 5
    v1_tv = [v1] * 5
    pot0_tv = [pot0] * 5
    pot1_tv = [pot1] * 5
    f0_tv = [f0] * 5
    f1_tv = [f1] * 5

    for ii in range(5):
        p0_tv[ii] = np.sin(np.pi * tv_vec[ii]) * p0_ref
        p1_tv[ii] = np.sin(np.pi * tv_vec[ii]) * p1_ref
        v0_tv[ii] = np.sin(np.pi * tv_vec[ii]) * v0_ref
        v1_tv[ii] = np.sin(np.pi * tv_vec[ii]) * v1_ref
        pot0_tv[ii] = np.sin(np.pi * tv_vec[ii]) * pot0_ref
        pot1_tv[ii] = np.sin(np.pi * tv_vec[ii]) * pot1_ref
        f0_tv[ii] = np.sin(np.pi * tv_vec[ii]) * f0_ref
        f1_tv[ii] = np.sin(np.pi * tv_vec[ii]) * f1_ref
    # END: create time-dependent data
    
    # write the data (in order!!!!)
    # this order needs to be the same as what was given above.
    # cell data first...
    ts_ugrid.append_data(p0_tv)
    ts_ugrid.append_data(p1_tv)
    ts_ugrid.append_data(v0_tv)
    ts_ugrid.append_data(v1_tv)
    
    # ... then point data.
    ts_ugrid.append_data(pot0_tv)
    ts_ugrid.append_data(pot1_tv)
    ts_ugrid.append_data(f0_tv)
    ts_ugrid.append_data(f1_tv)
    
    # close the file
    ts_ugrid.close_unstructuredGridToVTK()
    
    # now write a "fake" .vtu file for a mesh with one cell, which defines a "TimeValue" variable for that one cell.
    ts_ugrid.write_fake_grid_for_time(FILE_PATH_TIME)
    # this .vtu file can be included in the Paraview pipeline.  You can display the actual time variable by
    # using a Python Annotation Filter:  "Time: %1.2f" %TimeValue[0]
    # this is useful for non-uniform time steps
    # make sure to hide the fake grid

if __name__ == "__main__":
    run()
