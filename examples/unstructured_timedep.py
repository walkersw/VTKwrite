"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to use the high level unstructuredGridToVTK function.
This example shows how to export *only* an unstructured grid given its
vertices and topology through connectivity and offset lists.
It also shows how to store a static grid and a "time-sequence" of data
all in one .vtu file, with TimeSeriesUnstructuredGrid.  In addition, the time
values are stored separately using a "fake" grid.
Check the VTK file format for details of the unstructured grid.

(This example used the class timeseries_unstructuredGrid, which is deprecated; see
"Migrating from timeseries_unstructuredGrid" in the README.  For more output options,
see unstructured_timeseries.py.)

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
from VTKwrite import unstructuredGridToVTK, TimeSeriesUnstructuredGrid
from VTKwrite import VtkTriangle, VtkQuad
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
    conn = np.zeros(10, dtype = np.int64)

    conn[0], conn[1], conn[2] = 0, 1, 3              # first triangle
    conn[3], conn[4], conn[5] = 1, 4, 3              # second triangle
    conn[6], conn[7], conn[8], conn[9] = 1, 2, 5, 4  # rectangle

    # Define offset of last vertex of each element
    offset = np.zeros(3, dtype = np.int64)
    offset[0] = 3
    offset[1] = 6
    offset[2] = 10

    # Define cell types

    ctype = np.zeros(3, dtype = np.uint8)
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
    
    # vector of time values
    tv_vec = np.array([0.0, 0.1, 0.2, 0.3, 0.4])

    # Write the grid and all time steps into one .vtu file, one step at a time.  The time
    # values are also written to a "fake" grid (FILE_PATH_TIME + ".vtu") with one cell,
    # whose cell variable "TimeValue" holds the time of each step.  That .vtu file can be
    # included in the Paraview pipeline: display the actual time with a Python Annotation
    # Filter, "Time: %1.2f" %TimeValue[0].  This is useful for non-uniform time steps;
    # make sure to hide the fake grid.
    with TimeSeriesUnstructuredGrid(FILE_PATH_ALL, x, y, z, connectivity = conn, offsets = offset,
                                    cell_types = ctype, times = tv_vec, comments = comments,
                                    time_values_file = FILE_PATH_TIME) as ts:
        for t in tv_vec:
            s = np.sin(np.pi * t)
            # data defined on cells (vectors: one row per cell) ...
            cell_data = {"pressure0": s * p0_ref, "pressure1": s * p1_ref,
                         "vel0": s * v0_ref.reshape(3, 3), "vel1": s * v1_ref.reshape(3, 3)}
            # ... and at points (vectors: one row per point)
            point_data = {"potential0": s * pot0_ref, "potential1": s * pot1_ref,
                          "flux0": s * f0_ref.reshape(6, 3), "flux1": s * f1_ref.reshape(6, 3)}
            ts.write_step(t, point_data = point_data, cell_data = cell_data)

if __name__ == "__main__":
    run()
