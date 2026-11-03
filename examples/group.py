"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to use VtkGroup to write a .pvd collection: a list of data
files, each with a time value, that ParaView opens as one time series.

Use VtkGroup when the mesh itself changes in time (moving or deforming
domains, shape optimization, adaptive refinement), so that each time step
needs its own mesh.  For a fixed mesh, TimeSeriesUnstructuredGrid (see
unstructured_timeseries.py) is simpler and writes the mesh only once.

Here a square domain is stretched over time and re-meshed with a finer grid
at every step; u(x, y, t) = exp(-t) sin(pi x / a(t)) sin(pi y / b(t)) is
written on each mesh.

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
import shutil
import numpy as np
from VTKwrite.interface import unstructuredGridToVTK
from VTKwrite.vtkbin import VtkGroup, VtkTriangle

FILE_PATH = "./group"       # -> group.pvd, and the step files in the folder group/

def clean():
    try:
        os.remove(FILE_PATH + ".pvd")
    except OSError:
        pass
    shutil.rmtree(FILE_PATH, ignore_errors = True)

def rectangle_mesh(a, b, n):
    """ Points and triangles of [0, a] x [0, b], split into n x n squares of 2 triangles. """
    X, Y = np.meshgrid(np.linspace(0.0, a, n + 1), np.linspace(0.0, b, n + 1), indexing = "ij")
    x, y = X.ravel(), Y.ravel()
    idx = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    p, q = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel()     # corners of each square
    r, s = idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
    triangles = np.concatenate((np.column_stack((p, q, r)), np.column_stack((p, r, s))))
    return x, y, triangles

def run():
    print("Running group...")
    os.makedirs(FILE_PATH, exist_ok = True)
    times = [0.0, 0.25, 0.5, 1.0]

    with VtkGroup(FILE_PATH) as group:            # group.pvd is written when the block ends
        for k, t in enumerate(times):
            a, b = 1.0 + 0.5 * t, 1.0 + 0.25 * np.sin(np.pi * t)    # the domain at time t
            x, y, triangles = rectangle_mesh(a, b, n = 4 + 2 * k)  # a new, finer mesh
            z = np.zeros_like(x)
            u = np.exp(-t) * np.sin(np.pi * x / a) * np.sin(np.pi * y / b)
            ntri = len(triangles)
            step_file = unstructuredGridToVTK(
                os.path.join(FILE_PATH, "group_%04d" % k), x, y, z,
                connectivity = triangles.ravel(),
                offsets = 3 * np.arange(1, ntri + 1),
                cell_types = np.full(ntri, VtkTriangle.tid, dtype = np.uint8),
                all_cell_data = [["cell_id", "scalars", np.arange(ntri)]],
                all_point_data = [["u", "scalars", u]])
            group.addFile(filepath = step_file, sim_time = t)   # stored relative to group.pvd

if __name__ == "__main__":
    run()
