"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to use TimeSeriesUnstructuredGrid to write a time-dependent
simulation one step at a time: into a single .vtu file, into a single .vtkhdf
(HDF5) file, or as a .pvd collection with one .vtu file per step.

The mesh is a structured triangulation of the unit square; the "solution" is
u(x, y, t) = exp(-t) sin(pi x) sin(pi y), with its gradient as a vector field,
and the time steps are not uniform.
"""

import os
import shutil
import numpy as np
from VTKwrite.interface import TimeSeriesUnstructuredGrid
from VTKwrite.vtkbin import VtkTriangle

FILE_PATH_VTU = "./timeseries_single"
FILE_PATH_PVD = "./timeseries_collection"
FILE_PATH_HDF = "./timeseries_single_hdf"

def clean():
    for f in (FILE_PATH_VTU + ".vtu", FILE_PATH_VTU + "_time_values.vtu", FILE_PATH_PVD + ".pvd",
              FILE_PATH_HDF + ".vtkhdf"):
        try:
            os.remove(f)
        except OSError:
            pass
    shutil.rmtree(FILE_PATH_PVD, ignore_errors = True)

def unit_square_mesh(n):
    """ Points and triangles of an n x n grid of squares, each split into 2 triangles. """
    s = np.linspace(0.0, 1.0, n + 1)
    X, Y = np.meshgrid(s, s, indexing = "ij")
    x, y = X.ravel(), Y.ravel()
    z = np.zeros_like(x)
    idx = np.arange((n + 1) ** 2).reshape(n + 1, n + 1)
    a, b = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel()     # corners of each square
    c, d = idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
    triangles = np.concatenate((np.column_stack((a, b, c)), np.column_stack((a, c, d))))
    connectivity = triangles.ravel()
    offsets = 3 * np.arange(1, len(triangles) + 1)
    cell_types = np.full(len(triangles), VtkTriangle.tid, dtype = np.uint8)
    return x, y, z, connectivity, offsets, cell_types

def solution(x, y, t):
    u = np.exp(-t) * np.sin(np.pi * x) * np.sin(np.pi * y)
    grad_u = np.exp(-t) * np.pi * np.column_stack((np.cos(np.pi * x) * np.sin(np.pi * y),
                                                   np.sin(np.pi * x) * np.cos(np.pi * y),
                                                   np.zeros_like(x)))     # one row per point
    return u, grad_u

def run():
    print("running unstructured_timeseries...")
    x, y, z, conn, offsets, cell_types = unit_square_mesh(16)
    ncells = len(cell_types)
    region = np.where(np.arange(ncells) < ncells // 2, 1, 2)     # does not change in time

    times = [0.0, 0.05, 0.1, 0.2, 0.4, 0.8]                       # non-uniform steps

    # 1. Everything in one .vtu file.  The time values need not be known in advance.
    with TimeSeriesUnstructuredGrid(FILE_PATH_VTU, x, y, z, conn, offsets, cell_types,
                                    static_cell_data = {"region": region},
                                    comments = ["u = exp(-t) sin(pi x) sin(pi y)"]) as ts:
        for t in times:
            u, grad_u = solution(x, y, t)
            ts.write_step(t, point_data = {"u": u, "grad_u": grad_u})
    # -> timeseries_single.vtu, and timeseries_single_time_values.vtu with the actual times
    # NOTE: VTK 9.6+ (ParaView 6.1) reads only the first step of such files (a VTK bug);
    #       use the .vtkhdf or .pvd output below with those versions.

    # 2. Everything in one VTKHDF (HDF5) file, with the actual time values (needs h5py).
    try:
        import h5py
    except ImportError:
        print("  h5py is not installed: skipping the .vtkhdf output")
    else:
        with TimeSeriesUnstructuredGrid(FILE_PATH_HDF, x, y, z, conn, offsets, cell_types,
                                        output = "vtkhdf",
                                        static_cell_data = {"region": region}) as ts:
            for t in times:
                u, grad_u = solution(x, y, t)
                ts.write_step(t, point_data = {"u": u, "grad_u": grad_u})
        # -> timeseries_single_hdf.vtkhdf

    # 3. A .pvd collection: timeseries_collection.pvd + timeseries_collection/*.vtu
    #    (readable by all ParaView versions, and usable even if the run stops early)
    with TimeSeriesUnstructuredGrid(FILE_PATH_PVD, x, y, z, conn, offsets, cell_types,
                                    output = "pvd", static_cell_data = {"region": region}) as ts:
        for t in times:
            u, grad_u = solution(x, y, t)
            ts.write_step(t, point_data = {"u": u, "grad_u": grad_u})

if __name__ == "__main__":
    run()
