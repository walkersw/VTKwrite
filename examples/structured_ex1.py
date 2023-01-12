"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A more detailed example of how to export a structured grid.

Copyright (c) 01-12-2023,  Shawn W. Walker
"""

import os
from VTKwrite.interface import structuredToVTK
import numpy as np

FILE_PATH = "./structured_ex1"
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

    # Coordinates of the points of the Cartesian grid
    X = np.arange(0, lx + 0.1*dx, dx, dtype='float64')
    Y = np.arange(0, ly + 0.1*dy, dy, dtype='float64')
    Z = np.arange(0, lz + 0.1*dz, dz, dtype='float64')

    # create the "plaid" grid (sort of like Matlab's meshgrid)
    # NOTE: this is stored in the .vts file
    x = np.zeros((nx + 1, ny + 1, nz + 1))
    y = np.zeros((nx + 1, ny + 1, nz + 1))
    z = np.zeros((nx + 1, ny + 1, nz + 1))
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                x[i,j,k] = X[i]
                y[i,j,k] = Y[j]
                z[i,j,k] = Z[k]

    # create the "plaid" grid for the MIDPTS of the CELLS
    xm = np.zeros((nx, ny, nz))
    ym = np.zeros((nx, ny, nz))
    zm = np.zeros((nx, ny, nz))
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                xm[i,j,k] = X[i] + 0.5*dx
                ym[i,j,k] = Y[j] + 0.5*dy
                zm[i,j,k] = Z[k] + 0.5*dz

    # Cell variables
    P_cell = np.zeros((nx, ny, nz))
    print(P_cell.shape)
    # make up a smooth function
    P_cell = (xm - 0.3)**2 + (ym - 0.6)**2 + (zm - 0.5)**2
    P_cell = P_cell.reshape(ncells, order='F')
    print(P_cell.shape)
    Vec_cell = np.zeros((3, nx, ny, nz))
    print(Vec_cell.shape)
    # make up a smooth function
    Vec_cell[0,:,:,:] = 1*xm[:,:,:] + 0*ym[:,:,:] + 0*zm[:,:,:]
    Vec_cell[1,:,:,:] = 0*xm[:,:,:] + 1*ym[:,:,:] + 0*zm[:,:,:]
    Vec_cell[2,:,:,:] = 0*xm[:,:,:] + 0*ym[:,:,:] + 1*zm[:,:,:]
    Vec_cell = Vec_cell.reshape(3*ncells, order='F')
    print(Vec_cell.shape)

    # Point variables
    P_pnt = np.zeros((nx+1, ny+1, nz+1))
    print(P_pnt.shape)
    # make up a smooth function
    P_pnt = (x - 0.3)**2 + (y - 0.6)**2 + (z - 0.5)**2
    P_pnt = P_pnt.reshape(npoints, order='F')
    print(P_pnt.shape)
    Vec_pnt = np.zeros((3, nx+1, ny+1, nz+1))
    print(Vec_pnt.shape)
    # make up a smooth function
    Vec_pnt[0,:,:,:] = 1*x[:,:,:] + 0*y[:,:,:] + 0*z[:,:,:]
    Vec_pnt[1,:,:,:] = 0*x[:,:,:] + 1*y[:,:,:] + 0*z[:,:,:]
    Vec_pnt[2,:,:,:] = 0*x[:,:,:] + 0*y[:,:,:] + 1*z[:,:,:]
    Vec_pnt = Vec_pnt.reshape(3*npoints, order='F')
    print(Vec_pnt.shape)

    # data defined on cells
    all_cell_data = []
    all_cell_data.append(["pressure_cell", "scalars", P_cell])
    all_cell_data.append(["vector_cell", "vectors", Vec_cell])

    # data defined at points
    all_point_data = []
    all_point_data.append(["pressure_pnt", "scalars", P_pnt])
    all_point_data.append(["vector_pnt", "vectors", Vec_pnt])

    comments = [ "comment 1", "comment 2" ]
    structuredToVTK(FILE_PATH, x, y, z, all_cell_data = all_cell_data, all_point_data = all_point_data, comments = comments)

if __name__ == "__main__":
    run()
