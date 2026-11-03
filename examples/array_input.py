"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of the array forms that the writers accept (VTKwrite 0.1.0 and later):

* point and cell data as dictionaries {name: array}; the type (scalars, vectors,
  tensors) follows from the number of components per point or cell;
* multi-component data with one row per point/cell, e.g. of shape (npoints, 3);
* the points as one array, points = P, of shape (npoints, 3) or (npoints, 2);
* on structured grids, points and vector fields indexed [i, j, k, component].

The domain is a quarter of a thick-walled cylinder, 1 <= r <= 2, under an inner
pressure (the plane-strain Lame solution).  It is written twice: as an unstructured
grid of quads, with the displacement and the stress tensor, and as a structured grid.
"""

import os
import numpy as np
from VTKwrite import unstructuredGridToVTK, structuredToVTK, VtkQuad

FILE_PATH_UG = "./array_input_unstructured"
FILE_PATH_SG = "./array_input_structured"

def clean():
    for f in (FILE_PATH_UG + ".vtu", FILE_PATH_SG + ".vts"):
        try:
            os.remove(f)
        except OSError:
            pass

def lame(x, y, p = 1.0, a = 1.0, b = 2.0, E = 1.0, nu = 0.3):
    """ Displacement (..., 3) and stress (..., 3, 3) of a cylinder a <= r <= b with inner
        pressure p and a traction-free outer surface (plane strain). """
    r = np.hypot(x, y)
    A = p * a**2 / (b**2 - a**2)
    B = A * b**2
    s_rr, s_tt = A - B / r**2, A + B / r**2
    s_zz = nu * (s_rr + s_tt)
    u_r = (1 + nu) / E * ((1 - 2 * nu) * A * r + B / r)
    zero = np.zeros_like(r)
    e_r = np.stack((x / r, y / r, zero), axis = -1)
    e_t = np.stack((-y / r, x / r, zero), axis = -1)
    e_z = np.stack((zero, zero, zero + 1), axis = -1)
    outer = lambda e: e[..., :, None] * e[..., None, :]
    stress = (s_rr[..., None, None] * outer(e_r) + s_tt[..., None, None] * outer(e_t)
              + s_zz[..., None, None] * outer(e_z))
    return u_r[..., None] * e_r, stress

def run():
    print("running array_input...")
    nr, nt = 5, 9                                   # points in the radial and angular directions
    r = np.linspace(1.0, 2.0, nr)
    theta = np.linspace(0.0, 0.5 * np.pi, nt)
    R, TH = np.meshgrid(r, theta, indexing = "ij")  # [i, j]: radius i, angle j
    X, Y = R * np.cos(TH), R * np.sin(TH)
    u, stress = lame(X, Y)                          # shapes (nr, nt, 3) and (nr, nt, 3, 3)
    ring = np.repeat(np.arange(nr - 1), nt - 1).reshape(nr - 1, nt - 1)   # cell (i, j) -> i

    # 1. Unstructured grid of quads.  Point number i * nt + j is (r_i, theta_j).
    P = np.column_stack((X.ravel(), Y.ravel()))     # (npoints, 2): z = 0
    idx = np.arange(nr * nt).reshape(nr, nt)
    quads = np.column_stack((idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(),
                             idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()))
    unstructuredGridToVTK(FILE_PATH_UG, points = P,
                          connectivity = quads.ravel(),
                          offsets = 4 * np.arange(1, len(quads) + 1),
                          cell_types = np.full(len(quads), VtkQuad.tid, dtype = np.uint8),
                          all_point_data = {"displacement": u.reshape(-1, 3),     # vectors
                                            "stress": stress.reshape(-1, 9),      # tensors
                                            "inner_wall": R.ravel() == 1.0},      # bool -> UInt8
                          all_cell_data = {"ring": ring.ravel()})                 # scalars

    # 2. The same as a structured grid with one layer in z: arrays indexed [i, j, k(, c)].
    points = np.stack((X, Y), axis = -1)[:, :, None, :]         # (nr, nt, 1, 2): z = 0
    structuredToVTK(FILE_PATH_SG, points = points,
                    all_point_data = {"r": R[:, :, None],                    # (nr, nt, 1)
                                      "displacement": u[:, :, None, :]},     # (nr, nt, 1, 3)
                    all_cell_data = {"ring": ring[:, :, None]})              # (nr-1, nt-1, 1)

if __name__ == "__main__":
    run()
