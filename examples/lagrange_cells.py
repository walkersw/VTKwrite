"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of higher-order (Lagrange) cells: curved cells and high-order fields.

The domain is a quarter annulus, 1 <= r <= 2, 0 <= theta <= pi/2.

1. Continuous: a mesh of Lagrange quadrilaterals of degree 3 with one global set of
   nodes, shared between neighboring cells; the field u = r^2 cos(2 theta).
2. Discontinuous: Lagrange triangles of degree 2, each with its own nodes
   (cellwise_mesh), and a field that jumps between cells.

In ParaView, raise "Nonlinear Subdivision Level" (in the Properties panel, advanced) to
see the curved edges and the high-order fields.
"""

import os
import numpy as np
from VTKwrite import (unstructuredGridToVTK, lagrange_nodes, cellwise_mesh,
                      VtkLagrangeQuadrilateral, VtkLagrangeTriangle)

FILE_PATH_CG = "./lagrange_quads"
FILE_PATH_DG = "./lagrange_dg_triangles"

def clean():
    for f in (FILE_PATH_CG + ".vtu", FILE_PATH_DG + ".vtu"):
        try:
            os.remove(f)
        except OSError:
            pass

def annulus(rt):
    """ (r, theta) -> (x, y) """
    return np.stack((rt[..., 0] * np.cos(rt[..., 1]), rt[..., 0] * np.sin(rt[..., 1])), axis = -1)

def run():
    print("running lagrange_cells...")
    nr, nt = 2, 4                                   # cells in the radial and angular directions
    dr, dt = 1.0 / nr, 0.5 * np.pi / nt

    # 1. Degree-3 quadrilaterals.  The nodes of the whole mesh form a (3 nr + 1) x (3 nt + 1)
    #    lattice; node (I, J) has the number I + (3 nr + 1) J.  The nodes of cell (a, b) in
    #    VTK's order are the lattice points (3 a, 3 b) + p * lagrange_nodes(...).
    p = 3
    local = np.rint(p * lagrange_nodes(VtkLagrangeQuadrilateral, p)).astype(int)    # (16, 2)
    NI, NJ = p * nr + 1, p * nt + 1
    cells = [(p * a + local[:, 0]) + NI * (p * b + local[:, 1]) for b in range(nt) for a in range(nr)]
    I, J = np.meshgrid(np.arange(NI), np.arange(NJ), indexing = "ij")
    rt = np.stack((1.0 + I.T.ravel() * dr / p, J.T.ravel() * dt / p), axis = -1)   # number I + NI J
    u = rt[:, 0] ** 2 * np.cos(2 * rt[:, 1])
    unstructuredGridToVTK(FILE_PATH_CG, points = annulus(rt),
                          connectivity = np.concatenate(cells),
                          offsets = len(local) * np.arange(1, len(cells) + 1),
                          cell_types = np.full(len(cells), VtkLagrangeQuadrilateral.tid),
                          all_point_data = {"u": u})

    # 2. Degree-2 triangles with their own nodes: each (r, theta) rectangle is split into
    #    two triangles, and the reference nodes are mapped into each of them.
    q = 2
    ref = lagrange_nodes(VtkLagrangeTriangle, q)                                  # (6, 2)
    cell_rt = []
    for b in range(nt):
        for a in range(nr):
            c0 = np.array([1.0 + a * dr, b * dt])
            e_r, e_t = np.array([dr, 0.0]), np.array([0.0, dt])
            for v0, v1, v2 in ((c0, c0 + e_r, c0 + e_r + e_t), (c0, c0 + e_r + e_t, c0 + e_t)):
                cell_rt.append(v0 + ref[:, :1] * (v1 - v0) + ref[:, 1:] * (v2 - v0))
    cell_rt = np.array(cell_rt)                                                   # (ncells, 6, 2)
    ncells = len(cell_rt)
    # a field that is quadratic in each cell and jumps between cells
    v = cell_rt[..., 0] ** 2 + 0.1 * np.arange(ncells)[:, None]                   # (ncells, 6)
    mesh = cellwise_mesh(annulus(cell_rt), VtkLagrangeTriangle)
    unstructuredGridToVTK(FILE_PATH_DG, **mesh, all_point_data = {"v": v.reshape(-1)},
                          all_cell_data = {"cell": np.arange(ncells)})

if __name__ == "__main__":
    run()
