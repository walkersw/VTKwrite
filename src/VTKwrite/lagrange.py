"""
VTKwrite.lagrange.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Helpers for VTK's arbitrary-order Lagrange cells (VtkLagrangeCurve, VtkLagrangeTriangle,
VtkLagrangeQuadrilateral, VtkLagrangeTetrahedron, VtkLagrangeHexahedron, VtkLagrangeWedge).

A Lagrange cell of degree p has its nodes at the EQUISPACED points of VTK's reference
cell, listed in VTK's own order (vertices, then edges, faces and the interior).  These
functions give that order, so that data from a finite element code can be put into it:

    lagrange_nodes(cell_type, degree)     reference coordinates of the nodes, in VTK order
    lagrange_npoints(cell_type, degree)   number of nodes
    node_permutation(cell_type, degree, nodes, vertices = None)
                                          the order in which to take your element's nodes
    cellwise_mesh(cell_points, cell_type) a mesh in which every cell has its own nodes
                                          (for discontinuous data, or any basis)

If the degrees of freedom of a finite element space are the values at these same
equispaced points, node_permutation turns its cell-to-node map into VTK connectivity:
    connectivity = cell_nodes[:, perm].ravel()
Otherwise (other nodes, e.g. Gauss-Lobatto points, or a hierarchical basis), evaluate the
geometry and the fields at lagrange_nodes(...) in every cell, and write the result with
cellwise_mesh.

VTK's reference cells (vertices in VTK's order):
    curve          0, 1
    triangle       (0,0), (1,0), (0,1)
    quadrilateral  (0,0), (1,0), (1,1), (0,1)
    tetrahedron    (0,0,0), (1,0,0), (0,1,0), (0,0,1)
    hexahedron     the unit cube; (0,0,0), (1,0,0), (1,1,0), (0,1,0), then the same at z = 1
    wedge          triangle x [0, 1]; (0,0,0), (1,0,0), (0,1,0), then the same at z = 1

Only NumPy is needed.  Anisotropic degrees (VTK's "HigherOrderDegrees" cell data) and
VTK's bubble-enriched 7-node triangle, 15-node tetrahedron and 21-node wedge are not
covered by these helpers (the writers accept such cells).

Copyright (c) 2026,  Shawn W. Walker
"""

import numpy as np

from .vtkbin import (VtkCellType, VtkLagrangeCurve, VtkLagrangeTriangle,
                     VtkLagrangeQuadrilateral, VtkLagrangeTetrahedron, VtkLagrangeHexahedron,
                     VtkLagrangeWedge)

__all__ = ["lagrange_nodes", "lagrange_npoints", "node_permutation", "cellwise_mesh"]

# cell type id -> (name, reference dimension)
_LAGRANGE = {VtkLagrangeCurve.tid : ("curve", 1), VtkLagrangeTriangle.tid : ("triangle", 2),
             VtkLagrangeQuadrilateral.tid : ("quadrilateral", 2),
             VtkLagrangeTetrahedron.tid : ("tetrahedron", 3),
             VtkLagrangeHexahedron.tid : ("hexahedron", 3), VtkLagrangeWedge.tid : ("wedge", 3)}


def _cell_id(cell_type):
    if isinstance(cell_type, VtkCellType):
        return cell_type.tid
    if isinstance(cell_type, (int, np.integer)) and not isinstance(cell_type, bool):
        return int(cell_type)
    raise TypeError("cell_type must be a VtkCellType (e.g. VtkLagrangeTriangle) or a VTK cell "
                    "type id, got %r" % (cell_type,))


def _lagrange_id(cell_type):
    tid = _cell_id(cell_type)
    if tid not in _LAGRANGE:
        raise ValueError("cell type %d is not a Lagrange cell; use one of %s" % (
            tid, ", ".join("%d (%s)" % (t, n) for t, (n, _) in sorted(_LAGRANGE.items()))))
    return tid


def _check_degree(degree):
    if isinstance(degree, bool) or not isinstance(degree, (int, np.integer)) or degree < 1:
        raise ValueError("degree must be an integer >= 1, got %r" % (degree,))
    return int(degree)


# ---------------------------------------------------------------------------
# VTK's node order, as integer lattice points (reference coordinates = lattice / p).
# Checked against VTK's own index functions (vtkLagrange*::PointIndexFromIJK, ::Index).
# ---------------------------------------------------------------------------
def _curve(p):
    return [(0,), (p,)] + [(i,) for i in range(1, p)]


def _triangle(p):
    """ Vertices, the edges 0-1, 1-2, 2-0 (each from its first vertex), then the interior,
        which is ordered in the same way as a triangle of degree p - 3. """
    if p < 0:
        return []
    if p == 0:
        return [(0, 0)]
    r = range(1, p)
    pts = [(0, 0), (p, 0), (0, p)]
    pts += [(i, 0) for i in r] + [(p - i, i) for i in r] + [(0, p - i) for i in r]
    pts += [(i + 1, j + 1) for i, j in _triangle(p - 3)]
    return pts


_TET_EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
_TET_FACES = ((0, 1, 3), (2, 3, 1), (0, 3, 2), (0, 2, 1))     # (first, second, third) vertex


def _tetrahedron(p):
    """ Vertices, the 6 edges, the 4 faces (each ordered as a triangle of degree p - 3 with
        the face's vertices in the order of _TET_FACES), then the interior, ordered as a
        tetrahedron of degree p - 4.  A lattice point is (w1, w2, w3): its barycentric
        weights (times p) for the vertices 1, 2, 3. """
    if p < 0:
        return []
    if p == 0:
        return [(0, 0, 0)]

    def lattice(weights):          # weights for the 4 vertices -> (w1, w2, w3)
        return (weights[1], weights[2], weights[3])

    pts = []
    for v in range(4):
        w = [0, 0, 0, 0]
        w[v] = p
        pts.append(lattice(w))
    for a, b in _TET_EDGES:
        for t in range(1, p):
            w = [0, 0, 0, 0]
            w[a], w[b] = p - t, t
            pts.append(lattice(w))
    for a, b, c in _TET_FACES:
        for u, v in _triangle(p - 3):
            w = [0, 0, 0, 0]
            w[a], w[b], w[c] = p - u - v - 2, u + 1, v + 1
            pts.append(lattice(w))
    pts += [(i + 1, j + 1, k + 1) for i, j, k in _tetrahedron(p - 4)]
    return pts


def _quadrilateral(p):
    """ Vertices; edges (0,0)-(p,0), (p,0)-(p,p), (0,p)-(p,p), (0,0)-(0,p), each with the
        free index increasing; then the interior with i fastest. """
    r = range(1, p)
    pts = [(0, 0), (p, 0), (p, p), (0, p)]
    pts += [(i, 0) for i in r] + [(p, j) for j in r] + [(i, p) for i in r] + [(0, j) for j in r]
    pts += [(i, j) for j in r for i in r]
    return pts


def _hexahedron(p):
    """ Vertices; the edges of the bottom (k = 0) and top (k = p) as for the quadrilateral;
        the vertical edges at the bottom vertices 0, 1, 2, 3; the faces i = 0, i = p,
        j = 0, j = p, k = 0, k = p (lowest free index fastest); then the interior. """
    r = range(1, p)
    pts = [(0, 0, 0), (p, 0, 0), (p, p, 0), (0, p, 0),
           (0, 0, p), (p, 0, p), (p, p, p), (0, p, p)]
    for k in (0, p):
        pts += ([(i, 0, k) for i in r] + [(p, j, k) for j in r] +
                [(i, p, k) for i in r] + [(0, j, k) for j in r])
    for i, j in ((0, 0), (p, 0), (p, p), (0, p)):
        pts += [(i, j, k) for k in r]
    for i in (0, p):
        pts += [(i, j, k) for k in r for j in r]
    for j in (0, p):
        pts += [(i, j, k) for k in r for i in r]
    for k in (0, p):
        pts += [(i, j, k) for j in r for i in r]
    pts += [(i, j, k) for k in r for j in r for i in r]
    return pts


def _wedge(p):
    """ Vertices (bottom triangle, then top); the triangle edges of the bottom and of the
        top (as for the triangle); the vertical edges; the bottom and top triangle faces
        (interior lattice points, i fastest); the quadrilateral faces j = 0 (i increasing),
        i + j = p (i decreasing) and i = 0 (j decreasing), each layer k after layer; then
        the interior, layer by layer. """
    r = range(1, p)
    tri = ((0, 0), (p, 0), (0, p))
    tri_inside = [(i, j) for j in r for i in r if i + j < p]
    pts = [(i, j, 0) for i, j in tri] + [(i, j, p) for i, j in tri]
    for k in (0, p):
        pts += [(i, 0, k) for i in r] + [(p - i, i, k) for i in r] + [(0, p - i, k) for i in r]
    for i, j in tri:
        pts += [(i, j, k) for k in r]
    for k in (0, p):
        pts += [(i, j, k) for i, j in tri_inside]
    pts += [(i, 0, k) for k in r for i in r]
    pts += [(p - i, i, k) for k in r for i in r]
    pts += [(0, p - j, k) for k in r for j in r]
    pts += [(i, j, k) for k in r for i, j in tri_inside]
    return pts


_ORDER = {VtkLagrangeCurve.tid : _curve, VtkLagrangeTriangle.tid : _triangle,
          VtkLagrangeQuadrilateral.tid : _quadrilateral,
          VtkLagrangeTetrahedron.tid : _tetrahedron, VtkLagrangeHexahedron.tid : _hexahedron,
          VtkLagrangeWedge.tid : _wedge}

_NPOINTS = {VtkLagrangeCurve.tid : lambda p: p + 1,
            VtkLagrangeTriangle.tid : lambda p: (p + 1) * (p + 2) // 2,
            VtkLagrangeQuadrilateral.tid : lambda p: (p + 1) ** 2,
            VtkLagrangeTetrahedron.tid : lambda p: (p + 1) * (p + 2) * (p + 3) // 6,
            VtkLagrangeHexahedron.tid : lambda p: (p + 1) ** 3,
            VtkLagrangeWedge.tid : lambda p: (p + 1) ** 2 * (p + 2) // 2}

# VTK's bubble-enriched quadratic cells: valid in files, not described by lagrange_nodes
_SPECIAL_NPOINTS = {VtkLagrangeTriangle.tid : 7, VtkLagrangeTetrahedron.tid : 15,
                    VtkLagrangeWedge.tid : 21}


def lagrange_npoints(cell_type, degree):
    """ Number of nodes of a Lagrange cell of the given degree. """
    return _NPOINTS[_lagrange_id(cell_type)](_check_degree(degree))


def _lagrange_degree(tid, npoints):
    """ The degree of a Lagrange cell with npoints nodes, or None if no degree fits
        (a valid special case such as the 7-node triangle returns 2). """
    if _SPECIAL_NPOINTS.get(tid) == npoints:
        return 2
    count = _NPOINTS[tid]
    p = 1
    while count(p) < npoints:
        p += 1
    return p if count(p) == npoints else None


def _reference_vertices(tid):
    return np.array(_ORDER[tid](1), dtype = float)


def lagrange_nodes(cell_type, degree, vertices = None):
    """ The nodes of a Lagrange cell, in VTK's order.

        PARAMETERS:
            cell_type: VtkLagrangeCurve, VtkLagrangeTriangle, VtkLagrangeQuadrilateral,
                       VtkLagrangeTetrahedron, VtkLagrangeHexahedron or VtkLagrangeWedge
                       (or its id, 68 to 73).
            degree:    the polynomial degree p >= 1.
            vertices:  optional: the vertices of another reference cell, listed in VTK's
                       vertex order (see the module docstring), with shape (nvertices, d).
                       The nodes are then mapped to that cell by the affine map that takes
                       VTK's reference vertices to these (it must exist: e.g. [-1, 1]^2 for
                       a quadrilateral, or any triangle).

        RETURNS:
            Array of shape (npoints, dim): the reference coordinates of the nodes, on VTK's
            reference cell (dim = 1, 2 or 3) or, with vertices, on the given cell (dim = d).
    """
    tid = _lagrange_id(cell_type)
    p = _check_degree(degree)
    nodes = np.array(_ORDER[tid](p), dtype = float) / p
    if vertices is None:
        return nodes
    A = _affine_map(tid, vertices)
    return np.column_stack((nodes, np.ones(len(nodes)))) @ A


def _affine_map(tid, vertices):
    """ The matrix A, shape (dim + 1, d), with [x, 1] @ A mapping VTK's reference cell to the
        cell with the given vertices. """
    name, dim = _LAGRANGE[tid]
    ref = _reference_vertices(tid)
    V = np.asarray(vertices, dtype = float)
    if V.ndim == 1 and dim == 1:
        V = V[:, None]
    if V.ndim != 2 or V.shape[0] != len(ref):
        raise ValueError("vertices of a %s must have shape (%d, d), got %s"
                         % (name, len(ref), np.shape(vertices)))
    M = np.column_stack((ref, np.ones(len(ref))))
    A, *_ = np.linalg.lstsq(M, V, rcond = None)
    scale = max(np.abs(V).max(), 1.0)
    if np.abs(M @ A - V).max() > 1e-10 * scale:
        raise ValueError("the given vertices are not an affine image of VTK's reference %s "
                         "(listed in VTK's vertex order?)" % name)
    if np.linalg.matrix_rank(A[:dim], tol = 1e-12 * scale) < dim:
        raise ValueError("the given vertices of the %s are degenerate" % name)
    return A


def node_permutation(cell_type, degree, nodes, vertices = None, tol = 1e-8):
    """ The order in which to take an element's nodes so that they are in VTK's order.

        PARAMETERS:
            cell_type, degree: as for lagrange_nodes.
            nodes:    array of shape (npoints, d): the reference coordinates of the nodes of
                      your element, in your order.
            vertices: the vertices of your reference cell in VTK's vertex order, shape
                      (nvertices, d) (default: VTK's reference cell; see lagrange_nodes).
            tol:      tolerance for matching coordinates, relative to the size of the cell.

        RETURNS:
            Integer array perm with nodes[perm] = lagrange_nodes(cell_type, degree, vertices),
            i.e. node perm[a] of your element is VTK's node a.  For a cell-to-node map
            cell_nodes of shape (ncells, npoints) in your order, the VTK connectivity is
            cell_nodes[:, perm].ravel().

        Raises ValueError if the nodes are not VTK's (equispaced) nodes.  Then evaluate the
        geometry and the fields at lagrange_nodes(...) instead.
    """
    tid = _lagrange_id(cell_type)
    target = lagrange_nodes(tid, degree, vertices)
    X = np.asarray(nodes, dtype = float)
    if X.ndim == 1 and target.shape[1] == 1:
        X = X[:, None]
    if X.ndim != 2 or X.shape != target.shape:
        raise ValueError("nodes must have shape %s for a degree-%d %s, got %s"
                         % (target.shape, degree, _LAGRANGE[tid][0], np.shape(nodes)))
    scale = max(np.ptp(target, axis = 0).max(), 1e-300)
    dist = np.linalg.norm(target[:, None, :] - X[None, :, :], axis = 2)
    perm = np.argmin(dist, axis = 1)
    missing = np.flatnonzero(dist[np.arange(len(perm)), perm] > tol * scale)
    if missing.size:
        a = int(missing[0])
        raise ValueError("VTK's node %d (at %s) is not among the given nodes (%d of %d nodes "
                         "unmatched).  VTK's Lagrange nodes are equispaced; for other nodes "
                         "(e.g. Gauss-Lobatto points) or a hierarchical basis, evaluate the "
                         "fields at lagrange_nodes(...) instead"
                         % (a, np.array2string(target[a], precision = 6), missing.size, len(perm)))
    if len(np.unique(perm)) != len(perm):
        raise ValueError("the given nodes contain duplicates (within tol)")
    return perm.astype(np.int64)


def cellwise_mesh(cell_points, cell_type):
    """ A mesh in which every cell has its own copies of its nodes, so that data may be
        discontinuous between cells (e.g. a DG field), or be given per cell for any basis.

        PARAMETERS:
            cell_points: array of shape (ncells, nnodes, 3), or (ncells, nnodes, 2) for z = 0:
                         the coordinates of the nodes of each cell, in VTK's order for the
                         cell type (for Lagrange cells: the order of lagrange_nodes).
            cell_type:   a VtkCellType (or its id); the same for all cells.

        RETURNS:
            A dictionary with the keys points, connectivity, offsets and cell_types, to be
            passed on as keyword arguments, e.g.
                mesh = cellwise_mesh(X, VtkLagrangeTriangle)
                unstructuredGridToVTK("dg", **mesh, all_point_data = {"u": U.reshape(-1)})
            Point data is given per cell node: an array of shape (ncells, nnodes) or
            (ncells, nnodes, k) becomes .reshape(ncells * nnodes) or (ncells * nnodes, k).
            Cell data is given per cell, as usual.
    """
    tid = _cell_id(cell_type)
    X = np.asarray(cell_points)
    if X.ndim != 3 or X.shape[2] not in (2, 3):
        raise ValueError("cell_points must have shape (ncells, nnodes, 3) or (ncells, nnodes, 2), "
                         "got %s" % (X.shape,))
    ncells, nnodes = X.shape[:2]
    return {"points" : X.reshape(ncells * nnodes, X.shape[2]),
            "connectivity" : np.arange(ncells * nnodes, dtype = np.int64),
            "offsets" : nnodes * np.arange(1, ncells + 1, dtype = np.int64),
            "cell_types" : np.full(ncells, tid, dtype = np.uint8)}
