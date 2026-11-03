"""Higher-order cells: the Lagrange node order, the helpers in VTKwrite.lagrange, the
fixed-order quadratic/cubic cells, the point-count checks, and the file version that VTK
needs for Lagrange hexahedra.  Degrees 1 to 3; the node order alone also for degrees 4 to 6,
where the orientation of the nodes inside a face first matters (tetrahedron faces from
p = 4, the triangle faces of a wedge from p = 5).

The node order is checked against VTK in two independent ways: VTK's own index functions,
and geometry: a polynomial map of degree p, written as an order-p cell, must be reproduced
exactly by VTK at random parametric points (one misplaced node breaks this)."""
import os
import re

import numpy as np
import pytest

from VTKwrite import (unstructuredGridToVTK, TimeSeriesUnstructuredGrid,
                      lagrange_nodes, lagrange_npoints,
                      node_permutation, cellwise_mesh,
                      VtkLagrangeCurve, VtkLagrangeTriangle, VtkLagrangeQuadrilateral,
                      VtkLagrangeTetrahedron, VtkLagrangeHexahedron, VtkLagrangeWedge,
                      VtkQuadraticEdge, VtkQuadraticTriangle, VtkQuadraticQuad,
                      VtkQuadraticTetra, VtkQuadraticHexahedron, VtkQuadraticWedge,
                      VtkQuadraticPyramid, VtkBiQuadraticQuad, VtkTriQuadraticHexahedron,
                      VtkBiQuadraticTriangle, VtkCubicLine, VtkTriQuadraticPyramid,
                      VtkTriangle, VtkQuad, VtkPolygon, VtkPolyLine)
from conftest import read_vtk, vtk_array, no_leftovers

LAGRANGE = [VtkLagrangeCurve, VtkLagrangeTriangle, VtkLagrangeQuadrilateral,
            VtkLagrangeTetrahedron, VtkLagrangeHexahedron, VtkLagrangeWedge]
DIM = {68: 1, 69: 2, 70: 2, 71: 3, 72: 3, 73: 3}
DEGREES = [1, 2, 3]
CASES = [(c, p) for c in LAGRANGE for p in DEGREES]
IDS = ["%s-p%d" % (c.name, p) for c, p in CASES]


def content(path):
    with open(path, "rb") as f:
        return f.read()


def file_version(path):
    return re.search(rb'<VTKFile[^>]*version="([^"]*)"', content(path)).group(1).decode()


def polynomial_map(tid, p):
    """ A map of the reference cell of degree p in the space of the cell (P_p on
        simplices, Q_p on quadrilaterals and hexahedra, P_p x P_p on wedges). """
    def F(q):
        q = np.atleast_2d(q)
        x = np.zeros((len(q), 3))
        x[:, :q.shape[1]] = q
        r, s, t = x[:, 0], x[:, 1], x[:, 2]
        if tid == 68:
            g = (0 * r, 0.1 * r ** p, 0.05 * r ** max(p - 1, 0) * (p > 1))
        elif tid in (69, 71):                                  # total degree <= p
            g = (0.1 * s ** p, 0.1 * (r ** (p - 1) * t if p > 1 else t),
                 0.05 * (r * s ** (p - 1) if tid == 71 else 0 * r))
        elif tid == 70:                                        # degree <= p in each of r, s
            g = (0.1 * (r * s) ** p, 0.1 * r ** p * s, 0 * r)
        elif tid == 72:
            g = (0.1 * (s * t) ** p, 0.1 * r ** p * t, 0.05 * (r * s * t) ** p)
        else:                                                  # wedge: P_p(r, s) x P_p(t)
            g = (0.1 * s ** p * t ** p, 0.1 * r ** (p - 1) * t if p > 1 else 0.1 * t, 0.05 * r ** p)
        return x + np.column_stack(g)
    return F


def random_parametric_points(tid, n, seed=0):
    rng = np.random.default_rng(seed)
    q = rng.random((n, 3))
    if tid in (69, 73):
        q[:, :2] *= (rng.random((n, 1)) / q[:, :2].sum(axis=1, keepdims=True))
    elif tid == 71:
        q *= rng.random((n, 1)) / q.sum(axis=1, keepdims=True)
    q[:, DIM.get(tid, 3):] = 0.0
    return q


def max_map_error(vtk, cell, F, params):
    """ |x(xi) - F(xi)| at the given parametric points, x evaluated by VTK. """
    weights = [0.0] * cell.GetNumberOfPoints()
    err = 0.0
    for q in params:
        x = [0.0, 0.0, 0.0]
        cell.EvaluateLocation(vtk.reference(0), list(q), x, weights)
        err = max(err, np.abs(np.array(x) - F(q)[0]).max())
    return err


def vtk_lattice(vtk, tid, p):
    """ VTK's node order as lattice points, from VTK's index functions (curves: from the
        parametric coordinates of a cell). """
    pts = {}
    rng = range(p + 1)
    if tid == 68:
        n = p + 1
        c = vtk.vtkLagrangeCurve()
        c.GetPointIds().SetNumberOfIds(n)
        c.GetPoints().SetNumberOfPoints(n)
        for a in range(n):
            c.GetPointIds().SetId(a, a)
        pc = c.GetParametricCoords()
        return [(int(round(pc[3 * a] * p)),) for a in range(n)]
    for i in rng:
        for j in rng:
            for k in rng:
                if tid == 70 and k == 0:
                    pts[vtk.vtkLagrangeQuadrilateral.PointIndexFromIJK(i, j, [p, p])] = (i, j)
                elif tid == 72:
                    pts[vtk.vtkLagrangeHexahedron.PointIndexFromIJK(i, j, k, [p, p, p])] = (i, j, k)
                elif tid == 73 and i + j <= p:
                    pts[vtk.vtkLagrangeWedge.PointIndexFromIJK(i, j, k, [p, p, p])] = (i, j, k)
                elif tid == 69 and k == 0 and i + j <= p:
                    pts[vtk.vtkLagrangeTriangle.Index([i, j, p - i - j], p)] = (i, j)
                elif tid == 71 and i + j + k <= p:
                    pts[vtk.vtkLagrangeTetra.Index([i, j, k, p - i - j - k], p)] = (i, j, k)
    assert sorted(pts) == list(range(len(pts)))
    return [pts[a] for a in range(len(pts))]


# ---------------------------------------------------------------------------
# node order
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("cell, p", CASES, ids=IDS)
def test_node_order_matches_vtk_index_functions(vtk, cell, p):
    nodes = lagrange_nodes(cell, p)
    assert nodes.shape == (lagrange_npoints(cell, p), DIM[cell.tid])
    lattice = [tuple(int(v) for v in row) for row in np.rint(nodes * p)]
    np.testing.assert_allclose(nodes * p, lattice, atol=1e-12)
    assert lattice == vtk_lattice(vtk, cell.tid, p)


HIGH_CASES = [(c, p) for c in LAGRANGE for p in (4, 5, 6)]


@pytest.mark.parametrize("cell, p", HIGH_CASES, ids=["%s-p%d" % (c.name, p) for c, p in HIGH_CASES])
def test_node_order_matches_vtk_higher_degrees(vtk, cell, p):
    """ Order only (no files): degrees 4 to 6. """
    lattice = [tuple(int(v) for v in row) for row in np.rint(lagrange_nodes(cell, p) * p)]
    assert len(lattice) == lagrange_npoints(cell, p)
    assert lattice == vtk_lattice(vtk, cell.tid, p)


@pytest.mark.parametrize("cell, p", CASES, ids=IDS)
def test_polynomial_geometry_reproduced_by_vtk(tmp_path, vtk, cell, p):
    F = polynomial_map(cell.tid, p)
    nodes = lagrange_nodes(cell, p)
    n = len(nodes)
    path = unstructuredGridToVTK(str(tmp_path / "cell"), points=F(nodes), connectivity=np.arange(n),
                                 offsets=[n], cell_types=[cell.tid])
    out = read_vtk(vtk, path)
    vtk_cell = out.GetCell(0)
    assert vtk_cell.GetCellType() == cell.tid and vtk_cell.GetNumberOfPoints() == n
    assert max_map_error(vtk, vtk_cell, F, random_parametric_points(cell.tid, 25)) < 1e-12
    # VTK's parametric coordinates of the nodes are lagrange_nodes
    pc = vtk_cell.GetParametricCoords()
    pc = np.array([pc[i] for i in range(3 * n)]).reshape(n, 3)[:, :DIM[cell.tid]]
    np.testing.assert_allclose(pc, nodes, atol=1e-12)


@pytest.mark.parametrize("output", ["vtu", "pvd", "vtkhdf"])
def test_lagrange_hexahedra_in_time_series(tmp_path, vtk, output):
    """ Lagrange hexahedra in all three time-series outputs (the .vtu files need
        version 2.1, the .vtkhdf file is not affected). """
    if output == "vtkhdf":
        pytest.importorskip("h5py")
    p = 3
    F = polynomial_map(72, p)
    nodes = lagrange_nodes(VtkLagrangeHexahedron, p)
    n = len(nodes)
    mesh = dict(points=F(nodes), connectivity=np.arange(n), offsets=[n], cell_types=[72])
    with TimeSeriesUnstructuredGrid(str(tmp_path / "s"), **mesh, output=output) as ts:
        for k in range(2):
            ts.write_step(float(k), point_data={"u": np.full(n, float(k))})
    if output == "vtkhdf":
        r = vtk.vtkHDFReader()
        r.SetFileName(str(tmp_path / "s.vtkhdf"))
        r.UpdateTimeStep(1.0)
        out = r.GetOutput()
    else:
        path = str(tmp_path / "s.vtu") if output == "vtu" else str(tmp_path / "s" / "s_00001.vtu")
        assert file_version(path) == "2.1"
        out = read_vtk(vtk, path)
    assert max_map_error(vtk, out.GetCell(0), F, random_parametric_points(72, 25)) < 1e-12
    assert no_leftovers(tmp_path)


def test_file_version(tmp_path):
    nodes = lagrange_nodes(VtkLagrangeHexahedron, 2)
    hexa = unstructuredGridToVTK(str(tmp_path / "h"), points=nodes, connectivity=np.arange(27),
                                 offsets=[27], cell_types=[72])
    tet_nodes = lagrange_nodes(VtkLagrangeTetrahedron, 2)
    tet = unstructuredGridToVTK(str(tmp_path / "t"), points=tet_nodes, connectivity=np.arange(10),
                                offsets=[10], cell_types=[71])
    assert file_version(hexa) == "2.1" and file_version(tet) == "1.0"


# ---------------------------------------------------------------------------
# fixed-order cells (26-29, 34, 35, 37 are new; 21-25 for comparison)
# ---------------------------------------------------------------------------
FIXED = [VtkQuadraticEdge, VtkQuadraticTriangle, VtkQuadraticQuad, VtkQuadraticTetra,
         VtkQuadraticHexahedron, VtkQuadraticWedge, VtkQuadraticPyramid, VtkBiQuadraticQuad,
         VtkTriQuadraticHexahedron, VtkBiQuadraticTriangle, VtkCubicLine, VtkTriQuadraticPyramid]


@pytest.mark.parametrize("cell", FIXED, ids=[c.name for c in FIXED])
def test_fixed_order_cells(tmp_path, vtk, cell):
    """ Read back as the intended cell type: VTK's shape functions at VTK's parametric node
        coordinates give back the nodes; and (except for the pyramids, whose parametric
        space is a collapsed cube) a polynomial map is reproduced. """
    generic = vtk.vtkGenericCell()
    generic.SetCellType(cell.tid)
    ref_cell = generic.GetRepresentativeCell()
    n = ref_cell.GetNumberOfPoints()
    pc = ref_cell.GetParametricCoords()
    ref = np.array([pc[i] for i in range(3 * n)]).reshape(n, 3)
    X = np.random.default_rng(cell.tid).random((n, 3))
    path = unstructuredGridToVTK(str(tmp_path / "nodal"), points=X, connectivity=np.arange(n),
                                 offsets=[n], cell_types=[cell.tid])
    out = read_vtk(vtk, path).GetCell(0)
    assert out.GetCellType() == cell.tid and out.GetNumberOfPoints() == n
    for a in range(n):
        x, w = [0.0, 0.0, 0.0], [0.0] * n
        out.EvaluateLocation(vtk.reference(0), list(ref[a]), x, w)
        np.testing.assert_allclose(x, X[a], atol=1e-12)
    if cell in (VtkQuadraticPyramid, VtkTriQuadraticPyramid):
        return
    if cell is VtkCubicLine:                  # parametric coordinate in [-1, 1]
        F = lambda q: np.column_stack((np.atleast_2d(q)[:, 0], 0.1 * np.atleast_2d(q)[:, 0] ** 3,
                                       0 * np.atleast_2d(q)[:, 0]))
    else:                                       # complete quadratic map
        F = lambda q: (np.atleast_2d(q) + 0.1 * np.column_stack(
            (np.atleast_2d(q)[:, 1] ** 2, np.atleast_2d(q)[:, 0] * np.atleast_2d(q)[:, 2],
             np.atleast_2d(q)[:, 0] * np.atleast_2d(q)[:, 1])))
    path = unstructuredGridToVTK(str(tmp_path / "c"), points=F(ref), connectivity=np.arange(n),
                                 offsets=[n], cell_types=[cell.tid])
    assert file_version(path) == "1.0"
    out = read_vtk(vtk, path).GetCell(0)
    assert out.GetCellType() == cell.tid
    rng = np.random.default_rng(3)
    params = [ref[0] + 0.3 * rng.random() * (ref[1] - ref[0]) + 0.3 * rng.random() * (ref[-1] - ref[0])
              for _ in range(10)]
    assert max_map_error(vtk, out, F, params) < 1e-12


# ---------------------------------------------------------------------------
# number of points per cell
# ---------------------------------------------------------------------------
def write_cells(path, sizes, types, **kw):
    n = int(sum(sizes))
    P = np.random.default_rng(0).random((n, 3))
    return unstructuredGridToVTK(path, points=P, connectivity=np.arange(n), offsets=np.cumsum(sizes),
                                 cell_types=types, **kw)


@pytest.mark.parametrize("sizes, types, msg", [
    ([3, 9], [VtkTriangle.tid, VtkQuadraticTetra.tid], "cell 1 is a Quadratic_Tetra (type 24) with 9 points, but that cell type has 10 points"),
    ([4, 4, 3], [9, 9, 9], "cell 2 is a Quad (type 9) with 3 points"),
    ([2], [VtkPolygon.tid], "at least 3 points"),
    ([8], [VtkLagrangeTriangle.tid], "3, 6, 10, 15 points for the degrees 1, 2, 3, 4, ..., or 7"),
    ([16], [VtkLagrangeTetrahedron.tid], "4, 10, 20, 35 points"),
    ([12], [VtkLagrangeHexahedron.tid], "HigherOrderDegrees"),
    ([24], [VtkLagrangeWedge.tid], "6, 18, 40, 75 points"),
    ([1], [VtkLagrangeCurve.tid], "2, 3, 4, 5 points"),
    ([5, 5], [VtkQuadraticEdge.tid, VtkQuadraticEdge.tid], "(2 cells of this type are wrong)"),
])
def test_wrong_point_counts(tmp_path, sizes, types, msg):
    with pytest.raises(ValueError) as info:
        write_cells(str(tmp_path / "c"), sizes, types)
    assert msg in str(info.value)
    assert os.listdir(tmp_path) == []


def test_valid_point_counts(tmp_path, vtk):
    # VTK's bubble-enriched 7-node triangle, 15-node tetrahedron, 21-node wedge; a
    # 1-point polyline; unchecked ids; anisotropic Lagrange cells with HigherOrderDegrees
    write_cells(str(tmp_path / "a"), [7, 15, 21, 1, 2, 12], [69, 71, 73, VtkPolyLine.tid, 68, 0])
    path = write_cells(str(tmp_path / "b"), [6, 12], [VtkLagrangeQuadrilateral.tid, VtkLagrangeHexahedron.tid],
                       all_cell_data={"HigherOrderDegrees": np.array([[2, 1, 0], [2, 1, 1]])})
    out = read_vtk(vtk, path)
    assert out.GetNumberOfCells() == 2 and out.GetCell(1).GetNumberOfPoints() == 12
    with TimeSeriesUnstructuredGrid(str(tmp_path / "ts"), points=np.zeros((6, 3)), connectivity=np.arange(6),
                                    offsets=[6], cell_types=[70],
                                    static_cell_data={"HigherOrderDegrees": np.array([[2, 1, 0]])}) as ts:
        ts.write_step(0.0, point_data={"u": np.zeros(6)})


def test_point_counts_checked_in_time_series(tmp_path):
    with pytest.raises(ValueError, match="Lagrange_Triangle"):
        TimeSeriesUnstructuredGrid(str(tmp_path / "s"), points=np.zeros((5, 3)),
                                   connectivity=np.arange(5), offsets=[5], cell_types=[69])
    with pytest.raises(ValueError, match="Quad"):
        TimeSeriesUnstructuredGrid(str(tmp_path / "t"), points=np.zeros((3, 3)),
                                   connectivity=np.arange(3), offsets=[3], cell_types=[9])
    assert os.listdir(tmp_path) == []


# ---------------------------------------------------------------------------
# node_permutation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("cell, p", CASES, ids=IDS)
def test_permutation_recovers_vtk_order(cell, p):
    nodes = lagrange_nodes(cell, p)
    shuffle = np.random.default_rng(p).permutation(len(nodes))
    mine = nodes[shuffle]                                   # "my" element: some other order
    perm = node_permutation(cell, p, mine)
    np.testing.assert_array_equal(mine[perm], nodes)
    np.testing.assert_array_equal(shuffle[perm], np.arange(len(nodes)))


def test_permutation_on_another_reference_cell():
    # a triangle whose VTK vertices 0, 1, 2 sit at (1, 0), (0, 1), (0, 0)
    tri = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
    nodes = lagrange_nodes(VtkLagrangeTriangle, 3, vertices=tri)
    np.testing.assert_allclose(nodes[:3], tri, atol=1e-12)
    mine = nodes[::-1]
    perm = node_permutation(VtkLagrangeTriangle, 3, mine, vertices=tri)
    np.testing.assert_array_equal(perm, np.arange(10)[::-1])
    # [-1, 1]^2 and [-1, 1]^3
    square = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
    q = lagrange_nodes(VtkLagrangeQuadrilateral, 2, vertices=square)
    np.testing.assert_allclose(q, 2 * lagrange_nodes(VtkLagrangeQuadrilateral, 2) - 1, atol=1e-12)
    cube = 2 * lagrange_nodes(VtkLagrangeHexahedron, 1) - 1
    h = lagrange_nodes(VtkLagrangeHexahedron, 3, vertices=cube)
    perm = node_permutation(VtkLagrangeHexahedron, 3, h[::-1], vertices=cube)
    np.testing.assert_array_equal(perm, np.arange(64)[::-1])


def test_permutation_gives_connectivity(tmp_path, vtk):
    """ The documented use: a continuous P2 mesh numbered in another local order. """
    # two P2 triangles on the unit square; local order: vertices, then the midpoint of the
    # edge opposite vertex 0, 1, 2 (a common convention, not VTK's)
    my_ref = np.array([[0, 0], [1, 0], [0, 1], [0.5, 0.5], [0, 0.5], [0.5, 0]])
    perm = node_permutation(VtkLagrangeTriangle, 2, my_ref)
    corners = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    tris = np.array([[0, 1, 2], [0, 2, 3]])
    nodes, cell_nodes = {}, []
    for t in tris:
        row = []
        for b in my_ref:                     # node at barycentric position b of this cell
            x = corners[t[0]] + b[0] * (corners[t[1]] - corners[t[0]]) + b[1] * (corners[t[2]] - corners[t[0]])
            row.append(nodes.setdefault(tuple(np.round(x, 12)), len(nodes)))
        cell_nodes.append(row)
    cell_nodes = np.array(cell_nodes)
    X = np.array(list(nodes))
    F = lambda x: np.column_stack((x[:, 0] + 0.2 * x[:, 1] ** 2, x[:, 1] + 0.1 * x[:, 0] * x[:, 1],
                                   0.3 * x[:, 0] ** 2))
    path = unstructuredGridToVTK(str(tmp_path / "p2"), points=F(X),
                                 connectivity=cell_nodes[:, perm].ravel(), offsets=[6, 12],
                                 cell_types=[69, 69])
    out = read_vtk(vtk, path)
    assert out.GetNumberOfPoints() == 9                      # shared nodes stay shared
    for c, t in enumerate(tris):
        A = np.column_stack((corners[t[1]] - corners[t[0]], corners[t[2]] - corners[t[0]]))
        G = lambda q: F(corners[t[0]] + np.atleast_2d(q)[:, :2] @ A.T)
        assert max_map_error(vtk, out.GetCell(c), G, random_parametric_points(69, 10, seed=c)) < 1e-12


@pytest.mark.parametrize("args, exc, msg", [
    ((VtkLagrangeCurve, 3, np.array([0.0, 1.0, 0.5 - 0.5 / 5 ** 0.5, 0.5 + 0.5 / 5 ** 0.5])),
     ValueError, "equispaced"),                              # Gauss-Lobatto points
    ((VtkLagrangeTriangle, 2, np.zeros((5, 2))), ValueError, "must have shape (6, 2)"),
    ((VtkLagrangeTriangle, 1, np.array([[0, 0], [1, 0], [1, 0]])), ValueError, "not among"),
    ((VtkQuad, 1, np.zeros((4, 2))), ValueError, "not a Lagrange cell"),
    ((VtkLagrangeTriangle, 0, np.zeros((1, 2))), ValueError, "degree must be an integer >= 1"),
    (("triangle", 2, np.zeros((6, 2))), TypeError, "cell_type must be"),
])
def test_permutation_errors(args, exc, msg):
    with pytest.raises(exc) as info:
        node_permutation(*args)
    assert msg in str(info.value)


def test_reference_cell_errors():
    trapezoid = np.array([[0.0, 0.0], [2.0, 0.0], [1.5, 1.0], [0.5, 1.0]])
    with pytest.raises(ValueError, match="not an affine image"):
        lagrange_nodes(VtkLagrangeQuadrilateral, 2, vertices=trapezoid)
    with pytest.raises(ValueError, match="degenerate"):
        lagrange_nodes(VtkLagrangeTriangle, 2, vertices=[[0, 0], [1, 1], [2, 2]])
    with pytest.raises(ValueError, match="shape"):
        lagrange_nodes(VtkLagrangeTetrahedron, 2, vertices=np.zeros((3, 3)))


# ---------------------------------------------------------------------------
# cellwise_mesh
# ---------------------------------------------------------------------------
def test_cellwise_mesh(tmp_path, vtk):
    p = 2
    ref = lagrange_nodes(VtkLagrangeTriangle, p)                  # (6, 2)
    corners = np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]])
    tris = np.array([[0, 1, 2], [0, 2, 3]])
    X = np.array([corners[t[0]] + ref @ np.column_stack((corners[t[1]] - corners[t[0]],
                                                           corners[t[2]] - corners[t[0]])).T
                  for t in tris])                                  # (2, 6, 2)
    U = np.array([[1.0] * 6, [2.0] * 6])                           # discontinuous across the diagonal
    mesh = cellwise_mesh(X, VtkLagrangeTriangle)
    path = unstructuredGridToVTK(str(tmp_path / "dg"), **mesh, all_point_data={"u": U.reshape(-1)},
                                 all_cell_data={"k": np.array([0, 1])})
    out = read_vtk(vtk, path)
    assert out.GetNumberOfPoints() == 12 and out.GetNumberOfCells() == 2
    np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "u"), U.ravel())
    # the same mesh works for time series
    with TimeSeriesUnstructuredGrid(str(tmp_path / "ts"), **mesh, output="pvd") as ts:
        ts.write_step(0.0, point_data={"u": U.reshape(-1)})
    with pytest.raises(ValueError, match="shape"):
        cellwise_mesh(X[0], VtkLagrangeTriangle)
    with pytest.raises(ValueError, match="Lagrange_Triangle"):         # 5 nodes: not a triangle
        unstructuredGridToVTK(str(tmp_path / "bad"), **cellwise_mesh(X[:, :5], VtkLagrangeTriangle))
