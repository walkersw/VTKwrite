"""Friendlier input: field dictionaries, (n, k) arrays, points=, (nx, ny, nz, k) grids.

The new forms must give files that are byte for byte the same as the older forms
(list of [name, type, data] with interleaved 1-D arrays, and x, y, z); the structured
grids are also checked entry by entry with VTK's own (i, j, k) -> id functions.
"""
import os

import numpy as np
import pytest

from VTKwrite import (unstructuredGridToVTK, pointsToVTK, pointsToVTKAsTIN, linesToVTK,
                      polyLinesToVTK, cylinderToVTK, imageToVTK, rectilinearToVTK,
                      structuredToVTK, TimeSeriesUnstructuredGrid,
                      VtkTriangle, VtkQuad)
from conftest import read_vtk, vtk_array, no_leftovers

X = np.array([0., 1., 2., 0., 1., 2.])
Y = np.array([0., 0., 0., 1., 1., 1.])
Z = np.array([0., 0., 0.5, 0., 0., 0.25])
P = np.column_stack((X, Y, Z))
TOPO = dict(connectivity=[0, 1, 3, 1, 2, 5, 4], offsets=[3, 7],
            cell_types=[VtkTriangle.tid, VtkQuad.tid])
U = np.linspace(0, 1, 6)                         # scalars at the points
Q = np.arange(18.0).reshape(6, 3) - 4            # vectors at the points, one row per point
T = np.arange(54.0).reshape(6, 9)                # tensors at the points
UV = np.arange(12.0).reshape(6, 2) / 7           # 2 texture coordinates
C = np.array([1.5, -2.0])                        # scalars on the cells
CV = np.array([[1., 2., 3.], [4., 5., 6.]])      # vectors on the cells
MARK = np.array([3, 7], dtype=np.int32)

# the same fields in the older form
POINT_LIST = [["u", "scalars", U], ["q", "vectors", Q.ravel()], ["T", "tensors", T.ravel()],
              ["uv", "tcoords", UV.ravel()], ["n", "normals", Q.ravel()]]
CELL_LIST = [["c", "scalars", C], ["cv", "vectors", CV.ravel()], ["mark", "scalars", MARK]]
POINT_DICT = {"u": U, "q": Q, "T": T, "uv": ("tcoords", UV), "n": ("normals", Q)}
CELL_DICT = {"c": C, "cv": CV, "mark": MARK}


def content(path):
    with open(path, "rb") as f:
        return f.read()


def same_files(*paths):
    first = content(paths[0])
    return all(content(p) == first for p in paths[1:])


# ---------------------------------------------------------------------------
# field dictionaries and (n, k) arrays
# ---------------------------------------------------------------------------
def test_unstructured_dict_and_rows_match_the_list_form(tmp_path, vtk):
    old = unstructuredGridToVTK(str(tmp_path / "old"), X, Y, Z, **TOPO,
                                all_point_data=POINT_LIST, all_cell_data=CELL_LIST)
    new = unstructuredGridToVTK(str(tmp_path / "new"), X, Y, Z, **TOPO,
                                all_point_data=POINT_DICT, all_cell_data=CELL_DICT)
    # (n, k) rows in the list form too
    rows = unstructuredGridToVTK(str(tmp_path / "rows"), X, Y, Z, **TOPO,
                                 all_point_data=[["u", "scalars", U], ["q", "vectors", Q],
                                                 ["T", "tensors", T], ["uv", "tcoords", UV],
                                                 ["n", "normals", Q]],
                                 all_cell_data=[["c", "scalars", C], ["cv", "vectors", CV],
                                                ["mark", "scalars", MARK]])
    assert same_files(old, new, rows)
    out = read_vtk(vtk, new)
    pd, cd = out.GetPointData(), out.GetCellData()
    np.testing.assert_array_equal(vtk_array(vtk, pd, "q"), Q)
    np.testing.assert_array_equal(vtk_array(vtk, pd, "T"), T)
    np.testing.assert_array_equal(vtk_array(vtk, pd, "uv"), UV)
    np.testing.assert_array_equal(vtk_array(vtk, cd, "cv"), CV)
    np.testing.assert_array_equal(vtk_array(vtk, cd, "mark"), MARK)
    assert pd.GetScalars().GetName() == "u" and pd.GetVectors().GetName() == "q"
    assert pd.GetTensors().GetName() == "T" and pd.GetNormals().GetName() == "n"
    assert pd.GetTCoords().GetName() == "uv"


def test_dict_order_is_the_file_order(tmp_path):
    a = unstructuredGridToVTK(str(tmp_path / "a"), X, Y, Z, **TOPO,
                              all_point_data={"u": U, "w": 2 * U})
    b = unstructuredGridToVTK(str(tmp_path / "b"), X, Y, Z, **TOPO,
                              all_point_data={"w": 2 * U, "u": U})
    assert b'Scalars="u"' in content(a) and b'Scalars="w"' in content(b)


def test_interleaved_1d_arrays_in_a_dict(tmp_path):
    """ Same rule as TimeSeriesUnstructuredGrid: 3 values per point are vectors. """
    a = unstructuredGridToVTK(str(tmp_path / "a"), X, Y, Z, **TOPO, all_point_data={"q": Q.ravel()})
    b = unstructuredGridToVTK(str(tmp_path / "b"), X, Y, Z, **TOPO,
                              all_point_data=[["q", "vectors", Q.ravel()]])
    assert same_files(a, b)


def test_caller_data_not_modified(tmp_path):
    fields = {"q": Q.copy(), "uv": ("tcoords", UV.copy())}
    unstructuredGridToVTK(str(tmp_path / "a"), X, Y, Z, **TOPO, all_point_data=fields)
    assert list(fields) == ["q", "uv"] and fields["uv"][0] == "tcoords"
    np.testing.assert_array_equal(fields["q"], Q)
    assert fields["q"].shape == (6, 3)


def test_fortran_ordered_and_strided_rows(tmp_path):
    Qf = np.asfortranarray(Q)
    Qs = np.repeat(Q, 2, axis=1)[:, ::2]                 # non-contiguous
    paths = [unstructuredGridToVTK(str(tmp_path / name), X, Y, Z, **TOPO,
                                   all_point_data={"q": q})
             for name, q in (("c", Q), ("f", Qf), ("s", Qs))]
    assert same_files(*paths)


@pytest.mark.parametrize("writer", ["points", "lines", "polylines", "tin"])
def test_other_unstructured_writers(tmp_path, writer):
    if writer == "tin":
        pytest.importorskip("scipy")
    point_list, point_dict = [["q", "vectors", Q.ravel()]], {"q": Q}
    cell_kw = {}
    if writer == "points":
        call = lambda path, pd, cd: pointsToVTK(path, X, Y, Z, all_point_data=pd)
    elif writer == "tin":
        call = lambda path, pd, cd: pointsToVTKAsTIN(path, X, Y, Z, data=pd)
    elif writer == "lines":                              # 3 lines
        call = lambda path, pd, cd: linesToVTK(path, X, Y, Z, all_point_data=pd, all_cell_data=cd)
        cell_kw = dict(old=[["d", "vectors", np.arange(9.0)]], new={"d": np.arange(9.0).reshape(3, 3)})
    else:                                                # polylines of 2 and 4 points
        call = lambda path, pd, cd: polyLinesToVTK(path, X, Y, Z, [2, 4], all_point_data=pd,
                                                   all_cell_data=cd)
        cell_kw = dict(old=[["d", "vectors", np.arange(6.0)]], new={"d": np.arange(6.0).reshape(2, 3)})
    old = call(str(tmp_path / "old"), point_list, cell_kw.get("old"))
    new = call(str(tmp_path / "new"), point_dict, cell_kw.get("new"))
    assert same_files(old, new)


def test_tin_name_clash_with_elevation(tmp_path):
    pytest.importorskip("scipy")
    with pytest.raises(ValueError, match="'Elevation' is used twice"):
        pointsToVTKAsTIN(str(tmp_path / "t"), X, Y, Z, data={"Elevation": U})
    assert os.listdir(tmp_path) == []


def test_cylinder_dicts(tmp_path):
    nl, npil = 2, 4
    c, h = np.arange(8.0), np.arange(12.0)
    old = cylinderToVTK(str(tmp_path / "old"), 0, 0, 0, 1, 1, nl, npil,
                        all_cell_data=[["c", "scalars", c]], pointData={"h": h})
    new = cylinderToVTK(str(tmp_path / "new"), 0, 0, 0, 1, 1, nl, npil,
                        all_cell_data={"c": c}, all_point_data={"h": h})
    positional = cylinderToVTK(str(tmp_path / "pos"), 0, 0, 0, 1, 1, nl, npil, {"c": c}, {"h": h})
    assert same_files(old, new, positional)
    vec = cylinderToVTK(str(tmp_path / "vec"), 0, 0, 0, 1, 1, nl, npil,
                        all_point_data={"q": np.ones((12, 3))})
    assert b'Vectors="q"' in content(vec)
    with pytest.raises(ValueError, match="'c' is used twice"):
        cylinderToVTK(str(tmp_path / "dup"), 0, 0, 0, 1, 1, nl, npil,
                      all_cell_data={"c": c}, cellData={"c": c})


@pytest.mark.parametrize("kw, exc, msg", [
    (dict(all_point_data={"q": (Q[:, 0], Q[:, 1], Q[:, 2])}), TypeError, "column_stack"),
    (dict(all_point_data=[["q", "vectors", [Q[:, 0], Q[:, 1], Q[:, 2]]]]), TypeError, "column_stack"),
    (dict(all_point_data={"q": [[0, 0, 1]] * 6}), TypeError, "column_stack"),
    (dict(all_point_data={"q": Q[:, :2]}), ValueError, "give the type explicitly"),
    (dict(all_point_data={"q": ("vectors", Q[:, :2])}), ValueError, "2 components per point, but vectors have 3"),
    (dict(all_point_data=[["q", "scalars", Q]]), ValueError, "3 components per point, but scalars have 1"),
    (dict(all_point_data={"q": Q.T}), ValueError, "one row per point (6 rows)"),
    (dict(all_point_data={"q": np.zeros(7)}), ValueError, "not a multiple of 6 points"),
    (dict(all_cell_data={"c": np.zeros((2, 4))}), ValueError, "give the type explicitly"),
    (dict(all_point_data={"q": ("vector", Q)}), ValueError, "type 'vector' is not allowed"),
    (dict(all_point_data={1: U}), TypeError, "names must be non-empty strings"),
    (dict(all_point_data={"q": Q.reshape(2, 3, 3)}), ValueError, "a 3-D array is read as one value"),
    (dict(all_point_data="u"), TypeError, "a dictionary {name: array} or a list"),
])
def test_invalid_fields(tmp_path, kw, exc, msg):
    with pytest.raises(exc) as info:
        unstructuredGridToVTK(str(tmp_path / "ug"), X, Y, Z, **TOPO, **kw)
    assert msg in str(info.value)
    assert os.listdir(tmp_path) == []


def test_3d_arrays_still_taken_in_fortran_order(tmp_path):
    """ As before, for the unstructured writers: a 3-D array of values per point. """
    u3 = U.reshape(1, 2, 3)
    a = pointsToVTK(str(tmp_path / "a"), X, Y, Z, all_point_data={"u": u3})
    b = pointsToVTK(str(tmp_path / "b"), X, Y, Z, all_point_data=[["u", "scalars", u3]])
    c = pointsToVTK(str(tmp_path / "c"), X, Y, Z,
                    all_point_data=[["u", "scalars", u3.ravel(order="F")]])
    assert same_files(a, b, c)


# ---------------------------------------------------------------------------
# points=
# ---------------------------------------------------------------------------
def test_points_keyword_matches_xyz(tmp_path):
    kw = dict(all_point_data=POINT_DICT, all_cell_data=CELL_DICT)
    a = unstructuredGridToVTK(str(tmp_path / "a"), X, Y, Z, **TOPO, **kw)
    b = unstructuredGridToVTK(str(tmp_path / "b"), points=P, **TOPO, **kw)
    c = unstructuredGridToVTK(str(tmp_path / "c"), points=[tuple(p) for p in P], **TOPO, **kw)
    d = unstructuredGridToVTK(str(tmp_path / "d"), points=np.asfortranarray(P), **TOPO, **kw)
    assert same_files(a, b, c, d)


def test_points_2d_means_z_zero(tmp_path):
    a = unstructuredGridToVTK(str(tmp_path / "a"), X, Y, np.zeros(6), **TOPO)
    b = unstructuredGridToVTK(str(tmp_path / "b"), points=P[:, :2], **TOPO)
    assert same_files(a, b)
    P32 = P[:, :2].astype(np.float32)
    c = pointsToVTK(str(tmp_path / "c"), points=P32)
    d = pointsToVTK(str(tmp_path / "d"), P32[:, 0], P32[:, 1], np.zeros(6, np.float32))
    assert same_files(c, d) and b'type="Float32"' in content(c)


@pytest.mark.parametrize("writer", ["points", "lines", "polylines", "tin2", "tin3"])
def test_points_keyword_other_writers(tmp_path, writer):
    if writer.startswith("tin"):
        pytest.importorskip("scipy")
    fn, extra = {"points": (pointsToVTK, {}), "lines": (linesToVTK, {}),
                 "polylines": (polyLinesToVTK, {"pointsPerLine": [2, 4]}),
                 "tin2": (pointsToVTKAsTIN, {"ndim": 2}),
                 "tin3": (pointsToVTKAsTIN, {"ndim": 3})}[writer]
    a = fn(str(tmp_path / "a"), X, Y, Z, **extra)
    b = fn(str(tmp_path / "b"), points=P, **extra)
    assert same_files(a, b)


@pytest.mark.parametrize("output", ["vtu", "pvd", "vtkhdf"])
def test_points_keyword_time_series(tmp_path, output):
    if output == "vtkhdf":
        pytest.importorskip("h5py")
    def write(name, **mesh):
        with TimeSeriesUnstructuredGrid(str(tmp_path / name), **mesh, **TOPO, output=output,
                                        static_cell_data={"mark": MARK}) as ts:
            for k in range(3):
                ts.write_step(0.5 * k, point_data={"u": U + k, "q": Q * k})
        return ts
    write("a", x=X, y=Y, z=Z)
    write("b", points=P)
    if output == "pvd":
        for k in range(3):
            assert same_files(tmp_path / "a" / ("a_%05d.vtu" % k), tmp_path / "b" / ("b_%05d.vtu" % k))
    else:
        assert same_files(tmp_path / ("a." + output), tmp_path / ("b." + output))
    assert no_leftovers(tmp_path)


def test_time_series_positional_call_unchanged(tmp_path):
    ts = TimeSeriesUnstructuredGrid(str(tmp_path / "s"), X, Y, Z, TOPO["connectivity"],
                                    TOPO["offsets"], TOPO["cell_types"], [0.0], "vtu")
    ts.write_step(point_data={"u": U})
    assert ts.close().endswith("s.vtu")


@pytest.mark.parametrize("call, exc, msg", [
    (lambda p: unstructuredGridToVTK(p, X, Y, Z, **TOPO, points=P), TypeError, "not both"),
    (lambda p: unstructuredGridToVTK(p, X, Y, **TOPO), TypeError, "missing: z"),
    (lambda p: unstructuredGridToVTK(p, **TOPO), TypeError, "missing: x, y, z"),
    (lambda p: unstructuredGridToVTK(p, points=P), TypeError,
     "missing the arguments 'connectivity', 'offsets', 'cell_types'"),
    (lambda p: unstructuredGridToVTK(p, X, Y, Z), TypeError, "missing the arguments"),
    (lambda p: unstructuredGridToVTK(p, points=P, connectivity=TOPO["connectivity"],
                                     offsets=TOPO["offsets"]), TypeError,
     "missing the argument 'cell_types'"),
    (lambda p: polyLinesToVTK(p, points=P), TypeError, "missing the argument 'pointsPerLine'"),
    (lambda p: pointsToVTK(p, points=X), ValueError, "(n, 3) or (n, 2)"),
    (lambda p: pointsToVTK(p, points=np.zeros((6, 4))), ValueError, "(n, 3) or (n, 2)"),
    (lambda p: pointsToVTK(p, points=P.astype(bool)), TypeError, "'bool'"),
    (lambda p: TimeSeriesUnstructuredGrid(p, points=P), TypeError, "missing the arguments"),
    (lambda p: structuredToVTK(p, points=np.zeros((2, 2, 2))), ValueError, "(nx, ny, nz, 3)"),
])
def test_points_keyword_errors(tmp_path, call, exc, msg):
    with pytest.raises(exc) as info:
        call(str(tmp_path / "f"))
    assert msg in str(info.value)
    assert os.listdir(tmp_path) == []


# ---------------------------------------------------------------------------
# structured, rectilinear and image grids: (nx, ny, nz) and (nx, ny, nz, k) arrays
# ---------------------------------------------------------------------------
NX, NY, NZ = 4, 3, 2                          # points per direction; all different on purpose
I, J, K = np.meshgrid(np.arange(NX), np.arange(NY), np.arange(NZ), indexing="ij")
GX, GY, GZ = I + 0.1 * J * J, 2.0 * J + 0.05 * K, 0.5 * K + 0.01 * I    # a curved grid
GP = np.stack((GX, GY, GZ), axis=-1)                                     # (NX, NY, NZ, 3)
PV = np.stack((I + 0.0, 10.0 * J, 100.0 * K), axis=-1)                   # point vectors
PT = 1000.0 * I + 100.0 * J + 10.0 * K                                   # point scalars
CI, CJ, CK = np.meshgrid(np.arange(NX - 1), np.arange(NY - 1), np.arange(NZ - 1), indexing="ij")
GCV = np.stack((-CI - 0.5, CJ + 0.25, CK + 0.125), axis=-1)               # cell vectors
GCT = 7.0 * CI - CJ + 0.5 * CK                                            # cell scalars


def hand_interleaved(A):
    """ What users had to do before (examples/structured_ex1.py). """
    return np.moveaxis(A, -1, 0).reshape(-1, order="F")


def dimensions(out):
    dims = [0, 0, 0]
    out.GetDimensions(dims)
    return tuple(dims)


def check_grid_values(vtk, out, point_dims):
    """ Compare every entry with the input arrays, using VTK's (i, j, k) -> id. """
    cell_dims = [max(n - 1, 1) for n in point_dims]
    pd, cd = out.GetPointData(), out.GetCellData()
    v, t = vtk_array(vtk, pd, "v"), vtk_array(vtk, pd, "t")
    cv, ct = vtk_array(vtk, cd, "cv"), vtk_array(vtk, cd, "ct")
    for i, j, k in np.ndindex(*point_dims):
        pid = vtk.vtkStructuredData.ComputePointId(list(point_dims), [i, j, k])
        np.testing.assert_array_equal(v[pid], PV[i, j, k])
        assert t[pid] == PT[i, j, k]
    for i, j, k in np.ndindex(*cell_dims):
        cid = vtk.vtkStructuredData.ComputeCellId(list(point_dims), [i, j, k])
        np.testing.assert_array_equal(cv[cid], GCV[i, j, k])
        assert ct[cid] == GCT[i, j, k]
    assert pd.GetVectors().GetName() == "v" and cd.GetVectors().GetName() == "cv"


def test_structured_grid_arrays(tmp_path, vtk):
    new = structuredToVTK(str(tmp_path / "new"), points=GP,
                          all_point_data={"t": PT, "v": PV}, all_cell_data={"ct": GCT, "cv": GCV})
    old = structuredToVTK(str(tmp_path / "old"), GX, GY, GZ,
                          all_point_data=[["t", "scalars", PT], ["v", "vectors", hand_interleaved(PV)]],
                          all_cell_data=[["ct", "scalars", GCT], ["cv", "vectors", hand_interleaved(GCV)]])
    listed = structuredToVTK(str(tmp_path / "listed"), GX, GY, GZ,
                             all_point_data=[["t", "scalars", PT], ["v", "vectors", PV]],
                             all_cell_data=[["ct", "scalars", GCT], ["cv", "vectors", GCV]])
    assert same_files(new, old, listed)
    out = read_vtk(vtk, new)
    assert dimensions(out) == (NX, NY, NZ)
    from vtk.util.numpy_support import vtk_to_numpy
    pts = vtk_to_numpy(out.GetPoints().GetData())
    for i, j, k in np.ndindex(NX, NY, NZ):
        pid = vtk.vtkStructuredData.ComputePointId([NX, NY, NZ], [i, j, k])
        np.testing.assert_array_equal(pts[pid], GP[i, j, k])
    check_grid_values(vtk, out, (NX, NY, NZ))


def test_structured_points_2_components(tmp_path):
    a = structuredToVTK(str(tmp_path / "a"), points=GP[..., :2])
    b = structuredToVTK(str(tmp_path / "b"), GX, GY, np.zeros_like(GX))
    assert same_files(a, b)


def test_rectilinear_grid_arrays(tmp_path, vtk):
    xs, ys, zs = np.arange(NX) ** 2.0, np.arange(NY) - 0.5, np.arange(NZ) * 3.0
    path = rectilinearToVTK(str(tmp_path / "r"), xs, ys, zs,
                            all_point_data={"t": PT, "v": PV}, all_cell_data={"ct": GCT, "cv": GCV})
    check_grid_values(vtk, read_vtk(vtk, path), (NX, NY, NZ))


def test_image_arrays(tmp_path, vtk):
    # the image size comes from the first field, here a 4-D (vector) array
    path = imageToVTK(str(tmp_path / "i"), all_point_data={"v": PV, "t": PT},
                      all_cell_data={"cv": GCV, "ct": GCT})
    out = read_vtk(vtk, path)
    assert dimensions(out) == (NX, NY, NZ)
    check_grid_values(vtk, out, (NX, NY, NZ))
    cells_only = imageToVTK(str(tmp_path / "c"), all_cell_data={"cv": GCV})
    assert dimensions(read_vtk(vtk, cells_only)) == (NX, NY, NZ)


def test_grid_rows_are_vtk_point_ids(tmp_path):
    """ An (npoints, 3) array on a grid: row = VTK point id (i fastest). """
    a = structuredToVTK(str(tmp_path / "a"), points=GP, all_point_data={"v": PV})
    b = structuredToVTK(str(tmp_path / "b"), points=GP,
                        all_point_data={"v": hand_interleaved(PV).reshape(-1, 3)})
    assert same_files(a, b)


@pytest.mark.parametrize("fn, kw, exc, msg", [
    (structuredToVTK, dict(all_point_data={"v": np.swapaxes(PV, 0, 1)}), ValueError,
     "or (4, 3, 2) + (k,) for k components"),
    (structuredToVTK, dict(all_point_data={"v": PV[..., :2]}), ValueError, "give the type explicitly"),
    (structuredToVTK, dict(all_point_data={"T": np.zeros((NX, NY, NZ, 9))}), ValueError,
     "(tensors), which are not allowed here"),
    (structuredToVTK, dict(all_point_data={"v": (GX, GY, GZ)}), TypeError, "axis=-1"),
    (structuredToVTK, dict(all_point_data=[["t", "scalars", PT.T]]), ValueError,
     "expected (4, 3, 2)"),
    (rectilinearToVTK, dict(all_cell_data={"cv": np.moveaxis(GCV, -1, 0)}), ValueError,
     "or (3, 2, 1) + (k,)"),
])
def test_grid_errors(tmp_path, fn, kw, exc, msg):
    if fn is structuredToVTK:
        call = lambda: fn(str(tmp_path / "g"), points=GP, **kw)
    else:
        call = lambda: fn(str(tmp_path / "g"), np.arange(NX), np.arange(NY), np.arange(NZ), **kw)
    with pytest.raises(exc) as info:
        call()
    assert msg in str(info.value)
    assert os.listdir(tmp_path) == []


def test_image_errors_with_dicts(tmp_path):
    with pytest.raises(TypeError, match="axis=-1"):
        imageToVTK(str(tmp_path / "i"), all_point_data={"v": (PT, PT, PT)})
    with pytest.raises(ValueError, match="3-D arrays"):
        imageToVTK(str(tmp_path / "i"), all_point_data={"t": PT.ravel()})
    with pytest.raises(ValueError, match="does not match point data shape"):
        imageToVTK(str(tmp_path / "i"), all_point_data={"v": PV}, all_cell_data={"c": PT})
    assert os.listdir(tmp_path) == []
