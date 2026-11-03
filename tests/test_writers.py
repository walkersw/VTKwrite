"""Round trips for each high-level writer, read back with VTK."""
import numpy as np
import pytest

from VTKwrite.interface import (imageToVTK, rectilinearToVTK, structuredToVTK, pointsToVTK,
                                linesToVTK, polyLinesToVTK, unstructuredGridToVTK)
from VTKwrite.vtkbin import VtkTriangle, VtkQuad, VtkFile, VtkRectilinearGrid, VtkGroup
from conftest import read_vtk, vtk_array, no_leftovers

vtk_to_numpy = pytest.importorskip("vtk.util.numpy_support").vtk_to_numpy


def two_cells():
    x = np.array([0., 1., 2., 0., 1., 2.])
    y = np.array([0., 0., 0., 1., 1., 1.])
    z = np.zeros(6)
    conn = np.array([0, 1, 3, 1, 2, 5, 4], dtype=np.int64)
    offsets = np.array([3, 7], dtype=np.int64)
    types = np.array([VtkTriangle.tid, VtkQuad.tid], dtype=np.uint8)
    return x, y, z, conn, offsets, types


def test_unstructured_roundtrip(tmp_path, vtk):
    x, y, z, conn, offsets, types = two_cells()
    p = np.array([1.5, -2.0])
    v = np.arange(6.0)                      # 2 cells x 3 components, interleaved
    u = np.linspace(0, 1, 6)
    path = unstructuredGridToVTK(str(tmp_path / "ug"), x, y, z, conn, offsets, types,
                                 all_cell_data=[["p", "scalars", p], ["v", "vectors", v]],
                                 all_point_data=[["u", "scalars", u]], comments=["a comment"])
    out = read_vtk(vtk, path)
    assert out.GetNumberOfCells() == 2 and out.GetNumberOfPoints() == 6
    pts = vtk_to_numpy(out.GetPoints().GetData())
    np.testing.assert_array_equal(pts, np.column_stack((x, y, z)))
    assert [out.GetCellType(i) for i in range(2)] == [VtkTriangle.tid, VtkQuad.tid]
    np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "p"), p)
    np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "v"), v.reshape(2, 3))
    np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "u"), u)
    assert no_leftovers(tmp_path)


def test_active_attributes_are_set(tmp_path, vtk):
    """VTK matches 'Scalars'/'Vectors' case-sensitively; lowercase names were ignored."""
    x, y, z, conn, offsets, types = two_cells()
    path = unstructuredGridToVTK(str(tmp_path / "ug"), x, y, z, conn, offsets, types,
                                 all_point_data=[["u", "scalars", np.zeros(6)],
                                                 ["q", "vectors", np.zeros(18)]])
    pd = read_vtk(vtk, path).GetPointData()
    assert pd.GetScalars().GetName() == "u"
    assert pd.GetVectors().GetName() == "q"


def test_lists_are_accepted(tmp_path, vtk):
    x, y, z, conn, offsets, types = two_cells()
    path = unstructuredGridToVTK(str(tmp_path / "ug"), list(x), list(y), list(z), list(conn),
                                 list(offsets), list(types),
                                 all_cell_data=[("p", "scalars", [1.0, 2.0])])
    out = read_vtk(vtk, path)
    np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "p"), [1.0, 2.0])


def test_caller_lists_not_modified(tmp_path):
    x, y, z, conn, offsets, types = two_cells()
    data = [["p", "scalars", [1.0, 2.0]]]
    unstructuredGridToVTK(str(tmp_path / "ug"), x, y, z, conn, offsets, types, all_cell_data=data)
    assert data == [["p", "scalars", [1.0, 2.0]]]


def test_strided_arrays(tmp_path, vtk):
    x, y, z, conn, offsets, types = two_cells()
    big = np.zeros(12)
    big[::2] = x
    u = np.repeat(np.arange(6.0), 2)[::2]       # non-contiguous view
    path = unstructuredGridToVTK(str(tmp_path / "ug"), big[::2], y, z, conn, offsets, types,
                                 all_point_data=[["u", "scalars", u]])
    out = read_vtk(vtk, path)
    np.testing.assert_array_equal(vtk_to_numpy(out.GetPoints().GetData())[:, 0], x)
    np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "u"), np.arange(6.0))


def test_tcoords(tmp_path, vtk):
    x, y, z, conn, offsets, types = two_cells()
    uv = np.arange(12.0)
    path = unstructuredGridToVTK(str(tmp_path / "ug"), x, y, z, conn, offsets, types,
                                 all_point_data=[["uv", "tcoords", uv]])
    pd = read_vtk(vtk, path).GetPointData()
    np.testing.assert_array_equal(vtk_array(vtk, pd, "uv"), uv.reshape(6, 2))
    assert pd.GetTCoords().GetName() == "uv"


def test_image_roundtrip(tmp_path, vtk):
    rng = np.random.default_rng(0)
    pressure = rng.random((3, 4, 2))
    temp = rng.random((4, 5, 3))
    path = imageToVTK(str(tmp_path / "img"), origin=(1.0, 2.0, 3.0), spacing=(0.5, 0.25, 1.0),
                      all_cell_data=[["p", "scalars", pressure]],
                      all_point_data=[["T", "scalars", temp]])
    out = read_vtk(vtk, path)
    assert out.GetDimensions() == (4, 5, 3)
    assert out.GetOrigin() == (1.0, 2.0, 3.0) and out.GetSpacing() == (0.5, 0.25, 1.0)
    np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "p"), pressure.ravel(order="F"))
    np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "T"), temp.ravel(order="F"))


@pytest.mark.parametrize("origin", [(np.float64(0.5), np.float32(1.25), 2),
                                    np.array([0.5, 1.25, 2.0])])
def test_image_numpy_origin(tmp_path, vtk, origin):
    """NumPy 2 repr() of np.float64 gave Origin="np.float64(0.5) ..." (read as 0);
    an ndarray origin raised 'truth value of an array is ambiguous'."""
    path = imageToVTK(str(tmp_path / "img"), origin=origin, spacing=np.array([1.0, 1.0, 1.0]),
                      all_point_data=[["T", "scalars", np.zeros((2, 2, 2))]])
    assert b"np." not in open(path, "rb").read()
    assert read_vtk(vtk, path).GetOrigin() == (0.5, 1.25, 2.0)


def test_rectilinear_roundtrip(tmp_path, vtk):
    x, y, z = np.linspace(0, 1, 4), np.linspace(0, 2, 3), np.array([0.0, 0.5])
    p = np.arange(6.0).reshape(3, 2, 1)
    path = rectilinearToVTK(str(tmp_path / "rect"), x, y, z, all_cell_data=[["p", "scalars", p]])
    out = read_vtk(vtk, path)
    np.testing.assert_array_equal(vtk_to_numpy(out.GetXCoordinates()), x)
    np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "p"), p.ravel(order="F"))


def test_structured_roundtrip(tmp_path, vtk):
    X, Y, Z = np.meshgrid(np.linspace(0, 1, 3), np.linspace(0, 1, 4), np.linspace(0, 1, 2),
                          indexing="ij")
    T = X + 10 * Y + 100 * Z
    path = structuredToVTK(str(tmp_path / "s"), X, Y, Z, all_point_data=[["T", "scalars", T]])
    out = read_vtk(vtk, path)
    assert out.GetExtent() == (0, 2, 0, 3, 0, 1)
    np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "T"), T.ravel(order="F"))


def test_points_lines_polylines(tmp_path, vtk):
    x, y, z = np.arange(6.0), np.zeros(6), np.zeros(6)
    out = read_vtk(vtk, pointsToVTK(str(tmp_path / "pts"), x, y, z,
                                    all_point_data=[["u", "scalars", x]]))
    assert out.GetNumberOfCells() == 6
    out = read_vtk(vtk, linesToVTK(str(tmp_path / "lines"), x, y, z,
                                   all_cell_data=[["c", "scalars", np.arange(3.0)]]))
    assert out.GetNumberOfCells() == 3
    np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "c"), np.arange(3.0))
    out = read_vtk(vtk, polyLinesToVTK(str(tmp_path / "pl"), x, y, z, [2, 4]))
    assert out.GetNumberOfCells() == 2


def test_special_characters_in_names_and_comments(tmp_path, vtk):
    x, y, z, conn, offsets, types = two_cells()
    name = 'p&q "1" <a>'
    path = unstructuredGridToVTK(str(tmp_path / "ug"), x, y, z, conn, offsets, types,
                                 all_cell_data=[[name, "scalars", np.zeros(2)],
                                                ["température", "scalars", np.ones(2)]],
                                 comments=["h = 1e-2 -- refined", "ends with -"])
    cd = read_vtk(vtk, path).GetCellData()
    assert {cd.GetArrayName(i) for i in range(cd.GetNumberOfArrays())} == {name, "température"}


def test_low_level_context_manager(tmp_path, vtk):
    nx = 2
    x = y = z = np.linspace(0, 1, nx + 1)
    with VtkFile(str(tmp_path / "ll"), VtkRectilinearGrid) as w:
        w.openGrid(start=(0, 0, 0), end=(nx, nx, nx))
        w.openPiece(start=(0, 0, 0), end=(nx, nx, nx))
        w.openElement("Coordinates")
        for name, c in (("x", x), ("y", y), ("z", z)):
            w.addData(name, c)
        w.closeElement("Coordinates")
        w.closePiece()
        w.closeGrid()
        w.appendData(x).appendData(y).appendData(z)
    assert read_vtk(vtk, str(tmp_path / "ll.vtr")).GetNumberOfPoints() == 27


def test_group(tmp_path):
    with VtkGroup(str(tmp_path / "g")) as g:
        g.addFile(str(tmp_path / "a&b.vtu"), sim_time=np.float64(0.5))
    text = open(tmp_path / "g.pvd").read()
    assert 'file="a&amp;b.vtu"' in text and 'timestep="0.5"' in text


def test_flat_grids(tmp_path, vtk):
    """A direction with a single point still has one layer of cells (as in VTK)."""
    x, y, z = np.linspace(0, 1, 4), np.linspace(0, 1, 3), np.array([0.0])
    p = np.arange(6.0).reshape(3, 2, 1)
    out = read_vtk(vtk, rectilinearToVTK(str(tmp_path / "r"), x, y, z,
                                         all_cell_data=[["p", "scalars", p]]))
    assert out.GetNumberOfCells() == 6
    np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "p"), p.ravel(order="F"))

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    out = read_vtk(vtk, structuredToVTK(str(tmp_path / "s"), X, Y, Z,
                                        all_cell_data=[["p", "scalars", p.ravel(order="F")]]))
    assert out.GetNumberOfCells() == 6

    out = read_vtk(vtk, imageToVTK(str(tmp_path / "i"),
                                   all_cell_data=[["p", "scalars", p]],
                                   all_point_data=[["T", "scalars", np.ones((4, 3, 1))],
                                                   ["U", "scalars", np.ones(12)]]))
    assert out.GetDimensions() == (4, 3, 1) and out.GetNumberOfCells() == 6


def test_multidimensional_coordinates_for_unstructured(tmp_path, vtk):
    X, Y, Z = np.meshgrid(np.arange(2.0), np.arange(2.0), np.arange(2.0), indexing="ij")
    conn = np.array([0, 1, 3, 2, 4, 5, 7, 6])       # Fortran-order point numbering
    out = read_vtk(vtk, unstructuredGridToVTK(str(tmp_path / "hex"), X, Y, Z, conn,
                                              np.array([8]), np.array([12])))
    assert out.GetNumberOfPoints() == 8 and out.GetCellType(0) == 12
    np.testing.assert_array_equal(vtk_to_numpy(out.GetPoints().GetData())[:, 0],
                                  X.ravel(order="F"))


def test_mixed_coordinate_dtypes(tmp_path, vtk):
    x, y, z, conn, offsets, types = two_cells()
    out = read_vtk(vtk, unstructuredGridToVTK(str(tmp_path / "ug"), x.astype(np.int64), y,
                                              z.astype(np.float32), conn, offsets, types))
    np.testing.assert_array_equal(vtk_to_numpy(out.GetPoints().GetData()),
                                  np.column_stack((x, y, z)))


def test_symlinked_target_is_written_through(tmp_path, vtk):
    x, y, z, conn, offsets, types = two_cells()
    real = tmp_path / "real"
    real.mkdir()
    (real / "ug.vtu").write_text("old")
    link = tmp_path / "ug.vtu"
    try:
        link.symlink_to(real / "ug.vtu")
    except (OSError, NotImplementedError):
        pytest.skip("symbolic links not available")
    unstructuredGridToVTK(str(tmp_path / "ug"), x, y, z, conn, offsets, types)
    assert link.is_symlink()
    assert read_vtk(vtk, str(real / "ug.vtu")).GetNumberOfCells() == 2


def test_long_file_name(tmp_path, vtk):
    x, y, z, conn, offsets, types = two_cells()
    name = "a" * 200
    path = unstructuredGridToVTK(str(tmp_path / name), x, y, z, conn, offsets, types)
    assert read_vtk(vtk, path).GetNumberOfCells() == 2
