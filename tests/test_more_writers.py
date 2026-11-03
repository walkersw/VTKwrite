"""Writers fixed after the first review pass: cylinder, TIN, lines/polylines, topology dtypes."""
import numpy as np
import pytest

from VTKwrite.interface import (cylinderToVTK, pointsToVTKAsTIN, linesToVTK, polyLinesToVTK,
                                unstructuredGridToVTK)
from VTKwrite.vtkbin import VtkQuad, VtkTriangle, VtkTetra, VtkLine, VtkPolyLine
from conftest import read_appended, read_vtk, vtk_array


def test_cylinder(tmp_path, vtk):
    nl, npil = 3, 8
    path = cylinderToVTK(str(tmp_path / "cyl"), 1.0, 2.0, 0.0, 3.0, 0.5, nl, npil,
                         all_cell_data=[["c", "scalars", np.arange(nl * npil, dtype=float)]],
                         pointData={"h": np.arange(npil * (nl + 1), dtype=float)})
    out = read_vtk(vtk, path)
    assert out.GetNumberOfCells() == nl * npil
    assert out.GetNumberOfPoints() == npil * (nl + 1)
    assert {out.GetCellType(i) for i in range(out.GetNumberOfCells())} == {VtkQuad.tid}
    from vtk.util.numpy_support import vtk_to_numpy
    pts = vtk_to_numpy(out.GetPoints().GetData())
    np.testing.assert_allclose(np.hypot(pts[:, 0] - 1.0, pts[:, 1] - 2.0), 0.5)
    np.testing.assert_allclose(sorted(set(np.round(pts[:, 2], 12))), [0.0, 1.0, 2.0, 3.0])
    # every quad is a proper rectangle on the surface: total area = lateral area of the prism
    area = vtk.vtkMeshQuality()
    area.SetInputData(out)
    area.SetQuadQualityMeasureToArea()
    area.Update()
    a = vtk_to_numpy(area.GetOutput().GetCellData().GetArray("Quality"))
    side = 2 * 0.5 * np.sin(np.pi / npil)
    np.testing.assert_allclose(a.sum(), npil * side * 3.0)
    np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "h"),
                                  np.arange(npil * (nl + 1), dtype=float))


def test_cylinder_older_positional_dicts(tmp_path):
    path = cylinderToVTK(str(tmp_path / "cyl"), 0, 0, 0, 1, 1, 2, 4,
                         {"c": np.zeros(8)}, {"h": np.zeros(12)})
    blocks = read_appended(path)
    assert ("c", None) in blocks and ("h", None) in blocks


@pytest.mark.parametrize("ndim, tid, nverts", [(2, VtkTriangle.tid, 3), (3, VtkTetra.tid, 4)])
def test_tin(tmp_path, vtk, ndim, tid, nverts):
    pytest.importorskip("scipy")
    rng = np.random.default_rng(1)
    x, y, z = rng.random(30), rng.random(30), rng.random(30)
    path = pointsToVTKAsTIN(str(tmp_path / "tin"), x, y, z, ndim=ndim,
                            data=[["T", "scalars", x + y]], comments=["tin"])
    assert path is not None and path.endswith(".vtu")
    blocks = read_appended(path)
    assert blocks[("connectivity", None)].dtype.kind == "i"
    assert blocks[("types", None)].dtype == np.uint8
    out = read_vtk(vtk, path)
    assert {out.GetCellType(i) for i in range(out.GetNumberOfCells())} == {tid}
    assert out.GetCell(0).GetNumberOfPoints() == nverts
    np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "Elevation"), z)
    np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "T"), x + y)
    assert b"<!-- tin -->" in open(path, "rb").read()


def test_one_cell_type_per_line(tmp_path):
    x = np.arange(6.0)
    blocks = read_appended(linesToVTK(str(tmp_path / "l"), x, x, x))
    np.testing.assert_array_equal(blocks[("types", None)], [VtkLine.tid] * 3)
    blocks = read_appended(polyLinesToVTK(str(tmp_path / "pl"), x, x, x, [2, 4]))
    np.testing.assert_array_equal(blocks[("types", None)], [VtkPolyLine.tid] * 2)


def test_float_topology_written_as_integers(tmp_path):
    x = np.array([0., 1., 0.])
    blocks = read_appended(unstructuredGridToVTK(
        str(tmp_path / "ug"), x, x[::-1], np.zeros(3),
        np.array([0., 1., 2.]), np.array([3.]), np.array([5.])))
    assert blocks[("connectivity", None)].dtype == np.dtype("<i8")
    assert blocks[("offsets", None)].dtype == np.dtype("<i8")
    np.testing.assert_array_equal(blocks[("types", None)], np.array([5], dtype=np.uint8))


def test_cell_type_range_checked(tmp_path):
    x = np.array([0., 1., 0.])
    with pytest.raises(ValueError, match="between 0 and 255"):
        unstructuredGridToVTK(str(tmp_path / "ug"), x, x, x, [0, 1, 2], [3], [300])
