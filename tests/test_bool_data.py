"""bool data arrays (e.g. boundary markers) are written as UInt8 0/1."""
import os

import numpy as np
import pytest

from VTKwrite import (unstructuredGridToVTK, imageToVTK, pointsToVTK, structuredToVTK,
                      TimeSeriesUnstructuredGrid, timeseries_unstructuredGrid,
                      VtkTriangle, VtkQuad)
from conftest import read_appended, read_vtk, vtk_array, no_leftovers

X = np.array([0., 1., 2., 0., 1., 2.])
Y = np.array([0., 0., 0., 1., 1., 1.])
Z = np.zeros(6)
MESH = (X, Y, Z, [0, 1, 3, 1, 2, 5, 4], [3, 7], [VtkTriangle.tid, VtkQuad.tid])
ON_BOUNDARY = np.array([True, True, True, True, False, True])
CELL_FLAG = np.array([False, True])


def header_type(path, name):
    """The type="..." attribute of the DataArray called name."""
    with open(path, "rb") as f:
        text = f.read().split(b"<AppendedData")[0].decode()
    tag = [t for t in text.split("<DataArray")[1:] if 'Name="%s"' % name in t][0]
    return tag.split('type="')[1].split('"')[0]


def test_unstructured(tmp_path, vtk):
    path = unstructuredGridToVTK(str(tmp_path / "ug"), *MESH,
                                 all_cell_data=[["flag", "scalars", CELL_FLAG]],
                                 all_point_data=[["bnd", "scalars", list(ON_BOUNDARY)]])
    assert header_type(path, "bnd") == "UInt8" and header_type(path, "flag") == "UInt8"
    out = read_vtk(vtk, path)
    bnd = out.GetPointData().GetArray("bnd")
    assert bnd.GetDataType() == vtk.VTK_UNSIGNED_CHAR
    np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "bnd"), ON_BOUNDARY)
    np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "flag"), CELL_FLAG)


def test_image_and_structured_3d_arrays(tmp_path, vtk):
    mask = np.zeros((3, 4, 2), dtype=bool)
    mask[1, 2, 0] = mask[2, 3, 1] = True
    path = imageToVTK(str(tmp_path / "img"), all_point_data=[["mask", "scalars", mask]])
    np.testing.assert_array_equal(vtk_array(vtk, read_vtk(vtk, path).GetPointData(), "mask"),
                                  mask.ravel(order="F"))
    x, y, z = np.meshgrid(np.arange(3.0), np.arange(4.0), np.arange(2.0), indexing="ij")
    path = structuredToVTK(str(tmp_path / "sg"), x, y, z,
                           all_point_data=[["mask", "scalars", mask]])
    np.testing.assert_array_equal(vtk_array(vtk, read_vtk(vtk, path).GetPointData(), "mask"),
                                  mask.ravel(order="F"))


def test_points(tmp_path):
    path = pointsToVTK(str(tmp_path / "pts"), X, Y, Z, all_point_data=[["b", "scalars", ON_BOUNDARY]])
    blocks = read_appended(path)
    assert blocks[("b", None)].dtype == np.uint8
    np.testing.assert_array_equal(blocks[("b", None)], ON_BOUNDARY.astype(np.uint8))


@pytest.mark.parametrize("output", ["vtu", "pvd", "vtkhdf"])
def test_time_series(tmp_path, vtk, output):
    if output == "vtkhdf":
        pytest.importorskip("h5py")
    path = tmp_path / "sol"
    with TimeSeriesUnstructuredGrid(str(path), *MESH, output=output,
                                    static_point_data={"bnd": ON_BOUNDARY}) as ts:
        for k in range(3):
            ts.write_step(0.5 * k, point_data={"wet": X > k - 0.5},
                          cell_data={"flag": CELL_FLAG ^ bool(k % 2)})
    wet = [(X > k - 0.5).astype(np.uint8) for k in range(3)]
    flag = [(CELL_FLAG ^ bool(k % 2)).astype(np.uint8) for k in range(3)]
    if output == "vtu":
        blocks = read_appended(str(path) + ".vtu")
        assert blocks[("bnd", None)].dtype == np.uint8
        for k in range(3):
            np.testing.assert_array_equal(blocks[("wet", str(k))], wet[k])
            np.testing.assert_array_equal(blocks[("flag", str(k))], flag[k])
        np.testing.assert_array_equal(blocks[("bnd", None)], ON_BOUNDARY)
    elif output == "pvd":
        for k in range(3):
            out = read_vtk(vtk, os.path.join(str(path), "sol_%05d.vtu" % k))
            np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "wet"), wet[k])
            np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "flag"), flag[k])
            np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "bnd"), ON_BOUNDARY)
    else:
        r = vtk.vtkHDFReader()
        r.SetFileName(str(path) + ".vtkhdf")
        for k, t in enumerate([0.0, 0.5, 1.0]):
            r.UpdateTimeStep(t)
            pd, cd = r.GetOutput().GetPointData(), r.GetOutput().GetCellData()
            assert pd.GetArray("wet").GetDataType() == vtk.VTK_UNSIGNED_CHAR
            np.testing.assert_array_equal(vtk_array(vtk, pd, "wet"), wet[k])
            np.testing.assert_array_equal(vtk_array(vtk, cd, "flag"), flag[k])
            np.testing.assert_array_equal(vtk_array(vtk, pd, "bnd"), ON_BOUNDARY)
    assert no_leftovers(tmp_path)


def test_bool_step_after_float_first_step(tmp_path):
    """bool -> UInt8 -> the dtype of the first step (a safe conversion)."""
    with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH, times=[0.0, 1.0]) as ts:
        ts.write_step(point_data={"u": np.full(6, 0.5)})
        ts.write_step(point_data={"u": ON_BOUNDARY})
    blocks = read_appended(str(tmp_path / "sol.vtu"))
    assert blocks[("u", "1")].dtype == np.float64
    np.testing.assert_array_equal(blocks[("u", "1")], ON_BOUNDARY.astype(float))


@pytest.mark.filterwarnings("ignore:timeseries_unstructuredGrid is deprecated:DeprecationWarning")
def test_legacy_time_series(tmp_path):
    ts = timeseries_unstructuredGrid(str(tmp_path / "ts"), [0.0, 1.0])
    ts.init_unstructuredGridToVTK(*MESH, all_point_data=[["wet", "scalars", ON_BOUNDARY]])
    ts.append_data([ON_BOUNDARY, ~ON_BOUNDARY])
    blocks = read_appended(ts.close_unstructuredGridToVTK())
    np.testing.assert_array_equal(blocks[("wet", "0")], ON_BOUNDARY)
    np.testing.assert_array_equal(blocks[("wet", "1")], ~ON_BOUNDARY)
    assert blocks[("wet", "1")].dtype == np.uint8


def test_bool_coordinates_and_topology_still_rejected(tmp_path):
    with pytest.raises(TypeError, match="x: .*'bool'"):
        unstructuredGridToVTK(str(tmp_path / "a"), ON_BOUNDARY, Y, Z, *MESH[3:])
    with pytest.raises(TypeError, match="connectivity: .*'bool'"):
        unstructuredGridToVTK(str(tmp_path / "b"), X, Y, Z, np.ones(7, dtype=bool), *MESH[4:])
    assert os.listdir(tmp_path) == []
