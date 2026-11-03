"""Single-file time series (timeseries_unstructuredGrid, deprecated; kept until it is
removed)."""
import os
import warnings

import numpy as np
import pytest

from VTKwrite.interface import timeseries_unstructuredGrid
from VTKwrite.vtkbin import VtkTriangle
from conftest import read_appended, read_vtk, vtk_array, no_leftovers

TIMES = np.array([0.0, 0.1, 0.25])

# the class is deprecated; its tests run without the warning (checked once below)
pytestmark = pytest.mark.filterwarnings(
    "ignore:timeseries_unstructuredGrid is deprecated:DeprecationWarning")


def test_deprecation_warning():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        timeseries_unstructuredGrid("unused", TIMES)          # writes nothing
    assert len(caught) == 1 and caught[0].category is DeprecationWarning
    assert "TimeSeriesUnstructuredGrid" in str(caught[0].message)
    assert caught[0].filename == __file__                    # points at the caller


def start(tmp_path, name="ts"):
    x, y, z = np.array([0., 1., 0.]), np.array([0., 0., 1.]), np.zeros(3)
    ts = timeseries_unstructuredGrid(str(tmp_path / name), TIMES)
    ts.init_unstructuredGridToVTK(x, y, z, np.array([0, 1, 2]), np.array([3]),
                                  np.array([VtkTriangle.tid], dtype=np.uint8),
                                  all_cell_data=[["p", "scalars", np.zeros(1)]],
                                  all_point_data=[["u", "scalars", np.zeros(3)],
                                                  ["q", "vectors", np.zeros(9)]])
    return ts


def p_at(k): return np.array([10.0 + k])
def u_at(k): return np.full(3, float(k))
def q_at(k): return np.arange(9.0) + 100 * k


def test_values_per_step(tmp_path):
    ts = start(tmp_path)
    ts.append_data([p_at(k) for k in range(3)])
    ts.append_data([u_at(k) for k in range(3)])
    ts.append_data([q_at(k) for k in range(3)])
    path = ts.close_unstructuredGridToVTK()
    blocks = read_appended(path)
    for k in range(3):
        np.testing.assert_array_equal(blocks[("p", str(k))], p_at(k))
        np.testing.assert_array_equal(blocks[("u", str(k))], u_at(k))
        np.testing.assert_array_equal(blocks[("q", str(k))], q_at(k))
    assert no_leftovers(tmp_path)


def test_values_per_step_vtk(tmp_path, vtk, request):
    ts = start(tmp_path)
    ts.append_data([p_at(k) for k in range(3)])
    ts.append_data([u_at(k) for k in range(3)])
    ts.append_data([q_at(k) for k in range(3)])
    path = ts.close_unstructuredGridToVTK()
    version = tuple(int(v) for v in vtk.vtkVersion.GetVTKVersion().split(".")[:2])
    # the checks below run in any case, so that a fixed VTK makes this test fail (strict)
    request.applymarker(pytest.mark.xfail(version >= (9, 6), strict=True, reason=(
        "VTK >= 9.6 reads only the first time step of a single-file .vtu time series "
               "(upstream regression, VTK issue #20227). If this test passes, VTK is fixed: "
               "turn the version check into an upper bound, and update the README section "
               "'Known issue'")))
    for k in range(3):
        out = read_vtk(vtk, path, time_step=k)
        np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "p"), p_at(k))
        np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "u"), u_at(k))


def test_lossless_dtype_conversion(tmp_path):
    ts = start(tmp_path)
    ts.append_data([p_at(0), p_at(1).astype(np.float32), [12.0]])
    ts.append_data([u_at(k) for k in range(3)])
    ts.append_data([np.arange(9) + 100 * k for k in range(3)])   # int64 -> float64
    blocks = read_appended(ts.close_unstructuredGridToVTK())
    assert blocks[("p", "1")].dtype == np.dtype("<f8")
    np.testing.assert_array_equal(blocks[("p", "2")], [12.0])
    np.testing.assert_array_equal(blocks[("q", "2")], q_at(2))


def test_lossy_dtype_rejected(tmp_path):
    x, y, z = np.array([0., 1., 0.]), np.array([0., 0., 1.]), np.zeros(3)
    ts = timeseries_unstructuredGrid(str(tmp_path / "ts"), TIMES)
    ts.init_unstructuredGridToVTK(x, y, z, [0, 1, 2], [3], [5],
                                  all_point_data=[["u", "scalars", np.zeros(3, np.float32)]])
    with pytest.raises(TypeError, match="cannot be safely converted"):
        ts.append_data([np.zeros(3), np.zeros(3), np.zeros(3)])   # float64 -> float32


@pytest.mark.parametrize("bad, msg", [
    ([p_at(0), p_at(1)], "expected one per time value"),
    ([p_at(0), p_at(1), np.zeros(2)], "2 values, expected 1"),
])
def test_wrong_step_data(tmp_path, bad, msg):
    ts = start(tmp_path)
    with pytest.raises(ValueError, match=msg):
        ts.append_data(bad)
    # nothing was written by the failed call: the right data can still be appended
    ts.append_data([p_at(k) for k in range(3)])
    ts.append_data([u_at(k) for k in range(3)])
    ts.append_data([q_at(k) for k in range(3)])
    np.testing.assert_array_equal(read_appended(ts.close_unstructuredGridToVTK())[("p", "2")],
                                  p_at(2))


def test_too_many_appends(tmp_path):
    ts = start(tmp_path)
    for f in (p_at, u_at, q_at):
        ts.append_data([f(k) for k in range(3)])
    with pytest.raises(ValueError, match="more times than there are variables"):
        ts.append_data([p_at(k) for k in range(3)])


def test_close_with_missing_variable(tmp_path):
    ts = start(tmp_path)
    ts.append_data([p_at(k) for k in range(3)])
    with pytest.raises(ValueError, match="point 'u', point 'q'"):
        ts.close_unstructuredGridToVTK()
    assert os.listdir(tmp_path) == []


def test_fake_time_grid(tmp_path):
    ts = start(tmp_path)
    for f in (p_at, u_at, q_at):
        ts.append_data([f(k) for k in range(3)])
    ts.close_unstructuredGridToVTK()
    path = ts.write_fake_grid_for_time(str(tmp_path / "times"))
    text = open(path).read()
    for k, t in enumerate(TIMES):
        assert 'TimeStep="%d">%r</DataArray>' % (k, float(t)) in text
    assert no_leftovers(tmp_path)
