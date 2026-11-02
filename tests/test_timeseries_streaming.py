"""TimeSeriesUnstructuredGrid: streaming time series, single .vtu or .pvd."""
import glob
import os
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from VTKwrite.interface import TimeSeriesUnstructuredGrid
from VTKwrite.vtkbin import VtkTriangle, VtkQuad
from conftest import read_appended, read_vtk, vtk_array

X = np.array([0., 1., 2., 0., 1., 2.])
Y = np.array([0., 0., 0., 1., 1., 1.])
Z = np.zeros(6)
MESH = (X, Y, Z, [0, 1, 3, 1, 2, 5, 4], [3, 7], [VtkTriangle.tid, VtkQuad.tid])
TIMES = [0.0, 0.1, 0.25, 0.7]
STATIC = dict(static_cell_data={"material": np.array([7, 8], dtype=np.int32)},
              static_point_data={"h": np.linspace(0, 1, 6)})


def u_at(k): return np.full(6, float(k))
def q_at(k): return np.arange(18.0).reshape(6, 3) + 100 * k      # one row per point
def p_at(k): return np.array([10.0 + k, 20.0 + k])


def write_all(path, steps=4, **kw):
    with TimeSeriesUnstructuredGrid(str(path), *MESH, **STATIC, **kw) as ts:
        for k in range(steps):
            ts.write_step(TIMES[k], point_data={"u": u_at(k), "q": q_at(k)},
                          cell_data={"p": p_at(k)})
    return ts


def leftovers(directory):
    return [f for f in os.listdir(directory) if f.endswith((".part", ".steps"))]


def test_single_vtu_values(tmp_path):
    write_all(tmp_path / "sol", times=TIMES)
    blocks = read_appended(str(tmp_path / "sol.vtu"))
    for k in range(4):
        np.testing.assert_array_equal(blocks[("u", str(k))], u_at(k))
        np.testing.assert_array_equal(blocks[("q", str(k))], q_at(k).ravel())
        np.testing.assert_array_equal(blocks[("p", str(k))], p_at(k))
    np.testing.assert_array_equal(blocks[("material", None)], [7, 8])
    np.testing.assert_array_equal(blocks[("h", None)], np.linspace(0, 1, 6))
    text = open(tmp_path / "sol.vtu", "rb").read()
    assert b'TimeValues="0.0 0.1 0.25 0.7"' in text
    assert not leftovers(tmp_path)


def test_known_and_unknown_times_give_the_same_file(tmp_path):
    write_all(tmp_path / "a", times=TIMES)
    write_all(tmp_path / "b")
    assert open(tmp_path / "a.vtu", "rb").read() == open(tmp_path / "b.vtu", "rb").read()
    assert not leftovers(tmp_path)


def test_single_vtu_read_by_vtk(tmp_path, vtk):
    write_all(tmp_path / "sol")
    version = tuple(int(v) for v in vtk.vtkVersion.GetVTKVersion().split(".")[:2])
    if version >= (9, 6):
        pytest.xfail("VTK >= 9.6 drops all but the first time step (upstream regression)")
    for k in range(4):
        out = read_vtk(vtk, str(tmp_path / "sol.vtu"), time_step=k)
        pd, cd = out.GetPointData(), out.GetCellData()
        np.testing.assert_array_equal(vtk_array(vtk, pd, "u"), u_at(k))
        np.testing.assert_array_equal(vtk_array(vtk, pd, "q"), q_at(k))
        np.testing.assert_array_equal(vtk_array(vtk, cd, "p"), p_at(k))
        np.testing.assert_array_equal(vtk_array(vtk, cd, "material"), [7, 8])
        np.testing.assert_array_equal(vtk_array(vtk, pd, "h"), np.linspace(0, 1, 6))


def test_time_values_file(tmp_path):
    write_all(tmp_path / "sol")
    text = open(tmp_path / "sol_time_values.vtu").read()
    for k, t in enumerate(TIMES):
        assert 'TimeStep="%d">%r</DataArray>' % (k, t) in text
    write_all(tmp_path / "nofile", time_values_file=False)
    assert not os.path.exists(tmp_path / "nofile_time_values.vtu")
    write_all(tmp_path / "named", time_values_file=str(tmp_path / "times"))
    assert os.path.exists(tmp_path / "times.vtu")


def test_pvd_output(tmp_path, vtk):
    ts = write_all(tmp_path / "sol", output="pvd")
    root = ET.parse(tmp_path / "sol.pvd").getroot()
    entries = [(float(d.get("timestep")), d.get("file")) for d in root.iter("DataSet")]
    assert [t for t, _ in entries] == TIMES
    assert [f for _, f in entries] == ["sol/sol_%05d.vtu" % k for k in range(4)]
    for k, (t, f) in enumerate(entries):
        out = read_vtk(vtk, str(tmp_path / f))
        np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "u"), u_at(k))
        np.testing.assert_array_equal(vtk_array(vtk, out.GetPointData(), "q"), q_at(k))
        np.testing.assert_array_equal(vtk_array(vtk, out.GetCellData(), "material"), [7, 8])
    assert ts.nsteps == 4
    assert not leftovers(tmp_path) and not leftovers(tmp_path / "sol")


def test_pvd_keeps_finished_steps_after_an_error(tmp_path):
    with pytest.raises(RuntimeError, match="solver failed"):
        with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH, output="pvd") as ts:
            for k in range(4):
                if k == 2:
                    raise RuntimeError("solver failed")
                ts.write_step(TIMES[k], point_data={"u": u_at(k)})
    root = ET.parse(tmp_path / "sol.pvd").getroot()
    assert [float(d.get("timestep")) for d in root.iter("DataSet")] == TIMES[:2]
    assert sorted(os.listdir(tmp_path / "sol")) == ["sol_00000.vtu", "sol_00001.vtu"]


def test_pvd_in_current_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with TimeSeriesUnstructuredGrid("sol", *MESH, output="pvd") as ts:
        ts.write_step(0.0, point_data={"u": u_at(0)})
    assert 'file="sol/sol_00000.vtu"' in open("sol.pvd").read()


def test_error_in_vtu_mode_leaves_nothing(tmp_path):
    for kw in (dict(times=TIMES), {}):
        with pytest.raises(RuntimeError):
            with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH, **kw) as ts:
                ts.write_step(TIMES[0], point_data={"u": u_at(0)})
                raise RuntimeError("stop")
        assert os.listdir(tmp_path) == []


def test_steps_reach_the_disk_as_they_are_written(tmp_path):
    n = 50_000
    x = np.random.default_rng(0).random(n)
    with TimeSeriesUnstructuredGrid(str(tmp_path / "big"), x, x, x, np.arange(n),
                                    np.arange(1, n + 1), np.ones(n, np.uint8),
                                    times=[0.0, 1.0, 2.0]) as ts:
        sizes = []
        for k in range(3):
            ts.write_step(point_data={"u": x + k})
            sizes.append(sum(os.path.getsize(f) for f in glob.glob(str(tmp_path / ".big*"))))
        assert sizes[1] - sizes[0] >= 8 * n - 2**16 and sizes[2] - sizes[1] >= 8 * n - 2**16


def test_2d_arrays_and_explicit_types(tmp_path):
    with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH) as ts:
        ts.write_step(0.0, point_data={"q": q_at(0), "n": ("normals", q_at(0)),
                                       "uv": ("tcoords", np.zeros((6, 2)))})
    text = open(tmp_path / "sol.vtu", "rb").read()
    assert b'Vectors="q"' in text and b'Normals="n"' in text and b'TCoords="uv"' in text
    blocks = read_appended(str(tmp_path / "sol.vtu"))
    np.testing.assert_array_equal(blocks[("q", "0")], q_at(0).ravel())


def test_lossless_dtype_conversion(tmp_path):
    with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH) as ts:
        ts.write_step(0.0, point_data={"u": u_at(0)})
        ts.write_step(1.0, point_data={"u": u_at(1).astype(np.float32)})
        ts.write_step(2.0, point_data={"u": np.full(6, 2)})                 # int64 -> float64
    blocks = read_appended(str(tmp_path / "sol.vtu"))
    assert blocks[("u", "2")].dtype == np.dtype("<f8")
    np.testing.assert_array_equal(blocks[("u", "2")], u_at(2))


@pytest.mark.parametrize("step, exc, msg", [
    (dict(point_data={"u": u_at(1), "w": u_at(1)}), ValueError, "unexpected: w"),
    (dict(point_data={}), ValueError, "missing: u"),
    (dict(point_data={"u": q_at(1)}), ValueError, "the first step had scalars"),
    (dict(point_data={"u": np.arange(6) + 0j}), TypeError, "not supported"),
    (dict(point_data={"u": np.zeros((6, 2))}), ValueError, "give the type explicitly"),
])
def test_step_errors(tmp_path, step, exc, msg):
    with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH) as ts:
        ts.write_step(0.0, point_data={"u": u_at(0)})
        with pytest.raises(exc, match=msg):
            ts.write_step(1.0, **step)
        ts.write_step(1.0, point_data={"u": u_at(1)})       # a failed step writes nothing
    blocks = read_appended(str(tmp_path / "sol.vtu"))
    assert set(blocks) >= {("u", "0"), ("u", "1")} and ("u", "2") not in blocks


def test_unsafe_dtype_rejected(tmp_path):
    with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH) as ts:
        ts.write_step(0.0, point_data={"u": u_at(0).astype(np.float32)})
        with pytest.raises(TypeError, match="cannot be safely converted"):
            ts.write_step(1.0, point_data={"u": u_at(1)})
        ts.write_step(1.0, point_data={"u": u_at(1).astype(np.float32)})


def test_time_checks(tmp_path):
    with TimeSeriesUnstructuredGrid(str(tmp_path / "a"), *MESH) as ts:
        with pytest.raises(ValueError, match="needs the time value"):
            ts.write_step(point_data={"u": u_at(0)})
        ts.write_step(1.0, point_data={"u": u_at(0)})
        with pytest.raises(ValueError, match="must increase"):
            ts.write_step(1.0, point_data={"u": u_at(0)})
    with pytest.raises(ValueError, match="strictly increasing"):
        TimeSeriesUnstructuredGrid(str(tmp_path / "b"), *MESH, times=[0.0, 0.0])
    with TimeSeriesUnstructuredGrid(str(tmp_path / "c"), *MESH, times=[0.0, 0.5]) as ts:
        with pytest.raises(ValueError, match="does not match"):
            ts.write_step(0.1, point_data={"u": u_at(0)})
        ts.write_step(0.0, point_data={"u": u_at(0)})
        ts.write_step(point_data={"u": u_at(1)})
        with pytest.raises(ValueError, match="already been written"):
            ts.write_step(point_data={"u": u_at(2)})


def test_missing_steps_at_close(tmp_path):
    ts = TimeSeriesUnstructuredGrid(str(tmp_path / "v"), *MESH, times=TIMES)
    ts.write_step(point_data={"u": u_at(0)})
    with pytest.raises(ValueError, match="only 1 of the 4"):
        ts.close()
    assert os.listdir(tmp_path) == []
    ts = TimeSeriesUnstructuredGrid(str(tmp_path / "p"), *MESH, times=TIMES, output="pvd")
    ts.write_step(point_data={"u": u_at(0)})
    with pytest.raises(ValueError, match="only 1 of the 4"):
        ts.close()
    assert os.path.exists(tmp_path / "p.pvd")       # the step that was written is kept


def test_other_errors(tmp_path):
    with pytest.raises(ValueError, match="no time steps"):
        with TimeSeriesUnstructuredGrid(str(tmp_path / "a"), *MESH):
            pass
    with pytest.raises(ValueError, match="both as static and time-dependent"):
        with TimeSeriesUnstructuredGrid(str(tmp_path / "b"), *MESH, **STATIC) as ts:
            ts.write_step(0.0, point_data={"h": u_at(0)})
    with pytest.raises(ValueError, match="output must be 'vtu', 'vtkhdf' or 'pvd'"):
        TimeSeriesUnstructuredGrid(str(tmp_path / "c"), *MESH, output="vtk")
    ts = write_all(tmp_path / "d", steps=1)
    with pytest.raises(RuntimeError, match="closed"):
        ts.write_step(5.0, point_data={"u": u_at(0)})
    assert not leftovers(tmp_path)


# ---------- failures part-way, retries, and other edge cases ----------
import VTKwrite.interface as interface_module


def test_failure_inside_a_step_discards_the_output(tmp_path, monkeypatch):
    real = interface_module.writeArrayToFile
    calls = {"n": 0}

    def failing(stream, data):
        calls["n"] += 1
        if calls["n"] == 3:                         # second array of the second step
            raise OSError(28, "No space left on device")
        real(stream, data)
    for kw in (dict(times=TIMES), {}):
        calls["n"] = 0
        monkeypatch.setattr(interface_module, "writeArrayToFile", failing)
        ts = TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH, **kw)
        ts.write_step(TIMES[0], point_data={"u": u_at(0), "v": u_at(0)})
        with pytest.raises(OSError):
            ts.write_step(TIMES[1], point_data={"u": u_at(1), "v": u_at(1)})
        monkeypatch.setattr(interface_module, "writeArrayToFile", real)
        with pytest.raises(RuntimeError, match="closed"):
            ts.write_step(TIMES[2], point_data={"u": u_at(2), "v": u_at(2)})
        with pytest.raises(RuntimeError, match="aborted"):
            ts.close()
        assert os.listdir(tmp_path) == []


def test_failed_pvd_update_can_be_retried(tmp_path, monkeypatch):
    from VTKwrite import vtkbin
    real_save = vtkbin.VtkGroup.save
    state = {"fail": False}

    def save(self):
        if state["fail"]:
            state["fail"] = False
            self.abort()
            raise PermissionError("pvd file is open in another program")
        real_save(self)
    monkeypatch.setattr(vtkbin.VtkGroup, "save", save)
    with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH, output="pvd") as ts:
        ts.write_step(0.0, point_data={"u": u_at(0)})
        state["fail"] = True
        with pytest.raises(PermissionError):
            ts.write_step(1.0, point_data={"u": u_at(1)})
        ts.write_step(1.0, point_data={"u": u_at(1)})
        ts.write_step(2.0, point_data={"u": u_at(2)})
    entries = [(float(d.get("timestep")), d.get("file"))
               for d in ET.parse(tmp_path / "sol.pvd").getroot().iter("DataSet")]
    assert entries == [(0.0, "sol/sol_00000.vtu"), (1.0, "sol/sol_00001.vtu"),
                       (2.0, "sol/sol_00002.vtu")]
    blocks = read_appended(str(tmp_path / "sol" / "sol_00001.vtu"))
    np.testing.assert_array_equal(blocks[("u", None)], u_at(1))


def test_first_step_can_be_retried_after_opening_fails(tmp_path):
    for i, kw in enumerate((dict(times=TIMES[:2]), {})):
        folder = tmp_path / ("later%d" % i)
        path = folder / "sol"
        ts = TimeSeriesUnstructuredGrid(str(path), *MESH, time_values_file=False, **kw)
        with pytest.raises(OSError):
            ts.write_step(TIMES[0], point_data={"u": u_at(0)})   # folder does not exist yet
        os.makedirs(folder)
        ts.write_step(TIMES[0], point_data={"u": u_at(0)})
        ts.write_step(TIMES[1], point_data={"u": u_at(1)})
        ts.close()
        np.testing.assert_array_equal(read_appended(str(path) + ".vtu")[("u", "1")], u_at(1))


def test_time_values_file_folder_checked_early(tmp_path):
    with pytest.raises(ValueError, match="does not exist"):
        TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH,
                                   time_values_file=str(tmp_path / "nodir" / "tv"))


def test_non_finite_times_rejected(tmp_path):
    with pytest.raises(ValueError, match="finite"):
        TimeSeriesUnstructuredGrid(str(tmp_path / "a"), *MESH, times=[0.0, np.nan, 1.0])
    with TimeSeriesUnstructuredGrid(str(tmp_path / "b"), *MESH) as ts:
        with pytest.raises(ValueError, match="finite"):
            ts.write_step(np.inf, point_data={"u": u_at(0)})
        ts.write_step(0.0, point_data={"u": u_at(0)})


def test_tuple_of_component_arrays_rejected(tmp_path):
    x3 = np.array([0., 1., 0.])
    with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), x3, x3[::-1], np.zeros(3),
                                    [0, 1, 2], [3], [5]) as ts:
        with pytest.raises(TypeError, match="column_stack"):
            ts.write_step(0.0, point_data={"v": (np.ones(3), 10 * np.ones(3), 100 * np.ones(3))})
        ts.write_step(0.0, point_data={"v": np.column_stack((np.ones(3), 10 * np.ones(3),
                                                              100 * np.ones(3)))})
    np.testing.assert_array_equal(read_appended(str(tmp_path / "sol.vtu"))[("v", "0")][:3],
                                  [1, 10, 100])


def test_close_twice_returns_the_same_path(tmp_path):
    ts = TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH)
    ts.write_step(0.0, point_data={"u": u_at(0)})
    assert ts.close() == ts.close() == str(tmp_path / "sol.vtu")


# ---------- output = "vtkhdf" ----------
def write_hdf(path, steps=4, **kw):
    pytest.importorskip("h5py")
    return write_all(path, steps=steps, output="vtkhdf", **kw)


def test_vtkhdf_structure(tmp_path):
    h5py = pytest.importorskip("h5py")
    write_hdf(tmp_path / "sol")
    assert not leftovers(tmp_path)
    with h5py.File(tmp_path / "sol.vtkhdf", "r") as f:
        root = f["VTKHDF"]
        assert list(root.attrs["Version"]) == [2, 0]
        assert root.attrs["Type"] == b"UnstructuredGrid"
        np.testing.assert_array_equal(root["Points"][:], np.column_stack((X, Y, Z)))
        np.testing.assert_array_equal(root["Offsets"][:], [0, 3, 7])
        steps = root["Steps"]
        assert steps.attrs["NSteps"] == 4
        np.testing.assert_array_equal(steps["Values"][:], TIMES)
        np.testing.assert_array_equal(steps["PointDataOffsets/u"][:], [0, 6, 12, 18])
        np.testing.assert_array_equal(steps["PointDataOffsets/h"][:], [0, 0, 0, 0])
        np.testing.assert_array_equal(steps["CellDataOffsets/p"][:], [0, 2, 4, 6])
        for k in range(4):
            np.testing.assert_array_equal(root["PointData/u"][6 * k:6 * k + 6], u_at(k))
            np.testing.assert_array_equal(root["PointData/q"][6 * k:6 * k + 6], q_at(k))
            np.testing.assert_array_equal(root["CellData/p"][2 * k:2 * k + 2], p_at(k))
        np.testing.assert_array_equal(root["CellData/material"][:], [7, 8])
        assert root["PointData/u"].attrs["Attribute"] == b"Scalars"
        assert root["PointData/q"].attrs["Attribute"] == b"Vectors"


def test_vtkhdf_read_by_vtk(tmp_path, vtk):
    write_hdf(tmp_path / "sol", times=TIMES)
    r = vtk.vtkHDFReader()
    r.SetFileName(str(tmp_path / "sol.vtkhdf"))
    r.UpdateInformation()
    info = r.GetOutputInformation(0)
    key = vtk.vtkStreamingDemandDrivenPipeline.TIME_STEPS()
    assert [info.Get(key, i) for i in range(info.Length(key))] == TIMES
    for k, t in enumerate(TIMES):
        r.UpdateTimeStep(t)
        out = r.GetOutput()
        assert out.GetNumberOfPoints() == 6 and out.GetNumberOfCells() == 2
        pd, cd = out.GetPointData(), out.GetCellData()
        np.testing.assert_array_equal(vtk_array(vtk, pd, "u"), u_at(k))
        np.testing.assert_array_equal(vtk_array(vtk, pd, "q"), q_at(k))
        np.testing.assert_array_equal(vtk_array(vtk, cd, "p"), p_at(k))
        np.testing.assert_array_equal(vtk_array(vtk, cd, "material"), [7, 8])
        np.testing.assert_array_equal(vtk_array(vtk, pd, "h"), np.linspace(0, 1, 6))


def test_vtkhdf_other_types(tmp_path, vtk):
    pytest.importorskip("h5py")
    with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH, output="vtkhdf") as ts:
        for k in range(2):
            ts.write_step(float(k), point_data={"T": np.arange(54.0).reshape(6, 9) + k,
                                                "uv": ("tcoords", np.zeros((6, 2)) + k),
                                                "id": np.arange(6, dtype=np.int32)})
    r = vtk.vtkHDFReader()
    r.SetFileName(str(tmp_path / "sol.vtkhdf"))
    r.UpdateTimeStep(1.0)
    pd = r.GetOutput().GetPointData()
    np.testing.assert_array_equal(vtk_array(vtk, pd, "T"), np.arange(54.0).reshape(6, 9) + 1)
    np.testing.assert_array_equal(vtk_array(vtk, pd, "uv"), np.ones((6, 2)))
    np.testing.assert_array_equal(vtk_array(vtk, pd, "id"), np.arange(6))


def test_vtkhdf_names(tmp_path):
    pytest.importorskip("h5py")
    with pytest.raises(ValueError, match="may not contain"):
        TimeSeriesUnstructuredGrid(str(tmp_path / "a"), *MESH, output="vtkhdf",
                                   static_point_data={"a.b": np.zeros(6)})
    with TimeSeriesUnstructuredGrid(str(tmp_path / "b"), *MESH, output="vtkhdf") as ts:
        with pytest.raises(ValueError, match="may not contain"):
            ts.write_step(0.0, point_data={"u/v": u_at(0)})
        ts.write_step(0.0, point_data={"u": u_at(0)})
    assert sorted(os.listdir(tmp_path)) == ["b.vtkhdf"]


def test_vtkhdf_error_leaves_nothing(tmp_path, monkeypatch):
    pytest.importorskip("h5py")
    with pytest.raises(RuntimeError):
        with TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH, output="vtkhdf") as ts:
            ts.write_step(0.0, point_data={"u": u_at(0)})
            raise RuntimeError("stop")
    assert os.listdir(tmp_path) == []
    ts = TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH, output="vtkhdf", times=TIMES)
    ts.write_step(point_data={"u": u_at(0)})
    with pytest.raises(ValueError, match="only 1 of the 4"):
        ts.close()
    assert os.listdir(tmp_path) == []


def test_vtkhdf_needs_h5py(tmp_path, monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules, "h5py", None)
    with pytest.raises(ImportError, match="pip install h5py"):
        TimeSeriesUnstructuredGrid(str(tmp_path / "sol"), *MESH, output="vtkhdf")
