"""Bad input must raise a clear error and must not leave a file behind."""
import os
import subprocess
import sys

import numpy as np
import pytest

from VTKwrite.interface import (imageToVTK, structuredToVTK, unstructuredGridToVTK,
                                polyLinesToVTK, linesToVTK)
from VTKwrite.vtkbin import VtkTriangle, VtkQuad, VtkFile, VtkUnstructuredGrid
from conftest import no_leftovers


def two_cells():
    x = np.array([0., 1., 2., 0., 1., 2.])
    y = np.array([0., 0., 0., 1., 1., 1.])
    z = np.zeros(6)
    conn = np.array([0, 1, 3, 1, 2, 5, 4])
    offsets = np.array([3, 7])
    types = np.array([VtkTriangle.tid, VtkQuad.tid], dtype=np.uint8)
    return x, y, z, conn, offsets, types


def write(tmp_path, **kw):
    x, y, z, conn, offsets, types = two_cells()
    args = dict(x=x, y=y, z=z, connectivity=conn, offsets=offsets, cell_types=types)
    args.update(kw)
    return unstructuredGridToVTK(str(tmp_path / "ug"), **args)


@pytest.mark.parametrize("kw, exc, msg", [
    (dict(all_cell_data=[["p", "scalar", np.zeros(2)]]), ValueError, "type 'scalar' is not allowed"),
    (dict(all_cell_data=[["p", "scalars", np.zeros(3)]]), ValueError, "has 3 values, expected 2"),
    (dict(all_point_data=[["v", "vectors", np.zeros(6)]]), ValueError, "expected 18"),
    (dict(all_cell_data=[["p", "scalars", np.zeros(2)], ["p", "scalars", np.ones(2)]]),
     ValueError, "used twice"),
    (dict(all_cell_data=[["m", "scalars", np.zeros(2, dtype=np.float16)]]), TypeError, "'float16'"),
    (dict(all_cell_data=[["p", "scalars"]]), TypeError, "[name, type, data]"),
    (dict(all_point_data=[["uv", "tcoords", np.zeros(7)]]), ValueError, "tcoords"),
    (dict(y=np.zeros(5)), ValueError, "same number of entries"),
    (dict(offsets=np.array([3])), ValueError, "offsets has 1 entries"),
    (dict(offsets=np.array([3, 6])), ValueError, "length of connectivity"),
    (dict(connectivity=np.array([0, 1, 3, 1, 2, 6, 4])), ValueError, "only 6 points"),
    (dict(connectivity=np.array([0, 1, 3, 1, 2, 5, 4.5])), ValueError, "whole numbers"),
])
def test_invalid_input(tmp_path, kw, exc, msg):
    with pytest.raises(exc, match=None) as info:
        write(tmp_path, **kw)
    assert msg in str(info.value)
    assert os.listdir(tmp_path) == []          # no file, no temporary


def test_existing_file_kept_on_error(tmp_path):
    good = write(tmp_path)
    before = open(good, "rb").read()
    with pytest.raises(ValueError):
        write(tmp_path, all_cell_data=[["p", "scalars", np.zeros(5)]])
    assert open(good, "rb").read() == before
    assert no_leftovers(tmp_path)


def test_error_inside_low_level_context_manager(tmp_path):
    with pytest.raises(RuntimeError, match="boom"):
        with VtkFile(str(tmp_path / "ll"), VtkUnstructuredGrid) as w:
            w.openGrid()
            raise RuntimeError("boom")
    assert os.listdir(tmp_path) == []


def test_unsaved_low_level_file_leaves_no_temporary(tmp_path):
    w = VtkFile(str(tmp_path / "ll"), VtkUnstructuredGrid)
    w.openGrid()
    del w
    import gc; gc.collect()
    assert os.listdir(tmp_path) == []


def test_image_errors(tmp_path):
    with pytest.raises(ValueError, match="does not match point data shape"):
        imageToVTK(str(tmp_path / "i"), all_cell_data=[["p", "scalars", np.zeros((2, 2, 2))]],
                   all_point_data=[["T", "scalars", np.zeros((4, 3, 3))]])
    with pytest.raises(ValueError, match="3-D arrays"):
        imageToVTK(str(tmp_path / "i"), all_point_data=[["T", "scalars", np.zeros(8)]])
    with pytest.raises(ValueError, match="all_cell_data or all_point_data"):
        imageToVTK(str(tmp_path / "i"))
    assert os.listdir(tmp_path) == []


def test_structured_shape_errors(tmp_path):
    X = np.zeros((3, 3, 2))
    with pytest.raises(ValueError, match="same shape"):
        structuredToVTK(str(tmp_path / "s"), X, X, np.zeros((3, 2, 2)))
    with pytest.raises(ValueError, match="expected \\(3, 3, 2\\)"):
        structuredToVTK(str(tmp_path / "s"), X, X, X,
                        all_point_data=[["T", "scalars", np.zeros((3, 2, 3))]])
    assert os.listdir(tmp_path) == []


def test_lines_and_polylines_errors(tmp_path):
    x = np.zeros(5)
    with pytest.raises(ValueError, match="even number"):
        linesToVTK(str(tmp_path / "l"), x, x, x)
    with pytest.raises(ValueError, match="sums to 4"):
        polyLinesToVTK(str(tmp_path / "pl"), x, x, x, [2, 2])
    assert os.listdir(tmp_path) == []


def test_checks_run_under_python_O(tmp_path):
    """Input checks must not be asserts (which python -O removes)."""
    code = ("import numpy as np, sys\n"
            "from VTKwrite.interface import unstructuredGridToVTK\n"
            "try:\n"
            "    unstructuredGridToVTK(sys.argv[1], np.zeros(3), np.zeros(3), np.zeros(3),\n"
            "        np.array([0, 1, 2]), np.array([3]), np.array([5]),\n"
            "        all_cell_data=[['p', 'scalars', np.zeros(4)]])\n"
            "except ValueError:\n"
            "    sys.exit(0)\n"
            "sys.exit(1)\n")
    r = subprocess.run([sys.executable, "-O", "-c", code, str(tmp_path / "o")])
    assert r.returncode == 0
    assert os.listdir(tmp_path) == []
