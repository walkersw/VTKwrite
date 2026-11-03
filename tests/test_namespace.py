"""The package namespace: top-level exports, __version__ and the older import paths."""
import os
import re
import subprocess
import sys

import pytest

import VTKwrite
from VTKwrite import interface, vtkbin, lagrange

ROOT = os.path.join(os.path.dirname(__file__), os.pardir)

# Public names (other than modules) that VTKwrite.interface had in 0.0.4, partly through
# "from .vtkbin import *".  They must stay importable from there.
INTERFACE_0_0_4 = """
    TimeSeriesUnstructuredGrid VtkCellType VtkDataType VtkFile VtkFileType VtkFloat32
    VtkFloat64 VtkGroup VtkHexahedron VtkImageData VtkInt16 VtkInt32 VtkInt64 VtkInt8 VtkLine
    VtkPixel VtkPolyData VtkPolyLine VtkPolyVertex VtkPolygon VtkPyramid VtkQuad
    VtkQuadraticEdge VtkQuadraticHexahedron VtkQuadraticQuad VtkQuadraticTetra
    VtkQuadraticTriangle VtkRectilinearGrid VtkStructuredGrid VtkTetra VtkTriangle
    VtkTriangleStrip VtkUInt16 VtkUInt32 VtkUInt64 VtkUInt8 VtkUnstructuredGrid VtkVertex
    VtkVoxel VtkWedge XmlWriter cylinderToVTK imageToVTK linesToVTK np_to_vtk pointsToVTK
    pointsToVTKAsTIN polyLinesToVTK rectilinearToVTK sanitize_comment structuredToVTK
    timeseries_unstructuredGrid unstructuredGridToVTK writeArrayToFile writeArraysToFile
    writeBlockSize""".split()

CELL_TYPES = {"VtkVertex": "VTK_VERTEX", "VtkPolyVertex": "VTK_POLY_VERTEX",
              "VtkLine": "VTK_LINE", "VtkPolyLine": "VTK_POLY_LINE",
              "VtkTriangle": "VTK_TRIANGLE", "VtkTriangleStrip": "VTK_TRIANGLE_STRIP",
              "VtkPolygon": "VTK_POLYGON", "VtkPixel": "VTK_PIXEL", "VtkQuad": "VTK_QUAD",
              "VtkTetra": "VTK_TETRA", "VtkVoxel": "VTK_VOXEL",
              "VtkHexahedron": "VTK_HEXAHEDRON", "VtkWedge": "VTK_WEDGE",
              "VtkPyramid": "VTK_PYRAMID", "VtkQuadraticEdge": "VTK_QUADRATIC_EDGE",
              "VtkQuadraticTriangle": "VTK_QUADRATIC_TRIANGLE",
              "VtkQuadraticQuad": "VTK_QUADRATIC_QUAD",
              "VtkQuadraticTetra": "VTK_QUADRATIC_TETRA",
              "VtkQuadraticHexahedron": "VTK_QUADRATIC_HEXAHEDRON",
              "VtkQuadraticWedge": "VTK_QUADRATIC_WEDGE",
              "VtkQuadraticPyramid": "VTK_QUADRATIC_PYRAMID",
              "VtkBiQuadraticQuad": "VTK_BIQUADRATIC_QUAD",
              "VtkTriQuadraticHexahedron": "VTK_TRIQUADRATIC_HEXAHEDRON",
              "VtkBiQuadraticTriangle": "VTK_BIQUADRATIC_TRIANGLE",
              "VtkCubicLine": "VTK_CUBIC_LINE",
              "VtkTriQuadraticPyramid": "VTK_TRIQUADRATIC_PYRAMID",
              "VtkLagrangeCurve": "VTK_LAGRANGE_CURVE",
              "VtkLagrangeTriangle": "VTK_LAGRANGE_TRIANGLE",
              "VtkLagrangeQuadrilateral": "VTK_LAGRANGE_QUADRILATERAL",
              "VtkLagrangeTetrahedron": "VTK_LAGRANGE_TETRAHEDRON",
              "VtkLagrangeHexahedron": "VTK_LAGRANGE_HEXAHEDRON",
              "VtkLagrangeWedge": "VTK_LAGRANGE_WEDGE"}


def test_version():
    assert re.fullmatch(r"\d+\.\d+\.\d+", VTKwrite.__version__)


def test_pyproject_takes_the_version_from_the_package():
    tomllib = pytest.importorskip("tomllib")          # Python >= 3.11
    path = os.path.join(ROOT, "pyproject.toml")
    if not os.path.isfile(path):
        pytest.skip("pyproject.toml not found")
    with open(path, "rb") as f:
        meta = tomllib.load(f)
    assert "version" not in meta["project"] and "version" in meta["project"]["dynamic"]
    assert meta["tool"]["setuptools"]["dynamic"]["version"] == {"attr": "VTKwrite.__version__"}


def test_all_names_exist_and_are_the_originals():
    assert len(VTKwrite.__all__) == len(set(VTKwrite.__all__))
    for name in VTKwrite.__all__:
        obj = getattr(VTKwrite, name)
        source = next(m for m in (interface, vtkbin, lagrange) if hasattr(m, name))
        assert obj is getattr(source, name), name


def test_star_import():
    namespace = {}
    exec("__version__ = 'mine'\nfrom VTKwrite import *", namespace)
    assert set(VTKwrite.__all__) <= set(namespace)
    assert namespace["__version__"] == "mine"


def test_exports_cover_the_writers_and_cell_types():
    writers = [n for n in dir(interface) if n.endswith("ToVTK") and not n.startswith("_")]
    cells = [n for n, v in vars(vtkbin).items() if isinstance(v, vtkbin.VtkCellType)]
    assert writers and cells
    assert set(writers + cells) <= set(VTKwrite.__all__)


def test_old_import_paths_still_work():
    missing = [n for n in INTERFACE_0_0_4 if not hasattr(interface, n)]
    assert not missing
    for name in INTERFACE_0_0_4:
        if hasattr(vtkbin, name):
            assert getattr(interface, name) is getattr(vtkbin, name), name


def test_cell_type_ids_match_vtk(vtk):
    for name, constant in CELL_TYPES.items():
        assert getattr(VTKwrite, name).tid == getattr(vtk, constant), name
    assert set(CELL_TYPES) == {n for n, v in vars(vtkbin).items()
                               if isinstance(v, vtkbin.VtkCellType)}


def test_missing_numpy_raises_import_error():
    """NumPy is required: importing without it must fail, not print and continue."""
    code = "import sys; sys.modules['numpy'] = None; import VTKwrite"
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode != 0
    assert "ImportError" in r.stderr or "ModuleNotFoundError" in r.stderr
    assert "Numpy is not installed" not in r.stdout
