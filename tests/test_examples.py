"""Every example runs and produces files that VTK can read."""
import importlib.util
import os
import xml.etree.ElementTree as ET

import pytest

EXAMPLES = os.path.join(os.path.dirname(__file__), os.pardir, "examples")
NAMES = sorted(f[:-3] for f in os.listdir(EXAMPLES) if f.endswith(".py") and f != "runall.py")


@pytest.mark.parametrize("name", NAMES)
def test_example_runs(tmp_path, monkeypatch, vtk, name):
    if name == "points":
        pytest.importorskip("scipy")
    monkeypatch.chdir(tmp_path)
    spec = importlib.util.spec_from_file_location("example_" + name,
                                                  os.path.join(EXAMPLES, name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.run()
    written = []
    for folder, _, files in os.walk(tmp_path):                 # includes .pvd step folders
        written += [os.path.relpath(os.path.join(folder, f), tmp_path) for f in files]
    assert written and not [f for f in written if f.endswith((".part", ".steps"))]
    readers = {".vtu": vtk.vtkXMLUnstructuredGridReader, ".vts": vtk.vtkXMLStructuredGridReader,
               ".vtr": vtk.vtkXMLRectilinearGridReader, ".vti": vtk.vtkXMLImageDataReader,
               ".vtkhdf": vtk.vtkHDFReader, ".pvd": None}
    for f in written:
        reader = readers[os.path.splitext(f)[1]]
        if reader is None:
            # a .pvd collection: every file it lists must exist (paths are relative to it)
            folder = os.path.dirname(str(tmp_path / f))
            listed = [d.get("file") for d in ET.parse(str(tmp_path / f)).getroot().iter("DataSet")]
            assert listed, f
            missing = [g for g in listed if not os.path.isfile(os.path.join(folder, g))]
            assert not missing, "%s lists missing files: %s" % (f, missing)
            continue
        r = reader()
        r.SetFileName(str(tmp_path / f))
        r.Update()
        assert r.GetErrorCode() == 0 and r.GetOutput().GetNumberOfPoints() > 0, f
