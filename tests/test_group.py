"""VtkGroup (.pvd collections): stored paths, several files per time step, block names."""
import ntpath
import os
import types
import xml.etree.ElementTree as ET

import numpy as np

from VTKwrite import VtkGroup, pointsToVTK
from VTKwrite import vtkbin


def datasets(path):
    return [d.attrib for d in ET.parse(str(path)).getroot().iter("DataSet")]


def test_paths_relative_to_the_pvd(tmp_path):
    os.makedirs(tmp_path / "out" / "steps")
    f = pointsToVTK(str(tmp_path / "out" / "steps" / "a"), [0.0], [0.0], [0.0])
    with VtkGroup(str(tmp_path / "out" / "series")) as g:
        g.addFile(f, sim_time=0.0)
    with VtkGroup(str(tmp_path / "other")) as g:
        g.addFile(f, sim_time=0.0)
    assert datasets(tmp_path / "out" / "series.pvd")[0]["file"] == "steps/a.vtu"
    assert datasets(tmp_path / "other.pvd")[0]["file"] == "out/steps/a.vtu"


def test_windows_paths(tmp_path, monkeypatch):
    """ Simulated with ntpath: backslashes become "/", and a file on another drive is
        stored with its absolute path (os.path.relpath raises ValueError there). """
    with VtkGroup(str(tmp_path / "series")) as g:
        g.root = "C:\\work"
        monkeypatch.setattr(vtkbin, "os", types.SimpleNamespace(path=ntpath, sep="\\"))
        g.addFile("C:\\work\\steps\\a.vtu", sim_time=0.0)
        g.addFile("D:\\data\\b.vtu", sim_time=1.0)
        monkeypatch.undo()
    files = [d["file"] for d in datasets(tmp_path / "series.pvd")]
    assert files == ["steps/a.vtu", "D:/data/b.vtu"]


def test_several_files_per_time_step(tmp_path):
    with VtkGroup(str(tmp_path / "s")) as g:
        for t in (0.0, 0.5):
            for k, body in enumerate(("fluid", "solid & wall")):
                g.addFile(str(tmp_path / ("%s_%s.vtu" % (k, t))), sim_time=t, part=k, name=body)
        g.addFile(str(tmp_path / "x.vtu"), sim_time=1.0)
    ds = datasets(tmp_path / "s.pvd")
    assert [(d["timestep"], d["part"], d.get("name")) for d in ds] == [
        ("0.0", "0", "fluid"), ("0.0", "1", "solid & wall"),
        ("0.5", "0", "fluid"), ("0.5", "1", "solid & wall"), ("1.0", "0", None)]


def test_attributes_unchanged_without_name(tmp_path):
    """ The attributes and their order are as before when no name is given. """
    with VtkGroup(str(tmp_path / "s")) as g:
        g.addFile(str(tmp_path / "a.vtu"), sim_time=np.float64(0.25))
    with open(tmp_path / "s.pvd", "rb") as f:
        text = f.read()
    assert b'<DataSet timestep="0.25" group="" part="0" file="a.vtu"/>' in text
