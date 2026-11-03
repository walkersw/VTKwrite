"""Shared helpers for the VTKwrite tests.

Two independent ways to check a written file:
  * read_appended(): a small parser for the raw appended format written by VTKwrite
    (does not need VTK; used for exact byte-level checks);
  * read_vtk(): VTK's own XML readers (tests using it are skipped if vtk is missing).
"""
import os
import re
import struct

import numpy as np
import pytest

_TYPES = {"Int8": "i1", "UInt8": "u1", "Int16": "i2", "UInt16": "u2", "Int32": "i4",
          "UInt32": "u4", "Int64": "i8", "UInt64": "u8", "Float32": "f4", "Float64": "f8"}


def read_appended(path):
    """Return {(name, TimeStep or None): numpy array} for every appended DataArray."""
    with open(path, "rb") as f:
        raw = f.read()
    order = "<" if b'byte_order="LittleEndian"' in raw[:400] else ">"
    i = raw.index(b"<AppendedData")
    j = raw.index(b"_", i) + 1
    header = raw[:i].decode("utf-8")
    out = {}
    for m in re.finditer(r"<DataArray ([^>]*)/>", header):
        attrs = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
        off = int(attrs["offset"])
        nbytes = struct.unpack(order + "Q", raw[j + off:j + off + 8])[0]
        arr = np.frombuffer(raw[j + off + 8:j + off + 8 + nbytes],
                            dtype=order + _TYPES[attrs["type"]])
        out[(attrs.get("Name"), attrs.get("TimeStep"))] = arr
    return out


def no_leftovers(directory):
    """True if no temporary '.part' files are left in directory."""
    return not [f for f in os.listdir(directory) if f.endswith(".part")]


@pytest.fixture
def vtk():
    vtk = pytest.importorskip("vtk")
    vtk.vtkObject.GlobalWarningDisplayOff()
    return vtk


def read_vtk(vtk, path, time_step=None):
    readers = {".vtu": vtk.vtkXMLUnstructuredGridReader, ".vts": vtk.vtkXMLStructuredGridReader,
               ".vtr": vtk.vtkXMLRectilinearGridReader, ".vti": vtk.vtkXMLImageDataReader}
    r = readers[os.path.splitext(path)[1]]()
    r.SetFileName(str(path))
    if time_step is None:
        r.Update()
    else:
        r.UpdateTimeStep(time_step)
    assert r.GetErrorCode() == 0
    return r.GetOutput()


def vtk_array(vtk, data, name):
    from vtk.util.numpy_support import vtk_to_numpy
    a = data.GetArray(name)
    assert a is not None, "array %r not found" % name
    return vtk_to_numpy(a)
