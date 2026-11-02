"""The NumPy block writer must produce exactly the bytes of the original struct-based writer."""
import io
import struct
import sys

import numpy as np
import pytest

from VTKwrite import pyvtk

ORDER = "<" if sys.byteorder == "little" else ">"


def reference_array(data):
    """The original writeArrayToFile: values in Fortran order, packed with struct."""
    dd = np.ravel(data, order="F")
    return struct.pack(ORDER + str(dd.size) + pyvtk.np_to_struct[data.dtype.name], *dd)


def reference_arrays(x, y, z):
    """The original writeArraysToFile: x, y, z interleaved, packed with the format of x."""
    fmt = ORDER + pyvtk.np_to_struct[x.dtype.name]
    xx, yy, zz = (np.ravel(a, order="F") for a in (x, y, z))
    return b"".join(struct.pack(fmt, a) + struct.pack(fmt, b) + struct.pack(fmt, c)
                    for a, b, c in zip(xx, yy, zz))


def write(fn, *args):
    stream = io.BytesIO()
    fn(stream, *args)
    return stream.getvalue()


rng = np.random.default_rng(3)
ARRAYS = {
    "float64 1-D": rng.random(17),
    "float32 1-D": rng.random(17).astype(np.float32),
    "int16 1-D": rng.integers(-100, 100, 17).astype(np.int16),
    "uint8 1-D": rng.integers(0, 255, 17).astype(np.uint8),
    "non-native byte order": rng.random(17).astype(">f8" if ORDER == "<" else "<f8"),
    "strided view": rng.random(40)[::3],
    "3-D C order": rng.random((3, 4, 5)),
    "3-D Fortran order": np.asfortranarray(rng.random((3, 4, 5))),
    "3-D strided": rng.random((6, 4, 5))[::2, :, ::-1],
    "empty": np.zeros(0),
}


@pytest.mark.parametrize("name", ARRAYS)
def test_write_array_matches_reference(name):
    a = ARRAYS[name]
    assert write(pyvtk.writeArrayToFile, a) == reference_array(a)


@pytest.mark.parametrize("dtypes", [("float64",) * 3, ("float32",) * 3, ("int64",) * 3,
                                    ("float64", "int64", "float64")])
@pytest.mark.parametrize("shape", [(11,), (3, 4, 2)])
def test_write_arrays_matches_reference(dtypes, shape):
    x, y, z = (rng.integers(-50, 50, shape).astype(d) for d in dtypes)
    if shape != (11,):
        y = np.asfortranarray(y)
    assert write(pyvtk.writeArraysToFile, x, y, z) == reference_arrays(x, y, z)


def test_write_arrays_in_chunks(monkeypatch):
    monkeypatch.setattr(pyvtk, "_CHUNK_POINTS", 4)       # forces several chunks + a partial one
    x, y, z = rng.random(10), rng.random(10), rng.random(10)
    assert write(pyvtk.writeArraysToFile, x, y, z) == reference_arrays(x, y, z)


def test_unsupported_dtype():
    with pytest.raises(TypeError, match="not supported"):
        write(pyvtk.writeArrayToFile, np.zeros(3, dtype=bool))
