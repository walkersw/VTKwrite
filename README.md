# VTKwrite

Python package to write VTK files in XML format

- This class can write a static unstructured grid, with time-dependent data, all in the same .vtu file!

- This package was forked from another package "PyEVTK" by Paulo A. Herrera.

https://github.com/paulo-herrera/PyEVTK

- Installation: `pip install .` (see [Installation](#installation) below).

- You need Python 3.10 or better, and NumPy.

- Compatibility: the generated files can be read by Paraview versions: 5.4, 5.5, 5.6, 5.7, AND 5.11.  I added a time-dependent .vtu file example to the unit tests of Paraview, so it should keep working in the future.

## Installation

Run one of these in the package folder (the one containing `pyproject.toml`):

```
pip install .              # regular install
pip install -e .           # development ("editable") install: changes to src/ take effect without reinstalling
pip install --user .       # install for your user account only (if you are not an administrator)
```

Optional extras, which can be combined with any of the above:

```
pip install ".[hdf5]"      # also installs h5py, needed for TimeSeriesUnstructuredGrid(..., output="vtkhdf")
pip install -e ".[test]"   # also installs what the tests need: pytest, vtk, h5py, scipy
```

`python setup.py install` is no longer supported. Inside a virtual environment, use `pip install .` without `--user`.

## Usage

Everything can be imported from the package itself, for example
```python
from VTKwrite import unstructuredGridToVTK, TimeSeriesUnstructuredGrid, VtkTriangle, VtkQuad
```
(`import VTKwrite; VTKwrite.__all__` lists the names). The older imports from `VTKwrite.interface` and `VTKwrite.vtkbin` still work.

After installing, go to the ./examples/ sub-directory and do this:
```
python runall.py run
```
to run the examples; you can also do
```
python runall.py clean
```
to delete the output.

The individual Python files should give you enough information to use this package; `examples/array_input.py` shows the array forms described next.

## Points and data arrays

Point and cell data can be given as dictionaries, and the points as one array:
```python
from VTKwrite import unstructuredGridToVTK

# P: (npoints, 3) array, or (npoints, 2) for z = 0;  u: (npoints,);  grad_u: (npoints, 3)
# sigma: (npoints, 9);  material: (ncells,)
unstructuredGridToVTK("out/mesh", points=P, connectivity=conn, offsets=offsets, cell_types=types,
                      all_point_data={"u": u, "grad_u": grad_u, "stress": sigma},
                      all_cell_data={"material": material})
```
- The type of each array follows from its number of components: 1 is scalars, 3 is vectors, 9 is tensors. To choose another type, give a pair: `{"n": ("normals", n), "uv": ("tcoords", uv)}`.
- Multi-component data has one row per point or cell, e.g. shape `(npoints, 3)`. Interleaved 1-D arrays (`x0 y0 z0 x1 y1 z1 ...`) also work, so a 1-D array with 3 values per point is read as vectors. A tuple of component arrays such as `(vx, vy, vz)` is rejected, because with 3 points it could not be told apart from 3 rows; use `np.column_stack((vx, vy, vz))`.
- Arrays of dtype `bool` (e.g. boundary markers) are written as `UInt8`.
- `structuredToVTK`, `rectilinearToVTK` and `imageToVTK` take arrays of the grid's shape: `(nx, ny, nz)` for scalars and `(nx, ny, nz, 3)` for vectors, indexed `[i, j, k, component]`. `structuredToVTK(path, points=P)` takes `P` of shape `(nx, ny, nz, 3)`. The reordering VTK needs is done internally.
- `points=` works for `unstructuredGridToVTK`, `pointsToVTK`, `pointsToVTKAsTIN`, `linesToVTK`, `polyLinesToVTK`, `structuredToVTK` and `TimeSeriesUnstructuredGrid`; the arguments after it are then given by keyword.
- The older forms still work and give the same files: `x, y, z` coordinates, and lists `[["u", "scalars", u], ["grad_u", "vectors", grad_u_interleaved]]`.

## Time-dependent data

`TimeSeriesUnstructuredGrid` writes a simulation on a fixed unstructured mesh one time step at a time, so the whole history never has to be kept in memory:
```python
from VTKwrite import TimeSeriesUnstructuredGrid

with TimeSeriesUnstructuredGrid("out/sol", x, y, z, connectivity, offsets, cell_types,
                                static_cell_data={"material": material}) as ts:
    for t in times:
        ...  # compute u (one value per point) and q (shape (npoints, 3))
        ts.write_step(t, point_data={"u": u, "q": q})
```
Three kinds of output (choose with `output=`):

| `output=` | Files | Notes |
|---|---|---|
| `"vtu"` (default) | `out/sol.vtu` + `out/sol_time_values.vtu` | One XML file, mesh stored once. **Not read correctly by VTK 9.6+ / ParaView 6.1+** (see below). |
| `"vtkhdf"` | `out/sol.vtkhdf` | One HDF5 file, mesh stored once, actual time values shown by ParaView. Needs `h5py` (`pip install h5py`). Works with ParaView 5.12 and later, including 6.1. |
| `"pvd"` | `out/sol.pvd` + `out/sol/sol_00000.vtu`, ... | One file per step (mesh repeated), actual time values. Works with all ParaView versions; steps already written can be opened if the run stops early. |

The points may also be given as one array: `TimeSeriesUnstructuredGrid("out/sol", points=P, connectivity=..., offsets=..., cell_types=...)`.

See `examples/unstructured_timeseries.py`.

## Testing

The tests write files and read them back (with VTK's own readers where available). From the package folder:
```
pip install -e ".[test]"
python -m pytest tests
```
Tests that need `vtk`, `h5py` or `scipy` are skipped if that package is missing.

## Known issue: single-file .vtu time series in VTK 9.6 / ParaView 6.1 and later

VTK 9.6 and later, and therefore ParaView 6.1 and later, read only the first time step of a time series stored in a single `.vtu` file. That covers files written by `timeseries_unstructuredGrid` and by `TimeSeriesUnstructuredGrid` with `output="vtu"`. When such a file is opened, ParaView prints warnings like

```
Warning: In vtkXMLDataReader.cxx, line 286
vtkXMLUnstructuredGridReader (...): Another point data array named u already exists. Ignoring.
```

and every later time step shows the values of the first step (or garbage). The files themselves are correct: VTK 9.5 / ParaView 6.0 and earlier read them properly. This is a bug in VTK's XML reader, introduced in VTK 9.6; it has been reported upstream together with a fix.

Until a fixed VTK/ParaView is released:
- use `TimeSeriesUnstructuredGrid` with `output="vtkhdf"` (single file) or `output="pvd"`, or
- open single-file `.vtu` time series with ParaView 6.0 or earlier.

Single-step `.vtu` files, `.pvd` collections and `.vtkhdf` files are not affected.
