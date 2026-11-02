# VTKwrite

Python package to write VTK files in XML format

- This class can write a static unstructured grid, with time-dependent data, all in the same .vtu file!

- This package was forked from another package "PyEVTK" by Paulo A. Herrera.

https://github.com/paulo-herrera/PyEVTK

- Installation is done in the standard way:  ```python setup.py install```

- Or if you want a development version: ```python -m pip install -e .```

- You need Python 3.6 or better.

- Compatibility: the generated files can be read by Paraview versions: 5.4, 5.5, 5.6, 5.7, AND 5.11.  I added a time-dependent .vtu file example to the unit tests of Paraview, so it should keep working in the future.

## Local Install Help

If you need to make a "local" install (because you are not an administrator), then the install command will look like:

```
python setup.py install --prefix ~/.local
```

where you need to modify the ```~/.local``` text to reflect your system setup.

Note: the command above will install it to this directory: ```~/.local/lib/python3.x/site-packages/```

## Usage

After installing, go to the ./examples/ sub-directory and do this:
```
python runall.py run
```
to run the examples; you can also do
```
python runall.py clean
```
to delete the output.

The individual Python files should give you enough information to use this package.

## Time-dependent data

`TimeSeriesUnstructuredGrid` writes a simulation on a fixed unstructured mesh one time step at a time, so the whole history never has to be kept in memory:
```python
from VTKwrite.interface import TimeSeriesUnstructuredGrid

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

See `examples/unstructured_timeseries.py`.

## Testing

This test writes files and reads them back (with VTK's own readers where available):
```
python -m pip install pytest vtk
python -m pytest tests
```

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
