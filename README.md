# VTKwrite

Python package to write VTK files in XML format

- This class can write a static unstructured grid, with time-dependent data, all in the same .vtu file!

- This package was forked from another package "PyEVTK" by Paulo A. Herrera.

https://github.com/paulo-herrera/PyEVTK

- Installation: `pip install .` (see [Installation](#installation) below).

- You need Python 3.10 or better, and NumPy.

- Compatibility: the tests write files and read them back with VTK's own readers. They pass with VTK 9.5.2 and 9.7.1, and with VTK 9.2.6 except for time-dependent `.vtkhdf` files, which that version cannot read (VTK 9.4 and later can). Earlier versions of this package were checked with ParaView 5.4, 5.5, 5.6, 5.7 and 5.11. Single-file `.vtu` time series are not read correctly by VTK 9.6 and later (ParaView 6.1 and later); see [the known issue](#known-issue-single-file-vtu-time-series-in-vtk-96--paraview-61-and-later) below.

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

## Higher-order (Lagrange) cells

Besides the linear cells, the cell types include VTK's quadratic and cubic cells (`VtkQuadraticTriangle`, `VtkQuadraticTetra`, ..., `VtkQuadraticWedge`, `VtkQuadraticPyramid`, `VtkBiQuadraticQuad`, `VtkTriQuadraticHexahedron`, `VtkBiQuadraticTriangle`, `VtkCubicLine`, `VtkTriQuadraticPyramid`) and its arbitrary-order Lagrange cells: `VtkLagrangeCurve`, `VtkLagrangeTriangle`, `VtkLagrangeQuadrilateral`, `VtkLagrangeTetrahedron`, `VtkLagrangeHexahedron` and `VtkLagrangeWedge`. The degree of a Lagrange cell follows from its number of nodes (e.g. 10 for a cubic triangle).

VTK expects the nodes of a Lagrange cell at the equispaced points of its reference cell, in VTK's own order. Three helpers (NumPy only) deal with that:
- `lagrange_nodes(cell_type, degree)` gives the reference coordinates of the nodes in VTK's order. With `vertices=` they are mapped to another reference cell, e.g. `[-1, 1]^2`.
- `node_permutation(cell_type, degree, my_nodes)` gives the order in which to take the nodes of your element. With a cell-to-node map `cell_nodes` (one row per cell, in your local order), the connectivity is `cell_nodes[:, perm].ravel()`. This works when the degrees of freedom are values at these equispaced points; otherwise (e.g. Gauss-Lobatto nodes, or a hierarchical basis) evaluate the fields at `lagrange_nodes(...)` in each cell.
- `cellwise_mesh(cell_points, cell_type)` builds a mesh in which every cell has its own nodes, for discontinuous fields or anything evaluated cell by cell:
  ```python
  mesh = cellwise_mesh(X, VtkLagrangeTriangle)       # X: (ncells, nnodes, 3) node coordinates
  unstructuredGridToVTK("dg", **mesh, all_point_data={"u": U.reshape(-1)})   # U: (ncells, nnodes)
  ```

ParaView draws these cells curved; raise "Nonlinear Subdivision Level" (Properties panel, advanced) to see more detail. See `examples/lagrange_cells.py`.

Notes:
- The number of points of every cell is checked against its type. Lagrange quadrilaterals, hexahedra and wedges with different degrees per direction need VTK's cell data array `"HigherOrderDegrees"`; when it is given, their point counts are not checked.
- Files containing Lagrange hexahedra declare `VTKFile` version 2.1, because VTK reorders the nodes of hexahedra in files declaring an older version. This was checked with VTK 9.2.6, 9.5.2 and 9.7.1; much older VTK or ParaView releases may read such hexahedra wrongly. All other files are written as before.
- VTK has no working Lagrange pyramid (type 74 is read as an empty cell), so none is provided.

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

### Collections of files (`VtkGroup`)

`VtkGroup` writes a `.pvd` file that lists data files with their time values, e.g. for a mesh that changes in time (see `examples/group.py`):
```python
from VTKwrite import VtkGroup

with VtkGroup("out/series") as group:       # out/series.pvd is written when the block ends
    for k, t in enumerate(times):
        f = unstructuredGridToVTK("out/step_%04d" % k, ...)
        group.addFile(f, sim_time=t)
```
- The paths are stored relative to the `.pvd` file, with `/` as the separator. On Windows, a file on another drive than the `.pvd` is stored with its absolute path.
- Several files may belong to one time step, e.g. one per body: call `addFile` for each with exactly the same `sim_time`. ParaView then shows that step as a multiblock dataset with one block per file, in the order of the `addFile` calls; `name="..."` names the block. ParaView's PVD reader selects the files by time only and does not use the `group` and `part` attributes. (This follows ParaView's source code, `vtkPVDReader` and `vtkXMLCollectionReader`; it was not tested in ParaView.)

### Migrating from `timeseries_unstructuredGrid`

`timeseries_unstructuredGrid` is deprecated: it issues a `DeprecationWarning` and will be removed in the next release. `TimeSeriesUnstructuredGrid` writes the same kind of single `.vtu` file, one time step at a time:

| `timeseries_unstructuredGrid` | `TimeSeriesUnstructuredGrid` |
|---|---|
| `ts = timeseries_unstructuredGrid(path, times)`, then `ts.init_unstructuredGridToVTK(x, y, z, conn, offsets, types, all_cell_data=..., all_point_data=..., comments=...)` | `ts = TimeSeriesUnstructuredGrid(path, x, y, z, conn, offsets, types, times=times, comments=...)`; the fields are given with each step |
| one `ts.append_data([value_at_t0, value_at_t1, ...])` per variable, cell data first, in the order of the lists | one `ts.write_step(t, point_data={...}, cell_data={...})` per time step |
| `ts.close_unstructuredGridToVTK()` | `ts.close()`, or a `with` block |
| `ts.write_fake_grid_for_time(path)` | `time_values_file=` in the constructor: on by default (`<path>_time_values.vtu`), or a file name, or `False` |

The new file stores the data time step by time step, and its `TimeValues` holds the actual times instead of 0, 1, 2, .... VTK reads both files to the same data at every step (checked with VTK 9.5). `examples/unstructured_timedep.py` shows the new class.

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

## Reading the files with meshio

meshio 5.3.5 (the current release) fails to read many `.vtu` files with raw appended binary data, as VTKwrite writes them, when it runs on Python 3.13 or later. It fails with `binascii.Error: Incorrect padding`, depending on the sizes of the arrays. This is a meshio problem: it also fails on raw files written by VTK itself. meshio converts the raw data to base64 block by block and then decodes from each block to the end of the section as one string, which works only with the base64 decoder of Python 3.12 and earlier, which stops at the first padding character. With Python 3.10, 3.11 and 3.12, meshio 5.3.5 read the same files correctly. VTK and ParaView are not affected.
