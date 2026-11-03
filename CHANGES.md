# Changes

## 0.1.0 (unreleased)

Fixes for files that were silently wrong or unreadable (review section 2):

- **Numbers in attributes**: `Origin`/`Spacing` (and `.pvd` time values) are written as plain numbers. With NumPy 2, `repr(np.float64(0.5))` produced `np.float64(0.5)`, which VTK read as 0.
- **NumPy arrays as `origin`/`spacing`** no longer raise "truth value of an array is ambiguous".
- **XML escaping**: array names and other attributes are escaped (`&`, `<`, `"`, ...); comments containing `--` are made legal. Names may contain non-ASCII characters (files are UTF-8).
- **Active attributes**: `PointData`/`CellData` now carry `Scalars=`/`Vectors=`/... (VTK matches these names case-sensitively, so the previous lowercase names were ignored and no active scalars/vectors were set).
- **No partial files**: files are written to a temporary file and renamed when complete. An error part-way through leaves no file behind and does not overwrite an existing file. `VtkFile` and `VtkGroup` can be used as context managers (`with VtkFile(...) as w:`), which save on success and discard on error.
- **Input checks**: all high-level functions check their inputs before opening the file and raise `ValueError`/`TypeError` with a message naming the problem. This covers unknown data types (e.g. `"scalar"`), wrong array lengths, duplicate names, unsupported dtypes (e.g. `float16`), mismatched `x`/`y`/`z`, and inconsistent `connectivity`/`offsets`/`cell_types`. These checks no longer use `assert`, so they also run under `python -O`.
- **Time series** (`timeseries_unstructuredGrid`): every `append_data` call is checked against the header written at initialization: the number of time steps, the number of values, and the dtype. A lossless dtype change (e.g. float32 to float64) is converted; anything else raises. Nothing is written by a failed call. `close_unstructuredGridToVTK` raises, and discards the file, if a variable was never appended.

Side effects of the checks (previously crashed or wrote bad files):

- Lists are accepted wherever the docstrings say so (`unstructuredGridToVTK`, `linesToVTK`, `polyLinesToVTK`, ...).
- `"tcoords"` data (1–3 components) works.
- Non-contiguous (strided) arrays are accepted.
- Coordinates with different dtypes are converted to a common dtype.
- The caller's data lists are no longer modified.

Added a `pytest` suite in `tests/`.

Remaining functional bugs (review section 3):

- `cylinderToVTK` works. It always raised `TypeError` before. It now takes `all_cell_data`/`all_point_data` like the other writers; the older `cellData`/`pointData` dictionaries are still accepted. Points and quads are built with NumPy, and the angles/levels use `linspace`, so floating-point round-off can no longer add an extra ring. Cells are written as Quad with integer topology.
- `pointsToVTKAsTIN`:
  - returns the file path (it returned `None`);
  - writes integer connectivity/offsets and UInt8 cell types (they were Float64);
  - uses the `data` argument (extra point data after "Elevation") and `comments`;
  - with `ndim=3`, writes the tetrahedra that `Delaunay` produces. Before, it kept the first 3 vertices of each tetrahedron and labelled them as triangles;
  - raises a clear `ImportError` if SciPy is missing.
- `linesToVTK` and `polyLinesToVTK` write one cell type per cell (they wrote one per point).
- Topology arrays given as floats (e.g. `np.zeros(n)`) are written as Int64 connectivity/offsets and UInt8 cell types. Cell type ids outside 0–255 are rejected.
- Removed `VtkFile.appendHeader`, which could only raise `AttributeError` (it referred to attributes and functions that do not exist).
- Examples `unstructured.py` and `unstructured_timedep.py` use integer topology arrays.

Added tests for these and a test that runs every example and reads its output back.

Write speed (review section 5):

- The binary data is written with NumPy instead of one `struct.pack` call per value, about 15–20× faster for unstructured grids and about 10× faster for structured grids. Writing 10^6 points + 10^6 tetrahedra + 2 point fields went from about 1.8 s to 0.1 s.
- Points are interleaved in chunks of 2^20, so the extra memory stays bounded for large meshes. Peak memory in the benchmark went from 310 MB to 203 MB.
- The output is byte-for-byte identical to before (checked on all examples, on 50 000-cell meshes with mixed dtypes, non-native byte order, C/Fortran-ordered and strided 3-D arrays, and in `tests/test_binary.py` against the original algorithm).
- `polyLinesToVTK` computes offsets with `np.cumsum` instead of a Python loop.

Streaming time series (review section 4):

- New class `TimeSeriesUnstructuredGrid` (in `VTKwrite.interface`). It writes a simulation on a fixed unstructured mesh one step at a time: `write_step(t, point_data={...}, cell_data={...})`.
  - Each step goes to disk as it is written, so the history does not have to be kept in memory.
  - `output="vtu"` (default) writes everything into one .vtu file, as before: the mesh once, one array per field per step. It also writes the actual time values to `TimeValues` and to `<name>_time_values.vtu`. The time values do not have to be known in advance; if they are (`times=`), the data goes straight into the file, otherwise it is assembled when the series is closed.
  - `output="pvd"` writes `<name>.pvd` and one .vtu per step in `<name>/`, with the actual time values. The .pvd is updated after every step, so a run that stops early can still be viewed. This also works with VTK/ParaView versions that have the time-series bug.
  - Fields are dictionaries `{name: array}`. The type is inferred from the number of components (1: scalars, 3: vectors, 9: tensors), or given as `{name: (type, array)}`. Vectors may be arrays of shape (n, 3). Fields that do not change can be given once, as `static_point_data` / `static_cell_data`.
  - Every step is checked against the first (same fields, sizes and types; dtypes converted when safe), and a failed check writes nothing. Can be used as a context manager: on an error, a single .vtu is discarded, while for .pvd the finished steps are kept.
- The existing `timeseries_unstructuredGrid` class is unchanged and still works.
- New example `examples/unstructured_timeseries.py` (both outputs, non-uniform time steps, a static cell field).
- `VtkGroup` works when the .pvd path has no folder part (e.g. `VtkGroup("sol")`).
- `VtkFile.addHeader` accepts an explicit `offset`.

VTKHDF output:

- `TimeSeriesUnstructuredGrid(..., output="vtkhdf")` writes the whole time series into one VTKHDF (HDF5) file, `<name>.vtkhdf`. It is not affected by the VTK 9.6 XML time-series bug.
  - The mesh and the static fields are stored once, and each time-dependent field is one dataset holding all steps.
  - The actual time values are stored in the file, and steps are written as they come.
  - Checked with VTK 9.4, 9.5, 9.6 (the version in ParaView 6.1) and 9.7.
  - Needs `h5py` (optional dependency: `pip install VTKwrite[hdf5]` or `pip install h5py`).
- README: a fuller note on the VTK 9.6 / ParaView 6.1 single-file `.vtu` problem, and a table of the three time-series outputs.

Examples:

- `examples/group.py` now writes the step files that `group.pvd` lists (before, the .pvd pointed at files that did not exist). It shows when to use `VtkGroup`: a mesh that changes in time (here a stretching domain, re-meshed more finely at each step).
- `tests/test_examples.py` checks that every file listed in a .pvd exists.

Packaging (review section 7):

- Version 0.1.0. It is set in one place, `VTKwrite.__version__` in `src/VTKwrite/__init__.py`, which `pyproject.toml` reads.
- License metadata in the current format: `license = "GPL-3.0-only AND MIT"` (an SPDX expression: GPL-3.0-only for VTKwrite, MIT for the code derived from PyEVTK) and `license-files`. The license classifier is removed: it was not a valid trove classifier, so PyPI would have rejected an upload. Building no longer prints setuptools deprecation warnings. Building needs setuptools 77 or later (pip installs it automatically).
- The MIT notice of PyEVTK, on which VTKwrite is based, is also shipped as `LICENSES/MIT-PyEVTK.txt`, so it is part of every wheel and source distribution. It is still in `__init__.py` too.
- Requires Python 3.10 or later. The metadata said 3.6, which current NumPy and setuptools do not support. Tested with Python 3.10 and 3.13.
- New extra `test`: `pip install -e ".[test]"` installs pytest, vtk, h5py and scipy. The `hdf5` extra is unchanged.
- Removed `setup.py` (empty) and `setup_py.BAK` (out of date). Install with `pip install .`; `python setup.py install` is no longer supported.
- The source distribution now contains the examples, `tests/conftest.py` and `CHANGES.md` (new `MANIFEST.in`), so its tests can be run. Before, it had the test files but not `conftest.py` or the examples they need.
- Added a `.gitignore`.
- README: a new installation section (`pip install .`, `pip install -e .`, `pip install --user .`, the extras), and the tests are installed through the `test` extra.

Package namespace (review section 6):

- Everything public can be imported from the package itself, e.g. `from VTKwrite import unstructuredGridToVTK, TimeSeriesUnstructuredGrid, VtkTriangle`. This covers the high-level writers, both time-series classes, `VtkFile`, `VtkGroup`, the file types (`VtkUnstructuredGrid`, ...), `VtkCellType` and the cell types. `VTKwrite.__all__` lists them. The older imports from `VTKwrite.interface` and `VTKwrite.vtkbin` still work.
- `VTKwrite.__version__` gives the version.
- `bool` data arrays (e.g. boundary markers) are accepted as point and cell data, in all writers and both time-series classes, and are written as `UInt8` (0 and 1). Before, they raised `TypeError`. Coordinates and topology arrays of dtype `bool` are still rejected, as are dtypes that VTK files cannot hold (e.g. `float16`).
- Importing VTKwrite without NumPy raises `ImportError`. Before, a message was printed and the code failed later with `NameError`.
- `VtkFile(..., largeFile=...)`: the parameter is documented as ignored, because the files always use 64-bit block sizes. It is kept so that existing calls still work.
- `VTKwrite.interface` imports the names it needs from `vtkbin` explicitly, instead of `from .vtkbin import *`. Every name that could be imported from `VTKwrite.interface` before still can, except the modules `sys` and `numbers`.
- New tests: `tests/test_namespace.py` (exports, version, old import paths, cell type ids against VTK's constants) and `tests/test_bool_data.py`.

Friendlier array input (review section 6):

- Every writer accepts point and cell data as a dictionary, as `TimeSeriesUnstructuredGrid` does: `all_point_data={"u": u, "grad_u": G}`. The type follows from the number of components per point or cell (1: scalars, 3: vectors, 9: tensors), or is given as `{name: (type, array)}`. The order of the dictionary is the order in the file. This covers `unstructuredGridToVTK`, `pointsToVTK`, `pointsToVTKAsTIN` (`data`), `linesToVTK`, `polyLinesToVTK`, `cylinderToVTK`, `imageToVTK`, `rectilinearToVTK` and `structuredToVTK`. The list form `[name, type, array]` still works.
- Multi-component data may have one row per point or cell, e.g. shape `(npoints, 3)`, in both forms. Before, it had to be interleaved by hand into a 1-D array.
- Structured, rectilinear and image grids accept arrays of shape `(nx, ny, nz, k)` for k components, indexed `[i, j, k, component]`. The order VTK expects (component fastest, then i, j, k) is produced internally; before, it had to be done by hand, as in `examples/structured_ex1.py`.
- `imageToVTK` and `rectilinearToVTK` accept vectors, like `structuredToVTK`. Before, they accepted scalars only.
- The points can be given as one array with the keyword `points=`: shape `(n, 3)`, or `(n, 2)` for z = 0, in `unstructuredGridToVTK`, `pointsToVTK`, `pointsToVTKAsTIN`, `linesToVTK`, `polyLinesToVTK` and `TimeSeriesUnstructuredGrid`; `(nx, ny, nz, 3)` or `(nx, ny, nz, 2)` in `structuredToVTK`. So that existing positional calls keep working, `x, y, z` and the arguments after them now default to `None`. A missing argument still raises `TypeError`, which names it.
- The new forms give byte-for-byte the same files as the older ones; this is tested for every writer and for all three time-series outputs. Writing is as fast as before and needs no more memory (10^6 points and 10^6 tetrahedra: 0.09 s and 23 MB extra in every form), because the XML writers no longer copy the columns of a points array.
- Multi-component data given as a list or tuple of arrays, e.g. `(vx, vy, vz)`, raises `TypeError` with a hint (`np.column_stack`) in all writers, as it already did in `TimeSeriesUnstructuredGrid`: with 3 points it could not be told apart from 3 rows. (In the list form it was a `ValueError` before.)
- In a dictionary, a 3-D array that does not have the grid's shape must hold one value per point or cell; it is read in Fortran order, as in the list form.
- `cylinderToVTK`: dictionaries given as `all_cell_data`/`all_point_data` are read like any field dictionary, which gives the same files as before for scalar data. `cellData`/`pointData` are unchanged.
- `timeseries_unstructuredGrid` is unchanged: it still takes only the list form with 1-D arrays (it is to be deprecated).
- New example `examples/array_input.py`: a thick-walled cylinder with its displacement and stress tensor, written as an unstructured and as a structured grid. The README has a new section on these forms, and the module docstring of `VTKwrite.interface` describes them in full.
- New tests in `tests/test_array_input.py`. The structured, rectilinear and image grids are checked entry by entry with VTK's own (i, j, k) to id functions.

FEM cell types (review section 6):

- New cell types: `VtkQuadraticWedge` (26), `VtkQuadraticPyramid` (27), `VtkBiQuadraticQuad` (28), `VtkTriQuadraticHexahedron` (29), `VtkBiQuadraticTriangle` (34), `VtkCubicLine` (35), `VtkTriQuadraticPyramid` (37), and the arbitrary-order Lagrange cells `VtkLagrangeCurve` (68), `VtkLagrangeTriangle` (69), `VtkLagrangeQuadrilateral` (70), `VtkLagrangeTetrahedron` (71), `VtkLagrangeHexahedron` (72) and `VtkLagrangeWedge` (73). There is no Lagrange pyramid: VTK defines type 74 but reads it as an empty cell.
- New module `VTKwrite.lagrange` (NumPy only; also importable from `VTKwrite`):
  - `lagrange_nodes(cell_type, degree, vertices=None)`: the reference coordinates of the nodes of a Lagrange cell in VTK's order, on VTK's reference cell or mapped to another one;
  - `lagrange_npoints(cell_type, degree)`;
  - `node_permutation(cell_type, degree, nodes, vertices=None)`: the order in which to take the nodes of another element so that they are in VTK's order. It raises `ValueError` if the nodes are not VTK's equispaced nodes (e.g. Gauss-Lobatto points), in which case the fields must be evaluated at `lagrange_nodes(...)` instead;
  - `cellwise_mesh(cell_points, cell_type)`: a mesh in which every cell has its own nodes (discontinuous fields), returned as keyword arguments for `unstructuredGridToVTK` and `TimeSeriesUnstructuredGrid`.
- **Lagrange hexahedra were read wrongly.** VTK's XML reader swaps the nodes of hexahedron edges 10 and 11 in files that declare a `VTKFile` version below 2.1, and VTKwrite declared 1.0. Files that contain Lagrange hexahedra now declare 2.1 (checked with VTK 9.2.6, 9.5.2 and 9.7.1; the `.vtkhdf` output is not affected). All other files still declare 1.0 and are unchanged. `VtkFile` has a new keyword argument `version`.
- **The number of points of every cell is checked against its type** (`ValueError`, naming the first wrong cell). This covers the cell types with a fixed number of points (checked against VTK), the minimum for poly-vertices, polylines, triangle strips and polygons, and the Lagrange cells (including VTK's bubble-enriched 7-node triangle, 15-node tetrahedron and 21-node wedge). Lagrange quadrilaterals, hexahedra and wedges are not checked if the cell data has VTK's `"HigherOrderDegrees"` array (different degrees per direction). Before, a wrong count gave a broken file without an error.
- New example `examples/lagrange_cells.py`: a curved quarter annulus as degree-3 Lagrange quadrilaterals with shared nodes, and as degree-2 Lagrange triangles with a discontinuous field. README: a section on higher-order cells.
- New tests in `tests/test_lagrange.py`, for degrees 1 to 3. The node order is compared with VTK's own index functions, and a polynomial map of degree p, written as one order-p cell, must be reproduced by VTK to round-off at random points. The node order alone is also compared for degrees 4 to 6, where the orientation of the nodes inside a face first matters.

Collections (`VtkGroup`):

- `addFile` stores the path of a data file with `/` as the separator, so a `.pvd` written on Windows can also be read elsewhere. On Linux and macOS nothing changes.
- `addFile` no longer fails on Windows when the data file is on another drive than the `.pvd` (`os.path.relpath` raised `ValueError`); the absolute path is stored instead.
- New optional argument `name`: written as the `name` attribute, which ParaView uses as the block name.
- Documented how several files per time step are shown (one block per file), and that ParaView's PVD reader ignores `group` and `part`. This is taken from ParaView's source code and was not tested in ParaView.

Deprecated:

- `timeseries_unstructuredGrid` issues a `DeprecationWarning` and will be removed in the next release. Use `TimeSeriesUnstructuredGrid`:
  - `timeseries_unstructuredGrid(path, times)` plus `init_unstructuredGridToVTK(x, y, z, ...)` becomes the constructor `TimeSeriesUnstructuredGrid(path, x, y, z, ..., times=times)`;
  - the `append_data([...])` calls, one per variable, become one `write_step(t, point_data={...}, cell_data={...})` call per time step;
  - `close_unstructuredGridToVTK()` becomes `close()`, or a `with` block;
  - `write_fake_grid_for_time(path)` is replaced by `time_values_file=` (on by default).
  The README has the same migration note.
- `examples/unstructured_timedep.py` now uses `TimeSeriesUnstructuredGrid`. Its file `unstructured_timedep_all.vtu` has a different layout (steps stored time step by time step, the actual times in `TimeValues`), but VTK 9.5 reads it to the same arrays at every step as before. In `unstructured_timedep_time_values.vtu`, only a comment changed (the associated file name, without `./`). `unstructured_timedep_grid_only.vtu` is unchanged.

Documentation:

- README: a section on `VtkGroup`, and a note on meshio. meshio 5.3.5 fails with `Incorrect padding` on many raw-appended `.vtu` files on Python 3.13 and later, including files written by VTK itself; it reads them with Python 3.12 and earlier.

Still open: removing `timeseries_unstructuredGrid` in the next release. Not planned: a .vtp writer and compression.
