# Changes

## Unreleased

Fixes for files that were silently wrong or unreadable (review section 2):

- **Numbers in attributes**: `Origin`/`Spacing` (and `.pvd` time values) are written as plain numbers. With NumPy 2, `repr(np.float64(0.5))` produced `np.float64(0.5)`, which VTK read as 0.
- **NumPy arrays as `origin`/`spacing`** no longer raise "truth value of an array is ambiguous".
- **XML escaping**: array names and other attributes are escaped (`&`, `<`, `"`, ...); comments containing `--` are made legal. Names may contain non-ASCII characters (files are UTF-8).
- **Active attributes**: `PointData`/`CellData` now carry `Scalars=`/`Vectors=`/... (VTK matches these names case-sensitively, so the previous lowercase names were ignored and no active scalars/vectors were set).
- **No partial files**: files are written to a temporary file and renamed when complete. An error part-way through leaves no file behind and does not overwrite an existing file. `VtkFile` and `VtkGroup` can be used as context managers (`with VtkFile(...) as w:`), which save on success and discard on error.
- **Input checks**: all high-level functions check their inputs before opening the file and raise `ValueError`/`TypeError` with a message naming the problem. This covers unknown data types (e.g. `"scalar"`), wrong array lengths, duplicate names, unsupported dtypes (e.g. `bool`), mismatched `x`/`y`/`z`, and inconsistent `connectivity`/`offsets`/`cell_types`. These checks no longer use `assert`, so they also run under `python -O`.
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

Still open (see the review): API conveniences (e.g. (n, 3) arrays and dictionaries for the other writers, FEM cell types), and packaging metadata.
