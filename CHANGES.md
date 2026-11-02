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

Still open (see the review): write speed, a streaming time-series interface, and packaging metadata.
