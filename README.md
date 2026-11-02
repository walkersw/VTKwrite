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

## Testing

This test writes files and reads them back (with VTK's own readers where available):
```
python -m pip install pytest vtk
python -m pytest tests
```

## Known issue: VTK 9.6 and later

VTK 9.6 and later (including ParaView releases built on it) read only the first time step of a single-file time series written by `timeseries_unstructuredGrid`; later steps show wrong values. This is a bug in VTK's XML reader (reported upstream with a fix). Until a fixed VTK is released, use ParaView/VTK 9.5 or earlier for these files.
