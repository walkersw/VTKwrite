# VTKwrite

Python package to write VTK files in XML format

- This class can write a static unstructured grid, with time-dependent data, all in the same .vtu file!

- This package was forked from another package "PyEVTK" by Paulo A. Herrera.

https://github.com/paulo-herrera/PyEVTK

- Installation is done in the standard way:  ```python setup.py install```

- You need Python 3.6 or better.

- Compatibility: the generated files can be read by Paraview versions: 5.4, 5.5, 5.6, 5.7.  Eventually, they can be read by Paraview version 5.10.11, because of a recent bug fix.

## Local Install Help

If you need to make a "local" install (because you are not an administrator), then the install command will look like:

```
python3 setup.py install --prefix ~/.local
```

where you need to modify the ```~/.local``` text to reflect your system setup.

Note: the command above will install it to this directory: ```~/.local/lib/python3.8/site-packages```

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
