"""
VTKwrite.interface.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
High level Python library to export data to binary VTK file.

All inputs are checked before a file is opened, and errors raise ValueError or
TypeError (not assert, so the checks also run under "python -O"). A file is only
created under its final name once it has been written completely.

Copyright (c) 05-07-2021,  Shawn W. Walker
"""

from .vtkbin import * # VtkFile, VtkUnstructuredGrid, etc.
from .vtkbin import _vtk_data_type, _number_to_string
from .xmlwrite import sanitize_comment, XmlWriter
from .pyvtk import writeBlockSize, writeArrayToFile
import os
import shutil
import uuid
try:
    import numpy as np
except:
    print("Numpy is not installed. Please install it before running VTKwrite again.")

# =================================
#       Input checking
# =================================
# number of components for each data type ("tcoords" may have 1, 2 or 3)
_NCOMP = {"scalars" : 1, "vectors" : 3, "normals" : 3, "tensors" : 9, "tcoords" : None}
_ALL_TYPES = ("scalars", "vectors", "normals", "tensors", "tcoords")

def _to_array(what, a):
    """ Convert a list/tuple/array to a numpy array with a dtype VTK supports. """
    if isinstance(a, np.ndarray):
        arr = a
    elif isinstance(a, (list, tuple)):
        arr = np.array(a)
    else:
        raise TypeError("%s must be a numpy array, list or tuple, got %s"
                        % (what, type(a).__name__))
    if arr.dtype.byteorder not in ("=", "|"):
        arr = arr.astype(arr.dtype.newbyteorder("="))
    try:
        _vtk_data_type(arr.dtype.name)
    except TypeError as e:
        raise TypeError("%s: %s" % (what, e)) from None
    return arr

def _to_vector(what, a):
    """ A 1-D array (lists and tuples are converted). """
    arr = _to_array(what, a)
    if arr.ndim != 1:
        raise ValueError("%s must be 1-D, got shape %s" % (what, arr.shape))
    return arr

def _check_xyz(x, y, z):
    """ Coordinates of points given as three arrays with the same number of entries.
        Multi-dimensional arrays are accepted and taken in Fortran order (as before).
        If the dtypes differ, all three are converted to a common dtype. """
    x, y, z = _to_array("x", x), _to_array("y", y), _to_array("z", z)
    if not (x.size == y.size == z.size):
        raise ValueError("x, y, z must have the same number of entries, got %d, %d, %d"
                         % (x.size, y.size, z.size))
    if not (x.dtype == y.dtype == z.dtype):
        common = np.result_type(x, y, z)
        x, y, z = x.astype(common), y.astype(common), z.astype(common)
    return x, y, z

def _cells_per_direction(npoints_per_direction):
    """ Cells of a structured grid in each direction; a direction with a single
        point (a flat grid) still has one layer of cells, as in VTK. """
    return tuple(max(int(n) - 1, 1) for n in npoints_per_direction)

def _check_fields(fields, where, ntuples, allowed = _ALL_TYPES, shape3d = None):
    """ Check a list of [name, type, data] entries and return a NEW list of
        [name, type, data_array, ncomp] (the caller's list is not modified).

        where:   "all_cell_data" or "all_point_data" (used in error messages).
        ntuples: number of cells or points the data must cover.
        shape3d: if given, 3-D scalar arrays are also accepted and must have this shape.
    """
    if fields is None:
        return []
    if not isinstance(fields, (list, tuple)):
        raise TypeError("%s must be a list of [name, type, data] entries" % where)
    checked = []
    names = set()
    for ii, entry in enumerate(fields):
        if not isinstance(entry, (list, tuple)) or len(entry) != 3:
            raise TypeError("%s[%d] must be a list [name, type, data]" % (where, ii))
        name, dtype_str, data = entry
        if not isinstance(name, str) or not name:
            raise TypeError("%s[%d]: name must be a non-empty string, got %r" % (where, ii, name))
        if name in names:
            raise ValueError("%s: the name '%s' is used twice" % (where, name))
        names.add(name)
        if dtype_str not in allowed:
            raise ValueError("%s['%s']: type %r is not allowed here; use one of %s"
                             % (where, name, dtype_str, ", ".join(repr(t) for t in allowed)))
        arr = _to_array("%s['%s']" % (where, name), data)
        if arr.ndim not in (1, 3):
            raise ValueError("%s['%s'] must be a 1-D array%s, got shape %s"
                             % (where, name,
                                " (or a 3-D array of shape %s for scalars)" % (tuple(shape3d),)
                                if shape3d is not None else "", arr.shape))
        if arr.ndim == 3 and shape3d is not None and dtype_str == "scalars":
            if tuple(arr.shape) != tuple(shape3d):
                raise ValueError("%s['%s'] has shape %s, expected %s"
                                 % (where, name, arr.shape, tuple(shape3d)))
        ncomp = _NCOMP[dtype_str]
        if ncomp is None: # tcoords: 1, 2 or 3 components
            ncomp = arr.size // ntuples if ntuples > 0 else 0
            if ncomp not in (1, 2, 3) or ncomp * ntuples != arr.size:
                raise ValueError("%s['%s'] (tcoords) has %d values; expected 1, 2 or 3 per tuple "
                                 "for %d tuples" % (where, name, arr.size, ntuples))
        elif arr.size != ncomp * ntuples:
            raise ValueError("%s['%s'] has %d values, expected %d (%d %s x %d component%s)"
                             % (where, name, arr.size, ncomp * ntuples, ntuples,
                                "cells" if "cell" in where else "points", ncomp,
                                "" if ncomp == 1 else "s"))
        checked.append([name, dtype_str, arr, ncomp])
    return checked

def _check_comments(comments):
    if comments is None:
        return None
    if isinstance(comments, str):
        comments = [comments]
    return [str(c) for c in comments]

def _check_topology(connectivity, offsets, cell_types, npoints):
    """ Check the unstructured-grid topology arrays against each other. """
    connectivity = _to_vector("connectivity", connectivity)
    offsets = _to_vector("offsets", offsets)
    cell_types = _to_vector("cell_types", cell_types)
    for what, arr in (("connectivity", connectivity), ("offsets", offsets),
                      ("cell_types", cell_types)):
        if arr.size and arr.dtype.kind == "f" and not np.all(arr == np.round(arr)):
            raise ValueError("%s must contain whole numbers" % what)
    ncells = cell_types.size
    if offsets.size != ncells:
        raise ValueError("offsets has %d entries but cell_types has %d (one per cell)"
                         % (offsets.size, ncells))
    if ncells:
        if np.any(np.diff(offsets) < 0) or offsets[0] < 0:
            raise ValueError("offsets must be non-negative and non-decreasing "
                             "(offsets[i] = index one past the last vertex of cell i)")
        if int(offsets[-1]) != connectivity.size:
            raise ValueError("offsets[-1] = %d must equal the length of connectivity (%d)"
                             % (int(offsets[-1]), connectivity.size))
    if connectivity.size and (connectivity.min() < 0 or connectivity.max() >= npoints):
        raise ValueError("connectivity refers to point %d, but there are only %d points"
                         % (int(connectivity.max() if connectivity.max() >= npoints
                                else connectivity.min()), npoints))
    if cell_types.size and (cell_types.min() < 0 or cell_types.max() > 255):
        raise ValueError("cell_types must be VTK cell type ids between 0 and 255")
    # VTK expects integer topology arrays and UInt8 cell types; floats (e.g. from
    # np.zeros(n)) used to be written as Float64, which some readers reject.
    if connectivity.dtype.kind not in "iu":
        connectivity = connectivity.astype(np.int64)
    if offsets.dtype.kind not in "iu":
        offsets = offsets.astype(np.int64)
    if cell_types.dtype != np.uint8:
        cell_types = cell_types.astype(np.uint8)
    return connectivity, offsets, cell_types

# =================================
#       Helper functions
# =================================
from io import StringIO as _TextBuffer

def _write_text_atomically(path, text):
    """ Write text to path via a temporary file, so a failure leaves no partial file. """
    with XmlWriter(path, addDeclaration = False) as xml:
        xml.addText(text)

def _addDataToFile(vtkFile, all_cell_data, all_point_data, time_steps = None):
    """ Write the DataArray headers. Entries are [name, type, data, ncomp] as
        returned by _check_fields. """
    dt_strs = ["scalars", "vectors", "normals", "tensors", "tcoords"]

    all_cell_data = all_cell_data or []
    all_point_data = all_point_data or []

    for assoc, fields in (("Cell", all_cell_data), ("Point", all_point_data)):
        if not fields:
            continue
        # default ("active") variable for each data type: the first of that type
        default_vars = {}
        for dt in dt_strs:
            for entry in fields:
                if entry[1] == dt:
                    default_vars[dt] = entry[0]
                    break
        vtkFile.openData(assoc, **default_vars)
        for name, dt, data, ncomp in fields:
            if time_steps is not None:
                for ti in time_steps:
                    vtkFile.internal_addData(name, data, ncomp, str(ti))
            else:
                vtkFile.internal_addData(name, data, ncomp)
        vtkFile.closeData(assoc)

def _appendDataToFile(vtkFile, all_cell_data, all_point_data):
    # append data to binary section, in the same order as the headers
    for entry in (all_cell_data or []):
        vtkFile.appendData(entry[2])
    for entry in (all_point_data or []):
        vtkFile.appendData(entry[2])

def __convertListToArray(list1d):
    ''' If data is a list and no a Numpy array, then it convert it
        to an array, otherwise return the same array '''
    if (list1d is not None) and (not isinstance(list1d, np.ndarray)):
        if not isinstance(list1d, (list, tuple)):
            raise TypeError("expected a list, tuple or numpy array, got %s" % type(list1d).__name__)
        return np.array(list1d)
    else:
        return list1d

def __convertDictListToArrays(data):
    ''' If data in dictironary are lists and no a Numpy array,
        then it creates a new dictionary and convert the list to arrays,
        otherwise return the same dictionary '''
    if data is not None:
        dict = {}
        for k, list1d in data.items():
            dict[k] = __convertListToArray(list1d)
        return dict
    else:
        return data # None

def _write_unstructured(path, x, y, z, connectivity, offsets, cell_types,
                        cell_fields, point_fields, comments):
    """ Write a .vtu file from already checked inputs. """
    npoints = x.size
    ncells = offsets.size   # one offset per cell
    with VtkFile(path, VtkUnstructuredGrid) as w:
        if comments: w.addComments(comments)
        w.openGrid()
        w.openPiece(ncells = ncells, npoints = npoints)

        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity)
        w.addData("offsets", offsets)
        w.addData("types", cell_types)
        w.closeElement("Cells")

        _addDataToFile(w, all_cell_data = cell_fields, all_point_data = point_fields)

        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )
        w.appendData(connectivity).appendData(offsets).appendData(cell_types)

        _appendDataToFile(w, all_cell_data = cell_fields, all_point_data = point_fields)
    return w.getFileName()

def _write_time_values_grid(filepath, time_values, associated_file):
    """ Write filepath + '.vtu': a one-line "fake" grid whose cell array "TimeValue"
        holds the actual time value of each time step (see write_fake_grid_for_time). """
    full_file = filepath + '.vtu'
    _write_text_atomically(full_file, _time_values_grid_text(time_values, associated_file))
    return full_file

def _time_values_grid_text(time_values, associated_file):
    """ Content of the "fake" grid file written by _write_time_values_grid. """
    time_file = _TextBuffer()
    time_file.write('<?xml version="1.0"?>\n')
    time_file.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
    time_file.write('<!-- This is a fake grid, that is only used for plotting actual Time Values -->\n')
    Comment_str_2 = '<!-- ' + sanitize_comment('See the associated file: ' + associated_file) + ' -->\n'
    time_file.write(Comment_str_2)
    time_file.write('<!-- This .vtu file can be included in the Paraview pipeline. -->\n')
    time_file.write('<!-- Then use a Python Annotation Filter with this Expression: "Time: %1.2f" %TimeValue[0] -->\n')
    time_file.write('<!-- Make sure the "Array Association" is set to "Cell Data" -->\n')
    time_file.write('<!-- This is useful when the time step spacing is not uniform. -->\n')
    time_file.write('<!-- Note: make sure to hide the fake grid. -->\n')

    num_time_indices = len(time_values)
    time_indices = np.arange(num_time_indices)
    time_ind_str = ' '.join(map(str, time_indices))

    UG_TV_str = '<UnstructuredGrid TimeValues="' + time_ind_str + '">\n'
    time_file.write(UG_TV_str)
    time_file.write('<Piece NumberOfPoints="2" NumberOfCells="1">\n')
    time_file.write('<Points>\n')
    time_file.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
    time_file.write('0 0 0\n')
    time_file.write('1 0 0\n')
    time_file.write('</DataArray>\n')
    time_file.write('</Points>\n')
    time_file.write('<Cells>\n')
    time_file.write('<DataArray type="Int32" Name="connectivity" format="ascii">\n')
    time_file.write('0 1\n')
    time_file.write('</DataArray>\n')
    time_file.write('<DataArray type="Int32" Name="offsets" format="ascii">\n')
    time_file.write('2\n')
    time_file.write('</DataArray>\n')
    time_file.write('<DataArray type="UInt8" Name="types" format="ascii">\n')
    time_file.write('3\n')
    time_file.write('</DataArray>\n')
    time_file.write('</Cells>\n')
    time_file.write('<CellData Scalars="TimeValue">\n')

    for ii in range(num_time_indices):
        TV_ii_str = '<DataArray Name="TimeValue" NumberOfComponents="1" type="Float64" format="ascii" TimeStep="' + \
                               str(ii) + '">' + _number_to_string(float(time_values[ii])) + '</DataArray>\n'
        time_file.write(TV_ii_str)

    time_file.write('</CellData>\n')
    time_file.write('</Piece>\n')
    time_file.write('</UnstructuredGrid>\n')
    time_file.write('</VTKFile>\n')
    return time_file.getvalue()

# =================================
#       High level functions
# =================================
def imageToVTK(path, origin = (0.0,0.0,0.0), spacing = (1.0,1.0,1.0), all_cell_data = None, all_point_data = None, comments = None ):
    """ Exports data values as a rectangular image.
        
        PARAMETERS:
            path: name of the file without extension where data should be saved.
            origin: grid origin (default = (0,0,0))
            spacing: grid spacing (default = (1,1,1))
            
            all_cell_data: A List of Lists.  It has this format:
                    all_cell_data[ii] = cellData, where
                    cellData is a List (of length 3) that defines a variable associated with the grid cells:
                        cellData[0] = the *name* of the variable stored;
                        cellData[1] = the *type* of the variable, only "scalars" is allowed here.
                        cellData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_cell_data" must have the same number of elements (number of cells).
                    Note: the length of "all_cell_data" is the number of variables (defined on cells).
            all_point_data: A List of Lists.  It has this format:
                    all_point_data[ii] = pointData, where
                    pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                        pointData[0] = the *name* of the variable stored;
                        pointData[1] = the *type* of the variable, only "scalars" is allowed here.
                        pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_point_data" must have the same number of elements (number of vertices).
                    Note: the length of "all_point_data" is the number of variables (defined on vertices).
            comments: list of comment strings, which will be added to the header section of the file.
         
         RETURNS:
            Full path to saved file.

        NOTE: At least, all_cell_data or all_point_data must be present to infer the dimensions of the image.
    """
    if not all_cell_data and not all_point_data:
        raise ValueError("imageToVTK needs all_cell_data or all_point_data to infer the image size")
    origin, spacing = list(origin), list(spacing)
    if len(origin) != 3 or len(spacing) != 3:
        raise ValueError("origin and spacing must have 3 entries, got %r and %r" % (origin, spacing))
    for v in origin + spacing:
        _number_to_string(v)  # raises for anything that is not a number

    # Extract dimensions (from 3-D arrays)
    start = (0,0,0)
    def _first_shape(fields, where):
        if not fields:
            return None
        entry = fields[0]
        if not isinstance(entry, (list, tuple)) or len(entry) != 3:
            raise TypeError("%s[0] must be a list [name, type, data]" % where)
        arr = _to_array("%s[0]" % where, entry[2])
        if arr.ndim != 3:
            raise ValueError("imageToVTK needs 3-D arrays; %s[0] has shape %s" % (where, arr.shape))
        return tuple(arr.shape)
    cshape = _first_shape(all_cell_data, "all_cell_data")
    pshape = _first_shape(all_point_data, "all_point_data")
    if pshape is not None:
        # the point data defines the grid (as before); a direction with one point is flat
        point_shape = pshape
        cell_shape = _cells_per_direction(pshape)
        if cshape is not None and cshape != cell_shape:
            raise ValueError("cell data shape %s does not match point data shape %s "
                             "(expected cell shape %s)" % (cshape, pshape, cell_shape))
    else:
        cell_shape = cshape
        point_shape = tuple(c + 1 for c in cshape)
    end = tuple(p - 1 for p in point_shape)
    if min(end) < 0:
        raise ValueError("invalid image size: point shape %s" % (point_shape,))
    ncells = int(np.prod(cell_shape)); npoints = int(np.prod(point_shape))

    cell_fields = _check_fields(all_cell_data, "all_cell_data", ncells, ("scalars",), cell_shape)
    point_fields = _check_fields(all_point_data, "all_point_data", npoints, ("scalars",), point_shape)
    comments = _check_comments(comments)

    # Write data to file
    with VtkFile(path, VtkImageData) as w:
        if comments: w.addComments(comments)
        w.openGrid(start = start, end = end, origin = origin, spacing = spacing)
        w.openPiece(start = start, end = end)
        _addDataToFile(w, cell_fields, point_fields)
        w.closePiece()
        w.closeGrid()
        _appendDataToFile(w, cell_fields, point_fields)
    return w.getFileName()

# ==============================================================================
def rectilinearToVTK(path, x, y, z, all_cell_data = None, all_point_data = None, comments = None):
    """
        Writes data values as a rectilinear or rectangular grid.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: coordinates of the nodes of the grid as 1D arrays.
                     The grid should be Cartesian, i.e. faces in all cells are orthogonal.
                     Arrays size should be equal to the number of nodes of the grid in each direction.

            all_cell_data: A List of Lists.  It has this format:
                    all_cell_data[ii] = cellData, where
                    cellData is a List (of length 3) that defines a variable associated with the grid cells:
                        cellData[0] = the *name* of the variable stored;
                        cellData[1] = the *type* of the variable, only "scalars" is allowed here.
                        cellData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_cell_data" must have the same number of elements (number of cells).
                    Note: the length of "all_cell_data" is the number of variables (defined on cells).
            all_point_data: A List of Lists.  It has this format:
                    all_point_data[ii] = pointData, where
                    pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                        pointData[0] = the *name* of the variable stored;
                        pointData[1] = the *type* of the variable, only "scalars" is allowed here.
                        pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_point_data" must have the same number of elements (number of vertices).
                    Note: the length of "all_point_data" is the number of variables (defined on vertices).
            comments: list of comment strings, which will be added to the header section of the file.
            
        RETURNS:
            Full path to saved file.

    """
    x, y, z = _to_vector("x", x), _to_vector("y", y), _to_vector("z", z)
    if min(x.size, y.size, z.size) < 1:
        raise ValueError("x, y, z must each contain at least one coordinate")
    ftype = VtkRectilinearGrid
    nx, ny, nz = x.size - 1, y.size - 1, z.size - 1
    # Extract dimensions
    start = (0,0,0)
    end   = (nx, ny, nz)

    cell_shape = _cells_per_direction((x.size, y.size, z.size))
    cell_fields = _check_fields(all_cell_data, "all_cell_data", int(np.prod(cell_shape)),
                                ("scalars",), cell_shape)
    point_fields = _check_fields(all_point_data, "all_point_data", x.size * y.size * z.size,
                                 ("scalars",), (x.size, y.size, z.size))
    comments = _check_comments(comments)

    with VtkFile(path, ftype) as w:
        if comments: w.addComments(comments)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)

        w.openElement("Coordinates")
        w.addData("x_coordinates", x)
        w.addData("y_coordinates", y)
        w.addData("z_coordinates", z)
        w.closeElement("Coordinates")

        _addDataToFile(w, all_cell_data = cell_fields, all_point_data = point_fields)
        w.closePiece()
        w.closeGrid()
        # Write coordinates
        w.appendData(x).appendData(y).appendData(z)
        # Write data
        _appendDataToFile(w, all_cell_data = cell_fields, all_point_data = point_fields)
    return w.getFileName()


def structuredToVTK(path, x, y, z, all_cell_data = None, all_point_data = None, comments = None):
    """
        Writes data values as a rectilinear or rectangular grid.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: coordinates of the nodes of the grid as 3D arrays.
                     The grid should be structured, i.e. all cells should have the same number of neighbors.
                     Arrays size in each dimension should be equal to the number of nodes of the grid in each direction.
            all_cell_data: A List of Lists.  It has this format:
                    all_cell_data[ii] = cellData, where
                    cellData is a List (of length 3) that defines a variable associated with the grid cells:
                        cellData[0] = the *name* of the variable stored;
                        cellData[1] = the *type* of the variable, only "scalars" and "vectors" allowed here.
                        cellData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_cell_data" must have the same number of elements (number of cells).
                    Note: the length of "all_cell_data" is the number of variables (defined on cells).
            all_point_data: A List of Lists.  It has this format:
                    all_point_data[ii] = pointData, where
                    pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                        pointData[0] = the *name* of the variable stored;
                        pointData[1] = the *type* of the variable, only "scalars" and "vectors" allowed here.
                        pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_point_data" must have the same number of elements (number of vertices).
                    Note: the length of "all_point_data" is the number of variables (defined on vertices).
            comments: list of comment strings, which will be added to the header section of the file.
            
        RETURNS:
            Full path to saved file.

    """
    x, y, z = (_to_array("x", x), _to_array("y", y), _to_array("z", z))
    if not (x.ndim == 3 and y.ndim == 3 and z.ndim == 3):
        raise ValueError("x, y, z must be 3-D arrays, got shapes %s, %s, %s"
                         % (x.shape, y.shape, z.shape))
    if not (x.shape == y.shape == z.shape):
        raise ValueError("x, y, z must have the same shape, got %s, %s, %s"
                         % (x.shape, y.shape, z.shape))
    if not (x.dtype == y.dtype == z.dtype):
        common = np.result_type(x, y, z)
        x, y, z = x.astype(common), y.astype(common), z.astype(common)

    ftype = VtkStructuredGrid
    s = x.shape
    nx, ny, nz = s[0] - 1, s[1] - 1, s[2] - 1
    start = (0,0,0)
    end = (nx, ny, nz)

    cell_shape = _cells_per_direction(s)
    cell_fields = _check_fields(all_cell_data, "all_cell_data", int(np.prod(cell_shape)),
                                ("scalars", "vectors"), cell_shape)
    point_fields = _check_fields(all_point_data, "all_point_data", x.size,
                                 ("scalars", "vectors"), s)
    comments = _check_comments(comments)

    with VtkFile(path, ftype) as w:
        if comments: w.addComments(comments)
        w.openGrid(start = start, end = end)
        w.openPiece(start = start, end = end)
        w.openElement("Points")
        w.addData("points", (x,y,z))
        w.closeElement("Points")

        _addDataToFile(w, all_cell_data = cell_fields, all_point_data = point_fields)
        w.closePiece()
        w.closeGrid()
        w.appendData( (x,y,z) )
        _appendDataToFile(w, all_cell_data = cell_fields, all_point_data = point_fields)
    return w.getFileName()

# def gridToVTK(path, x, y, z, cellData = None, pointData = None, comments = None):
    # """
        # Writes data values as a rectilinear or rectangular grid.

        # PARAMETERS:
            # path: name of the file without extension where data should be saved.
            # x, y, z: coordinates of the nodes of the grid. They can be 1D or 3D depending if
                     # the grid should be saved as a rectilinear or logically structured grid, respectively.
                     # Arrays should contain coordinates of the nodes of the grid.
                     # If arrays are 1D, then the grid should be Cartesian, i.e. faces in all cells are orthogonal.
                     # If arrays are 3D, then the grid should be logically structured with hexahedral cells.
                     # In both cases the arrays dimensions should be equal to the number of nodes of the grid.
            # cellData: dictionary containing arrays with cell centered data.
                      # Keys should be the names of the data arrays.
                      # Arrays must have the same dimensions in all directions and must contain 
                      # only scalar data.
            # pointData: dictionary containing arrays with node centered data.
                       # Keys should be the names of the data arrays.
                       # Arrays must have same dimension in each direction and 
                       # they should be equal to the dimensions of the cell data plus one and
                       # must contain only scalar data.
            # comments: list of comment strings, which will be added to the header section of the file.
            
        # RETURNS:
            # Full path to saved file.

    # """
    # # Extract dimensions
    # start = (0,0,0)
    # nx = ny = nz = 0

    # if (x.ndim == 1 and y.ndim == 1 and z.ndim == 1):
        # nx, ny, nz = x.size - 1, y.size - 1, z.size - 1
        # isRect = True
        # ftype = VtkRectilinearGrid
    # elif (x.ndim == 3 and y.ndim == 3 and z.ndim == 3):
        # s = x.shape
        # nx, ny, nz = s[0] - 1, s[1] - 1, s[2] - 1
        # isRect = False
        # ftype = VtkStructuredGrid
    # else:
        # assert(False)
    # end = (nx, ny, nz)


    # w =  VtkFile(path, ftype)
    # if comments: w.addComments(comments)
    # w.openGrid(start = start, end = end)
    # w.openPiece(start = start, end = end)

    # if isRect:
        # w.openElement("Coordinates")
        # w.addData("x_coordinates", x)
        # w.addData("y_coordinates", y)
        # w.addData("z_coordinates", z)
        # w.closeElement("Coordinates")
    # else:
        # w.openElement("Points")
        # w.addData("points", (x,y,z))
        # w.closeElement("Points")

    # _addDataToFile(w, cellData, pointData)
    # w.closePiece()
    # w.closeGrid()
    # # Write coordinates
    # if isRect:
        # w.appendData(x).appendData(y).appendData(z)
    # else:
        # w.appendData( (x,y,z) )
    # # Write data
    # _appendDataToFile(w, cellData, pointData)
    # w.save()
    # return w.getFileName()

# ==============================================================================
def pointsToVTK(path, x, y, z, all_point_data = None, comments = None ):
    """
        Export points and associated data as an unstructured grid.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: 1D list-type object (list, tuple or numpy) with coordinates of the points.
            all_point_data: A List of Lists.  It has this format:
                    all_point_data[ii] = pointData, where
                    pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                        pointData[0] = the *name* of the variable stored;
                        pointData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                        pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_point_data" must have the same number of elements (number of vertices).
                    Note: the length of "all_point_data" is the number of variables (defined on vertices).
            comments: list of comment strings, which will be added to the header section of the file.
            
        RETURNS:
            Full path to saved file.

    """
    x, y, z = _check_xyz(x, y, z)
    npoints = x.size
    point_fields = _check_fields(all_point_data, "all_point_data", npoints)
    comments = _check_comments(comments)

    # create some temporary arrays to write grid topology
    offsets = np.arange(start = 1, stop = npoints + 1, dtype = 'int32') # index of last node in each cell
    connectivity = np.arange(npoints, dtype = 'int32')                 # each point is only connected to itself
    cell_types = np.empty(npoints, dtype = 'uint8')

    cell_types[:] = VtkVertex.tid

    return _write_unstructured(path, x, y, z, connectivity, offsets, cell_types,
                               [], point_fields, comments)

# ==============================================================================
def pointsToVTKAsTIN(path, x, y, z, data = None, comments = None, ndim = 2):
    """
        Export points and associated data as a triangular irregular grid.
        It builds a triangular grid that has the input points as nodes
        using the Delaunay triangulation function in Scipy, which requires
        a convex set of points (check the documentation for further details
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html).

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: 1D list-type object (list, tuple or numpy) with coordinates of the points.
            data: additional point data (optional), written after the "Elevation" array
                  (which holds z).  A List of Lists.  It has this format:
                    data[ii] = pointData, where
                    pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                        pointData[0] = the *name* of the variable stored;
                        pointData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                        pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "data" must have the same number of elements (number of vertices).
                    Note: the length of "data" is the number of variables (defined on vertices).
            comments: list of comment strings, which will be added to the header section of the file.
            ndim: is the number of dimensions considered when calling Delaunay.
                  If ndim = 2, then only coordinates x and y are passed, and the cells are triangles.
                  If ndim = 3, then x, y and z coordinates are passed, and the cells are tetrahedra.
            
        RETURNS:
            Full path to saved file.
        
        REQUIRES: Scipy > 1.2.
    """
    try:
        from scipy.spatial import Delaunay
    except ImportError:
        raise ImportError("pointsToVTKAsTIN needs scipy (scipy.spatial.Delaunay); "
                          "install it with 'pip install scipy'") from None

    if ndim not in (2, 3):
        raise ValueError("ndim must be 2 or 3, got %r" % (ndim,))
    x, y, z = _check_xyz(x, y, z)
    npts = x.size
    point_fields = _check_fields([["Elevation", "scalars", z]] + list(data or []),
                                 "data", npts)
    comments = _check_comments(comments)

    # Delaunay triangulation (ndim = 2) or tetrahedralization (ndim = 3)
    coords = (x, y) if ndim == 2 else (x, y, z)
    tri = Delaunay(np.column_stack([np.ravel(c, order = 'F') for c in coords]))
    simplices = tri.simplices                    # (ncells, ndim + 1)
    ncells, nverts = simplices.shape
    conn = simplices.astype(np.int64).ravel()
    offset = np.arange(1, ncells + 1, dtype = np.int64) * nverts
    tid = VtkTriangle.tid if ndim == 2 else VtkTetra.tid
    cell_type = np.full(ncells, tid, dtype = np.uint8)

    return _write_unstructured(path, x, y, z, conn, offset, cell_type,
                               [], point_fields, comments)

# ==============================================================================
def linesToVTK(path, x, y, z, all_cell_data = None, all_point_data = None, comments = None ):
    """
        Export line segments that join 2 points and associated data.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: 1D list-type object (list, tuple or numpy) with coordinates of the vertex of the lines. It is assumed that each line.
                     is defined by two points, then the length of the arrays should be equal to 2 * number of lines.
                     So each consecutive pair of points is a line.
            all_cell_data: A List of Lists.  It has this format:
                    all_cell_data[ii] = cellData, where
                    cellData is a List (of length 3) that defines a variable associated with the grid cells:
                        cellData[0] = the *name* of the variable stored;
                        cellData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                        cellData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_cell_data" must have the same number of elements (number of cells).
                    Note: the length of "all_cell_data" is the number of variables (defined on cells).
            all_point_data: A List of Lists.  It has this format:
                    all_point_data[ii] = pointData, where
                    pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                        pointData[0] = the *name* of the variable stored;
                        pointData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                        pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_point_data" must have the same number of elements (number of vertices).
                    Note: the length of "all_point_data" is the number of variables (defined on vertices).
            comments: list of comment strings, which will be added to the header section of the file.
                  
        RETURNS:
            Full path to saved file.

    """
    x, y, z = _check_xyz(x, y, z)
    if x.size % 2 != 0:
        raise ValueError("linesToVTK needs an even number of points (2 per line), got %d" % x.size)

    npoints = len(x)
    ncells = npoints // 2

    cell_fields = _check_fields(all_cell_data, "all_cell_data", ncells)
    point_fields = _check_fields(all_point_data, "all_point_data", npoints)
    comments = _check_comments(comments)

    # create some temporary arrays to write grid topology
    offsets = np.arange(start = 2, step = 2, stop = npoints + 1, dtype = 'int32') # index of last node in each cell
    connectivity = np.arange(npoints, dtype = 'int32')                          # each point is only connected to itself
    cell_types = np.full(ncells, VtkLine.tid, dtype = 'uint8')                  # one type per line

    return _write_unstructured(path, x, y, z, connectivity, offsets, cell_types,
                               cell_fields, point_fields, comments)

# ==============================================================================
def polyLinesToVTK(path, x, y, z, pointsPerLine, all_cell_data = None, all_point_data = None, comments = None ):
    """
        Export line segments that joint 2 points and associated data.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: 1D list-type object (list, tuple or numpy) arrays with coordinates of the vertices of the lines. It is assumed that each line.
                     has diffent number of points.
            pointsPerLine: 1D list-type object (list, tuple or numpy) array that defines the number of points associated to each line. Thus, 
                           the length of this array defines the number of lines. It also implicitly 
                           defines the connectivity or topology of the set of lines. It is assumed 
                           that points that define a line are consecutive in the x, y and z arrays.
            all_cell_data: A List of Lists.  It has this format:
                    all_cell_data[ii] = cellData, where
                    cellData is a List (of length 3) that defines a variable associated with the grid cells:
                        cellData[0] = the *name* of the variable stored;
                        cellData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                        cellData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_cell_data" must have the same number of elements (number of cells).
                    Note: the length of "all_cell_data" is the number of variables (defined on cells).
            all_point_data: A List of Lists.  It has this format:
                    all_point_data[ii] = pointData, where
                    pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                        pointData[0] = the *name* of the variable stored;
                        pointData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                        pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_point_data" must have the same number of elements (number of vertices).
                    Note: the length of "all_point_data" is the number of variables (defined on vertices).
            comments: list of comment strings, which will be added to the header section of the file.
            
        RETURNS:
            Full path to saved file.

    """
    x, y, z = _check_xyz(x, y, z)
    pointsPerLine = _to_vector("pointsPerLine", pointsPerLine)

    npoints = len(x)
    ncells = pointsPerLine.size
    if ncells and (pointsPerLine.min() < 1):
        raise ValueError("every entry of pointsPerLine must be at least 1")
    if int(np.sum(pointsPerLine)) != npoints:
        raise ValueError("pointsPerLine sums to %d, but there are %d points"
                         % (int(np.sum(pointsPerLine)), npoints))

    cell_fields = _check_fields(all_cell_data, "all_cell_data", ncells)
    point_fields = _check_fields(all_point_data, "all_point_data", npoints)
    comments = _check_comments(comments)

    # create some temporary arrays to write grid topology
    offsets = np.cumsum(pointsPerLine).astype('int32')  # index of last node in each cell

    connectivity = np.arange(npoints, dtype = 'int32')      # each line connects points that are consecutive

    cell_types = np.full(ncells, VtkPolyLine.tid, dtype = 'uint8')  # one type per polyline

    return _write_unstructured(path, x, y, z, connectivity, offsets, cell_types,
                               cell_fields, point_fields, comments)

def unstructuredGridToVTK(path, x, y, z, connectivity, offsets, cell_types, all_cell_data = None, all_point_data = None, \
                      comments = None):
    """
        Export unstructured grid and associated data.

        PARAMETERS:
            path: name of the file without extension where data should be saved.
            x, y, z: 1D list-type object (list, tuple or numpy) with coordinates of the vertices of cells. It is assumed that each element
                     has diffent number of vertices.
            connectivity: 1D list-type object (list, tuple or numpy) that defines the vertices associated to each element. 
                          Together with offset define the connectivity or topology of the grid. 
                          It is assumed that vertices in an element are listed consecutively.
            offsets: 1D list-type object (list, tuple or numpy) with the index of the last vertex of each element in the connectivity array.
                     It should have length nelem, where nelem is the number of cells or elements in the grid.
            cell_types: 1D list-type object (list, tuple or numpy) with an integer that defines the cell type of each element in the grid.
                        It should have size nelem. This should be assigned from pyvtk.vtkbin.VtkXXXX.tid, where XXXX represents
                        the type of cell. Please check the VTK file format specification for allowed cell types.                       
            all_cell_data: A List of Lists.  It has this format:
                    all_cell_data[ii] = cellData, where
                    cellData is a List (of length 3) that defines a variable associated with the grid cells:
                        cellData[0] = the *name* of the variable stored;
                        cellData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                        cellData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_cell_data" must have the same number of elements (number of cells).
                    Note: the length of "all_cell_data" is the number of variables (defined on cells).
            all_point_data: A List of Lists.  It has this format:
                    all_point_data[ii] = pointData, where
                    pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                        pointData[0] = the *name* of the variable stored;
                        pointData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                        pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "all_point_data" must have the same number of elements (number of vertices).
                    Note: the length of "all_point_data" is the number of variables (defined on vertices).
            comments: list of comment strings, which will be added to the header section of the file.
            
        RETURNS:
            Full path to saved file.

    """
    x, y, z = _check_xyz(x, y, z)
    npoints = x.size
    connectivity, offsets, cell_types = _check_topology(connectivity, offsets, cell_types, npoints)
    ncells = cell_types.size

    cell_fields = _check_fields(all_cell_data, "all_cell_data", ncells)
    point_fields = _check_fields(all_point_data, "all_point_data", npoints)
    comments = _check_comments(comments)

    return _write_unstructured(path, x, y, z, connectivity, offsets, cell_types,
                               cell_fields, point_fields, comments)

# ==============================================================================
def cylinderToVTK(path, x0, y0, z0, z1, radius, nlayers, npilars = 16, all_cell_data = None,
                  all_point_data = None, comments = None, cellData = None, pointData = None):
    """
        Export the lateral surface of a vertical cylinder as an unstructured grid of quads.

      PARAMETERS:
        path: path to file without extension.
        x0, y0: center of cylinder.
        z0, z1: lower and top elevation of the cylinder.
        radius: radius of cylinder.
        nlayers: Number of layers in z direction to divide the cylinder.
        npilars: Number of points around the circumference of the cylinder.
                 Higher value gives higher resolution to represent the curved shape.
        all_cell_data: A List of Lists, as in unstructuredGridToVTK, for the
                  ncells = npilars * nlayers quads.  Cell k * npilars + p is the quad between
                  pillars p and p+1 (counted counter-clockwise from the x-axis) in layer k.
        all_point_data: A List of Lists, as in unstructuredGridToVTK, for the
                  npoints = npilars * (nlayers + 1) points.  Point k * npilars + p is at
                  angle 2*pi*p/npilars and height level k.
        comments: list of comment strings, which will be added to the header section of the file.
        cellData, pointData: (older form) dictionaries {name: 1D array} of scalar data;
                  they are added to all_cell_data / all_point_data.

      RETURNS:
            Full path to saved file.

        NOTE: This function only export vertical shapes for now. However, it should be easy to
              rotate the cylinder to represent other orientations.
    """
    nlayers, npilars = int(nlayers), int(npilars)
    if nlayers < 1 or npilars < 3:
        raise ValueError("cylinderToVTK needs nlayers >= 1 and npilars >= 3")

    # Points: npilars around the circle on each of the nlayers + 1 levels
    angles = np.linspace(0.0, 2.0 * np.pi, npilars, endpoint = False)
    levels = np.linspace(z0, z1, nlayers + 1)
    xx = np.tile(x0 + radius * np.cos(angles), nlayers + 1)
    yy = np.tile(y0 + radius * np.sin(angles), nlayers + 1)
    zz = np.repeat(levels, npilars)

    # Quads (n0, n1, n3, n2), counter-clockwise seen from outside
    p0 = np.arange(npilars)
    p1 = (p0 + 1) % npilars
    layer = np.arange(nlayers)[:, None] * npilars
    n0, n1 = p0 + layer, p1 + layer
    n2, n3 = n0 + npilars, n1 + npilars
    conn = np.stack((n0, n1, n3, n2), axis = -1).astype(np.int64).ravel()
    ncells = npilars * nlayers
    offsets = np.arange(1, ncells + 1, dtype = np.int64) * 4
    ctype = np.full(ncells, VtkQuad.tid, dtype = np.uint8)

    # the older positional form passed dictionaries in these two places
    if isinstance(all_cell_data, dict) and cellData is None:
        all_cell_data, cellData = None, all_cell_data
    if isinstance(all_point_data, dict) and pointData is None:
        all_point_data, pointData = None, all_point_data

    def _dict_to_fields(d, where):
        if d is None:
            return []
        if not isinstance(d, dict):
            raise TypeError("%s must be a dictionary {name: array}" % where)
        return [[name, "scalars", arr] for name, arr in d.items()]
    all_cell_data = list(all_cell_data or []) + _dict_to_fields(cellData, "cellData")
    all_point_data = list(all_point_data or []) + _dict_to_fields(pointData, "pointData")

    return unstructuredGridToVTK(path, xx, yy, zz, connectivity = conn, offsets = offsets,
                                 cell_types = ctype, all_cell_data = all_cell_data or None,
                                 all_point_data = all_point_data or None, comments = comments)

# =================================
#  time-series on unstructuredGrid
# =================================
#unstructuredGridToVTK
class timeseries_unstructuredGrid:

    def __ts_convertListToArray(self, list1d):
        ''' If data is a list and no a Numpy array, then it convert it
            to an array, otherwise return the same array '''
        if (list1d is not None) and (not isinstance(list1d, np.ndarray)):
            if not isinstance(list1d, (list, tuple)):
                raise TypeError("expected a list, tuple or numpy array, got %s"
                                % type(list1d).__name__)
            return np.array(list1d)
        else:
            return list1d

    def __ts_convertDictListToArrays(self, data):
        ''' If data in dictironary are lists and no a Numpy array,
            then it creates a new dictionary and convert the list to arrays,
            otherwise return the same dictionary '''
        if data is not None:
            dict = {}
            for k, list1d in data.items():
                dict[k] = self.__ts_convertListToArray(list1d)
            return dict
        else:
            return data # None


    def __init__(self, filepath, time_values):
        """
            PARAMETERS:
                filepath: filename without extension (will be a .vtu file).
                time_values: numpy array of time values.
        """
        time_values = _to_vector("time_values", time_values)
        if time_values.size == 0:
            raise ValueError("time_values must contain at least one time value")
        self.ftype = VtkUnstructuredGrid
        self.filename = filepath
        self.VtkFile_obj = []
        self.time_values = time_values
        self.data_order = []
        # one entry per variable, in the order append_data must be called:
        # (association, name, numpy dtype, number of values per time step)
        self._expected = []
        self._next = 0

    def init_unstructuredGridToVTK(self, x, y, z, connectivity, offsets, cell_types, all_cell_data = None, all_point_data = None, comments = None):
        """
            INITIAL Export of unstructured grid and associated data (header info only).

            PARAMETERS:
                x, y, z: 1D list-type object (list, tuple or numpy) with coordinates of the vertices of cells. It is assumed that each element
                         has diffent number of vertices.
                connectivity: 1D list-type object (list, tuple or numpy) that defines the vertices associated to each element. 
                              Together with offset define the connectivity or topology of the grid. 
                              It is assumed that vertices in an element are listed consecutively.
                offsets: 1D list-type object (list, tuple or numpy) with the index of the last vertex of each element in the connectivity array.
                         It should have length nelem, where nelem is the number of cells or elements in the grid.
                cell_types: 1D list-type object (list, tuple or numpy) with an integer that defines the cell type of each element in the grid.
                            It should have size nelem. This should be assigned from pyvtk.vtkbin.VtkXXXX.tid, where XXXX represents
                            the type of cell. Please check the VTK file format specification for allowed cell types.                       
                all_cell_data: A List of Lists.  It has this format:
                        all_cell_data[ii] = cellData, where
                        cellData is a List (of length 3) that defines a variable associated with the grid cells:
                            cellData[0] = the *name* of the variable stored;
                            cellData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                            cellData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                        Note: all data arrays inside "all_cell_data" must have the same number of elements (number of cells).
                        Note: the length of "all_cell_data" is the number of variables (defined on cells).
                all_point_data: A List of Lists.  It has this format:
                        all_point_data[ii] = pointData, where
                        pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                            pointData[0] = the *name* of the variable stored;
                            pointData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                            pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                        Note: all data arrays inside "all_point_data" must have the same number of elements (number of vertices).
                        Note: the length of "all_point_data" is the number of variables (defined on vertices).
                comments: list of comment strings, which will be added to the header section of the file.
                
            RETURNS:
                XXX

        """
        if self.VtkFile_obj:
            raise RuntimeError("init_unstructuredGridToVTK was already called")
        x, y, z = _check_xyz(x, y, z)
        npoints = x.size
        connectivity, offsets, cell_types = _check_topology(connectivity, offsets, cell_types, npoints)
        ncells = cell_types.size

        cell_fields = _check_fields(all_cell_data, "all_cell_data", ncells)
        point_fields = _check_fields(all_point_data, "all_point_data", npoints)
        comments = _check_comments(comments)

        self._expected = ([("cell", f[0], f[2].dtype, f[2].size) for f in cell_fields] +
                          [("point", f[0], f[2].dtype, f[2].size) for f in point_fields])
        self._next = 0

        self.VtkFile_obj = VtkFile(self.filename, VtkUnstructuredGrid)
        try:
            if comments: self.VtkFile_obj.addComments(comments)

            num_time_indices = len(self.time_values)
            time_indices = np.arange(num_time_indices)
            time_values_str = ' '.join(map(str, time_indices))
            self.VtkFile_obj.openGrid(time_values = time_values_str)

            self.VtkFile_obj.openPiece(ncells = ncells, npoints = npoints)
            self.VtkFile_obj.openElement("Points")
            self.VtkFile_obj.addData("points", (x,y,z))
            self.VtkFile_obj.closeElement("Points")
            self.VtkFile_obj.openElement("Cells")
            self.VtkFile_obj.addData("connectivity", connectivity)
            self.VtkFile_obj.addData("offsets", offsets)
            self.VtkFile_obj.addData("types", cell_types)
            self.VtkFile_obj.closeElement("Cells")

            #self.data_order
            # create time-step indices
            time_steps = np.arange(len(self.time_values))
            _addDataToFile(self.VtkFile_obj, all_cell_data = cell_fields, all_point_data = point_fields, time_steps = time_steps)

            self.VtkFile_obj.closePiece()
            self.VtkFile_obj.closeGrid()

            self.VtkFile_obj.appendData( (x,y,z) )
            self.VtkFile_obj.appendData(connectivity).appendData(offsets).appendData(cell_types)
        except BaseException:
            self.VtkFile_obj.abort()
            raise

    def append_data(self, varData):
        """
        Append data array (for a variable) at a specific time value.

        PARAMETERS:
            varData is a List of Lists where
                    varData[ii] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).

        Variables must be appended in the order given to init_unstructuredGridToVTK:
        first all cell data, then all point data.  Each call must provide one array per
        time value, each with the same number of values as the array given at
        initialization.  Arrays of a different dtype are converted when NumPy's "safe"
        casting rules allow it (e.g. float32 -> float64, int -> float64); otherwise a
        TypeError is raised.  Nothing is written
        if a check fails.
        """
        if not self.VtkFile_obj:
            raise RuntimeError("call init_unstructuredGridToVTK before append_data")
        if self.VtkFile_obj.xml.closed:
            raise RuntimeError("the file has already been closed")
        if varData is None:
            raise ValueError("varData must be a list with one array per time value")
        if self._next >= len(self._expected):
            raise ValueError("append_data was called more times than there are variables (%d)"
                             % len(self._expected))
        assoc, name, dtype, size = self._expected[self._next]
        where = "append_data for %s variable '%s'" % (assoc, name)
        ntimes = len(self.time_values)
        if len(varData) != ntimes:
            raise ValueError("%s: got %d arrays, expected one per time value (%d)"
                             % (where, len(varData), ntimes))

        # check everything before writing anything
        arrays = []
        for ii in range(ntimes):
            arr = _to_array("%s, time index %d" % (where, ii), varData[ii])
            if arr.ndim not in (1, 3):
                raise ValueError("%s, time index %d: must be a 1-D array, got shape %s"
                                 % (where, ii, arr.shape))
            if arr.size != size:
                raise ValueError("%s, time index %d: %d values, expected %d"
                                 % (where, ii, arr.size, size))
            if arr.dtype != dtype:
                if not np.can_cast(arr.dtype, dtype, casting = "safe"):
                    raise TypeError("%s, time index %d: dtype %s cannot be safely converted to %s "
                                    "(the dtype given at initialization)"
                                    % (where, ii, arr.dtype, dtype))
                arr = arr.astype(dtype)
            arrays.append(arr)

        # append data to binary section
        for arr in arrays:
            self.VtkFile_obj.appendData(arr)
        self._next += 1

    def close_unstructuredGridToVTK(self):
        """
        Close the file.

        RETURNS:
            Full path to saved file.

        Raises ValueError (and discards the file) if append_data was not called for
        every variable.
        """
        if not self.VtkFile_obj:
            raise RuntimeError("call init_unstructuredGridToVTK first")
        if self._next != len(self._expected):
            missing = ", ".join("%s '%s'" % (a, n) for a, n, _, _ in self._expected[self._next:])
            self.VtkFile_obj.abort()
            raise ValueError("file discarded: no data was appended for %s" % missing)

        self.VtkFile_obj.save()
        return self.VtkFile_obj.getFileName()

    def write_fake_grid_for_time(self, filepath):
        """
        This writes a "self-contained" .vtu file that contains an unstructured grid consisting of one line element,
        and cellData consisting of a "TimeValue" variable that is indexed with the time-step.
        The purpose of this is to have a way of storing actual time values that can be non-uniformly spaced, without
        using a .pvd file.

        You can then use this in Paraview by including the .vtu file in your pipeline.  Then just create a
        Python Annotation Filter, and put this in the Expression field:  "Time: %1.2f" %TimeValue[0]
        Note: you will need to hide the "grid" of this .vtu, which is no big deal.

        RETURNS:
            Full path to saved file.
        """

        return _write_time_values_grid(filepath, self.time_values, self.filename + self.ftype.ext)

# =================================
#  Streaming time series: one .vtu file, or a .pvd file with one .vtu per step
# =================================
_INFERRED_TYPE = {1 : "scalars", 3 : "vectors", 9 : "tensors"}

def _step_field(where, name, value, ntuples):
    """ Check one field given as an array, or as (type, array).  A 2-D array with one
        row per point/cell, e.g. shape (n, 3), is accepted for multi-component data.
        Returns [name, type, flat array, ncomp] (the format used by _check_fields). """
    if not isinstance(name, str) or not name:
        raise TypeError("%s: field names must be non-empty strings, got %r" % (where, name))
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[0], str):
        dtype_str, data = value
    else:
        dtype_str, data = None, value
    if isinstance(data, (list, tuple)) and data and isinstance(data[0], (np.ndarray, list, tuple)):
        raise TypeError("%s['%s']: give multi-component data as one array with one row per "
                        "point/cell, e.g. np.column_stack((vx, vy, vz)), not as a %s of arrays"
                        % (where, name, type(data).__name__))
    arr = _to_array("%s['%s']" % (where, name), data)
    if arr.ndim == 2 and arr.shape[0] == ntuples:
        arr = np.ascontiguousarray(arr).reshape(-1)     # rows are tuples: x0 y0 z0 x1 ...
    elif arr.ndim != 1:
        raise ValueError("%s['%s'] must be 1-D, or 2-D with one row per %s; got shape %s"
                         % (where, name, "cell" if "cell" in where else "point", arr.shape))
    if dtype_str is None:
        ncomp = arr.size // ntuples if ntuples else 1
        if ntuples and arr.size % ntuples:
            raise ValueError("%s['%s'] has %d values, which is not a multiple of %d %s"
                             % (where, name, arr.size, ntuples,
                                "cells" if "cell" in where else "points"))
        dtype_str = _INFERRED_TYPE.get(ncomp)
        if dtype_str is None:
            raise ValueError("%s['%s'] has %d components per tuple; give the type explicitly, "
                             "e.g. ('tcoords', array), or pad 2-D vectors to 3 components"
                             % (where, name, ncomp))
    return _check_fields([[name, dtype_str, arr]], where, ntuples)[0]

def _block_size(arr):
    """ Bytes taken by an array in the appended section (8-byte size header + data). """
    return 8 + arr.size * arr.dtype.itemsize


class TimeSeriesUnstructuredGrid:
    """
    Write a time-dependent simulation on a fixed unstructured grid, one time step at a time.

    Two kinds of output:

    * output = "vtu" (default): everything in ONE .vtu file.  The mesh is stored once, and
      every time step is stored as a separate DataArray tagged with its TimeStep.  Each
      step is written to disk as soon as write_step() is called, so the whole history never
      has to be kept in memory.  If the time values are known in advance (argument
      "times"), the data goes straight into the file; otherwise it is collected in a
      temporary file next to the target and assembled when close() is called.
      ParaView/VTK use the step index (0, 1, 2, ...) as "time" for such files; the actual
      time values are written to TimeValues and, if time_values_file is set, to a small
      companion file (see timeseries_unstructuredGrid.write_fake_grid_for_time).
      The .vtu only exists once close() has succeeded.
      NOTE: VTK 9.6 and later (e.g. ParaView 6.1) have a bug that reads only the first
      time step of such files; use output = "vtkhdf" or "pvd" for those versions.

    * output = "vtkhdf": everything in ONE VTKHDF (HDF5) file, <path>.vtkhdf.  The mesh is
      stored once and each field is one HDF5 dataset holding all time steps; the actual
      time values are stored in the file, so ParaView shows them directly.  Each step is
      written to disk as soon as write_step() is called.  Needs the h5py package.  Read
      by ParaView 5.12 and later, including 6.1 (checked with VTK 9.4, 9.5, 9.6 and 9.7).
      Array names may not contain "/" or ".".  The file only exists once close() has
      succeeded.

    * output = "pvd": a ParaView collection.  Each step is a complete .vtu file in the
      folder <path>/, named <name>_00000.vtu, <name>_00001.vtu, ..., and <path>.pvd lists
      them with their actual time values.  The .pvd is updated after every step, so the
      steps written so far can be opened even if the simulation stops early.  Works in
      all ParaView versions; the cost is that the mesh is stored in every step file.

    Fields are given to write_step() as dictionaries {name: array}.  The type is inferred
    from the number of components (1: scalars, 3: vectors, 9: tensors); give
    {name: (type, array)} to choose another type (e.g. "normals", or "tcoords" with 1-3
    components).  Multi-component data may be interleaved 1-D arrays (x0 y0 z0 x1 ...)
    or 2-D arrays with one row per point/cell, e.g. shape (npoints, 3).  The first
    step fixes the set of fields, their types and dtypes; later steps must give the same
    fields (arrays are converted to the first step's dtype when NumPy's "safe" casting
    rules allow it, e.g. float32 -> float64).
    Fields that do not change in time can be given once, as static_point_data /
    static_cell_data.

    Example:

        with TimeSeriesUnstructuredGrid("out/sol", x, y, z, conn, offsets, types,
                                        static_cell_data={"material": mat}) as ts:
            for t in times:
                ...  # compute u (one value per point) and q (shape (npoints, 3))
                ts.write_step(t, point_data={"u": u, "q": q})
        # -> out/sol.vtu (and out/sol_time_values.vtu)

    PARAMETERS:
        path: file name without extension.
        x, y, z, connectivity, offsets, cell_types: the mesh, as in unstructuredGridToVTK.
        times: the time values (strictly increasing), if known in advance.  Then exactly
               this many steps must be written, and write_step(time) may omit the time.
        output: "vtu" (one .vtu file), "vtkhdf" (one .vtkhdf file) or "pvd" (one .vtu file
               per step plus a .pvd).
        static_point_data, static_cell_data: dictionaries {name: array} of data that do
               not depend on time; written once.
        comments: list of comment strings for the file header(s).
        time_values_file: for output = "vtu": True to also write
               <path>_time_values.vtu (the "fake grid" holding the time values), False for
               no such file, or a file name without extension.
    """

    def __init__(self, path, x, y, z, connectivity, offsets, cell_types, times = None,
                 output = "vtu", static_point_data = None, static_cell_data = None,
                 comments = None, time_values_file = True):
        if output not in ("vtu", "vtkhdf", "pvd"):
            raise ValueError("output must be 'vtu', 'vtkhdf' or 'pvd', got %r" % (output,))
        self.path = str(path)
        self.output = output
        if output == "vtkhdf":
            try:
                import h5py
            except ImportError:
                raise ImportError("output = 'vtkhdf' needs the h5py package; "
                                  "install it with 'pip install h5py'") from None
            self._h5py = h5py
        x, y, z = _check_xyz(x, y, z)
        connectivity, offsets, cell_types = _check_topology(connectivity, offsets, cell_types,
                                                            x.size)
        self._mesh = (x, y, z, connectivity, offsets, cell_types)
        self._npoints, self._ncells = x.size, offsets.size
        self._static = {"point" : self._fields(static_point_data, "static_point_data", "point"),
                        "cell" : self._fields(static_cell_data, "static_cell_data", "cell")}
        self._comments = _check_comments(comments)
        self._time_values_file = time_values_file
        if output == "vtkhdf":
            for assoc in ("point", "cell"):
                self._check_hdf_names(self._static[assoc], "static_%s_data" % assoc)

        if times is not None:
            times = _to_vector("times", times).astype(float)
            if times.size == 0:
                raise ValueError("times must contain at least one time value")
            if not np.all(np.isfinite(times)):
                raise ValueError("times must be finite numbers")
            if np.any(np.diff(times) <= 0):
                raise ValueError("times must be strictly increasing")
        if output == "vtu" and isinstance(time_values_file, str):
            tv_dir = os.path.dirname(os.path.abspath(time_values_file))
            if not os.path.isdir(tv_dir):
                raise ValueError("folder for time_values_file does not exist: %s" % tv_dir)
        self._planned_times = times
        self._times = []            # time values of the steps written so far
        self._pvd_files = []        # output = "pvd": (step file, time) of the steps written
        self._result = None         # path returned by close()
        self._layout = None         # {"point": [(name, type, dtype, size, ncomp)], "cell": [...]}
        self._vtk = None            # output = "vtu": the VtkFile being written
        self._scratch = None        # output = "vtu", times unknown: step data until close()
        self._scratch_path = None
        self._stream = None         # where step data is appended
        self._h5 = None             # output = "vtkhdf": the open h5py.File
        self._h5_tmppath = None
        self._closed = False

    # ---------- context manager ----------
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is not None:
            self.abort()
        elif not self._closed:
            self.close()
        return False

    def __del__(self):
        if not getattr(self, "_closed", True):
            try:
                self.abort()
            except Exception:
                pass

    # ---------- public methods ----------
    @property
    def nsteps(self):
        """ Number of time steps written so far. """
        return len(self._times)

    def write_step(self, time = None, point_data = None, cell_data = None):
        """
        Write one time step.

        PARAMETERS:
            time: the time value of this step.  May be omitted if "times" was given to
                  the constructor (if given, it must match).
            point_data: dictionary {name: array} of data at the points.
            cell_data: dictionary {name: array} of data on the cells.

        RETURNS:
            The index of this step (0, 1, 2, ...).

        Nothing is written if a check fails.
        """
        if self._closed:
            raise RuntimeError("this time series has been closed")
        k = len(self._times)
        if self._planned_times is not None:
            if k >= self._planned_times.size:
                raise ValueError("all %d time steps given in 'times' have already been written"
                                 % self._planned_times.size)
            t = float(self._planned_times[k])
            if time is not None and not np.isclose(float(time), t, rtol = 1e-9, atol = 1e-12):
                raise ValueError("time %r does not match times[%d] = %r" % (time, k, t))
        else:
            if time is None:
                raise ValueError("write_step needs the time value (no 'times' were given)")
            t = float(time)
            if not np.isfinite(t):
                raise ValueError("time must be a finite number, got %r" % (time,))
            if self._times and t <= self._times[-1]:
                raise ValueError("time values must increase: %r after %r" % (t, self._times[-1]))

        fields = {"point" : self._fields(point_data, "point_data", "point"),
                  "cell" : self._fields(cell_data, "cell_data", "cell")}
        if self._layout is None:
            self._check_names(fields)
            layout = {assoc : [(f[0], f[1], f[2].dtype, f[2].size, f[3]) for f in fields[assoc]]
                      for assoc in ("point", "cell")}
        else:
            layout = self._layout
            fields = {assoc : self._match_layout(assoc, fields[assoc]) for assoc in fields}

        if self.output == "vtu":
            if self._layout is None:
                self._layout = layout
                try:
                    self._start_vtu()
                except BaseException:
                    # e.g. the folder does not exist: nothing was written, a retry may work
                    self._layout = None
                    if self._vtk is not None:
                        self._vtk.abort()
                        self._vtk = None
                    self._remove_scratch()
                    self._stream = None
                    raise
            try:
                for assoc in ("point", "cell"):
                    for f in fields[assoc]:
                        writeBlockSize(self._stream, f[2].size * f[2].dtype.itemsize)
                        writeArrayToFile(self._stream, f[2])
            except BaseException:
                # a partly written step cannot be undone: discard the output
                self.abort()
                raise
        elif self.output == "vtkhdf":
            if self._layout is None:
                for assoc in ("point", "cell"):
                    self._check_hdf_names(fields[assoc], "%s_data" % assoc)
                self._layout = layout
                try:
                    self._start_hdf()
                except BaseException:
                    self._layout = None
                    self._remove_hdf()
                    raise
            try:
                self._write_hdf_step(k, t, fields)
            except BaseException:
                self.abort()
                raise
        else:
            self._layout = layout
            self._write_pvd_step(k, t, fields)
        self._times.append(t)
        return k

    def close(self):
        """
        Finish the output.

        RETURNS:
            Full path of the .vtu, .vtkhdf or .pvd file.
        """
        if self._closed:
            if self._result is not None:
                return self._result
            raise RuntimeError("this time series was aborted; nothing more can be saved")
        nsteps = len(self._times)
        if nsteps == 0:
            self.abort()
            raise ValueError("no time steps were written; nothing was saved")
        if self._planned_times is not None and nsteps != self._planned_times.size:
            missing = self._planned_times.size - nsteps
            if self.output != "pvd":
                self.abort()
                raise ValueError("only %d of the %d time steps in 'times' were written; "
                                 "the file was discarded" % (nsteps, self._planned_times.size))
            self._closed = True
            raise ValueError("only %d of the %d time steps in 'times' were written (%d missing); "
                             "the .pvd lists the steps that were written"
                             % (nsteps, self._planned_times.size, missing))
        if self.output == "pvd":
            self._closed = True
            self._result = os.path.abspath(self.path + ".pvd")
            return self._result
        if self.output == "vtkhdf":
            try:
                self._h5.close()
                self._h5 = None
                target = os.path.realpath(os.path.abspath(self.path + ".vtkhdf"))
                os.replace(self._h5_tmppath, target)
                self._h5_tmppath = None
            except BaseException:
                self.abort()
                raise
            self._closed = True
            self._result = os.path.abspath(self.path + ".vtkhdf")
            return self._result

        tv_writer = None
        try:
            if self._scratch is not None:
                # times were not known in advance: header now, then the collected step data
                self._open_vtu(self._times)
                self._scratch.flush()
                self._scratch.seek(0)
                shutil.copyfileobj(self._scratch, self._vtk.xml.stream, 16 * 2**20)
                self._remove_scratch()
            filename = self._vtk.getFileName()
            if self._time_values_file:
                # written to a temporary file first, so that a failure here leaves no .vtu
                tv_path = (self.path + "_time_values" if self._time_values_file is True
                           else str(self._time_values_file))
                tv_writer = XmlWriter(tv_path + ".vtu", addDeclaration = False)
                tv_writer.addText(_time_values_grid_text(self._times, os.path.basename(filename)))
            self._vtk.save()
            if tv_writer is not None:
                tv_writer.close()
        except BaseException:
            if tv_writer is not None:
                tv_writer.abort()
            self.abort()
            raise
        self._closed = True
        self._result = filename
        return filename

    def abort(self):
        """ Stop without saving.  For output = "vtu" or "vtkhdf" nothing is written; for
            output = "pvd" the steps written so far (and the .pvd listing them) are kept. """
        self._closed = True
        if self._vtk is not None:
            self._vtk.abort()
        self._remove_scratch()
        self._remove_hdf()

    # ---------- helpers ----------
    def _fields(self, data, where, assoc):
        if data is None:
            return []
        if not isinstance(data, dict):
            raise TypeError("%s must be a dictionary {name: array}" % where)
        n = self._npoints if assoc == "point" else self._ncells
        return [_step_field(where, name, value, n) for name, value in data.items()]

    def _check_names(self, fields):
        for assoc in ("point", "cell"):
            static = {f[0] for f in self._static[assoc]}
            clash = static & {f[0] for f in fields[assoc]}
            if clash:
                raise ValueError("%s data name(s) used both as static and time-dependent data: %s"
                                 % (assoc, ", ".join(sorted(clash))))

    def _match_layout(self, assoc, fields):
        """ Check one step's fields against the first step; return them in layout order. """
        given = {f[0] : f for f in fields}
        expected = [entry[0] for entry in self._layout[assoc]]
        missing = [n for n in expected if n not in given]
        extra = [n for n in given if n not in expected]
        if missing or extra:
            raise ValueError("step %d: %s data must have the same fields as the first step%s%s"
                             % (len(self._times), assoc,
                                "; missing: " + ", ".join(missing) if missing else "",
                                "; unexpected: " + ", ".join(extra) if extra else ""))
        ordered = []
        for name, dtype_str, dtype, size, ncomp in self._layout[assoc]:
            f = given[name]
            if f[1] != dtype_str or f[2].size != size:
                raise ValueError("step %d: %s data '%s' is %s with %d values; the first step had "
                                 "%s with %d values" % (len(self._times), assoc, name, f[1],
                                                        f[2].size, dtype_str, size))
            arr = f[2]
            if arr.dtype != dtype:
                if not np.can_cast(arr.dtype, dtype, casting = "safe"):
                    raise TypeError("step %d: %s data '%s' has dtype %s, which cannot be safely "
                                    "converted to %s (the dtype of the first step)"
                                    % (len(self._times), assoc, name, arr.dtype, dtype))
                arr = arr.astype(dtype)
            ordered.append([name, dtype_str, arr, ncomp])
        return ordered

    def _start_vtu(self):
        if self._planned_times is not None:
            self._open_vtu(self._planned_times)
            self._stream = self._vtk.xml.stream
        else:
            directory, base = os.path.split(os.path.realpath(os.path.abspath(self.path + ".vtu")))
            self._scratch_path = os.path.join(directory, ".%s.%s.steps" % (base[:40],
                                                                          uuid.uuid4().hex[:8]))
            self._scratch = open(self._scratch_path, "x+b")
            self._stream = self._scratch

    def _remove_scratch(self):
        if self._scratch is not None:
            try:
                self._scratch.close()
            finally:
                self._scratch = None
                try:
                    os.remove(self._scratch_path)
                except OSError:
                    pass

    def _open_vtu(self, times):
        """ Write the header for len(times) steps, then the mesh and the static data.
            Step data is laid out time step by time step after that:
                offset(field f, step k) = G + k * S + (sizes of the fields before f),
            where G is the size of the mesh and static blocks and S that of one step. """
        x, y, z, connectivity, offsets, cell_types = self._mesh
        nsteps = len(times)
        w = VtkFile(self.path, VtkUnstructuredGrid)
        self._vtk = w
        if self._comments: w.addComments(self._comments)
        w.openGrid(time_values = " ".join(_number_to_string(float(t)) for t in times))
        w.openPiece(ncells = self._ncells, npoints = self._npoints)
        w.openElement("Points")
        w.addData("points", (x, y, z))
        w.closeElement("Points")
        w.openElement("Cells")
        w.addData("connectivity", connectivity)
        w.addData("offsets", offsets)
        w.addData("types", cell_types)
        w.closeElement("Cells")

        # byte offsets
        pos = w.offset
        static_offsets = {}
        for assoc in ("point", "cell"):
            static_offsets[assoc] = []
            for f in self._static[assoc]:
                static_offsets[assoc].append(pos)
                pos += _block_size(f[2])
        G = pos
        step_offsets = {}
        pos = 0
        for assoc in ("point", "cell"):
            step_offsets[assoc] = []
            for name, dtype_str, dtype, size, ncomp in self._layout[assoc]:
                step_offsets[assoc].append(pos)
                pos += 8 + size * dtype.itemsize
        S = pos

        for assoc, tag in (("point", "Point"), ("cell", "Cell")):
            layout, static = self._layout[assoc], self._static[assoc]
            if not layout and not static:
                continue
            defaults = {}
            for dtype_str, name in ([(e[1], e[0]) for e in layout] + [(f[1], f[0]) for f in static]):
                defaults.setdefault(dtype_str, name)
            w.openData(tag, **defaults)
            # Time-dependent arrays first, step by step, then the static arrays: with this
            # order the arrays are also matched correctly by VTK readers that pair data
            # elements with output arrays by position (VTK 9.5 and earlier).
            for k in range(nsteps):
                for (name, dtype_str, dtype, size, ncomp), off in zip(layout, step_offsets[assoc]):
                    w.addHeader(name, dtype.name, size // ncomp, ncomp, time_value = str(k),
                                offset = G + k * S + off)
            for f, off in zip(static, static_offsets[assoc]):
                w.addHeader(f[0], f[2].dtype.name, f[2].size // f[3], f[3], offset = off)
            w.closeData(tag)
        w.closePiece()
        w.closeGrid()

        w.appendData((x, y, z))
        w.appendData(connectivity).appendData(offsets).appendData(cell_types)
        for assoc in ("point", "cell"):
            for f in self._static[assoc]:
                w.appendData(f[2])

    def _write_pvd_step(self, k, t, fields):
        directory = self.path
        base = os.path.basename(os.path.normpath(self.path))
        os.makedirs(directory, exist_ok = True)
        x, y, z, connectivity, offsets, cell_types = self._mesh
        step_file = _write_unstructured(os.path.join(directory, "%s_%05d" % (base, k)),
                                        x, y, z, connectivity, offsets, cell_types,
                                        fields["cell"] + self._static["cell"],
                                        fields["point"] + self._static["point"],
                                        self._comments)
        entries = self._pvd_files + [(step_file, t)]
        # the .pvd is rewritten completely (and atomically); if that fails, the step is not
        # counted, and writing it again later reuses the same step file name
        with VtkGroup(self.path) as group:
            for filename, time in entries:
                group.addFile(filepath = filename, sim_time = _number_to_string(time))
        self._pvd_files = entries

    # ---------- output = "vtkhdf" ----------
    # Layout (VTKHDF 2.0, temporal UnstructuredGrid with a static mesh, one part):
    #   /VTKHDF                  attributes Version = [2, 0], Type = "UnstructuredGrid"
    #     Points (n, 3), Connectivity, Offsets (ncells + 1, starting with 0), Types
    #     NumberOfPoints, NumberOfCells, NumberOfConnectivityIds: one entry (the one part)
    #     PointData/<name>, CellData/<name>: all steps stacked, (nsteps * n[, ncomp]);
    #                          static arrays hold a single copy, (n[, ncomp])
    #     Steps                attribute NSteps; per step: Values (the time), PartOffsets = 0,
    #                          NumberOfParts = 1, PointOffsets = 0, CellOffsets = 0,
    #                          ConnectivityIdOffsets = 0 (the mesh is reused), and
    #                          PointDataOffsets/<name>, CellDataOffsets/<name> = k * n
    #                          (0 for static arrays)
    # Arrays carry the "Attribute" string (Scalars, Vectors, ...) for the first array of
    # each type; readers older than VTKHDF 2.6 ignore it.
    _HDF_ATTRIBUTE = {"scalars" : "Scalars", "vectors" : "Vectors", "normals" : "Normals",
                      "tensors" : "Tensors", "tcoords" : "TCoords"}

    @staticmethod
    def _check_hdf_names(fields, where):
        for f in fields:
            if "/" in f[0] or "." in f[0]:
                raise ValueError("%s: VTKHDF array names may not contain '/' or '.', got '%s'"
                                 % (where, f[0]))

    def _start_hdf(self):
        h5py = self._h5py
        target = os.path.realpath(os.path.abspath(self.path + ".vtkhdf"))
        directory, base = os.path.split(target)
        self._h5_tmppath = os.path.join(directory, ".%s.%s.part" % (base[:40], uuid.uuid4().hex[:8]))
        f = h5py.File(self._h5_tmppath, "x")
        self._h5 = f
        x, y, z, connectivity, offsets, cell_types = self._mesh

        def fixed_string(text):
            data = text.encode("ascii")
            return data, h5py.string_dtype("ascii", len(data))

        root = f.create_group("VTKHDF")
        root.attrs["Version"] = np.array([2, 0], dtype = np.int64)
        data, dtype = fixed_string("UnstructuredGrid")
        root.attrs.create("Type", data, dtype = dtype)
        points = np.empty((self._npoints, 3), dtype = x.dtype)
        points[:, 0], points[:, 1], points[:, 2] = (np.ravel(c, order = 'F') for c in (x, y, z))
        root.create_dataset("Points", data = points)
        root.create_dataset("Connectivity", data = connectivity.astype(np.int64))
        root.create_dataset("Offsets", data = np.concatenate(([0], offsets)).astype(np.int64))
        root.create_dataset("Types", data = cell_types.astype(np.uint8))
        for name, value in (("NumberOfPoints", self._npoints), ("NumberOfCells", self._ncells),
                            ("NumberOfConnectivityIds", connectivity.size)):
            root.create_dataset(name, data = np.array([value], dtype = np.int64))

        steps = root.create_group("Steps")
        steps.attrs["NSteps"] = np.int64(0)
        for name, dtype in (("Values", np.float64), ("PartOffsets", np.int64),
                            ("NumberOfParts", np.int64), ("PointOffsets", np.int64)):
            steps.create_dataset(name, shape = (0,), maxshape = (None,), dtype = dtype,
                                 chunks = (128,))
        for name in ("CellOffsets", "ConnectivityIdOffsets"):
            steps.create_dataset(name, shape = (0, 1), maxshape = (None, 1), dtype = np.int64,
                                 chunks = (128, 1))

        for assoc, group_name in (("point", "PointData"), ("cell", "CellData")):
            n = self._npoints if assoc == "point" else self._ncells
            group = root.create_group(group_name)
            offsets_group = steps.create_group(group_name + "Offsets")
            marked = set()
            entries = ([(e[0], e[1], e[2], e[4], False) for e in self._layout[assoc]] +
                       [(f[0], f[1], f[2].dtype, f[3], True) for f in self._static[assoc]])
            for name, dtype_str, dtype, ncomp, static in entries:
                shape = (n,) if ncomp == 1 else (n, ncomp)
                if static:
                    arr = [f[2] for f in self._static[assoc] if f[0] == name][0]
                    ds = group.create_dataset(name, data = arr.reshape(shape))
                else:
                    # chunks of about 1 MB, at most one step long
                    rows = max(1, min(n, 2**20 // max(1, dtype.itemsize * ncomp)))
                    ds = group.create_dataset(name, shape = (0,) + shape[1:],
                                              maxshape = (None,) + shape[1:], dtype = dtype,
                                              chunks = (rows,) + shape[1:])
                if dtype_str not in marked:
                    marked.add(dtype_str)
                    data, string_dtype = fixed_string(self._HDF_ATTRIBUTE[dtype_str])
                    ds.attrs.create("Attribute", data, dtype = string_dtype)
                offsets_group.create_dataset(name, shape = (0,), maxshape = (None,),
                                             dtype = np.int64, chunks = (128,))

    def _write_hdf_step(self, k, t, fields):
        root = self._h5["VTKHDF"]
        steps = root["Steps"]

        def append(ds, value):
            ds.resize(ds.shape[0] + 1, axis = 0)
            ds[-1] = value

        for assoc, group_name in (("point", "PointData"), ("cell", "CellData")):
            n = self._npoints if assoc == "point" else self._ncells
            group, offsets_group = root[group_name], steps[group_name + "Offsets"]
            for name, dtype_str, arr, ncomp in fields[assoc]:
                ds = group[name]
                ds.resize((k + 1) * n, axis = 0)
                ds[k * n:(k + 1) * n] = arr if ncomp == 1 else arr.reshape(n, ncomp)
                append(offsets_group[name], k * n)
            for f in self._static[assoc]:
                append(offsets_group[f[0]], 0)
        append(steps["Values"], t)
        append(steps["PartOffsets"], 0)
        append(steps["NumberOfParts"], 1)
        append(steps["PointOffsets"], 0)
        append(steps["CellOffsets"], [0])
        append(steps["ConnectivityIdOffsets"], [0])
        steps.attrs["NSteps"] = np.int64(k + 1)
        self._h5.flush()

    def _remove_hdf(self):
        if self._h5 is not None:
            try:
                self._h5.close()
            except Exception:
                pass
            self._h5 = None
        if self._h5_tmppath is not None:
            try:
                os.remove(self._h5_tmppath)
            except OSError:
                pass
            self._h5_tmppath = None
