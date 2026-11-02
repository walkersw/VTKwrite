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
import os
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
            data: (THIS IS NOT ACTUALLY USED!) A List of Lists.  It has this format:
                    data[ii] = pointData, where
                    pointData is a List (of length 3) that defines a variable associated with the grid vertices:
                        pointData[0] = the *name* of the variable stored;
                        pointData[1] = the *type* of the variable, i.e. "scalars", "vectors", "normals", "tensors", or "tcoords".
                        pointData[2] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
                    Note: all data arrays inside "data" must have the same number of elements (number of vertices).
                    Note: the length of "data" is the number of variables (defined on vertices).
            comments: list of comment strings, which will be added to the header section of the file.
            ndim: is the number of dimensions considered when calling Delaunay.
                  If ndim = 2, then only coordinates x and y are passed.
                  If ndim = 3, then x, y and z coordinates are passed.
            
        RETURNS:
            Full path to saved file.
        
        REQUIRES: Scipy > 1.2.
    """
    # TODO: Check if it makes and it would be possible to add cellData.
    try:
        from scipy.spatial import Delaunay
    except:
        print("Failed to import scipy.spatial. Please install it if it is not installed.")

    if ndim not in (2, 3):
        raise ValueError("ndim must be 2 or 3, got %r" % (ndim,))
    x, y, z = _check_xyz(x, y, z)
    #XXX = __convertDictListToArrays(data)

    npts = len(x)

    points = np.zeros( (npts, ndim) ) # needs to create the 2D or 3D temporary array to call Delaunay
    for i in range(npts):
        points[i,0] = x[i]
        points[i,1] = y[i]
        if ndim > 2: points[i,2] = z[i]

    tri = Delaunay(points)

    # list of triangles that form the tesselation
    ncells, npoints_per_cell = tri.simplices.shape[0], tri.simplices.shape[1]
    conn =  np.zeros(ncells * 3)
    for i in range(ncells):
        ii = i * 3
        conn[ii]     = tri.simplices[i,0]
        conn[ii + 1] = tri.simplices[i,1]
        conn[ii + 2] = tri.simplices[i,2]

    offset = np.zeros(ncells)
    for i in range(ncells): offset[i] = (i + 1) * 3

    # initialize the data structure
    all_point_data = [["Elevation", "scalars", z]]

    cell_type = np.ones(ncells) * VtkTriangle.tid
    unstructuredGridToVTK(path, x, y, z, connectivity = conn, offsets = offset, cell_types = cell_type, all_cell_data = None, all_point_data = all_point_data, comments = None)

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
    cell_types = np.empty(npoints, dtype = 'uint8')

    cell_types[:] = VtkLine.tid

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
    offsets = np.zeros(ncells, dtype = 'int32')         # index of last node in each cell
    ii = 0
    for i in range(ncells):
        ii += pointsPerLine[i]
        offsets[i] = ii

    connectivity = np.arange(npoints, dtype = 'int32')      # each line connects points that are consecutive

    cell_types = np.empty(npoints, dtype = 'uint8')
    cell_types[:] = VtkPolyLine.tid

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
def cylinderToVTK(path, x0, y0, z0, z1, radius, nlayers, npilars = 16, cellData=None, pointData=None, comments = None ):
    """
        Export cylinder as VTK unstructured grid.  SWW: needs to be updated!!!!
    
      PARAMETERS:
        path: path to file without extension.
        x0, yo: center of cylinder.
        z0, z1: lower and top elevation of the cylinder.
        radius: radius of cylinder.
        nlayers: Number of layers in z direction to divide the cylinder.
        npilars: Number of points around the diameter of the cylinder. 
                 Higher value gives higher resolution to represent the curved shape.
        cellData: dictionary with 1D arrays that store cell data. 
                  Arrays should have number of elements equal to ncells = npilars * nlayers.
        pointData: dictionary with 1D arrays that store point data.
                  Arrays should have number of elements equal to npoints = npilars * (nlayers + 1).
        comments: list of comment strings, which will be added to the header section of the file.
                    
      RETURNS: 
            Full path to saved file.
        
        NOTE: This function only export vertical shapes for now. However, it should be easy to 
              rotate the cylinder to represent other orientations.
    """
    import math as m
    
    # Define x, y coordinates from polar coordinates.
    dpi = 2.0 * m.pi / npilars
    angles = np.arange(0.0, 2.0 * m.pi, dpi)

    x = radius * np.cos(angles) + x0
    y = radius * np.sin(angles) + y0

    dz = (z1 - z0) / nlayers
    z = np.arange(z0, z1+dz, step = dz)

    npoints = npilars * (nlayers + 1)
    ncells  = npilars * nlayers

    xx = np.zeros(npoints)
    yy = np.zeros(npoints)
    zz = np.zeros(npoints)

    ii = 0
    for k in range(nlayers + 1):
        for p in range(npilars):
            xx[ii] = x[p]
            yy[ii] = y[p]
            zz[ii] = z[k]
            ii = ii + 1

    # Define connectivity
    conn = np.zeros(4 * ncells, dtype = np.int64)
    ii = 0
    for l in range(nlayers):
        for p in range(npilars):
            p0 = p
            if(p + 1 == npilars):
                p1 = 0 
            else: 
                p1 = p + 1 # circular loop
       
            n0 = p0 + l * npilars
            n1 = p1 + l * npilars 
            n2 = n0 + npilars
            n3 = n1 + npilars
        
            conn[ii + 0] = n0
            conn[ii + 1] = n1
            conn[ii + 2] = n3
            conn[ii + 3] = n2 
            ii = ii + 4
 
    # Define offsets 
    offsets = np.zeros(ncells, dtype = np.int64)
    for i in range(ncells):
        offsets[i] = (i + 1) * 4

    # Define cell types
    ctype = np.ones(ncells) + VtkPixel.tid
    
    return unstructuredGridToVTK(path, xx, yy, zz, connectivity = conn, offsets = offsets, cell_types = ctype, cellData = cellData, pointData = pointData, comments = comments)

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
        initialization.  Arrays of a different dtype are converted when this is lossless
        (e.g. float32 -> float64); otherwise a TypeError is raised.  Nothing is written
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

        full_file = filepath + '.vtu'
        time_file = _TextBuffer()
        time_file.write('<?xml version="1.0"?>\n')
        time_file.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        time_file.write('<!-- This is a fake grid, that is only used for plotting actual Time Values -->\n')
        Comment_str_2 = '<!-- ' + sanitize_comment('See the associated file: ' + self.filename + self.ftype.ext) + ' -->\n'
        time_file.write(Comment_str_2)
        time_file.write('<!-- This .vtu file can be included in the Paraview pipeline. -->\n')
        time_file.write('<!-- Then use a Python Annotation Filter with this Expression: "Time: %1.2f" %TimeValue[0] -->\n')
        time_file.write('<!-- Make sure the "Array Association" is set to "Cell Data" -->\n')
        time_file.write('<!-- This is useful when the time step spacing is not uniform. -->\n')
        time_file.write('<!-- Note: make sure to hide the fake grid. -->\n')

        num_time_indices = len(self.time_values)
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
                                   str(ii) + '">' + _number_to_string(float(self.time_values[ii])) + '</DataArray>\n'
            time_file.write(TV_ii_str)

        time_file.write('</CellData>\n')
        time_file.write('</Piece>\n')
        time_file.write('</UnstructuredGrid>\n')
        time_file.write('</VTKFile>\n')
        _write_text_atomically(full_file, time_file.getvalue())

        return full_file
