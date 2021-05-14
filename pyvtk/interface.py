"""
pyvtk.interface.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
High level Python library to export data to binary VTK file.

Copyright (c) 05-07-2021,  Shawn W. Walker
"""

from .vtkbin import * # VtkFile, VtkUnstructuredGrid, etc.
try:
    import numpy as np
except:
    print("Numpy is not installed. Please install it before running VTKwrite again.")

# =================================
#       Helper functions
# =================================
def _addDataToFile(vtkFile, all_cell_data, all_point_data, time_steps = None):
    dt_strs = ["scalars", "vectors", "normals", "tensors", "tcoords"]
    ncomp_dict = {"scalars" : 1, "vectors" : 3, "normals" : 3, "tensors" : 9, "tcoords" : 0}
    # SWW: tcoords not supported properly...
    
    # initialize dicts for setting default scalars, vectors, etc...
    default_cell_vars  = {}
    default_point_vars = {}
    
    if all_cell_data is not None:
        len_all_cell_data = len(all_cell_data)
    else:
        len_all_cell_data = 0

    if all_point_data is not None:
        len_all_point_data = len(all_point_data)
    else:
        len_all_point_data = 0
    
    # get default vars for each data-type
    for dt in dt_strs:
        for ii in range(len_all_cell_data):
            if dt == all_cell_data[ii][1]:
                default_cell_vars[dt] = all_cell_data[ii][0]
                break

    # get default vars for each data-type
    for dt in dt_strs:
        for ii in range(len_all_point_data):
            if dt == all_point_data[ii][1]:
                default_point_vars[dt] = all_point_data[ii][0]
                break

    # Cell data
    if default_cell_vars:
        vtkFile.openData("Cell", **default_cell_vars)
        for ii in range(len_all_cell_data):
            if time_steps is not None:
                for ti in time_steps:
                    vtkFile.internal_addData(all_cell_data[ii][0], all_cell_data[ii][2], ncomp_dict[all_cell_data[ii][1]], str(ti))
            else:
                vtkFile.internal_addData(all_cell_data[ii][0], all_cell_data[ii][2], ncomp_dict[all_cell_data[ii][1]])
        vtkFile.closeData("Cell")

    # Point data
    if default_point_vars:
        vtkFile.openData("Point", **default_point_vars)
        for ii in range(len_all_point_data):
            if time_steps is not None:
                for ti in time_steps:
                    vtkFile.internal_addData(all_point_data[ii][0], all_point_data[ii][2], ncomp_dict[all_point_data[ii][1]], str(ti))
            else:
                vtkFile.internal_addData(all_point_data[ii][0], all_point_data[ii][2], ncomp_dict[all_point_data[ii][1]])
        vtkFile.closeData("Point")

def _appendDataToFile(vtkFile, all_cell_data, all_point_data):
    #dt_strs = ["scalars", "vectors", "normals", "tensors", "tcoords"]
    #ncomp_dict = {"scalars" : 1, "vectors" : 3, "normals" : 3, "tensors" : 9, "tcoords" : 0}
    # SWW: tcoords not supported properly...
    
    if all_cell_data is not None:
        len_all_cell_data = len(all_cell_data)
    else:
        len_all_cell_data = 0

    if all_point_data is not None:
        len_all_point_data = len(all_point_data)
    else:
        len_all_point_data = 0
    
    # append data to binary section
    
    # cell based data
    for ii in range(len_all_cell_data):
        data = all_cell_data[ii][2]
        vtkFile.appendData(data)

    # point based data
    for ii in range(len_all_point_data):
        data = all_point_data[ii][2]
        vtkFile.appendData(data)

def __convertListToArray(list1d):
    ''' If data is a list and no a Numpy array, then it convert it
        to an array, otherwise return the same array '''
    if (list1d is not None) and (not type(list1d).__name__ == "ndarray"):
        assert isinstance(list1d, (list, tuple))
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
    assert (all_cell_data != None or all_point_data != None)
    
    # Extract dimensions
    start = (0,0,0)
    end = None
    
    if all_cell_data != None:
        for ii in range(len(all_cell_data)):
            cellData = all_cell_data[ii]
            assert (cellData[1] == "scalars")
            data = cellData[2]
            end = data.shape
            break

    if all_point_data != None:
        for ii in range(len(all_point_data)):
            pointData = all_point_data[ii]
            assert (pointData[1] == "scalars")
            data = pointData[2]
            end = data.shape
            end = (end[0] - 1, end[1] - 1, end[2] - 1)
            break
    
    # Write data to file
    w = VtkFile(path, VtkImageData)
    if comments: w.addComments(comments)
    w.openGrid(start = start, end = end, origin = origin, spacing = spacing)
    w.openPiece(start = start, end = end)
    _addDataToFile(w, all_cell_data, all_point_data)
    w.closePiece()
    w.closeGrid()
    _appendDataToFile(w, all_cell_data, all_point_data)
    w.save()
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
    assert (x.ndim == 1 and y.ndim == 1 and z.ndim == 1), "Wrong array dimension"
    ftype = VtkRectilinearGrid
    nx, ny, nz = x.size - 1, y.size - 1, z.size - 1
    # Extract dimensions
    start = (0,0,0)
    end   = (nx, ny, nz)

    if all_cell_data != None:
        for ii in range(len(all_cell_data)):
            cellData = all_cell_data[ii]
            assert (cellData[1] == "scalars")

    if all_point_data != None:
        for ii in range(len(all_point_data)):
            pointData = all_point_data[ii]
            assert (pointData[1] == "scalars")
    
    w =  VtkFile(path, ftype)
    if comments: w.addComments(comments)
    w.openGrid(start = start, end = end)
    w.openPiece(start = start, end = end)

    w.openElement("Coordinates")
    w.addData("x_coordinates", x)
    w.addData("y_coordinates", y)
    w.addData("z_coordinates", z)
    w.closeElement("Coordinates")

    _addDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)
    w.closePiece()
    w.closeGrid()
    # Write coordinates
    w.appendData(x).appendData(y).appendData(z)
    # Write data
    _appendDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)
    w.save()
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
    assert (x.ndim == 3 and y.ndim == 3 and z.ndim == 3), "Wrong arrays dimensions"
    
    ftype = VtkStructuredGrid
    s = x.shape
    nx, ny, nz = s[0] - 1, s[1] - 1, s[2] - 1
    start = (0,0,0)
    end = (nx, ny, nz)

    if all_cell_data != None:
        for ii in range(len(all_cell_data)):
            cellData = all_cell_data[ii]
            assert (cellData[1] == "scalars")

    if all_point_data != None:
        for ii in range(len(all_point_data)):
            pointData = all_point_data[ii]
            assert (pointData[1] == "scalars")

    w =  VtkFile(path, ftype)
    if comments: w.addComments(comments)
    w.openGrid(start = start, end = end)
    w.openPiece(start = start, end = end)
    w.openElement("Points")
    w.addData("points", (x,y,z))
    w.closeElement("Points")

    _addDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)
    w.closePiece()
    w.closeGrid()
    w.appendData( (x,y,z) )
    _appendDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)
    w.save()
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
    assert ( len(x) == len(y) == len(z) )
    x = __convertListToArray(x)
    y = __convertListToArray(y)
    z = __convertListToArray(z)

    if all_point_data is not None:
        len_all_point_data = len(all_point_data)
    else:
        len_all_point_data = 0
    for ii in range(len_all_point_data):
        all_point_data[ii][2] = __convertListToArray(all_point_data[ii][2])
    
    npoints = len(x)
    
    # create some temporary arrays to write grid topology
    offsets = np.arange(start = 1, stop = npoints + 1, dtype = 'int32') # index of last node in each cell
    connectivity = np.arange(npoints, dtype = 'int32')                 # each point is only connected to itself
    cell_types = np.empty(npoints, dtype = 'uint8') 
   
    cell_types[:] = VtkVertex.tid

    w = VtkFile(path, VtkUnstructuredGrid)
    if comments: w.addComments(comments)
    w.openGrid()
    w.openPiece(ncells = npoints, npoints = npoints)
    
    w.openElement("Points")
    w.addData("points", (x,y,z))
    w.closeElement("Points")
    w.openElement("Cells")
    w.addData("connectivity", connectivity)
    w.addData("offsets", offsets)
    w.addData("types", cell_types)
    w.closeElement("Cells")
    
    _addDataToFile(w, all_cell_data = None, all_point_data = all_point_data)

    w.closePiece()
    w.closeGrid()
    w.appendData( (x,y,z) )
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)

    _appendDataToFile(w, all_cell_data = None, all_point_data = all_point_data)

    w.save()
    return w.getFileName()
    
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

    assert len(x) == len(y) and len(x) == len(z)
    assert (ndim == 2) or (ndim == 3)
    x = __convertListToArray(x)
    y = __convertListToArray(y)
    z = __convertListToArray(z)
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
    assert (x.size == y.size == z.size)
    assert (x.size % 2 == 0)
    
    x = __convertListToArray(x)
    y = __convertListToArray(y)
    z = __convertListToArray(z)

    if all_cell_data is not None:
        len_all_cell_data = len(all_cell_data)
    else:
        len_all_cell_data = 0
    for ii in range(len_all_cell_data):
        all_cell_data[ii][2] = __convertListToArray(all_cell_data[ii][2])

    if all_point_data is not None:
        len_all_point_data = len(all_point_data)
    else:
        len_all_point_data = 0
    for ii in range(len_all_point_data):
        all_point_data[ii][2] = __convertListToArray(all_point_data[ii][2])
    
    npoints = len(x)
    ncells = int(len(x) / 2.0)
    
    # Check all_cell_data has the same size that the number of cells
    
    # create some temporary arrays to write grid topology
    offsets = np.arange(start = 2, step = 2, stop = npoints + 1, dtype = 'int32') # index of last node in each cell
    connectivity = np.arange(npoints, dtype = 'int32')                          # each point is only connected to itself
    cell_types = np.empty(npoints, dtype = 'uint8') 
   
    cell_types[:] = VtkLine.tid

    w = VtkFile(path, VtkUnstructuredGrid)
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
    
    _addDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)

    w.closePiece()
    w.closeGrid()
    w.appendData( (x,y,z) )
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)

    _appendDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)

    w.save()
    return w.getFileName()

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
    assert (x.size == y.size == z.size)
    
    x = __convertListToArray(x)
    y = __convertListToArray(y)
    z = __convertListToArray(z)
    
    if all_cell_data is not None:
        len_all_cell_data = len(all_cell_data)
    else:
        len_all_cell_data = 0
    for ii in range(len_all_cell_data):
        all_cell_data[ii][2] = __convertListToArray(all_cell_data[ii][2])

    if all_point_data is not None:
        len_all_point_data = len(all_point_data)
    else:
        len_all_point_data = 0
    for ii in range(len_all_point_data):
        all_point_data[ii][2] = __convertListToArray(all_point_data[ii][2])

    npoints = len(x)
    ncells = pointsPerLine.size
    
    # create some temporary arrays to write grid topology
    offsets = np.zeros(ncells, dtype = 'int32')         # index of last node in each cell
    ii = 0
    for i in range(ncells):
        ii += pointsPerLine[i]
        offsets[i] = ii
    
    connectivity = np.arange(npoints, dtype = 'int32')      # each line connects points that are consecutive
   
    cell_types = np.empty(npoints, dtype = 'uint8') 
    cell_types[:] = VtkPolyLine.tid

    w = VtkFile(path, VtkUnstructuredGrid)
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
    
    _addDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)

    w.closePiece()
    w.closeGrid()
    w.appendData( (x,y,z) )
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)

    _appendDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)

    w.save()
    return w.getFileName()

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
    assert (x.size == y.size == z.size)
    x = __convertListToArray(x)
    y = __convertListToArray(y)
    z = __convertListToArray(z)
    connectivity = __convertListToArray(connectivity)
    offsets = __convertListToArray(offsets)
    cell_types = __convertListToArray(cell_types)
    
    if all_cell_data is not None:
        len_all_cell_data = len(all_cell_data)
    else:
        len_all_cell_data = 0
    for ii in range(len_all_cell_data):
        all_cell_data[ii][2] = __convertListToArray(all_cell_data[ii][2])

    if all_point_data is not None:
        len_all_point_data = len(all_point_data)
    else:
        len_all_point_data = 0
    for ii in range(len_all_point_data):
        all_point_data[ii][2] = __convertListToArray(all_point_data[ii][2])

    npoints = x.size
    ncells = cell_types.size
    assert (offsets.size == ncells)
    
    w = VtkFile(path, VtkUnstructuredGrid)
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
    
    _addDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)
    
    w.closePiece()
    w.closeGrid()

    w.appendData( (x,y,z) )
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)

    _appendDataToFile(w, all_cell_data = all_cell_data, all_point_data = all_point_data)

    w.save()
    return w.getFileName()

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
        if (list1d is not None) and (not type(list1d).__name__ == "ndarray"):
            assert isinstance(list1d, (list, tuple))
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
        self.ftype = VtkUnstructuredGrid
        self.filename = filepath
        self.VtkFile_obj = []
        self.time_values = time_values
        self.data_order = []

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
        assert (x.size == y.size == z.size)
        x = self.__ts_convertListToArray(x)
        y = self.__ts_convertListToArray(y)
        z = self.__ts_convertListToArray(z)
        connectivity = self.__ts_convertListToArray(connectivity)
        offsets = self.__ts_convertListToArray(offsets)
        cell_types = self.__ts_convertListToArray(cell_types)

        for ii in range(len(all_cell_data)):
            all_cell_data[ii][2] = self.__ts_convertListToArray(all_cell_data[ii][2])

        for ii in range(len(all_point_data)):
            all_point_data[ii][2] = self.__ts_convertListToArray(all_point_data[ii][2])

        npoints = x.size
        ncells = cell_types.size
        assert (offsets.size == ncells)

        self.VtkFile_obj = VtkFile(self.filename, VtkUnstructuredGrid)
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
        _addDataToFile(self.VtkFile_obj, all_cell_data = all_cell_data, all_point_data = all_point_data, time_steps = time_steps)

        self.VtkFile_obj.closePiece()
        self.VtkFile_obj.closeGrid()

        self.VtkFile_obj.appendData( (x,y,z) )
        self.VtkFile_obj.appendData(connectivity).appendData(offsets).appendData(cell_types)

    def append_data(self, varData):
        """
        Append data array (for a variable) at a specific time value.

        PARAMETERS:
            varData is a List of Lists where
                    varData[ii] = the data array itself, i.e. a 1D list-type object (list, tuple or numpy).
        """

        assert ( varData is not None )
        
        # append data to binary section

        # loop through the time sequence
        for ii in range(len(varData)):
            data = varData[ii]
            self.VtkFile_obj.appendData(data)

    def close_unstructuredGridToVTK(self):
        """
        Close the file.

        RETURNS:
            Full path to saved file.
        """

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
        time_file = open(full_file, 'w')
        time_file.write('<?xml version="1.0"?>\n')
        time_file.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        time_file.write('<!-- This is a fake grid, that is only used for plotting actual Time Values -->\n')
        Comment_str_2 = '<!-- See the associated file: ' + self.filename + self.ftype.ext + ' -->\n'
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
        time_file.write('<CellData scalars="TimeValue">\n')
        
        for ii in range(num_time_indices):
            TV_ii_str = '<DataArray Name="TimeValue" NumberOfComponents="1" type="Float64" format="ascii" TimeStep="' + \
                                   str(ii) + '">' + str(self.time_values[ii]) + '</DataArray>\n'
            time_file.write(TV_ii_str)
        
        time_file.write('</CellData>\n')
        time_file.write('</Piece>\n')
        time_file.write('</UnstructuredGrid>\n')
        time_file.write('</VTKFile>\n')
        time_file.close()
        
        return full_file
