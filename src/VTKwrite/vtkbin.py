"""
VTKwrite.vtkbin.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Low level Python library to export data to binary VTK file.

Copyright (c) 05-07-2021,  Shawn W. Walker
"""

from .pyvtk import writeBlockSize, writeArrayToFile, writeArraysToFile
from .xmlwrite import XmlWriter
import numpy as np
import numbers
import sys
import os

# ================================
#            VTK Types
# ================================

#     FILE TYPES
class VtkFileType:

    def __init__(self, name, ext):
        self.name = name
        self.ext  = ext

    def __str__(self):
        return "Name: %s  Ext: %s \n" % (self.name, self.ext)

VtkImageData       = VtkFileType("ImageData", ".vti")
VtkPolyData        = VtkFileType("PolyData", ".vtp")
VtkRectilinearGrid  = VtkFileType("RectilinearGrid", ".vtr")
VtkStructuredGrid   = VtkFileType("StructuredGrid", ".vts")
VtkUnstructuredGrid = VtkFileType("UnstructuredGrid", ".vtu")

#    DATA TYPES
class VtkDataType:

    def __init__(self, size, name):
        self.size = size
        self.name = name

    def __str__(self):
        return "Type: %s  Size: %d \n" % (self.name, self.size)

VtkInt8    = VtkDataType(1, "Int8")
VtkUInt8   = VtkDataType(1, "UInt8")
VtkInt16   = VtkDataType(2, "Int16")
VtkUInt16  = VtkDataType(2, "UInt16")
VtkInt32   = VtkDataType(4, "Int32")
VtkUInt32  = VtkDataType(4, "UInt32")
VtkInt64   = VtkDataType(8, "Int64")
VtkUInt64  = VtkDataType(8, "UInt64")
VtkFloat32 = VtkDataType(4, "Float32")
VtkFloat64 = VtkDataType(8, "Float64")

# Map numpy to VTK data types
np_to_vtk = { 'int8'    : VtkInt8,
            'uint8'   : VtkUInt8,
            'int16'   : VtkInt16,
            'uint16'  : VtkUInt16,
            'int32'   : VtkInt32,
            'uint32'  : VtkUInt32,
            'int64'   : VtkInt64,
            'uint64'  : VtkUInt64,
            'float32' : VtkFloat32,
            'float64' : VtkFloat64 }

#    CELL TYPES
class VtkCellType:

    def __init__(self, tid, name):
        self.tid = tid
        self.name = name

    def __str__(self):
        return "VtkCellType( %s ) \n" % ( self.name )

VtkVertex = VtkCellType(1, "Vertex")
VtkPolyVertex = VtkCellType(2, "PolyVertex")
VtkLine = VtkCellType(3, "Line")
VtkPolyLine = VtkCellType(4, "PolyLine")
VtkTriangle = VtkCellType(5, "Triangle")
VtkTriangleStrip = VtkCellType(6, "TriangleStrip")
VtkPolygon = VtkCellType(7, "Polygon")
VtkPixel = VtkCellType(8, "Pixel")
VtkQuad = VtkCellType(9, "Quad")
VtkTetra = VtkCellType(10, "Tetra")
VtkVoxel = VtkCellType(11, "Voxel")
VtkHexahedron = VtkCellType(12, "Hexahedron")
VtkWedge = VtkCellType(13, "Wedge")
VtkPyramid = VtkCellType(14, "Pyramid")
VtkQuadraticEdge = VtkCellType(21, "Quadratic_Edge")
VtkQuadraticTriangle = VtkCellType(22, "Quadratic_Triangle")
VtkQuadraticQuad = VtkCellType(23, "Quadratic_Quad")
VtkQuadraticTetra = VtkCellType(24, "Quadratic_Tetra")
VtkQuadraticHexahedron = VtkCellType(25, "Quadratic_Hexahedron")

def _vtk_data_type(dtype_name):
    """ Return the VtkDataType for a numpy dtype name, or raise TypeError. """
    try:
        return np_to_vtk[str(dtype_name)]
    except KeyError:
        raise TypeError("numpy dtype '%s' is not supported by VTK XML files; supported: %s"
                        % (dtype_name, ", ".join(np_to_vtk))) from None

# ==============================
#       Helper functions
# ==============================
def _mix_extents(start, end):
    if len(start) != 3 or len(end) != 3:
        raise ValueError("start and end must have 3 entries, got %r and %r" % (start, end))
    string = "%d %d %d %d %d %d" % (start[0], end[0], start[1], end[1], start[2], end[2])
    return string

def _number_to_string(num):
    """ Plain text for a number. repr() is not used directly because, since
        NumPy 2, repr(np.float64(0.5)) is 'np.float64(0.5)', which VTK cannot parse. """
    if isinstance(num, bool):
        raise TypeError("expected a number, got a bool")
    if isinstance(num, numbers.Integral):
        return str(int(num))
    return repr(float(num))

def _array_to_string(a):
    s = "".join([_number_to_string(num) + " " for num in a])
    return s

def _get_byte_order():
    if sys.byteorder == "little":
        return "LittleEndian"
    else:
        return "BigEndian"

# ================================
#        VtkGroup class
# ================================
class VtkGroup:
    
    def __init__(self, filepath):
        """ Creates a VtkGroup file that is stored in filepath.
            
            PARAMETERS:
                filepath: filename without extension.
        """
        self.xml = XmlWriter(filepath + ".pvd")
        self.xml.openElement("VTKFile")
        self.xml.addAttributes(type = "Collection", version = "0.1",  byte_order = _get_byte_order())
        self.xml.openElement("Collection")
        self.root = os.path.dirname(filepath)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is not None:
            self.xml.abort()
        elif not self.xml.closed:
            self.save()
        return False

    def save(self):
        """ Closes this VtkGroup. """
        self.xml.closeElement("Collection")
        self.xml.closeElement("VTKFile")
        self.xml.close()

    def abort(self):
        """ Discard this VtkGroup without writing it. """
        self.xml.abort()
    
    def addFile(self, filepath, sim_time, group = "", part = "0"):
        """ Adds file to this VTK group.

            PARAMETERS:
                filepath: full path to VTK file.
                sim_time: simulated time.
                group: This attribute is not required; it is only for informational purposes.
                part: It is an integer value greater than or equal to 0.
                
            See: http://www.paraview.org/Wiki/ParaView/Data_formats#PVD_File_Format for details.
        """
        # TODO: Check what the other attributes are for.
        filename = os.path.relpath(filepath, start = self.root)
        self.xml.openElement("DataSet")
        self.xml.addAttributes(timestep = sim_time, group = group, part = part, file = filename)
        self.xml.closeElement()
        


# ================================
#        VtkFile class         
# ================================
class VtkFile:
    
    def __init__(self, filepath, ftype, largeFile = False):
        """
            PARAMETERS:
                filepath: filename without extension.
                ftype: file type, e.g. VtkImageData, etc.
                largeFile: If size of the stored data cannot be represented by a UInt32.
        """
        self.ftype = ftype
        self.filename = filepath + ftype.ext
        self.xml = XmlWriter(self.filename)
        self.offset = 0  # offset in bytes after beginning of binary section
        self.appendedDataIsOpen = False
#        self.largeFile = largeFile

#       if largeFile == False:
#            self.xml.openElement("VTKFile").addAttributes(type = ftype.name,
#                                                          version = "0.1",
#                                                          byte_order = _get_byte_order())
#        else:
#           print "WARNING: output file only compatible with VTK 6.0 and later."
        self.xml.openElement("VTKFile").addAttributes(type = ftype.name,
                                                          version = "1.0",
                                                          byte_order = _get_byte_order(),
                                                          header_type = "UInt64")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        """ On error, discard the partial file; otherwise save it (if not saved yet). """
        if exc_type is not None:
            self.abort()
        elif not self.xml.closed:
            self.save()
        return False

    def abort(self):
        """ Discard the file being written. An existing file with the same name is kept. """
        self.xml.abort()

    def addComments(self, comments):
        """ Insert strings stored in comments list as comments into the xml header section of the file. 
            This method does not make any check for forbidden or special characters in XML comments. 
            This method should be called right after creating the file to avoid including comments in the 
            wrong place within the header. BE AWARE!!!
            Illegal sequences in comments ("--") are made safe automatically.
        """
        if self.appendedDataIsOpen:
            raise RuntimeError("comments must be added before the appended data section")
        for c in comments:
            self.xml.addComment(c)
        
    def getFileName(self):
        """ Returns absolute path to this file. """
        return  os.path.abspath(self.filename)

    def openPiece(self, start = None, end = None,
                        npoints = None, ncells = None,
                        nverts = None, nlines = None, nstrips = None, npolys = None): 
        """ Open piece section.
            
            PARAMETERS:
                Next two parameters must be given together.
                start: array or list with start indexes in each direction.
                end:   array or list with end indexes in each direction.

                npoints: number of points in piece (int).
                ncells: number of cells in piece (int). If present,
                        npoints must also be given.

                All the following parameters must be given together with npoints.
                They should all be integer values.
                nverts: number of vertices.
                nlines: number of lines.
                nstrips: number of strips.
                npolys: number of .

            RETURNS:
                this VtkFile to allow chained calls.
        """
        # TODO: Check what are the requirements for each type of grid.

        if (start is not None and end is not None):
            ext = _mix_extents(start, end)
            self.xml.openElement("Piece")
            self.xml.addAttributes( Extent = ext)

        elif (ncells is not None and npoints is not None):
            self.xml.openElement("Piece")
            self.xml.addAttributes(NumberOfPoints = int(npoints), NumberOfCells = int(ncells))

        elif any(v is not None for v in (npoints, nverts, nlines, nstrips, npolys)):
            if npoints is None: npoints = 0
            if nverts is None: nverts = 0
            if nlines is None: nlines = 0
            if nstrips is None: nstrips = 0
            if npolys is None: npolys = 0
            self.xml.openElement("Piece")
            self.xml.addAttributes(NumberOfPoints = int(npoints), NumberOfVerts = int(nverts),
                                   NumberOfLines = int(nlines), NumberOfStrips = int(nstrips),
                                   NumberOfPolys = int(npolys))
        else:
            raise ValueError("openPiece needs start and end, or npoints and ncells, "
                             "or npoints with nverts/nlines/nstrips/npolys")

        return self

    def closePiece(self):
        self.xml.closeElement("Piece")

    def openData(self, nodeType, scalars=None, vectors=None, normals=None, tensors=None, tcoords=None):
        """ Open data section.

            PARAMETERS:
                nodeType: Point or Cell.
                scalars: default data array name for scalar data.
                vectors: default data array name for vector data.
                normals: default data array name for normals data.
                tensors: default data array name for tensors data.
                tcoords: dafault data array name for tcoords data.

            RETURNS:
                this VtkFile to allow chained calls.
        """
        self.xml.openElement(nodeType + "Data")
        # VTK matches these attribute names case-sensitively ("Scalars", not "scalars");
        # with lowercase names the active scalars/vectors/... were silently not set.
        if scalars:
            self.xml.addAttributes(Scalars = scalars)
        if vectors:
            self.xml.addAttributes(Vectors = vectors)
        if normals:
            self.xml.addAttributes(Normals = normals)
        if tensors:
            self.xml.addAttributes(Tensors = tensors)
        if tcoords:
            self.xml.addAttributes(TCoords = tcoords)

        return self

    def closeData(self, nodeType):
        """ Close data section.

            PARAMETERS:
                nodeType: Point or Cell.
 
            RETURNS:
                this VtkFile to allow chained calls.
        """
        self.xml.closeElement(nodeType + "Data")


    def openGrid(self, start = None, end = None, origin = None, spacing = None, time_values = None):
        """ Open grid section.

            PARAMETERS:
                start: array or list of start indexes. Required for Structured, Rectilinear and ImageData grids.
                end: array or list of end indexes. Required for Structured, Rectilinear and ImageData grids.
                origin: 3D array or list with grid origin. Only required for ImageData grids.
                spacing: 3D array or list with grid spacing. Only required for ImageData grids.
                time_values: single string giving the time values (with spaces in-between). Only required for Unstructured grids.

            RETURNS:
                this VtkFile to allow chained calls.
        """
        gType = self.ftype.name
        if (gType == VtkImageData.name):
            if (start is None or end is None or origin is None or spacing is None):
                raise ValueError("ImageData needs start, end, origin and spacing")
            if len(origin) != 3 or len(spacing) != 3:
                raise ValueError("origin and spacing must have 3 entries, got %r and %r"
                                 % (origin, spacing))
            ext = _mix_extents(start, end)
            self.xml.openElement(gType)
            self.xml.addAttributes(WholeExtent = ext,
                                   Origin = _array_to_string(origin),
                                   Spacing = _array_to_string(spacing))

        elif (gType == VtkStructuredGrid.name or gType == VtkRectilinearGrid.name):
            if (start is None or end is None):
                raise ValueError("%s needs start and end" % gType)
            ext = _mix_extents(start, end)
            self.xml.openElement(gType)
            self.xml.addAttributes(WholeExtent = ext)
        else:
            self.xml.openElement(gType)
            if (gType == VtkUnstructuredGrid.name and time_values is not None):
                self.xml.addAttributes(TimeValues = time_values)
                
        return self

    def closeGrid(self):
        """ Closes grid element.

            RETURNS:
                this VtkFile to allow chained calls.
        """
        self.xml.closeElement(self.ftype.name)

    
    def addHeader(self, name, dtype, nelem, ncomp, time_value = None):
        """ Adds data array description to xml header section.

            PARAMETERS:
                name: data array name.
                dtype: string describing type of the data.
                       Format is the same as used by numpy, e.g. 'float64'.
                nelem: number of elements in the array.
                ncomp: number of components, 1 (=scalar) and 3 (=vector).
                time_value: string representing time value.

            RETURNS:
                This VtkFile to allow chained calls.
            
            NOTE: This is a low level function. Use addData if you want
                  to add a numpy array.
        """
        dtype = _vtk_data_type(dtype)

        self.xml.openElement("DataArray")
        if time_value is not None:
            self.xml.addAttributes( Name = name,
                                    NumberOfComponents = ncomp,
                                    type = dtype.name,
                                    format = "appended",
                                    offset = self.offset,
                                    TimeStep = time_value)
        else:
            self.xml.addAttributes( Name = name,
                                    NumberOfComponents = ncomp,
                                    type = dtype.name,
                                    format = "appended",
                                    offset = self.offset)
        self.xml.closeElement()

        #TODO: Check if 4/8 is platform independent
        #if self.largeFile == False:
        #    self.offset += nelem * ncomp * dtype.size + 4 # add 4 to indicate array size
        #else:
        self.offset += nelem * ncomp * dtype.size + 8 # add 8 to indicate array size
        return self

    # SWW: keep this for now...
    def addData(self, name, data):
        """ Adds array description to xml header section.
            
             PARAMETERS:
                name: data array name.
                data: one numpy array or a tuple with 3 numpy arrays. If a tuple, the individual
                      arrays must represent the components of a vector field.
                      All arrays must be one dimensional or three-dimensional.
        """
        if isinstance(data, tuple): # vector data
            if len(data) != 3:
                raise ValueError("vector data must be a tuple of 3 arrays, got %d" % len(data))
            x = data[0]
            self.addHeader(name, x.dtype.name, x.size, 3)
        elif isinstance(data, np.ndarray):
            if data.ndim == 1 or data.ndim == 3:
                self.addHeader(name, data.dtype.name, data.size, 1)
            else:
                raise ValueError("Bad array shape for '%s': %s (must be 1-D or 3-D)"
                                 % (name, data.shape))
        else:
            raise TypeError("data for '%s' must be a numpy array or a tuple of 3 arrays" % name)

    def internal_addData(self, name, data, ncomp, time_value = None):
        """ Adds array description to xml header section.
            
             PARAMETERS:
                name: data array name.
                data: one numpy array (a 1D array).
                ncomp: number of components to the data.
                time_value: string representing time value.
        """
        if isinstance(data, np.ndarray):
            if data.ndim == 1 or data.ndim == 3:
                if ncomp < 1 or data.size % ncomp != 0:
                    raise ValueError("size of '%s' (%d) is not a multiple of its %d components"
                                     % (name, data.size, ncomp))
                nelem = data.size // ncomp
                self.addHeader(name, data.dtype.name, nelem, ncomp, time_value)
            else:
                raise ValueError("Bad array shape for '%s': %s (must be 1-D or 3-D)"
                                 % (name, data.shape))
        else:
            raise TypeError("data for '%s' must be a numpy array" % name)

    def appendHeader(self, dtype, nelem, ncomp):
        """ This function only writes the size of the data block that will be appended.
            The data itself must be written immediately after calling this function.
            
            PARAMETERS:
                dtype: string with data type representation (same as numpy). For example, 'float64'
                nelem: number of elements.
                ncomp: number of components, 1 (=scalar) or 3 (=vector).
        """
        self.openAppendedData()
        dsize = np_to_vtk[dtype].size
        block_size = dsize * ncomp * nelem
        if self.largeFile == False:
            writeBlockSize(self.xml.stream, block_size)
        else:
            writeBlockSize64Bit(self.xml.stream, block_size)

            
    def appendData(self, data):
        """ Append data to binary section.
            This function writes the header section and the data to the binary file.

            PARAMETERS:
                data: one numpy array or a tuple with 3 numpy arrays. If a tuple, the individual
                      arrays must represent the components of a vector field.
                      All arrays must be one dimensional or three-dimensional.
                      The order of the arrays must coincide with the numbering scheme of the grid.
            
            RETURNS:
                This VtkFile to allow chained calls

            TODO: Extend this function to accept contiguous C order arrays.
        """
        self.openAppendedData()

        if isinstance(data, tuple): # 3 numpy arrays
            ncomp = len(data)
            if ncomp != 3:
                raise ValueError("vector data must be a tuple of 3 arrays, got %d" % ncomp)
            dsize = data[0].dtype.itemsize
            nelem = data[0].size
            block_size = ncomp * nelem * dsize
            #if self.largeFile == False:
            writeBlockSize(self.xml.stream, block_size)
            #else:
            #    writeBlockSize64Bit(self.xml.stream, block_size)
            x, y, z = data[0], data[1], data[2]
            writeArraysToFile(self.xml.stream, x, y, z)
            
        elif isinstance(data, np.ndarray) and (data.ndim == 1 or data.ndim == 3): # single numpy array
            ncomp = 1 
            dsize = data.dtype.itemsize
            nelem = data.size
            block_size = ncomp * nelem * dsize
            #if self.largeFile == False:
            writeBlockSize(self.xml.stream, block_size)
            #else:
            #    writeBlockSize64Bit(self.xml.stream, block_size)
            writeArrayToFile(self.xml.stream, data)

        else:
            raise TypeError("appendData expects a 1-D or 3-D numpy array or a tuple of 3 arrays, "
                            "got %s" % (getattr(data, "shape", type(data).__name__),))

        return self

    def openAppendedData(self):
        """ Opens binary section.

            It is not necessary to explicitly call this function from an external library.
        """
        if not self.appendedDataIsOpen:
            self.xml.openElement("AppendedData").addAttributes(encoding = "raw").addText("_")
            self.appendedDataIsOpen = True

    def closeAppendedData(self):
        """ Closes binary section.

            It is not necessary to explicitly call this function from an external library.
        """
        self.xml.closeElement("AppendedData")

    def openElement(self, tagName):
        """ Useful to add elements such as: Coordinates, Points, Verts, etc. """
        self.xml.openElement(tagName)

    def closeElement(self, tagName):
        self.xml.closeElement(tagName)

    def save(self):
        """ Closes file """
        if self.appendedDataIsOpen:
            self.xml.closeElement("AppendedData")
        self.xml.closeElement("VTKFile")
        self.xml.close()

#
