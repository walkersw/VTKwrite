"""
VTKwrite.pyvtk.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the Python interface for writing data.

Copyright (c) 05-07-2021,  Shawn W. Walker
"""

import struct
import sys
try:
    import numpy as np
except:
    print("Numpy is not installed. Please install it before running VTKwrite again.")

# Map numpy dtype to struct format
np_to_struct = { 'int8'    : 'b',
               'uint8'   : 'B',
               'int16'   : 'h',
               'uint16'  : 'H',
               'int32'   : 'i',
               'uint32'  : 'I',
               'int64'   : 'q',
               'uint64'  : 'Q',
               'float32' : 'f',
               'float64' : 'd' }
              
def _struct_code(dtype):
    try:
        return np_to_struct[dtype.name]
    except KeyError:
        raise TypeError("numpy dtype '%s' is not supported; supported: %s"
                        % (dtype.name, ", ".join(np_to_struct))) from None

def _get_byte_order_char():
# Check format in https://docs.python.org/3.5/library/struct.html
    if sys.byteorder == "little":
        return '<'
    else:
        return '>'

# ================================
#        write functions
# ================================  
def writeBlockSize(stream, block_size):
    fmt = _get_byte_order_char() + 'Q' # Write size as unsigned long long == 64 bits unsigned integer
    stream.write(struct.pack(fmt, block_size))

def writeArrayToFile(stream, data):
    if data.ndim not in (1, 3):
        raise ValueError("array must be 1-D or 3-D, got shape %s" % (data.shape,))
    fmt = _get_byte_order_char() + str(data.size) + _struct_code(data.dtype)  # > for big endian

    # NOTE: VTK expects data in FORTRAN order
    # This is only needed when a multidimensional array has C-layout.
    # np.ravel also handles non-contiguous views (it copies when needed).
    dd = np.ravel(data, order='F')

    bin = struct.pack(fmt, *dd)
    stream.write(bin)
    
# ==============================================================================
def writeArraysToFile(stream, x, y, z):
    # Check if arrays have same shape and data type
    if not (x.size == y.size == z.size):
        raise ValueError("x, y, z have different sizes: %d, %d, %d" % (x.size, y.size, z.size))
    if not (x.dtype.itemsize == y.dtype.itemsize == z.dtype.itemsize):
        raise ValueError("x, y, z have different item sizes (dtypes %s, %s, %s)"
                         % (x.dtype, y.dtype, z.dtype))

    nitems = x.size

    fmt = _get_byte_order_char() + str(1) + _struct_code(x.dtype)  # > for big endian

    # NOTE: VTK expects data in FORTRAN order
    # This is only needed when a multidimensional array has C-layout
    xx = np.ravel(x, order='F')
    yy = np.ravel(y, order='F')
    zz = np.ravel(z, order='F')    
        
    # eliminate this loop by creating a composed array.
    for i in range(nitems):
        bx = struct.pack(fmt, xx[i])
        by = struct.pack(fmt, yy[i])
        bz = struct.pack(fmt, zz[i])
        stream.write(bx)
        stream.write(by)
        stream.write(bz)

#