"""
VTKwrite.pyvtk.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This is the Python interface for writing data.

Copyright (c) 05-07-2021,  Shawn W. Walker
"""

import struct
import sys

import numpy as np

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

# Points are interleaved and written in chunks of this many points, which bounds the
# extra memory needed for large meshes (3 * 2**20 values = 24 MB for float64).
_CHUNK_POINTS = 2**20

def _native(a):
    """ Same values with the machine's byte order (the file declares the native order). """
    if not a.dtype.isnative:
        a = a.astype(a.dtype.newbyteorder("="))
    return a

# ================================
#        write functions
# ================================
def writeBlockSize(stream, block_size):
    fmt = _get_byte_order_char() + 'Q' # Write size as unsigned long long == 64 bits unsigned integer
    stream.write(struct.pack(fmt, block_size))

def writeArrayToFile(stream, data):
    """ Write the values of a 1-D or 3-D array in Fortran order (the order VTK expects
        for 3-D arrays), in the machine's byte order. """
    if data.ndim not in (1, 3):
        raise ValueError("array must be 1-D or 3-D, got shape %s" % (data.shape,))
    _struct_code(data.dtype)    # raises TypeError for unsupported dtypes
    # ravel(order='F') is a view for 1-D and Fortran-ordered arrays and a contiguous
    # copy otherwise; either way the result can be written as one buffer.
    stream.write(np.ravel(_native(data), order='F'))

# ==============================================================================
def writeArraysToFile(stream, x, y, z):
    """ Write three coordinate (or component) arrays interleaved as x0 y0 z0 x1 y1 z1 ...,
        each taken in Fortran order, using the dtype of x. """
    # Check if arrays have same shape and data type
    if not (x.size == y.size == z.size):
        raise ValueError("x, y, z have different sizes: %d, %d, %d" % (x.size, y.size, z.size))
    if not (x.dtype.itemsize == y.dtype.itemsize == z.dtype.itemsize):
        raise ValueError("x, y, z have different item sizes (dtypes %s, %s, %s)"
                         % (x.dtype, y.dtype, z.dtype))
    _struct_code(x.dtype)       # raises TypeError for unsupported dtypes

    # NOTE: VTK expects data in FORTRAN order
    # This is only needed when a multidimensional array has C-layout.  1-D arrays are used
    # as they are, also when strided (e.g. the columns of an (n, 3) points array), so that
    # they are not copied as a whole: the chunks below copy them piece by piece.
    xx, yy, zz = (c if c.ndim == 1 else np.ravel(c, order='F') for c in (x, y, z))

    nitems = xx.size
    dtype = x.dtype.newbyteorder("=")
    buf = np.empty((min(nitems, _CHUNK_POINTS), 3), dtype = dtype)
    for start in range(0, nitems, _CHUNK_POINTS):
        stop = min(start + _CHUNK_POINTS, nitems)
        block = buf[:stop - start]
        block[:, 0] = xx[start:stop]   # assignment converts y, z to the dtype of x
        block[:, 1] = yy[start:stop]
        block[:, 2] = zz[start:stop]
        stream.write(block)

#
