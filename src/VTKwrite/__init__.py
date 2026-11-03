"""
VTKwrite
~~~~~~
The VTKwrite package - a Python package for writing VTK data in XML format.

The writers, the time-series classes, VtkFile/VtkGroup and the cell types can be
imported from the package itself, e.g.

    from VTKwrite import unstructuredGridToVTK, TimeSeriesUnstructuredGrid, VtkTriangle

(__all__ lists them).  The older imports from VTKwrite.interface and VTKwrite.vtkbin
still work.

Copyright (c) 02-23-2025,  Shawn W. Walker
"""

######################################################################################
# GNU GENERAL PUBLIC LICENSE, Version 3
# 
# Copyright (c) 2021 Shawn W. Walker
#
# See the LICENSE file.
#
# This package is based on another package "PyEVTK" by Paulo A. Herrera.
# His copyright notice is here:
#
# MIT License
# 
# Copyright (c) 2010-2021 Paulo A. Herrera
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
######################################################################################

# The single place where the version is set (pyproject.toml reads it from here).
__version__ = "0.2.0"

from .vtkbin import (
    VtkFile, VtkGroup,
    # file types, for VtkFile
    VtkImageData, VtkPolyData, VtkRectilinearGrid, VtkStructuredGrid, VtkUnstructuredGrid,
    # cell types (VtkTriangle.tid is the VTK id of a triangle, etc.)
    VtkCellType,
    VtkVertex, VtkPolyVertex, VtkLine, VtkPolyLine, VtkTriangle, VtkTriangleStrip,
    VtkPolygon, VtkPixel, VtkQuad, VtkTetra, VtkVoxel, VtkHexahedron, VtkWedge, VtkPyramid,
    VtkQuadraticEdge, VtkQuadraticTriangle, VtkQuadraticQuad, VtkQuadraticTetra,
    VtkQuadraticHexahedron, VtkQuadraticWedge, VtkQuadraticPyramid, VtkBiQuadraticQuad,
    VtkTriQuadraticHexahedron, VtkBiQuadraticTriangle, VtkCubicLine, VtkTriQuadraticPyramid,
    # arbitrary-order Lagrange cells (node order: lagrange_nodes)
    VtkLagrangeCurve, VtkLagrangeTriangle, VtkLagrangeQuadrilateral, VtkLagrangeTetrahedron,
    VtkLagrangeHexahedron, VtkLagrangeWedge,
)
from .lagrange import lagrange_nodes, lagrange_npoints, node_permutation, cellwise_mesh
from .interface import (
    imageToVTK, rectilinearToVTK, structuredToVTK, pointsToVTK, pointsToVTKAsTIN,
    linesToVTK, polyLinesToVTK, unstructuredGridToVTK, cylinderToVTK,
    TimeSeriesUnstructuredGrid,
)

# (__version__ is deliberately not in __all__: "from VTKwrite import *" should not
# replace the __version__ of the importing module)
__all__ = [
    # high-level writers
    "imageToVTK", "rectilinearToVTK", "structuredToVTK", "pointsToVTK", "pointsToVTKAsTIN",
    "linesToVTK", "polyLinesToVTK", "unstructuredGridToVTK", "cylinderToVTK",
    # time series
    "TimeSeriesUnstructuredGrid",
    # low level
    "VtkFile", "VtkGroup",
    "VtkImageData", "VtkPolyData", "VtkRectilinearGrid", "VtkStructuredGrid", "VtkUnstructuredGrid",
    # cell types
    "VtkCellType",
    "VtkVertex", "VtkPolyVertex", "VtkLine", "VtkPolyLine", "VtkTriangle", "VtkTriangleStrip",
    "VtkPolygon", "VtkPixel", "VtkQuad", "VtkTetra", "VtkVoxel", "VtkHexahedron", "VtkWedge",
    "VtkPyramid", "VtkQuadraticEdge", "VtkQuadraticTriangle", "VtkQuadraticQuad",
    "VtkQuadraticTetra", "VtkQuadraticHexahedron", "VtkQuadraticWedge", "VtkQuadraticPyramid",
    "VtkBiQuadraticQuad", "VtkTriQuadraticHexahedron", "VtkBiQuadraticTriangle", "VtkCubicLine",
    "VtkTriQuadraticPyramid",
    "VtkLagrangeCurve", "VtkLagrangeTriangle", "VtkLagrangeQuadrilateral",
    "VtkLagrangeTetrahedron", "VtkLagrangeHexahedron", "VtkLagrangeWedge",
    # Lagrange cells
    "lagrange_nodes", "lagrange_npoints", "node_permutation", "cellwise_mesh",
]
