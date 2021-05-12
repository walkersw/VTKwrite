#! /usr/bin/env python

######################################################################################
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

# ************************************************************************
# * Example of how to use the high level unstructuredGridToVTK function. *
# * This example shows how to export a unstructured grid give its        *
# * nodes and topology through a connectivity and offset lists.          *
# * Check the VTK file format for details of the unstructured grid.      *
# ************************************************************************
import os
# from VTKwrite.interface import unstructuredGridToVTK
from VTKwrite.interface import my_unstructuredGridToVTK
from VTKwrite.vtkbin import VtkTriangle, VtkQuad
import numpy as np

FILE_PATH = "./unstructured"
def clean():
    try:
        os.remove(FILE_PATH + ".vtu")
    except:
        pass

def run():
    print("MY running unstructured...")

    # Define vertices
    x = np.zeros(6)
    y = np.zeros(6)
    z = np.zeros(6)

    x[0], y[0], z[0] = 0.0, 0.0, 0.0
    x[1], y[1], z[1] = 1.0, 0.0, 0.0
    x[2], y[2], z[2] = 2.0, 0.0, 0.0
    x[3], y[3], z[3] = 0.0, 1.0, 0.0
    x[4], y[4], z[4] = 1.0, 1.0, 0.0
    x[5], y[5], z[5] = 2.0, 1.0, 0.0

    # Define connectivity or vertices that belongs to each element
    conn = np.zeros(10)

    conn[0], conn[1], conn[2] = 0, 1, 3              # first triangle
    conn[3], conn[4], conn[5] = 1, 4, 3              # second triangle
    conn[6], conn[7], conn[8], conn[9] = 1, 2, 5, 4  # rectangle

    # Define offset of last vertex of each element
    offset = np.zeros(3)
    offset[0] = 3
    offset[1] = 6
    offset[2] = 10

    # Define cell types

    ctype = np.zeros(3)
    ctype[0], ctype[1] = VtkTriangle.tid, VtkTriangle.tid
    ctype[2] = VtkQuad.tid
    
    # initialize the data structure
    all_cell_data = {"scalars" : None, "vectors" : None}
    # scalars
    cellData_sc = {"pressure0" : np.random.rand(3)}
    cellData_sc["pressure1"] = np.random.rand(3)
    all_cell_data["scalars"] = cellData_sc
    # vectors
    cellData_vc = {"vel0" : np.random.rand(3*3)}
    cellData_vc["vel1"] = np.random.rand(3*3)
    all_cell_data["vectors"] = cellData_vc
    
    # initialize the data structure
    all_point_data = {"scalars" : None, "vectors" : None}
    # scalars
    pointData_sc = {"potential0" : np.random.rand(6)}
    pointData_sc["potential1"] = np.random.rand(6)
    all_point_data["scalars"] = pointData_sc
    # vectors
    pointData_vc = {"flux0" : np.random.rand(6*3)}
    pointData_vc["flux1"] = np.random.rand(6*3)
    all_point_data["vectors"] = pointData_vc
    
    comments = [ "comment 1", "comment 2" ]
    my_unstructuredGridToVTK(FILE_PATH, x, y, z, connectivity = conn, offsets = offset, cell_types = ctype, all_cell_data = all_cell_data, all_point_data = all_point_data, comments = comments)


if __name__ == "__main__":
    run()
