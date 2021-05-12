"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to use the high level pointsToVTK function.

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
from VTKwrite.interface import pointsToVTK, pointsToVTKAsTIN
import numpy as np

FILE_PATH1 = "./rnd_points"
FILE_PATH2 = "./rnd_points_TIN"
FILE_PATH3 = "./line_points"
FILE_PATH4 = "./points_as_lists"
def clean():
    try:
        os.remove(FILE_PATH1 + ".vtu")
        os.remove(FILE_PATH2 + ".vtu")
        os.remove(FILE_PATH3 + ".vtu")
        os.remove(FILE_PATH4 + ".vtu")
    except:
        pass
    
def run():
    print("Running points...")
    
    # Example 1: Random points
    npoints = 100
    x = np.random.rand(npoints)
    y = np.random.rand(npoints)
    z = np.random.rand(npoints)
    
    pressure = np.random.rand(npoints)
    temp = np.random.rand(npoints)
    
    # initialize the data structure
    all_point_data = {"scalars" : None, "vectors" : None}
    # scalars
    pointData_sc = {"1_temp" : temp, "2_pressure" : pressure}
    all_point_data["scalars"] = pointData_sc

    # keys are sorted before exporting, hence it is useful to prefix a number to determine an order
    comments = [ "comment 1", "comment 2" ]
    pointsToVTK(FILE_PATH1, x, y, z, all_point_data = all_point_data, comments = comments) 

    # Example 2: Export as TIN
    ndim = 2 #only consider x, y coordinates to create the triangulation
    pointsToVTKAsTIN(FILE_PATH2, x, y, z, ndim = ndim, data = {"1_temp" : temp, "2_pressure" : pressure}, comments = comments)

    # Example 3: Regular point set
    x = np.arange(1.0,10.0,0.1)
    y = np.arange(1.0,10.0,0.1)
    z = np.arange(1.0,10.0,0.1)

    # initialize the data structure
    del all_point_data
    all_point_data = {"scalars" : None, "vectors" : None}
    # scalars
    pointData_sc = {"elev" : z}
    all_point_data["scalars"] = pointData_sc

    comments = [ "comment 1", "comment 2" ]
    pointsToVTK(FILE_PATH3, x, y, z, all_point_data = all_point_data, comments = comments)

    # Example 4: Point set of 5 points
    x = [0.0, 1.0, 0.5, 0.368, 0.4]
    y = [0.3, 2.0, 0.7, 0.1, 0.6]
    z = [1.0, 1.0, 0.3, 0.75, 0.9]
    
    pressure = [1.0, 2.0, 3.0, 4.0, 5.0]
    temp = [1.0, 2.0, 3.0, 4.0, 5.0]

    # initialize the data structure
    del all_point_data
    all_point_data = {"scalars" : None, "vectors" : None}
    # scalars
    pointData_sc = {"1_temp" : temp, "2_pressure" : pressure}
    all_point_data["scalars"] = pointData_sc

    # keys are sorted before exporting, hence it is useful to prefix a number to determine an order
    comments = [ "comment 1", "comment 2" ]
    pointsToVTK(FILE_PATH4, x, y, z, all_point_data = all_point_data, comments = comments) 

if __name__ == "__main__":
    run()
