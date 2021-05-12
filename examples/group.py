"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Example of how to create a VTK group to visualize time-dependent data. 

Copyright (c) 05-11-2021,  Shawn W. Walker
"""

import os
from VTKwrite.vtkbin import VtkGroup

FILE_PATH = "./group"
def clean():
    try:
        os.remove(FILE_PATH + ".pvd")
    except:
        pass

def run():
    print("Running group...")
    g = VtkGroup(FILE_PATH)
    g.addFile(filepath = "sim0000.vtu", sim_time = 0.0)
    g.addFile(filepath = "sim0001.vtu", sim_time = 1.0)
    g.addFile(filepath = "sim0002.vtu", sim_time = 2.0)
    g.addFile(filepath = "sim0003.vtu", sim_time = 3.0)
    g.save()

if __name__ == "__main__":
	run()
