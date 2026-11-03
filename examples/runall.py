import shutil

import array_input
import group
import image
import lagrange_cells
import lines
import points
import poly_lines
import rectilinear 
import structured 
import structured_ex1 
import unstructured 
import unstructured_timedep 
import unstructured_timeseries
import low_level

def testit(test):
    try:
        test()
    except:
        print("  FAILED")

def clean_all():
    array_input.clean()
    group.clean()
    image.clean()
    lagrange_cells.clean()
    lines.clean()
    points.clean()
    poly_lines.clean()
    rectilinear.clean()
    structured.clean()
    structured_ex1.clean()
    unstructured.clean()
    unstructured_timedep.clean()
    unstructured_timeseries.clean()
    low_level.clean()
    try:
        shutil.rmtree("__pycache__")
    except:
        pass
    
def test_all():
    testit(array_input.run)
    testit(group.run)
    testit(image.run)
    testit(lagrange_cells.run)
    testit(lines.run)
    testit(points.run)
    testit(poly_lines.run)
    testit(rectilinear.run)
    testit(structured.run)
    testit(structured_ex1.run)
    testit(unstructured.run)
    testit(unstructured_timedep.run)
    testit(unstructured_timeseries.run)
    testit(low_level.run)

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        opt = sys.argv[1]
    else:
        opt = "-"
    
    if (opt == "run"):
        test_all()
    elif (opt == "clean"):
        clean_all()
    else:
        print("UNKNOWN OPTION: " + opt)
        print("USE: python runall.py [run|clean]")
        
    print("*** ALL DONE ***")

