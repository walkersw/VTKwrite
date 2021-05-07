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

import os
import setuptools

# get key package details from pyvtk/__version__.py
about = {}  # type: ignore
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'pyvtk', '__version__.py')) as f:
    exec(f.read(), about)

# load the README file and use it as the long_description for PyPI
with open('README.md', 'r') as f:
    readme = f.read()

# package configuration - for reference see:
# https://setuptools.readthedocs.io/en/latest/setuptools.html#id9
setuptools.setup(
    name='VTKwrite',
    description=about['__description__'],
    long_description=readme,
    long_description_content_type='text/markdown',
    version=about['__version__'],
    author=about['__author__'],
    author_email=about['__author_email__'],
    url=about['__url__'],
    packages=['VTKwrite'],
    package_dir = {'VTKwrite' : 'pyvtk'},
    package_data = {'VTKwrite' :  ['../LICENSE', '../examples/*.py']},
    project_urls={
        "Bug Tracker": "https://github.com/walkersw/VTKwrite",
        "Documentation": "https://github.com/walkersw/VTKwrite",
        "Source Code": "https://github.com/walkersw/VTKwrite",
    },
    #include_package_data=True,
    python_requires=">=3.6",
    install_requires=['numpy'],
    license=about['__license__'],
    #zip_safe=False,
    #entry_points={
    #    'console_scripts': ['py-package-template=py_pkg.entry_points:main'],
    #},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GENERAL PUBLIC LICENSE",
        "Operating System :: OS Independent",
    ],
    keywords='vtk'
)
