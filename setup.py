import re
import sys
from distutils.core import setup, Extension


# Following the recommendations of PEP 396 we parse the version number 
# out of the module.
def parse_version(module_file):
    """
    Parses the version string from the specified file.
    
    This implementation is ugly, but there doesn't seem to be a good way
    to do this in general at the moment.
    """ 
    f = open(module_file)
    s = f.read()
    f.close()
    match = re.findall("__version__ = '([^']+)'", s)
    return match[0]



f = open("README.txt")
discsim_readme = f.read()
f.close()
discsim_version = parse_version("discsim.py") 
d = "lib/"
 
_discsim_module = Extension('_discsim', 
    sources = ["_discsimmodule.c", d + "util.c", d + "sim.c", d + "avl.c", 
            d + "nystrom.c"],
    libraries = ["gsl", "gslcblas"],
    include_dirs = [d])

requirements = ["ercs"]

setup(
    name = "discsim",
    version = discsim_version, 
    description = "Efficient coalescent simulation in continuous space",
    author = "Jerome Kelleher",
    author_email = "jerome.kelleher@ed.ac.uk",
    url = "http://pypi.python.org/pypi/discsim", 
    keywords = ["Coalescent simulation", "Continuous space"],
    license = "GNU GPLv3",
    platforms = ["POSIX"], 
    classifiers = [
        "Programming Language :: C",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    requires = requirements, 
    long_description = discsim_readme,
    ext_modules = [_discsim_module],
    py_modules = ['discsim'],
)

    
