

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8221491.svg)](https://doi.org/10.5281/zenodo.8221491)


# TexGen - MacOS Fork 
TexGen is a geometric textile modelling software package to be used for obtaining engineering properties of woven textiles and textile composites. More information about the project can be found here http://texgen.sourceforge.net.
This repo is a modification of the original project, ready to build for MacOS on Apple Silicon. It has been tested on MacOS Sequioa 15.1. Functionality is limited to 
the python interface due to missing dependencies. Octree Refinements have also been deactivated. A number of functions have not been tested on this version, making it highly experimental.
The modifications are minimal, most of the complexity lies in finding settings that will comply to make a usable program. 



[![Download TexGen](https://img.shields.io/sourceforge/dt/texgen.svg)](https://sourceforge.net/projects/texgen/files/latest/download)

## Building Code


### Python Installation
This section is written under the assumption, that python 3.9 will be used and installed via homebrew. The most recent TexGen Version (which I forked here) can use newer 
versions of python, so this is not a requirement anymore. 
If not the case already, install python via homebrew: 

```
brew install python@3.9
```

You can verify the path of the installation by using 
```
ls -l /opt/homebrew/bin/python3.9
```
If it produces a result, you should be good. Adjust path is necessary with different versions of python.

### Boost Installation
The compatibility layer for python and C++ is required to succesfully run TexGen  (at least in my experience, on mac). If you use python 3.9, a custom build of boost might be required. Boost tarball can be downloaded on boost.org.
After downloading and unpacking, add your python installation path to the ```project-config.jam```:
```
using python
    : 3.9
    : /opt/homebrew/bin/python3.9
    : /opt/homebrew/opt/python@3.9/Frameworks/Python.framework/Versions/3.9/Headers
    : /opt/homebrew/opt/python@3.9/Frameworks/Python.framework/Versions/3.9/lib/libpython3.9.dylib
    ;
```
Change your paths as required and insert at the end of the file. 
In the boost directory, execute
```
./bootstrap.sh --with-libraries=python
```
Afterwards, run 
```
sudo ./b2 \
  --with-python               \
  link=shared threading=multi \
  --prefix=/opt/homebrew/opt/boost \
  install
```
to install boost for the specified python version. 

### Actual Build process
Clone the repository into a local folder and navigate into the ```TexGen```-Directory. 
Then, execute the following commands: 
```
cmake .. -G "Unix Makefiles" -DTRACE_CORE=ON\
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_OSX_ARCHITECTURES=arm64 \
  -DBUILD_SHARED_LIBS=ON \
  -DBUILD_GUI=OFF            -DBUILD_RENDERER=OFF \
  -DBUILD_EXPORT=OFF         -DBUILD_UNIT_TESTS=OFF \
  -DBUILD_DOCUMENTATION=OFF  -DBUILD_EXAMPLES=OFF \
  -DBUILD_CASCADE_EXPORT=OFF -DBUILD_PYTHON_INTERFACE=ON \
  \
  -DSWIG_EXECUTABLE=/opt/homebrew/bin/swig \
  \
  -DPYTHON_EXECUTABLE="/opt/homebrew/bin/python3.9" \
  -DPYTHON_INCLUDE_DIRS="/opt/homebrew/opt/python@3.9/Frameworks/Python.framework/Versions/3.9/include/python3.9" \
  -DPYTHON_LIBRARIES="/opt/homebrew/opt/python@3.9/Frameworks/Python.framework/Versions/3.9/Python" \
  -DPYTHON_SITEPACKAGES_DIR="/opt/homebrew/opt/python@3.9/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages" \
  -DCMAKE_INSTALL_RPATH="/opt/homebrew/lib" \
  \
  -DCMAKE_DISABLE_FIND_PACKAGE_Qt6=ON \
  -DCMAKE_DISABLE_FIND_PACKAGE_Qt6Widgets=ON
  ```
There is some debugging in here, just ignore that. The build instructions show 
the degree to which modules have been disabled. If this executes succesfully, 
use 
```
make  -j"$(sysctl -n hw.ncpu)"
sudo make install
```
Sudo is required here to be able to access the local python libraries, use at your own risk, or dont install python where I have it if you arent comfortable with it. 
At this point, the software should be building. 

## Using TexGen via the Python Interface
I am sure there is another, easier way to do this. This is the one that works for me: 
```
python_module_parent_dir = "/opt/homebrew/lib/python3.9/site-packages"
library_dir = "/opt/homebrew/lib"

if python_module_parent_dir not in sys.path:
    sys.path.insert(0, python_module_parent_dir)

if library_dir not in sys.path:
    sys.path.insert(0, library_dir)

os.environ["DYLD_LIBRARY_PATH"] = library_dir + os.pathsep + os.environ.get("DYLD_LIBRARY_PATH", "")
from TexGen.Core import * 

```
Once again, make sure to adjust local python installation directories if version or type of installation differ. 


## Citing TexGen
This is a fork of a really powerful open source project, 
make sure to cite the original creators when publishing work that relies on it! These are the original Citations given by the original authors of the software: 

L P Brown and A C Long. "Modelling the geometry of textile reinforcements for composites: TexGen", Chapter 8 in "Composite reinforcements for optimum performance (Second Edition)", ed. P Boisse, Woodhead Publishing Ltd, 2021, ISBN: 978-0-12-819005-0. https://doi.org/10.1016/B978-0-12-819005-0.00008-3

Lin, H., Brown, L. P. & Long, A. C. 2011. Modelling and Simulating Textile Structures using TexGen. Advanced Materials Research, 331, 44-47

Louise Brown, mike-matveev, & georgespackman. (2023). louisepb/TexGen: TexGen v3.13.1 (v3.13.1). Zenodo. https://doi.org/10.5281/zenodo.8221491
