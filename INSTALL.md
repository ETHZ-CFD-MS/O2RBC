Installation instructions                           {#installpage}
=========================

This file contains the installation instructions for O2RBC, including the Python
library used for pre- and postprocessing of simulations.

The library requires

  - [OpenFOAM 4.1](https://openfoam.org/release/4-1/)
  - the development version of swak4foam
  - boost
  - Python 2.7 with igraph, NumPy, SciPy and Matplotlib
  - [PyFoam](https://openfoamwiki.net/index.php/Contrib/PyFoam)

OpenFOAM
--------

OpenFOAM 4.1 can be downloaded [here](https://openfoam.org/version/4-1/).
Detailed compilation instructions are available for the every supported platform.
The code was compiled and tested on CentOS 7.3 with gcc 4.8.5 and openmpi 1.6.5.

As usual with OpenFOAM, its use requires sourcing the appropriate configuration
file. With bash, this is done with

    . $WM_PROJECT_DIR/etc/bashrc

and with tcsh

    source $WM_PROJECT_DIR/etc/tcshrc

swak4foam
---------

The library [swak4foam](http://openfoamwiki.net/index.php/Contrib/swak4Foam) 
(SWiss Army knife for Foam) provides function objects and expressions for source
terms which are used by O2RBC. At time of writing, only the development version
of swak4foam supports OpenFOAM 4.1. The original compilation instructions are
available 
[here](https://openfoamwiki.net/index.php/Installation/swak4Foam/Downloading#OpenFOAM_4.0_to_4.x)
and 
[here](https://openfoamwiki.net/index.php/Installation/swak4Foam/Downloading#swak4Foam_development_version).

In a nutshell, it can compiled as follows:

    git clone https://github.com/Unofficial-Extend-Project-Mirror/openfoam-extend-swak4Foam-dev.git swak4Foam
    cd swak4Foam
    git checkout branches/develop
    ./Allwmake

Depending on the installed version of Bison, it might be necessary to run

    ./maintainanceScripts/compileRequirements.sh

prior to running ./Allwmake.

The compilation of the solvers hbPOAxisymmetricFoam and hbPOGraphFoam requires the 
environment variable SWAK4FOAMDIR. In this case, it should be set to /path/to/swak4Foam.

boost
-----

The [boost](http://www.boost.org/) library is used for the graph structure
underlying the capillary networks. Code compilation requires the environment
variable BOOST_INCLUDEDIR to be set to the path that contains the folder 'boost'
with the library header (\*.hpp) files. In a bash shell, this is done with

   export BOOST_INCLUDEDIR="/path/to/folder"

Versions 1.53 to 1.64 have been successfully employed.

O2RBC
-----

The compilation of the OpenFOAM library extension requires the following
environment variables that describe commonly used paths:

  - OF_CBF for the O2RBC library root folder
  - OF_CBF_SRC for the src folder
  - OF_SCRIPTS for the folder that contains bash and python scripts

When using bash, it is most convenient to put the following lines into the
bash profile (typically, ~/.bashrc):

    export OF_CBF="/path/to/O2RBC"
    export OF_CBF_SRC=$OF_CBF/src
    export OF_SCRIPTS=$OF_CBF/scripts

The execution of test cases and various scripts also requires the PATH
environment variable to be completed. With bash:

    export PATH=$PATH:$OF_SCRIPTS:$OF_SCRIPTS/bash:$OF_SCRIPTS/python/bin

Once this has been set up, the library can be compiled by running

    ./Allwmake

in the library's root directory. If the library has previously been compiled
with another version of OpenFOAM or different compiler flags, it is advised to
run

    wclean all

before compiling the library.

Python library
--------------

Some tutorials that execute parameter studies, as well as many postprocessing
utilities are based on Python. Therefore

    export PYTHONPATH=$PYTHONPATH:$OF_CBF/scripts:$OF_CBF/scripts/python:$OF_CBF/scripts/python/bin

Python 2.7.5 was used. The code makes heavy use of NumPy, SciPy and Matplotlib. 
These packages are available [here](https://scipy.org/install.html).
Some functionalities require Scipy >= 0.17.0. For postprocessing,
[igraph](http://igraph.org/python/) >= 0.7 is required.

PyFoam
------

The library [PyFoam](https://openfoamwiki.net/index.php/Contrib/PyFoam) is used
to read and write OpenFOAM dictionaries using Python. Installation instructions
are available
[here](https://openfoamwiki.net/index.php/Contrib/PyFoam#Installation). For
these functionalities to work, the path where the PyFoam folder is located needs
to be in a default path used by Python or in the environment variable
PYTHONPATH.

PyFoam 0.6.6 was used.

Documentation
-------------

O2RBC is documented using [doxygen](www.doxygen.org). The Doxygen-generated 
HTML documentation provides the best way to navigate through the code. Doxygen
is available on multiple platforms and can be downloaded
[here](http://www.stack.nl/~dimitri/doxygen/download.html).

The documentation of O2RBC is built as follows

    cd $OF_CBF/doc
    ./Allwmake

After compilation, the HTML documentation will be available at

    $OF_CBF/doc/Doxygen/html/index.html

which can be opened using your favourite web browser.

Doxygen >= 1.8.0 should be used since Markdown support is required.

Notes for macOS
--------------

Currently, packaged versions of OpenFOAM for macOS are provided using Docker for
Mac (see [here](https://openfoam.org/download/4-1-macos)). While our OpenFOAM
library extension works well in this self-contained environment, some required Python
packages are absent. While they can probably be installed manually, it is
probably more convenient to run pre- and processing Python scripts outside the
Docker environment.

