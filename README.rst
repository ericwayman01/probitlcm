probitlcm
=========

This software implements the model described in the `manuscript <http://arxiv.org/abs/2408.13143>`_

  Eric Alan Wayman, Steven Andrew Culpepper, Jeff Douglas, and Jesse Bowers. "A Restricted Latent Class Model with Polytomous Ordinal Correlated Attributes and Respondent-Level Covariates." arXiv preprint arXiv:2408.13143, 2024.

If you make use of this software, please cite the above.

Requirements
------------

- Python (a version greater than or equal to 3.9), as well as the Python development headers (package ``python3-dev`` on Ubuntu)
- The Ninja and CMake build systems (packages ``ninja`` and ``cmake`` on Ubuntu)
- LAPACK and BLAS (packages ``libblas-dev`` and ``liblapack-dev`` on Ubuntu)

For Windows package names, see the Windows subsection below.

Installation
------------

The software has been tested on Linux Mint (a derivative of Ubuntu) and Red Hat Enterprise Linux. See the following subsection for Windows instructions.

To install, first download the source package and unzip it. Then, activate an existing virtual environment or create a new one and then activate it:

::
   
   python3.X -m venv myenvironment
   source myenvironment/bin/activate

where X is replaced with an appropriate number for the desired version of Python. Then run

::
   
   python3.X -m pip install pathtosourcedir/

Alternatively, if for some reason you do not wish to use a virtual environment while running the software, run

::
   
   python3.X -m pip install --user pathtosourcedir/

The results in the manuscript were produced using Python 3.9, the latest version available on the computing cluster that was used.
   
Windows
^^^^^^^

On Windows, one must use MSYS2, and virtual environments are not supported.

First, download the source package and unzip it.

Then go to https://www.msys2.org/ . Follow the instructions to download and install MSYS2.

Launch an MSYS2 terminal (please use MINGW64).

Enter the following commands:

::
   
   pacman -S mingw-w64-x86_64-python mingw-w64-x86_64-python-pip
   pacman -S mingw-w64-x86_64-ninja mingw-w64-x86_64-cmake
   pacman -S base-devel mingw-w64-x86_64-gcc
   pacman -S mingw-w64-x86_64-blas mingw-w64-x86_64-lapack
   pacman -S mingw-w64-x86_64-python3-numpy mingw-w64-x86_64-pybind11 mingw-w64-x86_64-python3-scipy mingw-w64-x86_64-python3-pandas mingw-w64-x86_64-python3-matplotlib mingw-w64-x86_64-python-scikit-learn

Then run

::
   
   python -m pip install --user pathtosourcedir/


Documentation
-------------

Documentation is accessible at https://ericwayman.net/software/docs/probitlcm/ .

To build the documentation locally, first run ``python3.X -m pip install sphinx furo``. Then from the base directory of the source package, run ``sphinx-build docs/ docs/build/``.

Notes on Armadillo
------------------

Note: the following lines were added to ``CMakeLists.txt`` in ``armadillo``:

::

   option(BUILD_SHARED_LIBS "build shared library" OFF)
   option(CMAKE_POSITION_INDEPENDENT_CODE "this is needed" TRUE)

Note: A further edit to Armadillo files before compiling: according to the section "Programming Environment", subsection "Libraries" of https://campuscluster.illinois.edu/resources/docs/user-guide/, "To use MKL with GCC, consult the Intel MKL link advisor for the link flags to include." This sounds like a bad idea, so we'll disable it as explained in the README file of Armadillo: If MKL is installed and it is persistently giving problems during linking, Support for MKL can be disabled by editing the ``CMakeLists.txt`` file, deleting ``CMakeCache.txt`` and re-running the cmake based installation. Comment out the line containing: ``INCLUDE(ARMA_FindMKL)``.

Current versions of included dependencies
-----------------------------------------
- Armadillo 12.8.4
- CARMA 0.8.0
- nlohmann's JSON library 3.11.3

License
-------

``probitlcm`` is proved under the GNU General Public License v3, a copy of which can be found in the LICENSE file. By using, distributing, or contributing to this project, you agree to the terms and conditions of this license.
