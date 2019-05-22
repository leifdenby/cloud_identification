# Cloud numbering software

This code can be used for identifying individual cloud or updraft objects
in 3D and 2D datasets. These regions are first defined using a mask (for
example a minimum constraint on vertical velocity) and are then split by
identifying local maxima in a provided scalar field (for example cloud
liquid water).

![Example of identification of updrafts in 2D from simulation of RICO test
case with UCLALES](doc/updraft-labels-rico-using-local-maxima.png)

![Example of identificatino of clouds in 3D from simulation of RICO test
case with UCLALES](doc/cloud-labels-3d.png)

This is very recent code, and there are no publications with it yet. 

Algorithm:

1. assign each point that fulfils a masking criterion (e.g. each cloudy
   point) to a local maximum. Assignment to maxima is done through
   a steepest gradient approach. 

2. merging algorithm is applied to get rid of the
   smaller local peaks.

*OBS*: The algorithm assumes that the domain is periodic and the outermost
edges in the x, y and z-direction are assumed to be ghost-cells and won't be
taken into account (which means that to e.g. remove objects that intersect with
the edge a mask that one index away from the edge should be used).

# Installation

Whether you'll be using the interface through python or via the command line
you will need the following dependencies:

**1. netcdf**

If this is available through `module` load with

```bash
module load netcdf4
```

which may require you to load hdf5 too

```bash
module load hdf5
```


**2. blitz**

Which can either be provided through scipy (version <=0.18)

```bash
pip install scipy
```

or through `weave`

```bash
pip install weave
```

**3. A compiler with c++11 support**

The python interface for the cloud identification library is written using
[pybind11](https://github.com/pybind/pybind11) which requires c++11
functionality. GCC 4.8 support this, on `cloud9@leeds.ac.uk` this may be loaded
with

```
module load gnu/4.8.1
```

Finally `pybind11` and `weave` are added as git submodules and must be fetched
from github:

```bash
git submodule init
git submodule update
```

Once these are set up you are ready to compile.

For the command line interface simply compile with cmake:

    mkdir build/
    cd build
    cmake ..
    make

The python interface may be compiled and installed with

```bash
python setup.py install
```

## Using from python

Simply import the `cloud_identification` module and use
`cloud_identification.number_objects`.

## Running through command line

    ./build/main

# Input

1. Read `mask_field.nc` which should contain 3D fields `fieldext` and `maskext`
which should be of the `short` datatype.

2. `maskext` is turned into boolean 3D array, from variable `maskshort` to
`maskext`

#

# Development

Test (defined in python in `tests/`) may be run with

    python setup.py test


# Other info

To get ale-vim to pick up the correct libraries when compiling for
error-checking the compiler arguments must be set. These can be found in
`build/compile_commands.json`, eg

    let g:ale_cpp_gcc_options="
    -I/nfs/see-fs-02_users/earlcd/git-repos/cloud_identification/src
    -I/apps/developers/compilers/canopy/1.7.4/1/bit-64/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/scipy/weave/blitz
    -I/nfs/see-fs-02_users/earlcd/git-repos/cloud_identification/lib/pybind11/include
    -I/apps/developers/compilers/canopy/1.7.4/1/bit-64/appdata/canopy-1.7.4.3348.rh5-x86_64/include/python2.7
    -Wall"
