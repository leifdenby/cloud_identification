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

# Installation

Whether you'll be using the interface through python or via the command line
you will need the following dependencies:

**1. netcdf**

If this is available through `module` load with

```bash
module load netcdf4
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
