# Cloud numbering software

The code to do cloud numbering, it is rather procedural (=faster?). Works on
both MONC and UM data given the right parameters.

This is very recent code, and there are no publications with it yet. 

Algorithm:

- first assign each point that fulfils a masking criterion (e.g. each cloudy
point) to a local maximum.

- assignment to maxima is done through a steepest gradient approach. 

- subsequently, a merging algorithm is applied to get rid of the smaller local
peaks.

# Input

1. Read `mask_field.nc` which should contain 3D fields `fieldext` and `maskext`
which should be of the `short` datatype.

2. `maskext` is turned into boolean 3D array, from variable `maskshort` to
`maskext`

# Installation

Compile is handled with cmake. The required dependencies are:

**1. netcdf**

If this is available through `module` load with

```bash
module load netcdf
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

Once these are installed compile with cmake:

    mkdir build/
    cd build
    cmake ..
    make
