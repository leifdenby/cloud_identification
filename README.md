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

# Compile on MONSOON

Compile on MONSOON

    g++ -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -fPIC
    -I/opt/python/gnu/2.7.9/lib/python2.7/site-packages/scipy/weave/blitz
    -I/opt/cray/netcdf-hdf5parallel/4.3.2/CRAY/83/include
    -I/opt/cray/hdf5/1.8.13/CRAY/83/include -c cpp_identify_seungbu.cpp -o
    cpp_identify_seungbu.o -pthread -O6 -march=native -mtune=native
    -funroll-all-loops -fomit-frame-pointer -march=native -mtune=native -msse4
    -ftree-vectorize -ftree-vectorizer-verbose=5 -ffast-math -funroll-loops
    -ftracer
    
    g++ -o cpp_identify_seungbu.exe cpp_identify_seungbu.o
    -L/opt/cray/netcdf-hdf5parallel/4.3.2/CRAY/83/lib -lnetcdf

# Compile on my local system: 

    g++ -pthread -fno-strict-aliasing -DNDEBUG -g -fwrapv -fPIC
    -I/usr/lib/python2.7/dist-packages/scipy/weave/blitz -I/usr/include/ -c
    cpp_identify_seungbu.cpp -o cpp_identify_seungbu.o -pthread -O6 -march=native
    -mtune=native -funroll-all-loops -fomit-frame-pointer -march=native
    -mtune=native -msse4 -ftree-vectorize -ftree-vectorizer-verbose=5 -ffast-math
    -funroll-loops -ftracer 
    
    g++ -o cpp_identify_seungbu.exe cpp_identify_seungbu.o -L/usr/lib/ -lnetcdf
