import StringIO
import warnings

import numpy as np

import minkowski
import cloud_identification

def parse_data_str(s, shape):
    return np.loadtxt(StringIO.StringIO(s)).reshape(*shape).astype(int)

FUNCS = "N3 N2 N1 N0"

def _assert_equal_floats(v1, v2):
    assert v1.dtype == v2.dtype
    eps = np.finfo(v1.dtype).eps
    max_error = 1000*eps
    print eps, max_error
    print v1 - v2
    print v1
    print v2
    if np.max(np.abs(v1 - v2) < max_error) > 100*eps:
        warnings.warn("Larger than 100 epsilon descrepancy between c++ and python result, should investigate")
    assert np.all(np.abs(v1 - v2) < max_error)

def _check_func_values(test_data):
    m = test_data.get('mask')

    for func_name in FUNCS.split():
        if not func_name in test_data:
            warnings.warn("`{}` not defined for test".format(func_name))
        else:
            func = getattr(minkowski, func_name)
            assert np.all(np.array(func(m)) == np.array(test_data.get(func_name)))

        cpp_func = getattr(cloud_identification, func_name, None)

        if cpp_func is None:
            warnings.warn("Function `{}` not defined in Minkowski cpp interface".format(func_name))
        else:
            assert np.all(cpp_func(m) == np.array(test_data.get(func_name)))

        scales_fn_py = minkowski.topological_scales
        scales_fn_cpp = cloud_identification.topological_scales

        dx = 25.0
        scales_py = scales_fn_py(m, dx)
        scales_cpp = scales_fn_cpp(m, dx)
        _assert_equal_floats(scales_py, scales_cpp)

def test_cube():
    test_data = dict(
        mask = parse_data_str("""
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 0 0 0
        0 1 0 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        """, shape=(3, 5, 4)),
        N3=1,
        N2=6,
        N1=12,
        N0=8,
    )
    _check_func_values(test_data)

def test_two_cubes():
    test_data = dict(
        mask = parse_data_str("""
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 0 0 0
        0 1 2 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        """, shape=(3, 5, 4)),
        N3=[1,1],
        N2=[6,6],
        N1=[12,12],
        N0=[8,8],
    )
    _check_func_values(test_data)

def test_three_cubes():
    test_data = dict(
        mask = parse_data_str("""
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 0 3 0
        0 1 2 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        """, shape=(3, 5, 4)),
        N3=[1,1,1],
        N2=[6,6,6],
        N1=[12,12,12],
        N0=[8,8,8],
    )
    _check_func_values(test_data)

def test_rectangle():
    test_data = dict(
        mask = parse_data_str("""
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 0 0 0
        0 1 1 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        """, shape=(3, 5, 4)),
        N3=2,
        N2=2*6-1,
        N1=7*2+6,
        N0=12,
    )
    _check_func_values(test_data)

def test_block():
    test_data = dict(
        mask = parse_data_str("""
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 1 1 0
        0 1 1 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        """, shape=(3, 5, 4)),
        N3=4,
        N2=4*6-4,
        N1=2*12+3*3,
        N0=18,
    )
    _check_func_values(test_data)

def test_t_shape():
    test_data = dict(
        mask = parse_data_str("""
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0

        0 0 0 0
        0 1 0 0
        0 1 1 0
        0 1 0 0
        0 0 0 0

        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        0 0 0 0
        """, shape=(3, 5, 4)),
        N3=4,
        N2=4*6-3,
        N1=13*2+10,
        N0=4*8-3*4,
    )
    _check_func_values(test_data)

def make_spheroid_object(r0, dx):
    if type(r0) == float:
        lx = ly = lz = 4*r0
        Nx = Ny = Nz = int(lx/dx)
        r0x = r0y = r0z = r0
    else:
        if not type(r0) == np.ndarray:
            r0 = np.array(r0)
        r0x, r0y, r0z = r0
        l = 4*r0
        lx, ly, lz = l
        Nx, Ny, Nz = (l/dx).astype(int)

    x_ = np.linspace(-lx/2., lx/2., Nx)
    y_ = np.linspace(-ly/2., ly/2., Ny)
    z_ = np.linspace(-lz/2., lz/2., Nz)

    x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

    r = np.sqrt((x/r0x)**2. + (y/r0y)**2. + (z/r0z)**2.)

    mask = r < 1.0

    return mask, r

def test_sphere():
    dx = 25.
    mask, s = make_spheroid_object(r0=100., dx=dx)

    labels = cloud_identification.number_objects(mask=mask, scalar_field=s)

    assert cloud_identification.N0(labels) == minkowski.N0(labels)
    assert cloud_identification.N1(labels) == minkowski.N1(labels)
    assert cloud_identification.N2(labels) == minkowski.N2(labels)
    assert cloud_identification.N3(labels) == minkowski.N3(labels)


    assert len(np.unique(labels)) == 2

    scales_cpp = cloud_identification.topological_scales(labels, dx=dx)
    scales_py = minkowski.topological_scales(labels, dx=dx).astype(scales_cpp.dtype)

    _assert_equal_floats(scales_py, scales_cpp)

def test_ellipsoid():
    dx = 25.
    mask, s = make_spheroid_object(r0=[100., 400., 400.], dx=dx)

    labels = cloud_identification.number_objects(mask=mask, scalar_field=s)

    assert cloud_identification.N0(labels) == minkowski.N0(labels)
    assert cloud_identification.N1(labels) == minkowski.N1(labels)
    assert cloud_identification.N2(labels) == minkowski.N2(labels)
    assert cloud_identification.N3(labels) == minkowski.N3(labels)


    assert len(np.unique(labels)) == 2

    scales_cpp = cloud_identification.topological_scales(labels, dx=dx)
    scales_py = minkowski.topological_scales(labels, dx=dx).astype(scales_cpp.dtype)

    _assert_equal_floats(scales_py, scales_cpp)

if __name__ == "__main__":
    l = 100.
    N = 300

    dx = l/N

    x_ = np.linspace(-l/2., l/2., N)
    y_ = np.linspace(-l/2., l/2., N)
    z_ = np.linspace(-l/2., l/2., N)

    x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

    r = np.sqrt(x**2. + y**2. + z**2.)

    lr = 0.3*r.max()

    mask = r < lr

    labels = cloud_identification.number_objects(mask=mask, scalar_field=r)
    scales_cpp = cloud_identification.topological_scales(labels, dx=dx)

    print lr
    print scales_cpp
