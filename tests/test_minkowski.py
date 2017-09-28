import StringIO
import warnings

import numpy as np

import minkowski
import cloud_identification

def parse_data_str(s, shape):
    return np.loadtxt(StringIO.StringIO(s)).reshape(*shape).astype(int)

FUNCS = "N3 N2 N1 N0"

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
        assert np.all(np.abs(scales_py - scales_cpp) < 10*np.finfo(scales_cpp.dtype).eps)

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

def test_sphere():
    l = 100.
    N = 20
    dx = l/N

    x_ = np.linspace(-l/2., l/2., N)
    y_ = np.linspace(-l/2., l/2., N)
    z_ = np.linspace(-l/2., l/2., N)

    x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

    r = np.sqrt(x**2. + y**2. + z**2.)

    mask = r < 0.3*r.max()

    labels = cloud_identification.number_objects(mask=mask, scalar_field=r)

    assert cloud_identification.N0(labels) == minkowski.N0(labels)
    assert cloud_identification.N1(labels) == minkowski.N1(labels)
    assert cloud_identification.N2(labels) == minkowski.N2(labels)
    assert cloud_identification.N3(labels) == minkowski.N3(labels)


    assert len(np.unique(labels)) == 2

    scales_cpp = cloud_identification.topological_scales(labels, dx=dx)
    scales_py = minkowski.topological_scales(labels, dx=dx)

    assert np.all(np.abs(scales_py - scales_cpp) < 1.0e2*np.finfo(scales_cpp.dtype).eps)

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
