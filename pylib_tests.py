# coding: utf-8
"""
Routines for creating synthetic datasets and checking the output from the
cloud-identification algorithm
"""

# add build path to PYTHONPATH so we can import built module
import os
import sys
build_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build')
sys.path.append(build_path)


import numpy as np
import pytest

import py_cloud_identification

lx, ly = 200, 200
nx, ny = 200, 200
nz = 1

DO_CLEANUP = True

def create_circular_mask(x, y, x_offset=0.0):
    return np.sqrt((x - x_offset)**2. + y*y) < 25


def get_grid():
    x_ = np.linspace(-lx/2., lx/2., nx)
    y_ = np.linspace(-ly/2., ly/2., ny)
    x, y = np.meshgrid(x_, y_, indexing='ij')

    return x, y


def test_one_circle_x_periodic_scalar_field():
    x, y = get_grid()

    d = 1. - np.cos(x/10)
    m = create_circular_mask(x, y)

    d_out = run_classifier(data=d, mask=m)
    print np.unique(d_out)
    assert len(np.unique(d_out)) == 3

def run_classifier(data, mask):
    assert data.shape == mask.shape

    if len(data.shape) == 2:
        d = np.expand_dims(data, axis=-1)
        m = np.expand_dims(mask, axis=-1)
        assert d.base is data
    else:
        d = data
        m = mask

    return py_cloud_identification.number_objects(scalar_field=d, mask=m)
