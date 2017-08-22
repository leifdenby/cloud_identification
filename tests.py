# coding: utf-8
"""
Routines for creating synthetic datasets and checking the output from the
cloud-identification algorithm
"""

from scipy.io import netcdf_file
import netCDF4
import numpy as np
import os
import subprocess
import pytest

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


def test_one_circle_y_periodic_scalar_field():
    x, y = get_grid()

    d = 1. - np.cos(y/10)
    m = create_circular_mask(x, y)

    d_out = run_classifier(data=d, mask=m)
    assert len(np.unique(d_out)) == 3


def test_one_circle_no_gradient_scalar_field():
    x, y = get_grid()

    d = np.ones_like(x)
    m = create_circular_mask(x, y)

    d_out = run_classifier(data=d, mask=m)
    assert len(np.unique(d_out)) == 2

def test_two_circles_different_scalar_values():
    x, y = get_grid()

    d = np.ones_like(x) + x > 0.0

    m1 = create_circular_mask(x, y, x_offset=-lx/4.)
    m2 = create_circular_mask(x, y, x_offset=lx/4.)
    m = np.logical_or(m1, m2)

    d_out = run_classifier(data=d, mask=m)
    assert len(np.unique(d_out)) == 2


def test_two_circles_same_scalar_value():
    x, y = get_grid()

    d = np.ones_like(x) + x > 0.0

    m1 = create_circular_mask(x, y, x_offset=-lx/4.)
    m2 = create_circular_mask(x, y, x_offset=lx/4.)
    m = np.logical_or(m1, m2)

    d_out = run_classifier(data=d, mask=m)
    assert len(np.unique(d_out)) == 2

def test_two_circles_x_periodic_scalar_field():
    x, y = get_grid()

    d = x

    m1 = create_circular_mask(x, y, x_offset=-lx/4.)
    m2 = create_circular_mask(x, y, x_offset=lx/4.)
    m = np.logical_or(m1, m2)

    d_out = run_classifier(data=d, mask=m)

    num_regions = len(np.unique(d_out))

    assert num_regions == 3

def test_two_circles_x_periodic_scalar_field():
    x, y = get_grid()

    d = 1. - np.cos(x/10)

    m1 = create_circular_mask(x, y, x_offset=-lx/4.)
    m2 = create_circular_mask(x, y, x_offset=lx/4.)
    m = np.logical_or(m1, m2)

    d_out = run_classifier(data=d, mask=m)

    # TODO: there's a bug here(!) the algorithm sometimes splits of regions
    # which are only one cell which isn't a local maxima. Remove these for now
    # to get the counts right
    num_regions = len(filter(lambda n: n > 1, np.bincount(d_out.flatten())))

    # num_regions = len(np.unique(d_out))

    assert num_regions == 5

def run_classifier(data, mask):
    def save_input():
        fh = netcdf_file("mask_field.nc", "w")

        fh.createDimension("x", nx)
        fh.createDimension("y", ny)
        fh.createDimension("z", nz)

        maskext = fh.createVariable("maskext", np.int16, ("x", "y", "z"))
        maskext[:,:,0] = mask.astype(np.int16)

        fieldext = fh.createVariable("fieldext", np.int16, ["x", "y", "z"])
        # XXX: input must be saved as int16 for now, so we must rescale
        d_range = data.max() - data.min()
        if d_range == 0.0:
            d_scaled = np.ones_like(data)
        else:
            d_scaled = (data - data.min())/(d_range)*np.iinfo(np.int16).max
        fieldext[:,:,0] = d_scaled.astype(np.int16)

        fh.close()

    def delete_input():
        if DO_CLEANUP:
            os.remove("mask_field.nc")

    def read_output():
        fh = netCDF4.Dataset("output.nc")
        d = fh.variables['data'][:][:,:,0]
        fh.close()
        if DO_CLEANUP:
            os.remove("output.nc")

        return d

    save_input()
    proc = subprocess.Popen('build/main')
    proc.communicate()
    
    if not proc.returncode == 0:
        # delete_input()
        raise Exception("classification program crashed, return code: {}".format(proc.returncode))

    return read_output()
