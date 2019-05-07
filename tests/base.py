# coding: utf-8
"""
Routines for creating synthetic datasets and checking the output from the
cloud-identification algorithm
"""

import numpy as np
import xarray as xr
import os
import pytest

def create_circular_mask(grid, x_offset=0.0):
    x, y = grid.x, grid.y
    r = np.sqrt((x - x_offset)**2. + y*y)

    return r < grid.lx/8.

def get_grid():
    lx, ly = 200, 200
    nx, ny = 200, 200
    nz = 1

    x_ = np.linspace(-lx/2., lx/2., nx)
    y_ = np.linspace(-ly/2., ly/2., ny)
    ds = xr.Dataset(coords=dict(x=x_, y=y_))
    ds.attrs['lx'] = lx
    ds.attrs['ly'] = ly

    return ds

class BaseTestClass(object):
    def test_one_circle_x_periodic_scalar_field(self):
        grid = get_grid()
        x, y = grid.x, grid.y

        d = 1. - np.cos(x/10) + 0.*y
        m = create_circular_mask(grid)

        d_out = self.run_classifier(data=d, mask=m)
        print(np.unique(d_out))
        assert len(np.unique(d_out)) == 3


    def test_circle_domain_edge_x_periodic_scalar_field(self):
        grid = get_grid()
        x, y = grid.x, grid.y

        d = 1. - np.cos(x/10) + 0.*y
        m = create_circular_mask(grid, x_offset=grid.lx/2.)

        d_out = self.run_classifier(data=d, mask=m)
        assert len(np.unique(d_out)) == 2


    def test_one_circle_y_periodic_scalar_field(self):
        grid = get_grid()
        x, y = grid.x, grid.y

        d = 1. - np.cos(y/10) + 0.*x
        m = create_circular_mask(grid)

        d_out = self.run_classifier(data=d, mask=m)
        assert len(np.unique(d_out)) == 3

    # XXX: disabled until bug is fixed for no-gradient regions
    def _test_one_circle_no_gradient_scalar_field(self):
        grid = get_grid()
        x, y = grid.x, grid.y

        d = np.ones_like(x) + 0.*y
        m = create_circular_mask(grid)

        d_out = self.run_classifier(data=d, mask=m)
        assert len(np.unique(d_out)) == 2

    # XXX: disabled until bug is fixed for no-gradient regions
    def _test_two_circles_different_scalar_values(self):
        grid = get_grid()
        x, y = grid.x, grid.y

        d = np.ones_like(x) + x > 0.0 + 0.*y

        m1 = create_circular_mask(grid, x_offset=-grid.lx/4.)
        m2 = create_circular_mask(grid, x_offset=grid.lx/4.)
        m = np.logical_or(m1, m2)

        d_out = self.run_classifier(data=d, mask=m)
        assert len(np.unique(d_out)) == 3


    # XXX: disabled until bug is fixed for no-gradient regions
    def _test_two_circles_same_scalar_value(self):
        grid = get_grid()
        x, y = grid.x, grid.y

        d = np.ones_like(x) + x > 0.0 + 0.*y

        m1 = create_circular_mask(grid, x_offset=-grid.lx/4.)
        m2 = create_circular_mask(grid, x_offset=grid.lx/4.)
        m = np.logical_or(m1, m2)

        d_out = self.run_classifier(data=d, mask=m)

        assert len(np.unique(d_out)) == 3

    def test_two_circles_x_periodic_scalar_field(self):
        grid = get_grid()
        x, y = grid.x, grid.y

        d = x + 0.*y

        m1 = create_circular_mask(grid, x_offset=-grid.lx/4.)
        m2 = create_circular_mask(grid, x_offset=grid.lx/4.)
        m = np.logical_or(m1, m2)

        d_out = self.run_classifier(data=d, mask=m)

        num_regions = len(np.unique(d_out))

        assert num_regions == 3

    def test_two_circles_x_periodic_scalar_field(self):
        grid = get_grid()
        x, y = grid.x, grid.y

        d = 1. - np.cos((x+y)/10)

        m1 = create_circular_mask(grid, x_offset=-grid.lx/4.)
        m2 = create_circular_mask(grid, x_offset=grid.lx/4.)
        m = np.logical_or(m1, m2)

        d_out = self.run_classifier(data=d, mask=m)

        num_regions = len(np.unique(d_out))

        assert num_regions == 5

    def run_classifier(data, mask):
        raise NotImplementedError
