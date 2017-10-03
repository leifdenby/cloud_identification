# coding: utf-8
"""
Routines for testing CLI interface where netcdf files are used
"""
from scipy.io import netcdf_file
import netCDF4
import subprocess
import numpy as np
import os
import warnings

from base import BaseTestClass

DO_CLEANUP = True

EXECUTABLE_PATH = os.path.join(os.environ.get("BUILD_TEMP_DIR", "build"), "main")
if not os.path.exists(EXECUTABLE_PATH):
    warnings.warn("Couldn't find main executable, skipping tests for it")
else:
    class TestNCFile(BaseTestClass):

        def run_classifier(self, data, mask):
            def save_input():
                nx, ny = data.shape
                nz = 1

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

            proc = subprocess.Popen(EXECUTABLE_PATH)
            proc.communicate()
            delete_input()

            if not proc.returncode == 0:
                raise Exception("classification program crashed, return code: {}".format(proc.returncode))

            return read_output()

