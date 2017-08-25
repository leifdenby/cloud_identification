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

import py_cloud_identification
from base import BaseTestClass


class TestPythonInterface(BaseTestClass):

    def run_classifier(self, data, mask):
        assert data.shape == mask.shape

        if len(data.shape) == 2:
            d = np.expand_dims(data, axis=-1)
            m = np.expand_dims(mask, axis=-1)
            assert d.base is data
        else:
            d = data
            m = mask

        return py_cloud_identification.number_objects(scalar_field=d, mask=m)
