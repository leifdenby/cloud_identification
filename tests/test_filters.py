"""
Contains test for routines which perform manipulations on 3D arrays of objects
labels
"""
import numpy as np

import base
import cloud_identification

def _get_twocircle_dataset():
    grid = base.get_grid()
    x, y = grid.x, grid.y

    mask1 = base.create_circular_mask(grid)
    x_offset = grid.lx/2.
    mask2 = base.create_circular_mask(grid, x_offset=x_offset)

    mask = np.logical_or(mask1, mask2)

    d1 = np.sqrt(x*x + y*y)
    d2 = np.sqrt((x - x_offset)**2. + y*y)
    d = (d1.max() - d1) + (d2.max() - d2)

    d = np.expand_dims(d, axis=-1)
    m = np.expand_dims(mask, axis=-1)

    return d, m

def test_intersection_filter_no_overlap():
    """
    Check that labels are unchanged if there is no intersection with mask
    provided for filter
    """
    d, m = _get_twocircle_dataset()

    labels = cloud_identification.number_objects(d, m)
    print(labels[0:2,:,0])
    assert len(np.unique(labels)) == 3

    mask_overlap = np.zeros_like(m).astype(bool)
    cloud_identification.remove_intersecting(labels, mask_overlap)
    assert len(np.unique(labels)) == 3

def test_intersection_filter_complete_overlap():
    """
    Check that all labels are removed if there is intersection with mask
    provided for filter
    """
    d, m = _get_twocircle_dataset()

    labels = cloud_identification.number_objects(d, m)
    print(labels[0:2,:,0])
    assert len(np.unique(labels)) == 3

    mask_overlap = np.ones_like(m).astype(bool)
    cloud_identification.remove_intersecting(labels, mask_overlap)
    assert len(np.unique(labels)) == 1


def test_intersection_filter_edge_overlap():
    """
    Check that labels are removed if they overlap with edge
    """
    d, m = _get_twocircle_dataset()

    labels = cloud_identification.number_objects(d, m)
    assert len(np.unique(labels)) == 3

    mask_overlap = np.zeros_like(m).astype(bool)
    mask_overlap[1,:] = True
    mask_overlap[-2,:] = True
    mask_overlap[:,1] = True
    mask_overlap[:,-2] = True

    cloud_identification.remove_intersecting(labels, mask_overlap)
    assert len(np.unique(labels)) == 2


def test_intersection_filter_no_overlap():
    """
    Check that labels are correctly removed when only aiming to keep one object
    """
    d, m = _get_twocircle_dataset()

    labels = cloud_identification.number_objects(d, m)
    assert len(np.unique(labels)) == 3

    idxs_keep = np.unique(labels)
    cloud_identification.filter_labels(labels, idxs_keep)
    assert len(np.unique(labels)) == 3

    idxs_keep = np.unique(labels)
    idxs_keep = idxs_keep[idxs_keep != 0]
    cloud_identification.filter_labels(labels, idxs_keep)
    # label == 0 is "outside object" so the number of unique labels should be
    # unchanged
    assert len(np.unique(labels)) == 3

    idxs_keep = np.unique(labels)
    idxs_keep = idxs_keep[idxs_keep != 1]
    cloud_identification.filter_labels(labels, idxs_keep)
    # now we should have one label (2) and environment
    assert np.array_equal(np.unique(labels), [0, 2])
