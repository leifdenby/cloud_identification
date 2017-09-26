
from scipy import ndimage
import numpy as np

import minkowski_test_data


def test_one():
    test_data = minkowski_test_data.DATA1
    m = test_data.get('mask')

    assert N3(m) == test_data.get('N3')
    assert N2(m) == test_data.get('N2')

def test_two():
    test_data = minkowski_test_data.DATA2
    m = test_data.get('mask')

    assert N3(m) == test_data.get('N3')
    assert N2(m) == test_data.get('N2')

def test_three():
    test_data = minkowski_test_data.DATA3
    m = test_data.get('mask')

    assert N3(m) == test_data.get('N3')
    assert N2(m) == test_data.get('N2')

def test_four():
    test_data = minkowski_test_data.DATA4
    m = test_data.get('mask')

    assert N3(m) == test_data.get('N3')
    assert N2(m) == test_data.get('N2')
    assert N1(m) == test_data.get('N1')

def N3(m):
    """
    Number of cubes per label
    """
    return ndimage.sum(np.ones_like(m), m, np.unique(m)[1:])

def N2(m):
    """
    Number of faces
    """
    n_faces = np.zeros_like(np.unique(m))

    nx, ny, nz = m.shape

    dirs = [
        np.array([0, 0, 1]),
        np.array([0, 0, -1]),
        np.array([0, 1, 0]),
        np.array([0, -1, 0]),
        np.array([1, 0, 0]),
        np.array([-1, 0, 0]),
    ]

    def get_item(v, pos):
        return m[pos[0]%nx, pos[1]%ny, pos[2]%nz]


    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                pos = np.array([i, j, k])
                m_ = m[tuple(pos)] 
                if m_ == 0:
                    # skip calculating total faces of unlabelled region
                    continue

                for d in dirs:
                    if m_ != get_item(m, pos + d):
                        n_faces[m_] += 1

    return n_faces[1:]


def N1(m):
    """
    Number of edges
    """
    n_edges = np.zeros_like(np.unique(m))

    nx, ny, nz = m.shape

    dirs = [
        np.array([0, 0, 1]),
        np.array([0, 0, -1]),
        np.array([0, 1, 0]),
        np.array([0, -1, 0]),
        np.array([1, 0, 0]),
        np.array([-1, 0, 0]),
    ]

    def get_item(v, pos):
        return m[pos[0]%nx, pos[1]%ny, pos[2]%nz]




    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                pos = np.array([i, j, k])
                m_ = m[tuple(pos)] 

                if m_ == 0:
                    # skip calculating total edges of unlabelled region
                    continue

                for d1 in [-1, 1]:
                    for d2 in [-1, 1]:
                        # edges in x-plane
                        if m


                for d in dirs:
                    if m_ != get_item(m, pos + d):
                        n_faces[m_] += 1

    return n_faces[1:]

if __name__ == "__main__":
    from mayavi import mlab

    def make_grid(v, dx):
        nx, ny, nz = v.shape
        x_ = np.arange(nx)*dx
        y_ = np.arange(ny)*dx
        z_ = np.arange(nz)*dx
        
        return np.meshgrid(x_, y_, z_, indexing='ij')

    def plot_labels(cloud_labels, dx):
        m_c = cloud_labels != 0
        fig = mlab.figure()
        x, y, z = make_grid(cloud_labels, dx=dx)
        pts = mlab.points3d(x[m_c], y[m_c], z[m_c], cloud_labels[m_c], mode='cube', scale_factor=25., scale_mode='none')
        pts.glyph.glyph.clamping = False
        return fig

    test_data = minkowski_test_data.DATA3
    m = test_data.get('mask')
    plot_labels(m, 25.)
    mlab.show()
