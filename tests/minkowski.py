
from scipy import ndimage
import numpy as np


def N3(m):
    """
    Number of cubes per label
    """
    return ndimage.sum(np.ones_like(m), m, np.unique(m)[1:])

def N2(m):
    """
    Number of faces
    """
    n_faces = np.zeros(np.unique(m).shape)

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
                        n_faces[m_] += 1.0
                    else:
                        n_faces[m_] += 0.5

    return n_faces[1:].astype(int)

def N1(m):
    """
    Number of edges
    """
    n_edges = np.zeros(np.unique(m).shape)

    nx, ny, nz = m.shape

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
                        jj = (j+d1)%ny
                        kk = (k+d2)%nz
                        x_n = 1.0
                        if m[i ,jj,k ] == m_:
                            x_n += 1.0
                        if m[i ,j ,kk] == m_:
                            x_n += 1.0
                        if m[i ,jj,kk] == m_:
                            x_n += 1.0
                        n_edges[m] += 1.0/x_n

                        # edges in y-plane
                        ii = (i+d1)%nx
                        kk = (k+d2)%nz
                        x_n = 1.0
                        if m[ii,j ,k ] == m_:
                            x_n += 1.0
                        if m[i ,j ,kk] == m_:
                            x_n += 1.0
                        if m[ii,j ,kk] == m_:
                            x_n += 1.0
                        n_edges[m] += 1.0/x_n

                        # edges in z-plane
                        ii = (i+d1)%nx
                        jj = (j+d2)%ny
                        x_n = 1.0
                        if m[ii,j ,k ] == m_:
                            x_n += 1.0
                        if m[i ,jj,k ] == m_:
                            x_n += 1.0
                        if m[ii,jj,k ] == m_:
                            x_n += 1.0
                        n_edges[m] += 1.0/x_n

    return n_edges[1:]

def N0(m):
    """
    Number of vertices
    """
    n_vertices = np.zeros(np.unique(m).shape)

    nx, ny, nz = m.shape

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
                        for d3 in [-1, 1]:
                            ii = (i+d1)%nx
                            jj = (j+d2)%ny
                            kk = (k+d3)%nz

                            x_n = 1.0
                            if m[i ,jj,k ] == m_:
                                x_n += 1.0
                            if m[i ,j ,kk] == m_:
                                x_n += 1.0
                            if m[i ,jj,kk] == m_:
                                x_n += 1.0
                            if m[ii,jj,k ] == m_:
                                x_n += 1.0
                            if m[ii,j ,kk] == m_:
                                x_n += 1.0
                            if m[ii,jj,kk] == m_:
                                x_n += 1.0
                            if m[ii,j ,k ] == m_:
                                x_n += 1.0
                            n_vertices[m] += 1.0/x_n

    # make sure we're close enough to an int, doing differences in
    # seven operations so can get 7 machine epsilon out

    assert np.all(np.abs(v - int(v) < 7*np.finfo(float).eps) for v in n_vertices)

    return n_vertices[1:].astype(int)
