
from scipy import ndimage
import numpy as np


def N3(labels):
    """
    Number of cubes per label
    """
    return ndimage.sum(np.ones_like(labels), labels, np.unique(labels)[1:])

def N2(labels):
    """
    Number of faces
    """
    n_faces = np.zeros(np.unique(labels).shape)

    nx, ny, nz = labels.shape

    dirs = [
        np.array([0, 0, 1]),
        np.array([0, 0, -1]),
        np.array([0, 1, 0]),
        np.array([0, -1, 0]),
        np.array([1, 0, 0]),
        np.array([-1, 0, 0]),
    ]

    def get_item(v, pos):
        return labels[pos[0]%nx, pos[1]%ny, pos[2]%nz]


    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                pos = np.array([i, j, k])
                l_ = labels[tuple(pos)] 
                if l_ == 0:
                    # skip calculating total faces of unlabelled region
                    continue

                for d in dirs:
                    if l_ != get_item(labels, pos + d):
                        n_faces[l_] += 1.0
                    else:
                        n_faces[l_] += 0.5

    return n_faces[1:].astype(int)

def N1(labels):
    """
    Number of edges
    """
    n_edges = np.zeros(np.unique(labels).shape)

    nx, ny, nz = labels.shape

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                l_ = labels[i,j,k] 

                if l_ == 0:
                    # skip calculating total edges of unlabelled region
                    continue

                for d1 in [-1, 1]:
                    for d2 in [-1, 1]:
                        # edges in x-plane
                        jj = (j+d1)%ny
                        kk = (k+d2)%nz
                        x_n = 1.0
                        if labels[i ,jj,k ] == l_:
                            x_n += 1.0
                        if labels[i ,j ,kk] == l_:
                            x_n += 1.0
                        if labels[i ,jj,kk] == l_:
                            x_n += 1.0
                        n_edges[l_] += 1.0/x_n

                        # edges in y-plane
                        ii = (i+d1)%nx
                        kk = (k+d2)%nz
                        x_n = 1.0
                        if labels[ii,j ,k ] == l_:
                            x_n += 1.0
                        if labels[i ,j ,kk] == l_:
                            x_n += 1.0
                        if labels[ii,j ,kk] == l_:
                            x_n += 1.0
                        n_edges[l_] += 1.0/x_n

                        # edges in z-plane
                        ii = (i+d1)%nx
                        jj = (j+d2)%ny
                        x_n = 1.0
                        if labels[ii,j ,k ] == l_:
                            x_n += 1.0
                        if labels[i ,jj,k ] == l_:
                            x_n += 1.0
                        if labels[ii,jj,k ] == l_:
                            x_n += 1.0
                        n_edges[l_] += 1.0/x_n

    return n_edges[1:]

def N0(labels):
    """
    Number of vertices
    """
    n_vertices = np.zeros(np.unique(labels).shape)

    nx, ny, nz = labels.shape

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                pos = np.array([i, j, k])
                l_ = labels[tuple(pos)] 

                if l_ == 0:
                    # skip calculating total edges of unlabelled region
                    continue

                for d1 in [-1, 1]:
                    for d2 in [-1, 1]:
                        for d3 in [-1, 1]:
                            ii = (i+d1)%nx
                            jj = (j+d2)%ny
                            kk = (k+d3)%nz

                            x_n = 1.0
                            if labels[i ,jj,k ] == l_:
                                x_n += 1.0
                            if labels[i ,j ,kk] == l_:
                                x_n += 1.0
                            if labels[i ,jj,kk] == l_:
                                x_n += 1.0
                            if labels[ii,jj,k ] == l_:
                                x_n += 1.0
                            if labels[ii,j ,kk] == l_:
                                x_n += 1.0
                            if labels[ii,jj,kk] == l_:
                                x_n += 1.0
                            if labels[ii,j ,k ] == l_:
                                x_n += 1.0
                            n_vertices[l_] += 1.0/x_n

    # make sure we're close enough to an int, doing differences in
    # seven operations so can get 7 machine epsilon out

    assert np.all(np.abs(v - int(v) < 7*np.finfo(float).eps) for v in n_vertices)

    return n_vertices[1:].astype(int)

def topological_scales(labels, dx):
    """
    Compute characteristic thickness, width, length and genus using Minkowski
    functionals in 3D of all objects labelled objects in `labels`
    """
    n0 = N0(labels).astype(float)
    n1 = N1(labels).astype(float)
    n2 = N2(labels).astype(float)
    n3 = N3(labels).astype(float)

    V0 = n3
    V1 = 2.0*(n2-3.0*n3)/(9.0*dx)
    V2 = 2.0*(n1-2*n2+3*n3)/(9.0*dx*dx)
    V3 = (n0-n1+n2-n3)/(dx*dx*dx)

    thickness = np.abs(V0/(2.0*V1))
    width = np.abs(2.0*V1/(3.14159*V2))
    length = np.abs(dx*dx*dx*3.0*V2/4.0)
    genus = 1.0-0.5*(n0-n1+n2-n3)

    return np.array([thickness, width, length, genus])
