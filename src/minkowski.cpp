#include <blitz/array.h>

namespace minkowski {

  /** Compute number of cells for each labelled object (i.e. the volume)
  */
  blitz::Array<int,1> N3(blitz::Array<int,3> &labels) {
    int label_max = max(labels);

    blitz::TinyVector<int,3> shape = labels.shape();

    int nx = shape[0];
    int ny = shape[1];
    int nz = shape[2];

    blitz::Array<int,1> n_cells = blitz::Array<int,1>(label_max+1);
    n_cells = 0;

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          n_cells(labels(i,j,k)) += 1;
        }
      }
    }

    return n_cells(blitz::Range(1, blitz::toEnd));
  }


  /** Compute number of faces for each labelled object
  */
  blitz::Array<int,1> N2(blitz::Array<int,3> &labels) {
    int label_max = max(labels);

    blitz::TinyVector<int,3> shape = labels.shape();

    int nx = shape[0];
    int ny = shape[1];
    int nz = shape[2];

    blitz::Array<double,1> n_faces_fractions = blitz::Array<double,1>(label_max+1);

    int l_ = -1;

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          l_ = labels(i,j,k);

          if (l_ == 0) {
            //n_vertices_fractions(l_) = 0.0;
            // skip calculating total vertices of unlabelled region
          }
          else {
            n_faces_fractions(l_) = n_faces_fractions(l_)
                      + ( labels(i  ,j+1,k  ) == l_ ? 0.5 : 1.0)
                      + ( labels(i  ,j-1,k  ) == l_ ? 0.5 : 1.0)
                      + ( labels(i+1,j  ,k  ) == l_ ? 0.5 : 1.0)
                      + ( labels(i-1,j  ,k  ) == l_ ? 0.5 : 1.0)
                      + ( labels(i  ,j  ,k+1) == l_ ? 0.5 : 1.0)
                      + ( labels(i  ,j  ,k-1) == l_ ? 0.5 : 1.0);
          }
        }
      }
    }

    blitz::Array<int,1> n_faces = blitz::Array<int,1>(label_max+1);
    for (int i=0; i<n_faces.size(); i++) {
      n_faces(i) = (int)n_faces_fractions(i);
    }

    return n_faces(blitz::Range(1, blitz::toEnd));
  }


  /** Compute number of edges for each labelled object
  */
  blitz::Array<int,1> N1(blitz::Array<int,3> &labels) {
    int label_max = max(labels);

    blitz::TinyVector<int,3> shape = labels.shape();

    int nx = shape[0];
    int ny = shape[1];
    int nz = shape[2];

    blitz::Array<double,1> n_edges_fractions = blitz::Array<double,1>(label_max+1);

    int ii = -1, jj = -1, kk = -1;
    int l_ = -1;
    double x_n = -1.0;

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          l_ = labels(i,j,k);

          if (l_ == 0) {
            //n_edges_fractions(l_) = 0.0;
            // skip calculating total vertices of unlabelled region
          }
          else {
            for (int d1=-1; d1<2; d1+=2) {
              for (int d2=-1; d2<2; d2+=2) {
                // edges in x-plane
                jj = (j+d1)%ny;
                kk = (k+d2)%nz;
                x_n = 1.0;
                if (labels(i ,jj,k ) == l_) { x_n += 1.0; }
                if (labels(i ,j ,kk) == l_) { x_n += 1.0; }
                if (labels(i ,jj,kk) == l_) { x_n += 1.0; }
                n_edges_fractions(l_) += 1.0/x_n;

                // edges in y-plane
                ii = (i+d1)%nx;
                kk = (k+d2)%nz;
                x_n = 1.0;
                if (labels(ii,j ,k ) == l_) { x_n += 1.0; }
                if (labels(i ,j ,kk) == l_) { x_n += 1.0; }
                if (labels(ii,j ,kk) == l_) { x_n += 1.0; }
                n_edges_fractions(l_) += 1.0/x_n;

                // edges in z-plane
                ii = (i+d1)%nx;
                jj = (j+d2)%ny;
                x_n = 1.0;
                if (labels(ii,j ,k ) == l_) { x_n += 1.0; }
                if (labels(i ,jj,k ) == l_) { x_n += 1.0; }
                if (labels(ii,jj,k ) == l_) { x_n += 1.0; }
                n_edges_fractions(l_) += 1.0/x_n;
              }
            }
          }
        }
      }
    }

    blitz::Array<int,1> n_edges = blitz::Array<int,1>(label_max+1);
    for (int i=0; i<n_edges.size(); i++) {
      n_edges(i) = (int)n_edges_fractions(i);
    }

    return n_edges(blitz::Range(1, blitz::toEnd));
  }


  /** Compute number of vertices for each labelled object
  */
  blitz::Array<int,1> N0(blitz::Array<int,3> &labels) {
    int label_max = max(labels);

    blitz::TinyVector<int,3> shape = labels.shape();

    int nx = shape[0];
    int ny = shape[1];
    int nz = shape[2];

    blitz::Array<double,1> n_vertices_fractions = blitz::Array<double,1>(label_max+1);

    int ii = -1, jj = -1, kk = -1;
    int l_ = -1;
    double x_n = -1.0;

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          l_ = labels(i,j,k);

          if (l_ == 0) {
            //n_vertices_fractions(l_) = 0.0;
            // skip calculating total vertices of unlabelled region
          }
          else {
            for (int d1=-1; d1<2; d1+=2) {
              for (int d2=-1; d2<2; d2+=2) {
                for (int d3=-1; d3<2; d3+=2) {
                  ii = (i+nx+d1)%nx;
                  jj = (j+ny+d2)%ny;
                  kk = (k+nz+d3)%nz;

                  x_n = 1.0;
                  if (labels(i ,jj,k ) == l_){ x_n += 1.0; }
                  if (labels(i ,j ,kk) == l_){ x_n += 1.0; }
                  if (labels(i ,jj,kk) == l_){ x_n += 1.0; }
                  if (labels(ii,jj,k ) == l_){ x_n += 1.0; }
                  if (labels(ii,j ,kk) == l_){ x_n += 1.0; }
                  if (labels(ii,jj,kk) == l_){ x_n += 1.0; }
                  if (labels(ii,j ,k ) == l_){ x_n += 1.0; }

                  n_vertices_fractions(l_) += 1.0/x_n;
                }
              }
            }
          }
        }
      }
    }

    //# make sure we're close enough to an int, doing differences in
    //# seven operations so can get 7 machine epsilon out

    //assert np.all(np.abs(v - int(v) < 7*np.finfo(float).eps) for v in n_vertices)
    
    blitz::Array<int,1> n_vertices = blitz::Array<int,1>(label_max+1);
    for (int i=0; i<n_vertices.size(); i++) {
      n_vertices(i) = (int)n_vertices_fractions(i);
    }

    return n_vertices(blitz::Range(1, blitz::toEnd));
  }


  //def N2(m):
  //"""
  //Number of faces
  //"""
  //n_faces = np.zeros_like(np.unique(m))

  //nx, ny, nz = m.shape

  //dirs = [
  //np.array([0, 0, 1]),
  //np.array([0, 0, -1]),
  //np.array([0, 1, 0]),
  //np.array([0, -1, 0]),
  //np.array([1, 0, 0]),
  //np.array([-1, 0, 0]),
  //]

  //def get_item(v, pos):
  //return m[pos[0]%nx, pos[1]%ny, pos[2]%nz]


  //for i in range(nx):
  //for j in range(ny):
  //for k in range(nz):
  //pos = np.array([i, j, k])
  //m_ = m[tuple(pos)] 
  //if m_ == 0:
  //# skip calculating total faces of unlabelled region
  //continue

  //for d in dirs:
  //if m_ != get_item(m, pos + d):
  //n_faces[m_] += 1

  //return n_faces[1:]


  //def N1(m):
  //"""
  //Number of edges
  //"""
  //n_edges = np.zeros(np.unique(m).shape)

  //nx, ny, nz = m.shape

  //for i in range(nx):
  //for j in range(ny):
  //for k in range(nz):
  //pos = np.array([i, j, k])
  //m_ = m[tuple(pos)] 

  //if m_ == 0:
  //# skip calculating total edges of unlabelled region
  //continue

  //for d1 in [-1, 1]:
  //for d2 in [-1, 1]:
  //# edges in x-plane
  //jj = (j+d1)%ny
  //kk = (k+d2)%nz
  //x_n = 1.0
  //if m[i ,jj,k ] == m_:
  //x_n += 1.0
  //if m[i ,j ,kk] == m_:
  //x_n += 1.0
  //if m[i ,jj,kk] == m_:
  //x_n += 1.0
  //n_edges[m] += 1.0/x_n

  //# edges in y-plane
  //ii = (i+d1)%nx
  //kk = (k+d2)%nz
  //x_n = 1.0
  //if m[ii,j ,k ] == m_:
  //x_n += 1.0
  //if m[i ,j ,kk] == m_:
  //x_n += 1.0
  //if m[ii,j ,kk] == m_:
  //x_n += 1.0
  //n_edges[m] += 1.0/x_n

  //# edges in z-plane
  //ii = (i+d1)%nx
  //jj = (j+d2)%ny
  //x_n = 1.0
  //if m[ii,j ,k ] == m_:
  //x_n += 1.0
  //if m[i ,jj,k ] == m_:
  //x_n += 1.0
  //if m[ii,jj,k ] == m_:
  //x_n += 1.0
  //n_edges[m] += 1.0/x_n

  //return n_edges[1:]

  //def N0(m):

}
