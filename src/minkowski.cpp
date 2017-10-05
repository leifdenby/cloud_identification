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
    n_faces_fractions = 0;

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
      n_faces(i) = (int)round(n_faces_fractions(i));
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
    n_edges_fractions = 0;

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
      n_edges(i) = (int)round(n_edges_fractions(i));
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
    n_vertices_fractions = 0;

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
      n_vertices(i) = (int)round(n_vertices_fractions(i));
    }

    return n_vertices(blitz::Range(1, blitz::toEnd));
  }



  /** Compute characteristic thickness, width, length and genus for each
   * labelled object in 3D using Minkowski functionals
   */
  blitz::Array<float,2> topological_scales(blitz::Array<int,3> &labels, float dx) {
    int N_labels = max(labels);

    blitz::Array<float,1> n0(N_labels), n1(N_labels), n2(N_labels), n3(N_labels);
    n0 = blitz::cast<float>(N0(labels));
    n1 = blitz::cast<float>(N1(labels));
    n2 = blitz::cast<float>(N2(labels));
    n3 = blitz::cast<float>(N3(labels));

    blitz::Array<float,1> V0(N_labels), V1(N_labels), V2(N_labels), V3(N_labels);
    V0 = n3;
    V1 = 2.0*(n2-3.0*n3)/(9.0*dx);
    V2 = 2.0*(n1-2*n2+3*n3)/(9.0*dx*dx);
    V3 = (n0-n1+n2-n3)/(dx*dx*dx);

    blitz::Array<float,1> thickness(N_labels), width(N_labels), length(N_labels), genus(N_labels);

    thickness = abs(V0/(2.0*V1));
    width =  abs(2.0*V1/(3.14159*V2));
    length =  abs(3.0*V2/(4.0*V3));
    genus = 1.0-0.5*(n0-n1+n2-n3);

    blitz::Array<float,2> topological_scales = blitz::Array<float,2>(4, N_labels);
    topological_scales(0, blitz::Range::all()) = thickness;
    topological_scales(1, blitz::Range::all()) = width;
    topological_scales(2, blitz::Range::all()) = length;
    topological_scales(3, blitz::Range::all()) = genus;

    return topological_scales;
  }
}
