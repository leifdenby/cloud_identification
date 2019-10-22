#include <blitz/array.h>

#include "common.h"

namespace filters {
  /*
   * Remove all objects (identified by non-zero integers in `labels`) which
   * intersect with `mask`, by setting all elements with same integer value to
   * zero.
   *
   * OBS: updates `labels` in-place to reduce memory footprint
   */
  void remove_intersecting(
    blitz::Array<indexint,3> &labels,
    blitz::Array<bool, 3> &maskext
  ) {
    int label_max = max(labels);
    blitz::TinyVector<int,3> shape = labels.shape();

    int nx = shape[0];
    int ny = shape[1];
    int nz = shape[2];

    int d_val = -1;

    // OBS: could use linked-list data-structure to speed things up a bit...
    blitz::Array<bool,1> n_remove = blitz::Array<bool,1>(label_max+1);
    n_remove = 0;

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          if (maskext(i,j,k)) {
            d_val = labels(i,j,k);
            n_remove(d_val) = true;
          }
        }
      }
    }

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          d_val = labels(i,j,k);
          if (n_remove(d_val)) {
            labels(i,j,k) = 0;
          }
        }
      }
    }
  }

  /*
   * Remove all objects (identified by non-zero integers in `labels`) which
   * intersect at (i,j) position in 2D `mask`, by setting all elements with
   * same integer value to zero.
   *
   * OBS: updates `labels` in-place to reduce memory footprint
   */
  void remove_intersecting_xy(
    blitz::Array<indexint,3> &labels,
    blitz::Array<bool, 2> &maskext
  ) {
    int label_max = max(labels);
    blitz::TinyVector<int,3> shape = labels.shape();

    int nx = shape[0];
    int ny = shape[1];
    int nz = shape[2];

    int d_val = -1;

    // OBS: could use linked-list data-structure to speed things up a bit...
    blitz::Array<bool,1> n_remove = blitz::Array<bool,1>(label_max+1);
    n_remove = 0;

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          if (maskext(i,j)) {
            d_val = labels(i,j,k);
            n_remove(d_val) = true;
          }
        }
      }
    }

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          d_val = labels(i,j,k);
          if (n_remove(d_val)) {
            labels(i,j,k) = 0;
          }
        }
      }
    }
  }

  /*
   * Remove all objects (identified by non-zero integers in `labels`) which
   * aren't in the `keep`
   *
   * OBS: updates `labels` in-place to reduce memory footprint
   */
  void filter_labels(
    blitz::Array<indexint,3> &labels,
    blitz::Array<indexint,1> &keep
  ) {
    int label_max = max(labels);
    blitz::TinyVector<int,3> shape = labels.shape();

    int nx = shape[0];
    int ny = shape[1];
    int nz = shape[2];

    int n_to_keep = keep.size();

    int d_val = -1;

    bool should_keep = false;

    for (int i=0; i<nx; i++) {
      for (int j=0; j<ny; j++) {
        for (int k=0; k<nz; k++) {
          d_val = labels(i,j,k);
          should_keep = false;

          for (int n=0; n<n_to_keep; n++) {
            if (d_val == keep(n)) {
              should_keep = true;
              break;
            }
          }

          if (!should_keep) {
            labels(i,j,k) = 0;
          }
        }
      }
    }
  }
}
