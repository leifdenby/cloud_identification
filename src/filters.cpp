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
            //if (maskext(i,j,k)) {
              printf("%d %d %d => %d %s\n", i, j, k, labels(i,j,k), maskext(i,j,k) ? "T" : "F");
            //}
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
    for (int m=0; m<label_max+1; m++) {
      printf("%d: %d\n", m, n_remove(m));
    }
  }
}
