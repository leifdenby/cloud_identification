#include "blitz/array.h"

namespace filters {
  void remove_intersecting(blitz::Array<indexint,3> &labels, blitz::Array<bool, 3> &maskext);
  void remove_intersecting_xy(blitz::Array<indexint,3> &labels, blitz::Array<bool, 2> &maskext);
}
