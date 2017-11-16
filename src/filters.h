#include "blitz/array.h"

namespace filters {
  void remove_intersecting(blitz::Array<indexint,3> &labels, blitz::Array<bool, 3> &maskext);
}
