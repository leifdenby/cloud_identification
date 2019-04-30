#include "blitz/array.h"

#include "common.h"

void find_objects(
  blitz::Array<short,3> &fieldext,
  blitz::Array<bool, 3> &maskext,
  blitz::Array<indexint,3> &dataext
);

void scale_field(
  blitz::Array<double,3> &field_input,
  blitz::Array<short,3> &field_scaled
);
