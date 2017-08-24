#include "common.h"
#include "blitz/array.h"

int write_netcdf(blitz::Array<indexint,3> dataext);
void load_mask(blitz::Array<bool,3> &maskext);
void load_field(blitz::Array<short,3> &fieldext);
