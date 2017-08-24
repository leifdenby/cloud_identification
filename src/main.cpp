#include "blitz/array.h"

#include "file_io.h"
#include "cloud_identification.h"

// MAIN PROGRAM    
int main() {
  // array that holds the mask as boolean (including halo cells)
  blitz::Array<bool,3> maskext;

  // array that holds field values
  blitz::Array<short,3> fieldext;

  // array that hold actual numbers
  blitz::Array<indexint,3> dataext;
  find_objects(fieldext, maskext, dataext);

  write_netcdf(dataext);
}
