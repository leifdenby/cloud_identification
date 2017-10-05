#include <blitz/array.h>

namespace minkowski {

  blitz::Array<int,1> N0(blitz::Array<int,3> &labels);
  blitz::Array<int,1> N1(blitz::Array<int,3> &labels);
  blitz::Array<int,1> N2(blitz::Array<int,3> &labels);
  blitz::Array<int,1> N3(blitz::Array<int,3> &labels);
  blitz::Array<float,2> topological_scales(blitz::Array<int,3> &labels, float dx);
  blitz::Array<float,1> planarity(blitz::Array<float,2> &mf);
  blitz::Array<float,1> filamentarity(blitz::Array<float,2> &mf);
}
