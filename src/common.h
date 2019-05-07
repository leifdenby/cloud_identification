#include <stdint.h>
#include <math.h>

#ifndef CLOUD_IDENTIFICATION_COMMON_H
#define CLOUD_IDENTIFICATION_COMMON_H 1

// when scaling and casting to short ints we want to keep the same resolution
// for all input so we define the effective resolution here
static const float max_scalar_value=1000.;
static const float min_scalar_value=-max_scalar_value;

// when scaling we transform the data by taking the square root, so that's
// needed for the scaling too
static const float scaling_parameter=((float)SHRT_MAX)/sqrt(max_scalar_value);
static const float inv_scaling_parameter=1.0/scaling_parameter;

static const bool lperiodic = true; // is the domain periodic?

// PARAMETERS FOR CLUSTERING ALGORITHM
static const float mincolratio = 0.7; // fractional height used for deciding whether to merge cols
                                      // 0 corresponds to merging everything
                                      // 1 corresponds to not merging (steepest gradient association
                                      // with local maximum

// PARAMETERS WHICH DETERMINE MEMORY FOOTPRINT
// By what factor far can we make arrays smaller? (reduce memory footprint)
static const int blobfac = 2; // number of blobs wrt number of points
static const int borderfac = 1; // number of border points wrt number of points
static const int colfac = 1; // number of cols wrt number of points

typedef uint32_t indexint; //use unsigned array indices

#endif
