#include <stdint.h>

#ifndef CLOUD_IDENTIFICATION_COMMON_H
#define CLOUD_IDENTIFICATION_COMMON_H 1

// PROPERTIES OF INPUT
static const float inv_scaling_parameter=1.0/3200.0; // data conversion parameter
                                                       // corresponding to scaling in pre-processing script
static const bool lperiodic = true; // is the domain periodic?

// PARAMETERS FOR CLUSTERING ALGORITHM
static const float mincolratio = 0.7; // fractional height used for deciding whether to merge cols
                                      // 0 corresponds to merging everything
                                      // 1 corresponds to not merging (steepest gradient association
                                      // with local maximum

// PARAMETERS WHICH DETERMINE MEMORY FOOTPRINT
// By what factor far can we make arrays smaller? (reduce memory footprint)
static const int blobfac = 2; // number of blobs wrt number of points
static const int borderfac = 2; // number of border points wrt number of points
static const int colfac = 2; // number of cols wrt number of points

typedef uint32_t indexint; //use unsigned array indices

#endif
