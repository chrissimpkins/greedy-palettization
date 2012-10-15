// Copyright 2011 Google Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef GREEDY_PALETTIZATION_COLOR_DISTANCE_INL_H_
#define GREEDY_PALETTIZATION_COLOR_DISTANCE_INL_H_

#include <stdint.h>
#include <stdlib.h>
#include <math.h>

// Fast computation of the distance of colors used by the palettization
// algorithm. The distance computation consists of two parts, first
// we compute a the weighted sums of squared color-differences, which is
// a monotone function of the actual distance and is enough to decide whether
// a color is closer to a palette color than an other color, and the second
// part is the scaling and square-rooting to the [0..255] range.

namespace greedy_palettization {

// weights for the intensity calculation
static const int ired = 2;
static const int igreen = 5;
static const int iblue = 1;
// normalization factor for the intensity calculation (2^ishift)
static const int ishift = 3;

// weights for the color similarity function
static const int wr = 1;
static const int wg = 1;
static const int wb = 1;
static const int wi = 1;
static const int wRGB = 1;
static const int wsum_RGB = wr + wg + wb + wi;
static const int wa = 1;
static const int wwr = wRGB * wr;
static const int wwg = wRGB * wg;
static const int wwb = wRGB * wb;
static const int wwi = wRGB * wi;
static const int wsum_RGBA = wwr + wwg + wwb + wwi + wa;

typedef uint8_t color_t;

inline int ColorIntManhDistanceRGB(color_t r1, color_t g1, color_t b1,
                                   color_t r2, color_t g2, color_t b2) {
  const int rd = r1 - r2;
  const int gd = g1 - g2;
  const int bd = b1 - b2;
  const int id = ired * rd + igreen * gd + iblue * bd;
  return wr * abs(rd) + wg * abs(gd) + wb * abs(bd) + wi * (abs(id) >> ishift);
}

inline int ScaleManhDistanceRGB(int d) {
  return d / wsum_RGB;
}

inline int ColorIntQuadDistanceRGB(color_t r1, color_t g1, color_t b1,
                                   color_t r2, color_t g2, color_t b2) {
  const int rd = r1 - r2;
  const int gd = g1 - g2;
  const int bd = b1 - b2;
  const int id = ired * rd + igreen * gd + iblue * bd;
  return wr * rd*rd + wg * gd*gd + wb * bd*bd + wi * ((id*id) >> (2 * ishift));
}

inline int ScaleQuadDistanceRGB(int d) {
  return static_cast<int>(sqrt(d * (1.0 / wsum_RGB)) + 0.5);
}

inline int ColorIntManhDistanceRGBA(color_t r1, color_t g1,
                                    color_t b1, color_t a1,
                                    color_t r2, color_t g2,
                                    color_t b2, color_t a2) {
  const int rd = r1 - r2;
  const int gd = g1 - g2;
  const int bd = b1 - b2;
  const int ad = a1 - a2;
  const int id = (ired * rd + igreen * gd + iblue * bd) >> ishift;
  return (wwr * abs(rd) + wwg * abs(gd) + wwb * abs(bd) + wwi * abs(id) +
          wa * abs(ad));
}

inline int ScaleManhDistanceRGBA(int d) {
  return d / wsum_RGBA;
}

inline int ColorIntQuadDistanceRGBA(color_t r1, color_t g1,
                                    color_t b1, color_t a1,
                                    color_t r2, color_t g2,
                                    color_t b2, color_t a2) {
  const int rd = r1 - r2;
  const int gd = g1 - g2;
  const int bd = b1 - b2;
  const int ad = a1 - a2;
  const int id = (ired * rd + igreen * gd + iblue * bd) >> ishift;
  return wwr * rd*rd + wwg * gd*gd + wwb * bd*bd + wwi * id*id + wa * ad*ad;
}

inline int ScaleQuadDistanceRGBA(int d) {
  return static_cast<int>(sqrt(d * (1.0 / wsum_RGBA)) + 0.5);
}

}  // namespace greedy_palettization

#endif  // GREEDY_PALETTIZATION_COLOR_DISTANCE_INL_H_
