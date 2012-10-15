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
//
// This file contains the implementation of the palettization algorithm.
// The input constists of:
//   - an RGB or RGBA image as a sequence of bytes, together with the dimensions
//     of the image
//   - an array of static colors that must be represented exactly in the palette
//   - the maximum size of the palette (this must be <= 256)
// The output consists of:
//   - the actual size of the palette
//   - the palette as a sequence of bytes in RGB(A) order
//   - the palettized image in indexed format
// The algorithm has the following phases:
//   - Counting and indexing the colors: each color of the input image gets
//     a unique integer id in order of appearance, and the population count
//     for each color is determined.
//   - Initializing the palette: at first, the colors that were specified in
//     the colors to keep array and those with a population count higher than
//     a certain threshold (chosen so that at most one quarter of the palette
//     can be filled this way) are inserted to the palette. If the palette is
//     still empty then the color with the highest population count is selected
//     first.
//   - Greedy selection of the colors: the next palette color is the one that
//     maximizes the population count times the minimal distance to the current
//     palette. For this, we maintain for each color the minimal distance to
//     the palette and the closest palette color and update these values after
//     the selection of each new palette color.
//   - The mapping of the colors to the palette does not require further
//     computation, because we periodically updated this mapping during the
//     selection process.

#include <stdint.h>

#include <algorithm>
#ifdef ADDITIONAL_DEBUG_OUTPUT
#include <cstdio>
#endif
#include <limits>
#include <map>
#include <set>
#include <vector>

#include "interlace_bits-inl.h"
#include "color_distance-inl.h"
#include "greedy-palettization/palettize.h"

namespace greedy_palettization {

// 2^13 priority levels for the PQ seems to be a good compromise between
// accuracy, running time and stack space usage.
static const int kMaxPriority = 1 << 13;
static const int kMaxLevel = 3;
// The maximum allowed average pixel error when the quality is half of
// kPalettizationQualityMax.
static const int kErrorScaleRGB = 1;
static const int kErrorScaleRGBA = 1;

typedef uint8_t color_t;

int AveragePixelErrorThresholdRGB(int quality) {
  quality = std::max(1, std::min(quality, kPalettizationQualityMax));
  return kErrorScaleRGB * kPalettizationQualityMax / quality - kErrorScaleRGB;
}

int AveragePixelErrorThresholdRGBA(int quality) {
  quality = std::max(1, std::min(quality, kPalettizationQualityMax));
  return kErrorScaleRGBA * kPalettizationQualityMax / quality - kErrorScaleRGBA;
}

// This function is used in the multi-resolution grid to be able to compute
// the keys for the different resolutions by just shifting the first key.
inline int InterlaceBitsRGB(color_t r, color_t g, color_t b) {
  int z = 0;
  for (int i = 0; i < 7; ++i) {
    z += (r >> 5) & 4;
    z += (g >> 6) & 2;
    z += (b >> 7);
    z <<= 3;
    r <<= 1;
    g <<= 1;
    b <<= 1;
  }
  z += (r >> 5) & 4;
  z += (g >> 6) & 2;
  z += (b >> 7);
  return z;
}

// This function will compute the actual priorities of the colors based on
// the current distance from the palette, the population count and the signals
// from the multi-resolution grid.
inline int Priority(int d, int n, const int* density, const int* radius) {
  int p = d * n;
  for (int level = 0; level < kMaxLevel; ++level) {
    if (d > radius[level]) {
      p += density[level] * (d - radius[level]);
    }
  }
  return std::min(kMaxPriority - 1, p >> 4);
}


// The function updates the minimal distances, the clustering and the
// quantization error after the insertion of the new color into the palette.
void AddToRGBPalette(const color_t* red,
                     const color_t* green,
                     const color_t* blue,
                     const int* count,  // histogram of colors
                     const int index,   // index of color to be added
                     const int k,       // size of current palette
                     const int n,       // number of colors
                     int* dist,         // array of distances from palette
                     int* cluster,      // mapping of color indices to palette
                     int* center,       // the inverse mapping
                     int64_t* error) {    // measure of the quantization error
  center[k] = index;
  cluster[index] = k;
  *error -= dist[index] * count[index];
  dist[index] = 0;
  for (int j = 0; j < n; ++j) {
    if (dist[j] > 0) {
      const int d = ColorIntQuadDistanceRGB(red[index],
                                            green[index],
                                            blue[index],
                                            red[j], green[j], blue[j]);
      if (d < dist[j]) {
        *error += (d - dist[j]) * count[j];
        dist[j] = d;
        cluster[j] = k;
      }
    }
  }
}

// The function updates the minimal distances, the clustering and the
// quantization error after the insertion of the new color into the palette.
void AddToRGBAPalette(const color_t* red,
                      const color_t* green,
                      const color_t* blue,
                      const color_t* alpha,
                      const int* count,  // histogram of colors
                      const int index,   // index of color to be added
                      const int k,       // size of current palette
                      const int n,       // number of colors
                      int* dist,         // array of distances from palette
                      int* cluster,      // mapping of color indices to palette
                      int* center,       // the inverse mapping
                      int64_t* error) {    // measure of the quantization error
  center[k] = index;
  cluster[index] = k;
  *error -= dist[index] * count[index];
  dist[index] = 0;
  for (int j = 0; j < n; ++j) {
    if (dist[j] > 0) {
      const int d = ColorIntQuadDistanceRGBA(red[index], green[index],
                                             blue[index], alpha[index],
                                             red[j], green[j],
                                             blue[j], alpha[j]);
      if (d < dist[j]) {
        *error += (d - dist[j]) * count[j];
        dist[j] = d;
        cluster[j] = k;
      }
    }
  }
}

bool PalettizeWithGreedySelectionRGB(const int xsize, const int ysize,
                                     const uint8* const image,
                                     const uint8* const colors_to_keep,
                                     const int num_colors_to_keep,
                                     const int max_palette_size,
                                     const int quality,
                                     int* const palette_size,
                                     uint8* const palette,
                                     uint8* const palettized_image) {
  const int num_pixels = xsize * ysize;
  int* indexed_image(new int[num_pixels]);
  color_t* red(new color_t[num_pixels]);
  color_t* green(new color_t[num_pixels]);
  color_t* blue(new color_t[num_pixels]);
  std::vector<int> count(num_pixels, 0);
  int n = 0;  // number of colors

  // Here we first build an index of all the different colors in the input
  // image. To do this we map the 24 bit RGB representation of the colors
  // to a unique integer index assigned to the different colors in order of
  // appearence in the image.
  // The indexing of the colors is based on std::map, but map tends to be
  // slow if it has many elements, so we distribute the colors bewteen a
  // couple of maps to make things faster.
  std::map<uint32_t, int> index_map_array[4096];
  std::map<uint32_t, int>::iterator index_map_lookup;
  const uint8* imagep = &image[0];
  int* indexp = &indexed_image[0];
  // to make sure that prev_pixel != first pixel
  uint32_t prev_pixel = 0x01000000;
  int index = 0;
  for (int y = 0; y < ysize; ++y) {
    for (int x = 0; x < xsize; ++x) {
      color_t r = *imagep++;
      color_t g = *imagep++;
      color_t b = *imagep++;
      uint32_t pixel = (r << 16) | (g << 8) | b;
      if (pixel != prev_pixel) {
        prev_pixel = pixel;
        std::map<uint32_t, int>& index_map = index_map_array[pixel >> 12];
        index_map_lookup = index_map.find(pixel);
        if (index_map_lookup != index_map.end()) {
          index = index_map_lookup->second;
        } else {
          index_map[pixel] = index = n++;
          red[index] = r;
          green[index] = g;
          blue[index] = b;
        }
      }
      ++count[index];
      *indexp++ = index;
    }
  }

#ifdef ADDITIONAL_DEBUG_OUTPUT
  printf("Number of colors in image: %d\n", n);
#endif

  std::vector<int> dist(n, std::numeric_limits<int>::max());
  std::vector<int> cluster(n);
  std::vector<bool> in_palette(n, false);
  int center[256];
  int k = 0;  // palette size
  const int count_threshold = (xsize * ysize * 4) / max_palette_size;
  const int64_t error_threshold =
      num_pixels * AveragePixelErrorThresholdRGB(quality);
  int64_t error = 0;  // quantization error

  if (colors_to_keep != NULL) {
    std::set<uint32_t> static_colors;
    const uint8* colorp = &colors_to_keep[0];
    for (int i = 0; i < num_colors_to_keep; ++i) {
      color_t r = *colorp++;
      color_t g = *colorp++;
      color_t b = *colorp++;
      static_colors.insert((r << 16) | (g << 8) | b);
    }
    for (int i = 0; i < n; ++i) {
      uint32_t color = (red[i] << 16) | (green[i] << 8) | blue[i];
      if (static_colors.find(color) != static_colors.end()) {
        AddToRGBPalette(&red[0], &green[0], &blue[0], &count[0], i, k++, n,
                        &dist[0], &cluster[0], &center[0], &error);
        in_palette[i] = true;
      }
    }
  }
  int max_count = 0;
  int winner = 0;
  for (int i = 0; i < n; ++i) {
    if (count[i] > max_count) {
      max_count = count[i];
      winner = i;
    }
    if (!in_palette[i] && count[i] > count_threshold) {
#ifdef ADDITIONAL_DEBUG_OUTPUT
      printf("%3d %5d    %02x%02x%02x\n", k, count[i],
             red[i], green[i], blue[i]);
#endif
      AddToRGBPalette(&red[0], &green[0], &blue[0], &count[0], i, k++, n,
                      &dist[0], &cluster[0], &center[0], &error);
      in_palette[i] = true;
    }
  }
  if (k == 0) {
#ifdef ADDITIONAL_DEBUG_OUTPUT
    printf("%3d %5d    %02x%02x%02x\n", k, count[winner],
           red[winner], green[winner], blue[winner]);
#endif
    AddToRGBPalette(&red[0], &green[0], &blue[0], &count[0], winner, k++, n,
                    &dist[0], &cluster[0], &center[0], &error);
    in_palette[winner] = true;
  }

  // Calculation of the multi-resolution density grid.
  std::vector<int> density(n * kMaxLevel);
  std::vector<int> radius(n * kMaxLevel);
  std::vector<int> histogram[kMaxLevel];
  for (int level = 0; level < kMaxLevel; ++level) {
    histogram[level].resize(1 << (18 - kMaxLevel * level), 0);
  }
  for (int i = 0; i < n; ++i) {
    if (!in_palette[i]) {
      const int key = InterlaceBitsRGB(red[i], green[i], blue[i]) >> 6;
      for (int level = 0; level < kMaxLevel; ++level) {
        histogram[level][key >> (3 * level)] += count[i];
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    if (!in_palette[i]) {
      for (int level = 0; level < kMaxLevel; ++level) {
        const int mask = (4 << level) - 1;
        const int rd = std::max(red[i] & mask, mask - (red[i] & mask));
        const int gd = std::max(green[i] & mask, mask - (green[i] & mask));
        const int bd = std::max(blue[i] & mask, mask - (blue[i] & mask));
        radius[i * kMaxLevel + level] =
            ScaleQuadDistanceRGB(ColorIntQuadDistanceRGB(0, 0, 0, rd, gd, bd));
      }
      const int key = InterlaceBitsRGB(red[i], green[i], blue[i]) >> 6;
      if (kMaxLevel > 0) {
        density[i * kMaxLevel] = histogram[0][key] - count[i];
      }
      for (int level = 1; level < kMaxLevel; ++level) {
        density[i * kMaxLevel + level] =
            (histogram[level][key >> (3*level)] -
             histogram[level-1][key >> (3*level-3)]);
      }
    }
  }

  // Calculate the initial error now that the palette has been initialized.
  error = 0;
  for (int i = 0; i < n; ++i) {
    error += dist[i] * count[i];
  }

  std::vector<int>* bucket_array(new std::vector<int>[kMaxPriority]);
  int top_priority = -1;
  for (int i = 0; i < n; ++i) {
    if (!in_palette[i]) {
      int priority = Priority(ScaleQuadDistanceRGB(dist[i]), count[i],
                              &density[i * kMaxLevel],
                              &radius[i * kMaxLevel]);
      bucket_array[priority].push_back(i);
      top_priority = std::max(priority, top_priority);
    }
  }
  double error_accum = 0;
  while (top_priority >= 0 && k < max_palette_size) {
    if (error < error_threshold) {
      error_accum += std::min(error_threshold, error_threshold - error);
      if (error_accum >= 10 * error_threshold) {
        break;
      }
    }
    int i = bucket_array[top_priority].back();
    int priority = Priority(ScaleQuadDistanceRGB(dist[i]), count[i],
                            &density[i * kMaxLevel],
                            &radius[i * kMaxLevel]);
    if (priority < top_priority) {
      bucket_array[priority].push_back(i);
    } else {
#ifdef ADDITIONAL_DEBUG_OUTPUT
      int d = ScaleQuadDistanceRGB(dist[i]);
      int objval = 0;
      for (int j = 0; j < n; ++j) {
        objval += count[j] * ScaleQuadDistanceRGB(dist[j]);
      }
      printf("%3d %5d %3d  ", k, count[i], d);
      for (int level = 0; level < kMaxLevel; ++level) {
        printf(" %5d %2d", density[i * kMaxLevel + level],
               radius[i * kMaxLevel + level]);
      }
      printf("  %5d    %02x%02x%02x    %7d  %10lld\n", priority,
             red[i], green[i], blue[i], objval, error);
#endif
      AddToRGBPalette(&red[0], &green[0], &blue[0], &count[0], i, k++, n,
                      &dist[0], &cluster[0], &center[0], &error);
    }
    bucket_array[top_priority].pop_back();
    while (top_priority >= 0 && bucket_array[top_priority].empty()) {
      --top_priority;
    }
  }

#ifdef ADDITIONAL_DEBUG_OUTPUT
  int objval = 0;
  for (int i = 0; i < n; ++i) {
    objval += count[i] * ScaleQuadDistanceRGB(dist[i]);
  }
  printf("                                                               "
         "%7d  %10lld\n", objval, error);
#endif

  uint8* palettep = &palette[0];
  for (int j = 0; j < k; ++j) {
    const int index = center[j];
    *palettep++ = red[index];
    *palettep++ = green[index];
    *palettep++ = blue[index];
  }
  for (int i = 0; i < xsize * ysize; ++i) {
    palettized_image[i] = cluster[indexed_image[i]];
  }
  if (palette_size != NULL) { *palette_size = k; }

  delete[] bucket_array;
  delete[] blue;
  delete[] green;
  delete[] red;
  delete[] indexed_image;

  return true;
}

bool PalettizeWithGreedySelectionRGBA(const int xsize, const int ysize,
                                      const uint8* const image,
                                      const uint8* const colors_to_keep,
                                      const int num_colors_to_keep,
                                      const int max_palette_size,
                                      const int quality,
                                      int* const palette_size,
                                      uint8* const palette,
                                      uint8* const palettized_image) {
  const int num_pixels = xsize * ysize;
  int* indexed_image(new int[num_pixels]);
  color_t* red(new color_t[num_pixels]);
  color_t* green(new color_t[num_pixels]);
  color_t* blue(new color_t[num_pixels]);
  color_t* red_scaled(new color_t[num_pixels]);
  color_t* green_scaled(new color_t[num_pixels]);
  color_t* blue_scaled(new color_t[num_pixels]);
  color_t* alpha(new color_t[num_pixels]);
  std::vector<int> count(num_pixels, 0);
  int n = 0;  // number of colors

  // Here we first build an index of all the different colors in the input
  // image. To do this we map the 32 bit RGBA representation of the colors
  // to a unique integer index assigned to the different colors in order of
  // appearence in the image.
  // The indexing of the colors is based on std::map, but map tends to be
  // slow if it has many elements, so we distribute the colors bewteen a
  // couple of maps to make things faster.
  std::map<uint32_t, int> index_map_array[4096];
  std::map<uint32_t, int>::iterator index_map_lookup;
  const uint8* imagep = &image[0];
  int* indexp = &indexed_image[0];
  color_t r = image[0];
  color_t g = image[1];
  color_t b = image[2];
  color_t a = image[3];
  // to make sure that prev_pixel != first pixel
  uint32_t prev_pixel = (r << 24) + (g << 16) + (b << 8) + a + 1;
  int index = 0;
  for (int y = 0; y < ysize; ++y) {
    for (int x = 0; x < xsize; ++x) {
      color_t r = *imagep++;
      color_t g = *imagep++;
      color_t b = *imagep++;
      color_t a = *imagep++;
      // If the alpha is zero, then the (r, g, b) value does not matter
      uint32_t pixel = a == 0 ? 0 : (r << 24) | (g << 16) | (b << 8) | a;
      if (pixel != prev_pixel) {
        prev_pixel = pixel;
        std::map<uint32_t, int>& index_map = index_map_array[pixel >> 20];
        index_map_lookup = index_map.find(pixel);
        if (index_map_lookup != index_map.end()) {
          index = index_map_lookup->second;
        } else {
          index_map[pixel] = index = n++;
          red[index] = r;
          green[index] = g;
          blue[index] = b;
          alpha[index] = a;
          red_scaled[index]   = (r * (a + 1)) >> 8;
          green_scaled[index] = (g * (a + 1)) >> 8;
          blue_scaled[index]  = (b * (a + 1)) >> 8;
        }
      }
      ++count[index];
      *indexp++ = index;
    }
  }

#ifdef ADDITIONAL_DEBUG_OUTPUT
  printf("Number of colors in image: %d\n", n);
#endif

  std::vector<int> dist(n, std::numeric_limits<int>::max());
  std::vector<int> cluster(n);
  std::vector<bool> in_palette(n, false);
  int center[256];
  int k = 0;  // palette size
  const int count_threshold = (xsize * ysize * 4) / max_palette_size;
  const int64_t error_threshold =
      num_pixels * AveragePixelErrorThresholdRGBA(quality);
  int64_t error = 0;  // quantization error

  if (colors_to_keep != NULL) {
    std::set<uint32_t> static_colors;
    const uint8* colorp = &colors_to_keep[0];
    for (int i = 0; i < num_colors_to_keep; ++i) {
      color_t r = *colorp++;
      color_t g = *colorp++;
      color_t b = *colorp++;
      color_t a = *colorp++;
      static_colors.insert((r << 24) | (g << 16) | (b << 8) | a);
    }
    for (int i = 0; i < n; ++i) {
      uint32_t color = ((red[i] << 24) | (green[i] << 16) |
                            (blue[i] << 8) | alpha[i]);
      if (static_colors.find(color) != static_colors.end()) {
        AddToRGBAPalette(&red_scaled[0], &green_scaled[0], &blue_scaled[0],
                         &alpha[0], &count[0], i, k++, n,
                         &dist[0], &cluster[0], center, &error);
        in_palette[i] = true;
      }
    }
  }
  int max_count = 0;
  int winner = 0;
  for (int i = 0; i < n; ++i) {
    if (count[i] > max_count) {
      max_count = count[i];
      winner = i;
    }
    if (!in_palette[i] && count[i] > count_threshold) {
      AddToRGBAPalette(&red_scaled[0], &green_scaled[0], &blue_scaled[0],
                       &alpha[0], &count[0], i, k++, n,
                       &dist[0], &cluster[0], center, &error);
      in_palette[i] = true;
    }
  }
  if (k == 0) {
#ifdef ADDITIONAL_DEBUG_OUTPUT
    printf("%3d %5d %3d %4d    %02x%02x%02x%02x\n", k, count[winner], -1, -1,
           red_scaled[winner], green_scaled[winner], blue_scaled[winner],
           alpha[winner]);
#endif
    AddToRGBAPalette(&red_scaled[0], &green_scaled[0], &blue_scaled[0],
                     &alpha[0], &count[0], winner, k++, n,
                     &dist[0], &cluster[0], center, &error);
    in_palette[winner] = true;
  }

  // Calculation of the multi-resolution density grid.
  std::vector<int> density(n * kMaxLevel);
  std::vector<int> radius(n * kMaxLevel);
  std::map<uint32_t, int> histogram[kMaxLevel];
  for (int i = 0; i < n; ++i) {
    if (!in_palette[i]) {
      const int key = Interlace4Bytes(red_scaled[i], green_scaled[i],
                                      blue_scaled[i], alpha[i]) >> 8;
      for (int level = 0; level < kMaxLevel; ++level) {
        histogram[level][key >> (4*level)] += count[i];
      }
    }
  }
  for (int i = 0; i < n; ++i) {
    if (!in_palette[i]) {
      for (int level = 0; level < kMaxLevel; ++level) {
        const int mask = (4 << level) - 1;
        const int rd = std::max(red_scaled[i] & mask,
                                mask - (red_scaled[i] & mask));
        const int gd = std::max(green_scaled[i] & mask,
                                mask - (green[i] & mask));
        const int bd = std::max(blue_scaled[i] & mask,
                                mask - (blue_scaled[i] & mask));
        const int ad = std::max(alpha[i] & mask, mask - (alpha[i] & mask));
        radius[i * kMaxLevel + level] =
            ScaleQuadDistanceRGBA(ColorIntQuadDistanceRGBA(0, 0, 0, 0,
                                                           rd, gd, bd, ad));
      }
      const int key = Interlace4Bytes(red_scaled[i], green_scaled[i],
                                      blue_scaled[i], alpha[i]) >> 8;
      if (kMaxLevel > 0) {
        density[i * kMaxLevel] = histogram[0][key] - count[i];
      }
      for (int level = 1; level < kMaxLevel; ++level) {
        density[i * kMaxLevel + level] =
            (histogram[level][key >> (4*level)] -
             histogram[level-1][key >> (4*level -4)]);
      }
    }
  }

  // Calculate the initial error now that the palette has been initialized.
  error = 0;
  for (int i = 0; i < n; ++i) {
    error += dist[i] * count[i];
  }

  std::vector<int> bucket_array[kMaxPriority];
  int top_priority = -1;
  for (int i = 0; i < n; ++i) {
    if (!in_palette[i]) {
      int priority = Priority(ScaleQuadDistanceRGBA(dist[i]), count[i],
                              &density[i * kMaxLevel],
                              &radius[i * kMaxLevel]);
      bucket_array[priority].push_back(i);
      top_priority = std::max(priority, top_priority);
    }
  }

  double error_accum = 0;
  while (top_priority >= 0 && k < max_palette_size) {
    if (error < error_threshold) {
      error_accum += std::min(error_threshold, error_threshold - error);
      if (error_accum >= 10 * error_threshold) {
        break;
      }
    }
    int i = bucket_array[top_priority].back();
    int priority = Priority(ScaleQuadDistanceRGBA(dist[i]), count[i],
                            &density[i * kMaxLevel],
                            &radius[i * kMaxLevel]);
    if (priority < top_priority) {
      bucket_array[priority].push_back(i);
    } else {
#ifdef ADDITIONAL_DEBUG_OUTPUT
      printf("%3d %5d %3d %4d    %02x%02x%02x%02x   %10lld\n", k, count[i],
             ScaleQuadDistanceRGBA(dist[i]), priority,
             red_scaled[i], green_scaled[i], blue_scaled[i], alpha[i], error);
#endif
      AddToRGBAPalette(&red_scaled[0], &green_scaled[0], &blue_scaled[0],
                       &alpha[0], &count[0], i, k++, n,
                       &dist[0], &cluster[0], center, &error);
    }
    bucket_array[top_priority].pop_back();
    while (top_priority >= 0 && bucket_array[top_priority].empty()) {
      --top_priority;
    }
  }

  uint8* palettep = &palette[0];
  for (int j = 0; j < k; ++j) {
    const int index = center[j];
    *palettep++ = red[index];
    *palettep++ = green[index];
    *palettep++ = blue[index];
    *palettep++ = alpha[index];
  }
  for (int i = 0; i < xsize * ysize; ++i) {
    palettized_image[i] = cluster[indexed_image[i]];
  }
  if (palette_size != NULL) { *palette_size = k; }

  delete[] alpha;
  delete[] blue_scaled;
  delete[] green_scaled;
  delete[] red_scaled;
  delete[] blue;
  delete[] green;
  delete[] red;
  delete[] indexed_image;

  return true;
}

}  // namespace greedy_palettization
