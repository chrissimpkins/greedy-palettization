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

#ifndef GREEDY_PALETTIZATION_PALETTIZE_H_
#define GREEDY_PALETTIZATION_PALETTIZE_H_

typedef unsigned char uint8;

namespace greedy_palettization {

// These functions are the main interfaces of the palettization algorithm.
// Their input is an RGB or RGBA image as a sequence of bytes in red, green,
// blue, alpha order, and their output is the palette and the palettized image.
// The palette is also a sequence of bytes in the same RGB(A) order.
// Optionally, an array of colors can be specified which will be preserved by
// the palettization.
// The return value is false if for some reason the image can not really well
// palettized.
//
// Quality is an integer between 1 and 1024, with 128 normalized to reasonable
// results. The quality controls only the number of the colors in the palette.
// Pass kPalettizationQualityMax into it if you want to use the all the
// max_palette_size colors.
// Pass QUALITY_HIGH into it if you don't care to figure out an optimal value
// for your application.

static const int kPalettizationQualityMax = 1024;

enum PaletteEstimationQuality {
  QUALITY_LOW = 32,
  QUALITY_MEDIUM = 128,
  QUALITY_HIGH = 512,
};

bool GREEDY_PALETTIZATION_DLL_DECL
PalettizeWithGreedySelectionRGB(const int xsize, const int ysize,
                                const uint8* const image,
                                const uint8* const colors_to_keep,
                                const int num_colors_to_keep,
                                const int max_palette_size,
                                const int quality,
                                int* const palette_size,
                                uint8* const palette,
                                uint8* const palettized_image);

bool GREEDY_PALETTIZATION_DLL_DECL
PalettizeWithGreedySelectionRGBA(const int xsize, const int ysize,
                                 const uint8* const image,
                                 const uint8* const colors_to_keep,
                                 const int num_colors_to_keep,
                                 const int max_palette_size,
                                 const int quality,
                                 int* const palette_size,
                                 uint8* const palette,
                                 uint8* const palettized_image);

}  // namespace greedy_palettization

#endif  // GREEDY_PALETTIZATION_PALETTIZE_H_
