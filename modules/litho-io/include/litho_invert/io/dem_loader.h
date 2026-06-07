#pragma once
#include <litho_invert/core/common.h>
#include <vector>
#include <string>

namespace litho_invert {

// Load a DEM (digital elevation model) and sample it at the given (x,y) positions.
// Returns one z value per query position. Uses bilinear interpolation for points
// within the DEM extent; nearest-neighbor for points outside.
//
// Supported formats:
//   - XYZ ASCII grid: "x y z" per line, whitespace-separated.
//     Detects unique x and y values to determine grid dimensions.
//
//   - CSV: "x,y,z" per line, comma-separated. Same detection logic.
//     Header rows (non-numeric) are automatically skipped.
//
//   - Surface mesh (.ts/.obj): loaded via SurfaceMesh and interpolated.
//     (Stub — falls back to 0.0 if mesh loading is unavailable.)
//
std::vector<double> loadDEM(const std::string& path,
                            const std::vector<double>& xPositions,
                            const std::vector<double>& yPositions);

} // namespace litho_invert
