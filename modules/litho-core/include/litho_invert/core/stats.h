#pragma once
#include <litho_invert/core/common.h>
#include <vector>

namespace litho_invert {

struct DurbinWatsonResult {
    double dw_x = 2.0;  // DW along x-direction (sorted by x, then y)
    double dw_y = 2.0;  // DW along y-direction (sorted by y, then x)
};

// Compute Durbin-Watson statistic for 2D gridded observation data.
// Positions and residuals must have the same length.
// DW near 2.0 = random residuals (good).
// DW near 0.0 = positive autocorrelation (model underfitting).
// DW near 4.0 = negative autocorrelation (oscillating residuals).
DurbinWatsonResult computeDurbinWatson(const std::vector<Vector3d>& positions,
                                        const VectorXd& residuals);

} // namespace litho_invert

