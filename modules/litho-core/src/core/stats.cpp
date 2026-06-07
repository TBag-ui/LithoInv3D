#include <litho_invert/core/stats.h>
#include <algorithm>
#include <cmath>

namespace litho_invert {

namespace {

// Compute 1D Durbin-Watson on pre-sorted residuals.
// DW = Σ(r_i - r_{i-1})² / Σ(r_i)²
double computeDW1D(const VectorXd& residuals) {
    Index n = residuals.size();
    if (n < 2) return 2.0;

    double denom = 0.0;
    double numer = 0.0;
    for (Index i = 0; i < n; ++i) {
        denom += residuals[i] * residuals[i];
    }
    if (denom < 1e-30) return 2.0;  // all residuals ~0

    for (Index i = 1; i < n; ++i) {
        double diff = residuals[i] - residuals[i - 1];
        numer += diff * diff;
    }
    return numer / denom;
}

} // anonymous namespace

DurbinWatsonResult computeDurbinWatson(const std::vector<Vector3d>& positions,
                                        const VectorXd& residuals) {
    Index n = static_cast<Index>(positions.size());
    if (n < 2 || residuals.size() != n) {
        return {2.0, 2.0};
    }

    // --- DW along x: sort by x, then y ---
    {
        std::vector<std::pair<double, double>> idx(n);
        for (Index i = 0; i < n; ++i) {
            idx[i] = {positions[static_cast<size_t>(i)].x(),
                      positions[static_cast<size_t>(i)].y()};
        }
        std::vector<Index> order(n);
        for (Index i = 0; i < n; ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](Index a, Index b) {
            if (std::abs(idx[a].first - idx[b].first) > 1e-9) {
                return idx[a].first < idx[b].first;
            }
            return idx[a].second < idx[b].second;
        });

        VectorXd sortedR(n);
        for (Index i = 0; i < n; ++i) {
            sortedR[i] = residuals[order[i]];
        }
        double dw_x = computeDW1D(sortedR);

        // --- DW along y: sort by y, then x ---
        std::sort(order.begin(), order.end(), [&](Index a, Index b) {
            if (std::abs(idx[a].second - idx[b].second) > 1e-9) {
                return idx[a].second < idx[b].second;
            }
            return idx[a].first < idx[b].first;
        });

        for (Index i = 0; i < n; ++i) {
            sortedR[i] = residuals[order[i]];
        }
        double dw_y = computeDW1D(sortedR);

        return {dw_x, dw_y};
    }
}

} // namespace litho_invert

