#include <litho_invert/io/dem_loader.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <cmath>
#include <limits>

namespace litho_invert {

namespace {

// Bilinear interpolation on a regular grid.
// gridX, gridY: sorted unique coordinates
// gridZ: row-major data, gridZ[r * nx + c]
double bilinearInterp(double x, double y,
                      const std::vector<double>& gridX,
                      const std::vector<double>& gridY,
                      const std::vector<double>& gridZ) {
    int nx = static_cast<int>(gridX.size());
    int ny = static_cast<int>(gridY.size());
    if (nx < 2 || ny < 2) return 0.0;

    // Find bracketing x indices
    int ix = -1;
    if (x <= gridX.front()) {
        ix = 0;
    } else if (x >= gridX.back()) {
        ix = nx - 2;
    } else {
        for (int i = 0; i < nx - 1; ++i) {
            if (x >= gridX[i] && x <= gridX[i + 1]) {
                ix = i;
                break;
            }
        }
    }
    if (ix < 0) ix = 0;

    // Find bracketing y indices
    int iy = -1;
    if (y <= gridY.front()) {
        iy = 0;
    } else if (y >= gridY.back()) {
        iy = ny - 2;
    } else {
        for (int i = 0; i < ny - 1; ++i) {
            if (y >= gridY[i] && y <= gridY[i + 1]) {
                iy = i;
                break;
            }
        }
    }
    if (iy < 0) iy = 0;

    double dx = (gridX[ix + 1] - gridX[ix]);
    double dy = (gridY[iy + 1] - gridY[iy]);
    double tx = (dx > 0.0) ? (x - gridX[ix]) / dx : 0.0;
    double ty = (dy > 0.0) ? (y - gridY[iy]) / dy : 0.0;

    double z00 = gridZ[iy * nx + ix];
    double z10 = gridZ[iy * nx + ix + 1];
    double z01 = gridZ[(iy + 1) * nx + ix];
    double z11 = gridZ[(iy + 1) * nx + ix + 1];

    return (1.0 - tx) * (1.0 - ty) * z00
         + tx * (1.0 - ty) * z10
         + (1.0 - tx) * ty * z01
         + tx * ty * z11;
}

} // namespace

std::vector<double> loadDEM(const std::string& path,
                            const std::vector<double>& xPositions,
                            const std::vector<double>& yPositions) {
    const size_t n = xPositions.size();
    std::vector<double> result(n, 0.0);

    if (path.empty()) return result;

    std::ifstream file(path);
    if (!file.is_open()) return result;

    // Read all x,y,z triples (supports space-separated .xyz and comma-separated .csv)
    std::vector<double> xs, ys, zs;
    std::string line;
    while (std::getline(file, line)) {
        // Skip blank lines and comments
        if (line.empty() || line[0] == '#' || line[0] == ';') continue;

        // Replace commas with spaces so we can parse either format
        for (char& c : line) {
            if (c == ',') c = ' ';
        }

        std::istringstream ss(line);
        double x, y, z;
        if (ss >> x >> y >> z) {
            xs.push_back(x);
            ys.push_back(y);
            zs.push_back(z);
        }
    }

    if (xs.empty()) return result;

    // Find unique sorted x and y values to determine grid
    std::set<double> uniqueX(xs.begin(), xs.end());
    std::set<double> uniqueY(ys.begin(), ys.end());
    std::vector<double> gridX(uniqueX.begin(), uniqueX.end());
    std::vector<double> gridY(uniqueY.begin(), uniqueY.end());

    int nx = static_cast<int>(gridX.size());
    int ny = static_cast<int>(gridY.size());

    // Build a lookup: (x_index, y_index) -> z
    std::vector<double> gridZ(static_cast<size_t>(nx * ny),
                               std::numeric_limits<double>::quiet_NaN());

    // Map x/y values to indices
    for (size_t k = 0; k < xs.size(); ++k) {
        auto ix = static_cast<int>(
            std::lower_bound(gridX.begin(), gridX.end(), xs[k]) - gridX.begin());
        auto iy = static_cast<int>(
            std::lower_bound(gridY.begin(), gridY.end(), ys[k]) - gridY.begin());
        if (ix >= 0 && ix < nx && iy >= 0 && iy < ny) {
            gridZ[static_cast<size_t>(iy * nx + ix)] = zs[k];
        }
    }

    // Interpolate at each query position
    for (size_t i = 0; i < n; ++i) {
        result[i] = bilinearInterp(xPositions[i], yPositions[i],
                                   gridX, gridY, gridZ);
    }

    return result;
}

} // namespace litho_invert
