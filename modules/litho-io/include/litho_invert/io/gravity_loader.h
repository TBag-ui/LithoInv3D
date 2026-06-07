#pragma once
#include <litho_invert/core/common.h>
#include <string>

namespace litho_invert {

class CSVGravityLoader {
public:
    // Read gravity data from CSV file.
    // Expected format: x,y,z,g_obs[,g_std]
    // Lines starting with # are ignored.
    static GravityData load(const std::string& filepath);
};

} // namespace litho_invert

