#pragma once
#include <litho_invert/core/common.h>
#include <string>

namespace litho_invert {

class CSVMagneticLoader {
public:
    // Read magnetic data from CSV file.
    // Expected format: x,y,z,t_obs[,t_std]
    // Lines starting with # are ignored.
    static MagneticData load(const std::string& filepath);
};

} // namespace litho_invert

