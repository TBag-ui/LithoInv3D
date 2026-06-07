#pragma once
#include <litho_invert/litho/litho_group.h>
#include <vector>
#include <string>

namespace litho_invert {

struct LithoConfig {
    std::vector<LithoGroup> groups;
    double bottomDepth = -5000.0;
};

class JSONLithoConfigLoader {
public:
    // Read lithology configuration from a minimal JSON file.
    // Expected format:
    // {
    //   "groups": [
    //     {"id": 0, "name": "...", "density": 2.70},
    //     ...
    //   ],
    //   "bottom_depth": 5000.0
    // }
    static LithoConfig load(const std::string& filepath);
};

} // namespace litho_invert

