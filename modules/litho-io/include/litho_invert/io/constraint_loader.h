#pragma once
#include <litho_invert/core/common.h>
#include <vector>
#include <string>

namespace litho_invert {

class CSVConstraintLoader {
public:
    // Read constraints from CSV file.
    // Expected format: x,y,z_top,z_bottom,litho_group_id
    // Lines starting with # are ignored.
    static std::vector<Constraint> load(const std::string& filepath);
};

} // namespace litho_invert

