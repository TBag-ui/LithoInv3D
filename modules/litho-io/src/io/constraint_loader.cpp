#include <litho_invert/io/constraint_loader.h>
#include <fstream>
#include <sstream>
#include <iostream>

namespace litho_invert {

std::vector<Constraint> CSVConstraintLoader::load(const std::string& filepath) {
    std::vector<Constraint> constraints;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "CSVConstraintLoader: Could not open file: " << filepath << std::endl;
        return constraints;
    }

    std::string line;
    int lineNum = 0;
    while (std::getline(file, line)) {
        ++lineNum;

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string token;
        std::vector<double> values;

        while (std::getline(ss, token, ',')) {
            // Trim whitespace
            size_t start = token.find_first_not_of(" \t\r");
            size_t end = token.find_last_not_of(" \t\r");
            if (start == std::string::npos) continue;
            token = token.substr(start, end - start + 1);

            try {
                values.push_back(std::stod(token));
            } catch (...) {
                std::cerr << "CSVConstraintLoader: Parse error at line " << lineNum
                          << ": '" << token << "' is not a number" << std::endl;
                values.clear();
                break;
            }
        }

        if (values.size() < 5) {
            std::cerr << "CSVConstraintLoader: Skipping line " << lineNum
                      << ": expected 5 values (x,y,z_top,z_bottom,litho_group_id), got "
                      << values.size() << std::endl;
            continue;
        }

        Constraint c;
        c.position = Vector3d(values[0], values[1], 0.0);
        c.z_top = values[2];
        c.z_bottom = values[3];
        c.litho_group_id = static_cast<int>(values[4]);
        constraints.push_back(c);
    }

    return constraints;
}

} // namespace litho_invert

