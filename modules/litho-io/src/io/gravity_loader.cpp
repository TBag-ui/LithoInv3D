#include <litho_invert/io/gravity_loader.h>
#include <fstream>
#include <sstream>
#include <iostream>

namespace litho_invert {

GravityData CSVGravityLoader::load(const std::string& filepath) {
    GravityData data;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "CSVGravityLoader: Could not open file: " << filepath << std::endl;
        return data;
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

            // Skip header lines (first token is non-numeric, e.g. "x" or "X")
            if (values.empty() && !token.empty() && std::isalpha(static_cast<unsigned char>(token[0]))) {
                values.clear();
                break;
            }

            try {
                values.push_back(std::stod(token));
            } catch (...) {
                std::cerr << "CSVGravityLoader: Parse error at line " << lineNum
                          << ": '" << token << "' is not a number" << std::endl;
                values.clear();
                break;
            }
        }

        if (values.size() < 4) {
            std::cerr << "CSVGravityLoader: Skipping line " << lineNum
                      << ": expected at least 4 values (x,y,z,g_obs), got "
                      << values.size() << std::endl;
            continue;
        }

        GravityPoint pt;
        pt.position = Vector3d(values[0], values[1], values[2]);
        pt.g_obs = values[3];
        if (values.size() >= 5) {
            pt.g_std = values[4];
        }
        data.push_back(pt);
    }

    return data;
}

} // namespace litho_invert

