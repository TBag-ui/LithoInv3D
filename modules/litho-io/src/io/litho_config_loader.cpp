#include <litho_invert/io/litho_config_loader.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>

namespace litho_invert {

namespace {

// Trim whitespace from both ends
std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    return s.substr(start, end - start + 1);
}

// Remove surrounding quotes from a string
std::string unquote(const std::string& s) {
    std::string t = trim(s);
    if (t.size() >= 2 && ((t.front() == '"' && t.back() == '"') ||
                           (t.front() == '\'' && t.back() == '\''))) {
        return t.substr(1, t.size() - 2);
    }
    return t;
}

// Extract a numeric value after a key string in a line
// Handles "key": value, "key" : value, "key": "value", etc.
bool extractDouble(const std::string& line, const std::string& key, double& value) {
    // Find key in the line
    size_t pos = line.find(key);
    if (pos == std::string::npos) return false;

    // Find colon after key
    size_t colon = line.find(':', pos + key.size());
    if (colon == std::string::npos) return false;

    // Extract text after colon
    std::string after = line.substr(colon + 1);

    // Check if it's a quoted string (for error handling, not used for doubles)
    // Find the start of the numeric value
    size_t numStart = after.find_first_of("-0123456789.");
    if (numStart == std::string::npos) return false;

    size_t numEnd = after.find_first_not_of("0123456789.eE+-", numStart);
    std::string numStr = after.substr(numStart, numEnd - numStart);

    try {
        value = std::stod(numStr);
        return true;
    } catch (...) {
        return false;
    }
}

// Extract a string value after a key
bool extractString(const std::string& line, const std::string& key, std::string& value) {
    size_t pos = line.find(key);
    if (pos == std::string::npos) return false;

    size_t colon = line.find(':', pos + key.size());
    if (colon == std::string::npos) return false;

    std::string after = line.substr(colon + 1);

    // Find quoted string
    size_t quoteStart = after.find('"');
    if (quoteStart == std::string::npos) return false;

    size_t quoteEnd = after.find('"', quoteStart + 1);
    if (quoteEnd == std::string::npos) return false;

    value = after.substr(quoteStart + 1, quoteEnd - quoteStart - 1);
    return true;
}

// Extract an integer value after a key
bool extractInt(const std::string& line, const std::string& key, int& value) {
    double dval;
    if (extractDouble(line, key, dval)) {
        value = static_cast<int>(dval);
        return true;
    }
    return false;
}

} // anonymous namespace

LithoConfig JSONLithoConfigLoader::load(const std::string& filepath) {
    LithoConfig config;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "JSONLithoConfigLoader: Could not open file: " << filepath << std::endl;
        return config;
    }

    std::string content((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());

    // Check for bottom_depth anywhere in the file
    double bottomDepth = 5000.0;
    if (extractDouble(content, "bottom_depth", bottomDepth)) {
        config.bottomDepth = -bottomDepth;  // Convert to positive-up convention
    }

    // Find the groups array
    size_t groupsStart = content.find("\"groups\"");
    if (groupsStart == std::string::npos) {
        std::cerr << "JSONLithoConfigLoader: No 'groups' array found in config" << std::endl;
        return config;
    }

    // Find opening bracket of the groups array
    size_t arrayStart = content.find('[', groupsStart);
    if (arrayStart == std::string::npos) return config;

    // Find closing bracket
    size_t arrayEnd = content.find(']', arrayStart);
    if (arrayEnd == std::string::npos) return config;

    // Extract the array content
    std::string arrayContent = content.substr(arrayStart + 1, arrayEnd - arrayStart - 1);

    // Split into individual group objects by finding "{" and "}"
    size_t objStart = 0;
    while (true) {
        size_t braceOpen = arrayContent.find('{', objStart);
        if (braceOpen == std::string::npos) break;

        size_t braceClose = arrayContent.find('}', braceOpen);
        if (braceClose == std::string::npos) break;

        std::string objStr = arrayContent.substr(braceOpen, braceClose - braceOpen + 1);

        // Extract id, name, density
        int id = -1;
        std::string name;
        double density = 0.0;

        extractInt(objStr, "id", id);
        extractString(objStr, "name", name);
        extractDouble(objStr, "density", density);

        if (id >= 0 && !name.empty() && density > 0.0) {
            config.groups.push_back({id, name, density});
        } else {
            std::cerr << "JSONLithoConfigLoader: Skipping group with id="
                      << id << ", name='" << name << "', density=" << density << std::endl;
        }

        objStart = braceClose + 1;
    }

    if (config.groups.empty()) {
        std::cerr << "JSONLithoConfigLoader: No valid groups found" << std::endl;
    }

    return config;
}

} // namespace litho_invert

