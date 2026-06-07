#include <litho_invert/io/ini_config.h>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

namespace litho_invert {

namespace {

std::string trim(const std::string& s) {
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])))
        ++start;
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1])))
        --end;
    return s.substr(start, end - start);
}

std::string toLower(const std::string& s) {
    std::string r = s;
    for (auto& c : r) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return r;
}

} // namespace

bool IniConfig::load(const std::string& filePath) {
    std::ifstream file(filePath);
    if (!file.is_open()) return false;

    m_data.clear();
    std::string currentSection;
    int lineNum = 0;

    std::string line;
    while (std::getline(file, line)) {
        ++lineNum;
        std::string trimmed = trim(line);

        // Skip blank lines and comments
        if (trimmed.empty() || trimmed[0] == '#' || trimmed[0] == ';')
            continue;

        // Section header: [name]
        if (trimmed[0] == '[') {
            size_t end = trimmed.find(']');
            if (end == std::string::npos) continue;
            currentSection = trim(trimmed.substr(1, end - 1));
            if (m_data.find(currentSection) == m_data.end()) {
                m_data[currentSection] = {};
            }
            continue;
        }

        // Key = value  or  Key: value
        size_t sep = trimmed.find('=');
        if (sep == std::string::npos) sep = trimmed.find(':');
        if (sep == std::string::npos) continue;

        std::string key = trim(trimmed.substr(0, sep));
        std::string value = trim(trimmed.substr(sep + 1));

        if (!key.empty() && !currentSection.empty()) {
            m_data[currentSection][key] = value;
        }
    }

    return true;
}

void IniConfig::save(const std::string& filePath) const {
    std::ofstream file(filePath);
    if (!file.is_open()) return;

    for (const auto& section : m_data) {
        file << "[" << section.first << "]\n";
        for (const auto& kv : section.second) {
            file << kv.first << " = " << kv.second << "\n";
        }
        file << "\n";
    }
}

bool IniConfig::hasSection(const std::string& section) const {
    return m_data.find(section) != m_data.end();
}

std::vector<std::string> IniConfig::sections() const {
    std::vector<std::string> result;
    for (const auto& s : m_data) result.push_back(s.first);
    return result;
}

std::string IniConfig::getString(const std::string& section, const std::string& key,
                                  const std::string& defaultValue) const {
    auto sit = m_data.find(section);
    if (sit == m_data.end()) return defaultValue;
    auto kit = sit->second.find(key);
    if (kit == sit->second.end()) return defaultValue;
    return kit->second;
}

double IniConfig::getDouble(const std::string& section, const std::string& key,
                             double defaultValue) const {
    std::string s = getString(section, key, "");
    if (s.empty()) return defaultValue;
    try {
        return std::stod(s);
    } catch (...) {
        return defaultValue;
    }
}

int IniConfig::getInt(const std::string& section, const std::string& key,
                       int defaultValue) const {
    std::string s = getString(section, key, "");
    if (s.empty()) return defaultValue;
    try {
        return std::stoi(s);
    } catch (...) {
        return defaultValue;
    }
}

bool IniConfig::getBool(const std::string& section, const std::string& key,
                         bool defaultValue) const {
    std::string s = toLower(getString(section, key, ""));
    if (s.empty()) return defaultValue;
    return (s == "true" || s == "yes" || s == "1" || s == "on");
}

void IniConfig::setString(const std::string& section, const std::string& key,
                           const std::string& value) {
    m_data[section][key] = value;
}

void IniConfig::setDouble(const std::string& section, const std::string& key, double value) {
    m_data[section][key] = std::to_string(value);
}

void IniConfig::setInt(const std::string& section, const std::string& key, int value) {
    m_data[section][key] = std::to_string(value);
}

void IniConfig::setBool(const std::string& section, const std::string& key, bool value) {
    m_data[section][key] = value ? "true" : "false";
}

} // namespace litho_invert

