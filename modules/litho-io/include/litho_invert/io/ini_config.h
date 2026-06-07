#pragma once

#include <map>
#include <string>
#include <vector>

namespace litho_invert {

class IniConfig {
public:
    bool load(const std::string& filePath);
    void save(const std::string& filePath) const;

    bool hasSection(const std::string& section) const;
    std::vector<std::string> sections() const;

    std::string getString(const std::string& section, const std::string& key,
                          const std::string& defaultValue = "") const;
    double getDouble(const std::string& section, const std::string& key,
                     double defaultValue = 0.0) const;
    int getInt(const std::string& section, const std::string& key,
               int defaultValue = 0) const;
    bool getBool(const std::string& section, const std::string& key,
                 bool defaultValue = false) const;

    void setString(const std::string& section, const std::string& key,
                   const std::string& value);
    void setDouble(const std::string& section, const std::string& key, double value);
    void setInt(const std::string& section, const std::string& key, int value);
    void setBool(const std::string& section, const std::string& key, bool value);

private:
    std::map<std::string, std::map<std::string, std::string>> m_data;
};

} // namespace litho_invert

