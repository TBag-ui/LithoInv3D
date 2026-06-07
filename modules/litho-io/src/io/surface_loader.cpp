#include <litho_invert/io/surface_loader.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>

namespace litho_invert {

std::shared_ptr<SurfaceMesh> OBJSurfaceLoader::load(const std::string& filepath) {
    auto mesh = std::make_shared<SurfaceMesh>();

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "OBJSurfaceLoader: Could not open file: " << filepath << std::endl;
        return nullptr;
    }

    std::string line;
    int lineNum = 0;
    while (std::getline(file, line)) {
        ++lineNum;

        // Trim leading whitespace
        size_t start = line.find_first_not_of(" \t\r");
        if (start == std::string::npos) continue;
        line = line.substr(start);

        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        std::string keyword;
        ss >> keyword;

        if (keyword == "v") {
            // Vertex: v x y z
            double x, y, z;
            if (!(ss >> x >> y >> z)) {
                std::cerr << "OBJSurfaceLoader: Invalid vertex at line " << lineNum << std::endl;
                continue;
            }
            mesh->addVertex(Vector3d(x, y, z));
        } else if (keyword == "f" || keyword == "fo") {
            // Face: f v1 v2 v3  (triangles only, 1-based indices)
            std::vector<int> indices;
            std::string token;
            while (ss >> token) {
                // Handle "v1/vt1/vn1" format - extract the vertex index
                size_t slash = token.find('/');
                if (slash != std::string::npos) {
                    token = token.substr(0, slash);
                }
                try {
                    int idx = std::stoi(token);
                    indices.push_back(idx);
                } catch (...) {
                    std::cerr << "OBJSurfaceLoader: Invalid face index at line "
                              << lineNum << std::endl;
                    indices.clear();
                    break;
                }
            }

            if (indices.size() < 3) {
                std::cerr << "OBJSurfaceLoader: Face at line " << lineNum
                          << " has fewer than 3 vertices" << std::endl;
                continue;
            }

            // OBJ indices are 1-based; convert to 0-based
            for (int& idx : indices) {
                idx -= 1;
            }

            // Triangulate: convert polygons to triangles using fan triangulation
            for (size_t i = 1; i + 1 < indices.size(); ++i) {
                mesh->addTriangle(
                    static_cast<uint32_t>(indices[0]),
                    static_cast<uint32_t>(indices[i]),
                    static_cast<uint32_t>(indices[i + 1]));
            }
        }
        // Ignore other keywords (vn, vt, etc.)
    }

    if (mesh->vertexCount() == 0) {
        std::cerr << "OBJSurfaceLoader: No vertices found in " << filepath << std::endl;
        return nullptr;
    }

    return mesh;
}

std::shared_ptr<SurfaceMesh> TSurfaceLoader::load(const std::string& filepath) {
    auto mesh = std::make_shared<SurfaceMesh>();

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "TSurfaceLoader: Could not open file: " << filepath << std::endl;
        return nullptr;
    }

    std::string line;
    int lineNum = 0;
    bool inHeader = false;

    while (std::getline(file, line)) {
        ++lineNum;

        size_t start = line.find_first_not_of(" \t\r");
        if (start == std::string::npos) continue;
        line = line.substr(start);

        if (line.empty() || line[0] == '#') continue;

        if (inHeader) {
            if (line[0] == '}') inHeader = false;
            continue;
        }

        if (line.find("HEADER") != std::string::npos && line.find('{') != std::string::npos) {
            inHeader = true;
            continue;
        }

        if (line == "END") break;

        if (line == "TFACE" || line == "TRIANGLES") continue;
        if (line.rfind("GOCAD", 0) == 0) continue;

        std::stringstream ss(line);
        std::string keyword;
        ss >> keyword;

        if (keyword == "VRTX" || keyword == "PVRTX") {
            int idx;
            double x, y, z;
            if (!(ss >> idx >> x >> y >> z)) {
                std::cerr << "TSurfaceLoader: Invalid vertex at line " << lineNum << std::endl;
                continue;
            }
            mesh->addVertex(Vector3d(x, y, z));
        } else if (keyword == "TRGL") {
            int v1, v2, v3;
            if (!(ss >> v1 >> v2 >> v3)) {
                std::cerr << "TSurfaceLoader: Invalid triangle at line " << lineNum << std::endl;
                continue;
            }
            mesh->addTriangle(
                static_cast<uint32_t>(v1 - 1),
                static_cast<uint32_t>(v2 - 1),
                static_cast<uint32_t>(v3 - 1));
        }
    }

    if (mesh->vertexCount() == 0) {
        std::cerr << "TSurfaceLoader: No vertices found in " << filepath << std::endl;
        return nullptr;
    }

    return mesh;
}

std::shared_ptr<SurfaceMesh> loadSurfaceMesh(const std::string& filepath) {
    std::string ext;
    size_t dot = filepath.rfind('.');
    if (dot != std::string::npos) {
        ext = filepath.substr(dot);
        for (auto& c : ext) c = static_cast<char>(std::tolower(c));
    }

    if (ext == ".ts") {
        return TSurfaceLoader::load(filepath);
    } else if (ext == ".obj") {
        return OBJSurfaceLoader::load(filepath);
    }

    std::cerr << "loadSurfaceMesh: Unrecognized extension '" << ext
              << "' for file: " << filepath << std::endl;
    return nullptr;
}

} // namespace litho_invert

