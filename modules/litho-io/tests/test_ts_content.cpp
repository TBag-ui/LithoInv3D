#include <catch2/catch2.hpp>
#include <litho_invert/io/exporters.h>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cstring>

namespace litho_invert {
namespace {

static std::shared_ptr<SurfaceMesh> makeSimpleMesh() {
    auto mesh = std::make_shared<SurfaceMesh>();
    mesh->addVertex(0.0, 0.0, -10.0);
    mesh->addVertex(1.0, 0.0, -12.0);
    mesh->addVertex(1.0, 1.0, -11.0);
    mesh->addVertex(0.0, 1.0, -13.0);
    mesh->addTriangle(0, 1, 2);
    mesh->addTriangle(0, 2, 3);
    return mesh;
}

static std::string readFile(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) return "";
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

static std::string makeTempDir() {
    auto dir = std::filesystem::temp_directory_path() / "litho_io_ts_tests";
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    return dir.string();
}

TEST_CASE("exportTS produces valid GOCAD TSurf header", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "header_test");

    std::string content = readFile(dir + "/test_header_test.ts");
    REQUIRE(!content.empty());

    // GOCAD TSurf header
    REQUIRE(content.find("GOCAD TSurf 1") != std::string::npos);
    REQUIRE(content.find("HEADER {") != std::string::npos);
    REQUIRE(content.find("name:") != std::string::npos);
    REQUIRE(content.find("TFACE") != std::string::npos);
    REQUIRE(content.find("END") != std::string::npos);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS vertex count matches mesh", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "vcount");

    std::string content = readFile(dir + "/test_vcount.ts");

    // Count VRTX lines (case-insensitive)
    size_t vrtxCount = 0;
    std::istringstream stream(content);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.size() >= 4 && line.substr(0, 4) == "VRTX") ++vrtxCount;
    }
    REQUIRE(vrtxCount == mesh->vertexCount());

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS triangle count matches mesh", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "tcount");

    std::string content = readFile(dir + "/test_tcount.ts");

    // Count TRGL lines
    size_t trglCount = 0;
    std::istringstream stream(content);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.size() >= 4 && line.substr(0, 4) == "TRGL") ++trglCount;
    }
    REQUIRE(trglCount == mesh->triangleCount());

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS uses 1-based indexing", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "indexing");

    std::string content = readFile(dir + "/test_indexing.ts");

    // First vertex should be VRTX 1 (1-based)
    REQUIRE(content.find("VRTX 1 ") != std::string::npos);
    // Should NOT have VRTX 0
    REQUIRE(content.find("VRTX 0 ") == std::string::npos);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS vertex positions match mesh data", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "positions");

    std::string content = readFile(dir + "/test_positions.ts");

    // Verify first vertex position — format: "VRTX 1 X Y Z" in scientific notation
    const auto& v0 = mesh->vertex(0).position;
    std::ostringstream expected;
    expected << "VRTX 1 "
             << std::scientific << std::setprecision(8)
             << v0.x() << " " << v0.y() << " " << v0.z();
    REQUIRE(content.find(expected.str()) != std::string::npos);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS scientific precision is 8 decimal places", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = std::make_shared<SurfaceMesh>();
    mesh->addVertex(123.456789012345, 0.0, -50.0);
    mesh->addVertex(0.0, 0.0, 0.0);
    mesh->addVertex(1.0, 1.0, 1.0);
    mesh->addTriangle(0, 1, 2);

    exp.exportTS(*mesh, "precision");
    std::string content = readFile(dir + "/test_precision.ts");

    // Scientific with 8 decimal places: the coordinate should appear with 8 digits
    // after decimal point in scientific notation, e.g. "1.23456789e+02"
    std::istringstream stream(content);
    std::string line;
    bool foundPreciseCoord = false;
    while (std::getline(stream, line)) {
        if (line.find("VRTX 1") != std::string::npos) {
            size_t dotPos = line.find('.');
            REQUIRE(dotPos != std::string::npos);
            std::string afterDot = line.substr(dotPos + 1);
            // Count digits before 'e'
            size_t ePos = afterDot.find('e');
            if (ePos == std::string::npos) ePos = afterDot.find('E');
            REQUIRE(ePos != std::string::npos);
            std::string decimalDigits = afterDot.substr(0, ePos);
            REQUIRE(decimalDigits.size() == 8);
            foundPreciseCoord = true;
        }
    }
    REQUIRE(foundPreciseCoord);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS has no embedded nulls", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "nonulls");

    std::string content = readFile(dir + "/test_nonulls.ts");
    REQUIRE(!content.empty());

    for (char c : content) {
        REQUIRE(c != '\0');
    }

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS with empty mesh", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = std::make_shared<SurfaceMesh>();
    exp.exportTS(*mesh, "empty");

    std::string content = readFile(dir + "/test_empty.ts");
    REQUIRE(!content.empty());

    // Should still have header and END
    REQUIRE(content.find("GOCAD TSurf 1") != std::string::npos);
    REQUIRE(content.find("TFACE") != std::string::npos);
    REQUIRE(content.find("END") != std::string::npos);

    // Should have zero VRTX and zero TRGL lines
    std::istringstream stream(content);
    std::string line;
    size_t vrtxCount = 0, trglCount = 0;
    while (std::getline(stream, line)) {
        if (line.size() >= 4 && line.substr(0, 4) == "VRTX") ++vrtxCount;
        if (line.size() >= 4 && line.substr(0, 4) == "TRGL") ++trglCount;
    }
    REQUIRE(vrtxCount == 0);
    REQUIRE(trglCount == 0);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS output is pure ASCII text", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "ascii");

    std::string content = readFile(dir + "/test_ascii.ts");
    REQUIRE(!content.empty());

    for (char c : content) {
        // Allow printable ASCII, newline, carriage return
        bool valid = (c >= 0x20 && c <= 0x7E) || c == '\n' || c == '\r';
        REQUIRE(valid);
    }

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS filename follows naming convention", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "mybase");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "mysuffix");

    // File should be at dir/mybase_mysuffix.ts
    bool exists = std::filesystem::exists(dir + "/mybase_mysuffix.ts");
    REQUIRE(exists);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS with subfolder", "[io][ts]") {
    auto dir = makeTempDir();
    // Create subfolder so the exporter can write into it
    std::filesystem::create_directory(dir + "/sub");

    InversionExporter exp(dir, "test");
    exp.setSubfolder("sub");
    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "subfolder");

    bool exists = std::filesystem::exists(dir + "/sub/test_subfolder.ts");
    REQUIRE(exists);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS idempotent — same mesh twice produces identical files", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "idem");

    std::string first = readFile(dir + "/test_idem.ts");

    // Write again with a new exporter but same mesh
    InversionExporter exp2(dir, "test");
    exp2.exportTS(*mesh, "idem");

    std::string second = readFile(dir + "/test_idem.ts");

    REQUIRE(first == second);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS all vertices are written exactly once", "[io][ts]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeSimpleMesh();
    exp.exportTS(*mesh, "allverts");

    std::string content = readFile(dir + "/test_allverts.ts");

    // Collect vertex indices referenced in TRGL lines
    std::set<size_t> referencedIndices;
    std::istringstream stream(content);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.size() >= 4 && line.substr(0, 4) == "TRGL") {
            std::istringstream trgl(line);
            std::string keyword;
            size_t v1, v2, v3;
            trgl >> keyword >> v1 >> v2 >> v3;
            referencedIndices.insert(v1);
            referencedIndices.insert(v2);
            referencedIndices.insert(v3);
        }
    }

    // Every referenced vertex should be in range [1, vertexCount]
    for (size_t idx : referencedIndices) {
        REQUIRE(idx >= 1);
        REQUIRE(idx <= mesh->vertexCount());
    }

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

} // namespace
} // namespace litho_invert
