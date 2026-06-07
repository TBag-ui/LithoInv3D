#include <catch2/catch2.hpp>
#include <litho_invert/inversion/runner.h>
#include <litho_invert/forward/gravity_forward.h>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <chrono>
#include <thread>

namespace litho_invert {
namespace {

static std::shared_ptr<SurfaceMesh> makeBoxMesh(
    double x0, double y0, double z0,
    double x1, double y1, double z1,
    VertexFreedom f = VertexFreedom::XYZ_FREE)
{
    auto mesh = std::make_shared<SurfaceMesh>();
    mesh->addVertex(x0, y0, z0, f); mesh->addVertex(x1, y0, z0, f);
    mesh->addVertex(x1, y1, z0, f); mesh->addVertex(x0, y1, z0, f);
    mesh->addVertex(x0, y0, z1, f); mesh->addVertex(x1, y0, z1, f);
    mesh->addVertex(x1, y1, z1, f); mesh->addVertex(x0, y1, z1, f);

    mesh->addTriangle(0, 1, 2); mesh->addTriangle(0, 2, 3);
    mesh->addTriangle(4, 7, 6); mesh->addTriangle(4, 6, 5);
    mesh->addTriangle(3, 2, 6); mesh->addTriangle(3, 6, 7);
    mesh->addTriangle(0, 4, 5); mesh->addTriangle(0, 5, 1);
    mesh->addTriangle(1, 5, 6); mesh->addTriangle(1, 6, 2);
    mesh->addTriangle(0, 3, 7); mesh->addTriangle(0, 7, 4);

    return mesh;
}

static GravityData makeGridObs(double x0, double x1, int nx,
                                double y0, double y1, int ny, double z) {
    GravityData data;
    double dx = (nx > 1) ? (x1 - x0) / (nx - 1) : 0.0;
    double dy = (ny > 1) ? (y1 - y0) / (ny - 1) : 0.0;
    for (int iy = 0; iy < ny; ++iy)
        for (int ix = 0; ix < nx; ++ix)
            data.push_back(GravityPoint(Vector3d(x0 + ix * dx, y0 + iy * dy, z), 0.0, 0.01));
    return data;
}

static std::string makeTempDir() {
    auto dir = std::filesystem::temp_directory_path() / "litho_inv_iter_ts_tests";
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    return dir.string();
}

static std::string readFile(const std::string& path) {
    std::ifstream in(path);
    if (!in.is_open()) return "";
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

TEST_CASE("iterationExportDir triggers per-iteration .ts export", "[inversion][iter_export]") {
    auto dir = makeTempDir();
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;
    config.iterationExportDir = dir;

    InversionRunner runner(config);
    runner.run();

    // Check that .ts files exist for each group
    bool hasGroup0 = std::filesystem::exists(dir + "/iter_group_0.ts");
    REQUIRE(hasGroup0);

    // Verify file content
    std::string ts0 = readFile(dir + "/iter_group_0.ts");
    REQUIRE(!ts0.empty());
    REQUIRE(ts0.find("GOCAD TSurf 1") != std::string::npos);
    REQUIRE(ts0.find("TFACE") != std::string::npos);
    REQUIRE(ts0.find("END") != std::string::npos);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("iteration .ts files have correct vertex count", "[inversion][iter_export]") {
    auto dir = makeTempDir();
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);
    uint32_t expectedVerts = mesh->vertexCount();
    uint32_t expectedTris = mesh->triangleCount();

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;
    config.iterationExportDir = dir;

    InversionRunner runner(config);
    runner.run();

    std::string content = readFile(dir + "/iter_group_0.ts");
    REQUIRE(!content.empty());

    // Count VRTX and TRGL lines
    size_t vrtxLines = 0, trglLines = 0;
    std::istringstream stream(content);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.size() >= 4 && line.substr(0, 4) == "VRTX") ++vrtxLines;
        if (line.size() >= 4 && line.substr(0, 4) == "TRGL") ++trglLines;
    }
    REQUIRE(vrtxLines == expectedVerts);
    REQUIRE(trglLines == expectedTris);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("iteration export uses group naming when provided", "[inversion][iter_export]") {
    auto dir = makeTempDir();
    auto model = std::make_shared<LithologyModel>();

    auto mesh0 = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);
    auto mesh1 = makeBoxMesh(0.0, 0.0, -20.0, 10.0, 10.0, -50.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "Granite"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "Sulfide"; g1.density = 3.00;
    model->addGroup(g0);
    model->addGroup(g1);
    model->setGroupMesh(0, mesh0);
    model->setGroupMesh(1, mesh1);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;
    config.iterationExportDir = dir;
    config.groupExportNames = {"Granite", "Sulfide"};

    InversionRunner runner(config);
    runner.run();

    // Files should use the provided names
    bool hasGranite = std::filesystem::exists(dir + "/iter_Granite.ts");
    bool hasSulfide = std::filesystem::exists(dir + "/iter_Sulfide.ts");
    REQUIRE(hasGranite);
    REQUIRE(hasSulfide);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("iteration export directory is created if missing", "[inversion][iter_export]") {
    auto baseDir = makeTempDir();
    auto dir = baseDir + "/nested/subdir/export";

    // Directory should not exist yet
    REQUIRE_FALSE(std::filesystem::exists(dir));

    auto model = std::make_shared<LithologyModel>();
    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);
    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;
    config.iterationExportDir = dir;

    InversionRunner runner(config);
    runner.run();

    // Directory should now exist
    REQUIRE(std::filesystem::exists(dir));
    REQUIRE(std::filesystem::exists(dir + "/iter_group_0.ts"));

    std::error_code ec;
    std::filesystem::remove_all(baseDir, ec);
}

TEST_CASE("empty iterationExportDir does not export", "[inversion][iter_export]") {
    auto dir = makeTempDir();
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);
    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;
    // iterationExportDir left empty — no export should happen
    config.iterationExportDir = "";

    REQUIRE_NOTHROW([](const InversionConfig& c) { InversionRunner runner(c); runner.run(); }(config));

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("iteration export file content is valid GOCAD TSurf", "[inversion][iter_export]") {
    auto dir = makeTempDir();
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);
    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;
    config.iterationExportDir = dir;

    InversionRunner runner(config);
    runner.run();

    std::string content = readFile(dir + "/iter_group_0.ts");

    // Full structural validation
    REQUIRE(content.find("GOCAD TSurf 1\n") != std::string::npos);
    REQUIRE(content.find("HEADER {\n") != std::string::npos);
    REQUIRE(content.find("name:") != std::string::npos);
    REQUIRE(content.find("}\n") != std::string::npos);
    REQUIRE(content.find("TFACE\n") != std::string::npos);
    REQUIRE(content.find("END\n") != std::string::npos);

    // No embedded nulls
    for (char c : content) {
        REQUIRE(c != '\0');
    }

    // All lines are printable ASCII or empty
    std::istringstream stream(content);
    std::string line;
    while (std::getline(stream, line)) {
        for (char c : line) {
            bool ok = (c >= 0x20 && c <= 0x7E) || c == '\r';
            REQUIRE(ok);
        }
    }

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

} // namespace
} // namespace litho_invert
