#include <catch2/catch2.hpp>
#include <litho_invert/io/exporters.h>
#include <litho_invert/litho/lithology_model.h>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <set>
#include <algorithm>

namespace litho_invert {
namespace {

static std::shared_ptr<SurfaceMesh> makeBoxMesh(
    double x0, double y0, double z0,
    double x1, double y1, double z1)
{
    auto mesh = std::make_shared<SurfaceMesh>();
    mesh->addVertex(x0, y0, z0); mesh->addVertex(x1, y0, z0);
    mesh->addVertex(x1, y1, z0); mesh->addVertex(x0, y1, z0);
    mesh->addVertex(x0, y0, z1); mesh->addVertex(x1, y0, z1);
    mesh->addVertex(x1, y1, z1); mesh->addVertex(x0, y1, z1);

    mesh->addTriangle(0, 1, 2); mesh->addTriangle(0, 2, 3);
    mesh->addTriangle(4, 7, 6); mesh->addTriangle(4, 6, 5);
    mesh->addTriangle(3, 2, 6); mesh->addTriangle(3, 6, 7);
    mesh->addTriangle(0, 4, 5); mesh->addTriangle(0, 5, 1);
    mesh->addTriangle(1, 5, 6); mesh->addTriangle(1, 6, 2);
    mesh->addTriangle(0, 3, 7); mesh->addTriangle(0, 7, 4);

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
    auto dir = std::filesystem::temp_directory_path() / "litho_io_dedup_tests";
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    return dir.string();
}

static size_t countFiles(const std::string& dir, const std::string& ext) {
    size_t n = 0;
    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (entry.path().extension() == ("." + ext)) ++n;
    }
    return n;
}

TEST_CASE("exportTS produces consistent GOCAD format", "[io][dedup]") {
    // All .ts output must use the same GOCAD TSurf format:
    // GOCAD TSurf 1 header, TFACE keyword, scientific precision(8), 1-based indexing, END footer
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeBoxMesh(0, 0, 0, 10, 10, -20);

    // exportTS via direct call
    exp.exportTS(*mesh, "direct");

    // exportTS via exportStartingModel
    LithologyModel model;
    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model.addGroup(g0);
    model.setGroupMesh(0, mesh);

    exp.setSubfolder("starting");
    exp.exportStartingModel(model);

    // Both files should have the same format markers
    std::string direct = readFile(dir + "/test_direct.ts");
    std::string starting = readFile(dir + "/starting/test_group_0_body.ts");

    REQUIRE(!direct.empty());
    REQUIRE(!starting.empty());

    // Same header
    REQUIRE(direct.find("GOCAD TSurf 1") != std::string::npos);
    REQUIRE(starting.find("GOCAD TSurf 1") != std::string::npos);

    // Same keyword (TFACE, not TRIANGLES)
    REQUIRE(direct.find("TFACE") != std::string::npos);
    REQUIRE(starting.find("TFACE") != std::string::npos);
    REQUIRE(direct.find("TRIANGLES") == std::string::npos);
    REQUIRE(starting.find("TRIANGLES") == std::string::npos);

    // Same footer
    REQUIRE(direct.find("END") != std::string::npos);
    REQUIRE(starting.find("END") != std::string::npos);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportStartingModel calls exportTS once per group", "[io][dedup]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    LithologyModel model;

    auto mesh0 = makeBoxMesh(0, 0, 0, 10, 10, -20);
    auto mesh1 = makeBoxMesh(0, 0, -20, 10, 10, -50);

    LithoGroup g0; g0.id = 0; g0.name = "upper"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "lower"; g1.density = 3.00;

    model.addGroup(g0);
    model.addGroup(g1);
    model.setGroupMesh(0, mesh0);
    model.setGroupMesh(1, mesh1);

    exp.setSubfolder("start");
    exp.exportStartingModel(model);

    // exportStartingModel exports TS + OBJ + closed volume for each group
    // So for 2 groups: 2 TS files, 2 OBJ files, 2 closed volume TS files = 6 total
    // But we care about .ts files being exactly 1 per group (the direct TS, plus closed volume)
    // Actually: exportTS + exportOBJ + exportClosedVolume per group
    // exportClosedVolume also writes .ts, so 2 .ts per group

    // Check that the direct TS files exist (one per group)
    bool hasUpperTS = std::filesystem::exists(dir + "/start/test_group_0_upper.ts");
    bool hasLowerTS = std::filesystem::exists(dir + "/start/test_group_1_lower.ts");
    REQUIRE(hasUpperTS);
    REQUIRE(hasLowerTS);

    // Verify each .ts file has valid format
    std::string upperTS = readFile(dir + "/start/test_group_0_upper.ts");
    std::string lowerTS = readFile(dir + "/start/test_group_1_lower.ts");

    REQUIRE(!upperTS.empty());
    REQUIRE(!lowerTS.empty());
    REQUIRE(upperTS.find("GOCAD TSurf 1") != std::string::npos);
    REQUIRE(lowerTS.find("GOCAD TSurf 1") != std::string::npos);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportStartingModel handles missing meshes gracefully", "[io][dedup]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    LithologyModel model;

    LithoGroup g0; g0.id = 0; g0.name = "empty"; g0.density = 2.67;
    model.addGroup(g0);
    // Don't set a mesh for group 0

    exp.setSubfolder("missing");
    REQUIRE_NOTHROW(exp.exportStartingModel(model));

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportAll produces expected .ts files", "[io][dedup]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeBoxMesh(0, 0, 0, 10, 10, -20);

    LithologyModel model;
    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model.addGroup(g0);
    model.setGroupMesh(0, mesh);

    GravityData observed;
    observed.push_back(GravityPoint(Vector3d(5, 5, 0), 0.5, 0.01));

    InversionResult result;
    result.finalModel = std::make_shared<LithologyModel>(model);
    result.predictedData = VectorXd::Ones(1);
    result.converged = true;
    result.totalIterations = 10;
    result.finalMisfit = 0.01;
    result.finalRMS = 0.5;

    exp.setSubfolder("final");
    exp.exportAll(result, observed, 0, 10, 0, 10, -20, 0, 50.0);

    // Check that group TS file exists
    bool hasGroupTS = std::filesystem::exists(dir + "/final/test_group_0_body.ts");
    REQUIRE(hasGroupTS);

    // Check format
    std::string ts = readFile(dir + "/final/test_group_0_body.ts");
    REQUIRE(!ts.empty());
    REQUIRE(ts.find("GOCAD TSurf 1") != std::string::npos);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportAll with group naming uses provided names", "[io][dedup]") {
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh0 = makeBoxMesh(0, 0, 0, 10, 10, -20);
    auto mesh1 = makeBoxMesh(0, 0, -20, 10, 10, -50);

    LithologyModel model;
    LithoGroup g0; g0.id = 0; g0.name = "Granite"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "Sulfide"; g1.density = 3.00;
    model.addGroup(g0);
    model.addGroup(g1);
    model.setGroupMesh(0, mesh0);
    model.setGroupMesh(1, mesh1);

    exp.setGroupNaming({"Granite", "Sulfide"});

    GravityData observed;
    observed.push_back(GravityPoint(Vector3d(5, 5, 0), 0.5, 0.01));

    InversionResult result;
    result.finalModel = std::make_shared<LithologyModel>(model);
    result.predictedData = VectorXd::Ones(1);

    exp.setSubfolder("named");
    exp.exportAll(result, observed, 0, 10, 0, 10, -50, 0, 50.0);

    // With group naming, files use provided names, not group_N_groupname
    bool hasNamed1 = std::filesystem::exists(dir + "/named/test_Granite.ts");
    bool hasNamed2 = std::filesystem::exists(dir + "/named/test_Sulfide.ts");
    REQUIRE(hasNamed1);
    REQUIRE(hasNamed2);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS format is consistent across all export paths", "[io][dedup]") {
    // Verify that no matter how .ts files are generated, they use the
    // same consistent format — no alternate implementations slipping through.
    auto dir = makeTempDir();

    auto mesh = makeBoxMesh(0, 0, 0, 10, 10, -20);

    // Path 1: Direct exportTS
    {
        InversionExporter exp(dir, "test");
        exp.exportTS(*mesh, "path1");
    }

    // Path 2: Via exportStartingModel
    {
        InversionExporter exp(dir, "test");
        LithologyModel model;
        LithoGroup g; g.id = 0; g.name = "g"; g.density = 2.67;
        model.addGroup(g);
        model.setGroupMesh(0, mesh);
        exp.setSubfolder("path2");
        exp.exportStartingModel(model);
    }

    // Path 3: Via exportClosedVolume
    {
        InversionExporter exp(dir, "test");
        exp.exportClosedVolume(*mesh, "path3");
    }

    // Check all paths produce files with same structural markers
    std::vector<std::string> files = {
        dir + "/test_path1.ts",
        dir + "/path2/test_group_0_g.ts",
        dir + "/test_path3.ts"
    };

    for (const auto& f : files) {
        std::string content = readFile(f);
        // All must be valid GOCAD TSurf
        REQUIRE(!content.empty());
        REQUIRE(content.find("GOCAD TSurf 1") != std::string::npos);
        REQUIRE(content.find("TFACE") != std::string::npos);
        REQUIRE(content.find("END") != std::string::npos);
        // No file should have null bytes
        for (char c : content) {
            REQUIRE(c != '\0');
        }
    }

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

TEST_CASE("exportTS is the canonical .ts writer — no TRIANGLES keyword", "[io][dedup]") {
    // The DUPLICATE_TRACKING.txt identifies 4 .ts implementations.
    // The authoritative one (InversionExporter::exportTS) uses "TFACE".
    // The deprecated copies in Inversion_Cluster_API/ use "TRIANGLES".
    // This test ensures the modules/ path always uses TFACE.
    auto dir = makeTempDir();
    InversionExporter exp(dir, "test");

    auto mesh = makeBoxMesh(0, 0, 0, 10, 10, -20);
    exp.exportTS(*mesh, "keyword_check");

    std::string content = readFile(dir + "/test_keyword_check.ts");
    REQUIRE(content.find("TFACE") != std::string::npos);
    REQUIRE(content.find("TRIANGLES") == std::string::npos);

    std::error_code ec;
    std::filesystem::remove_all(dir, ec);
}

} // namespace
} // namespace litho_invert
