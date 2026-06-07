#include <catch2/catch2.hpp>
#include <litho_invert/inversion/runner.h>
#include <litho_invert/forward/gravity_forward.h>

using Catch::Approx;

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

static std::vector<Vector3d> capturePositions(const LithologyModel& model) {
    std::vector<Vector3d> positions;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            positions.push_back(mesh->vertex(vi).position);
        }
    }
    return positions;
}

static bool anyVertexMoved(const std::vector<Vector3d>& before,
                           const std::vector<Vector3d>& after,
                           double threshold = 1e-6) {
    for (size_t i = 0; i < before.size(); ++i) {
        if ((before[i] - after[i]).norm() > threshold) return true;
    }
    return false;
}

TEST_CASE("Vertices move during inversion", "[inversion][movement]") {
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    // Capture positions before inversion
    auto beforePositions = capturePositions(*model);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 3;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    runner.run();

    // Capture positions after inversion
    auto afterPositions = capturePositions(*model);

    REQUIRE(beforePositions.size() == afterPositions.size());
    REQUIRE(afterPositions.size() > 0);

    // At least one vertex should have changed position
    // (unless the model is already a perfect fit for zero observations)
    // With zero observations, the inversion should adjust vertices to reduce misfit
    REQUIRE(anyVertexMoved(beforePositions, afterPositions));
}

TEST_CASE("FIXED vertices do not move", "[inversion][movement]") {
    auto model = std::make_shared<LithologyModel>();

    // Create mesh with a mix of FIXED and Z_ONLY vertices
    auto mesh = std::make_shared<SurfaceMesh>();
    // First 4 vertices FIXED, last 4 Z_ONLY
    mesh->addVertex(0.0, 0.0, 0.0, VertexFreedom::FIXED);
    mesh->addVertex(10.0, 0.0, 0.0, VertexFreedom::FIXED);
    mesh->addVertex(10.0, 10.0, 0.0, VertexFreedom::FIXED);
    mesh->addVertex(0.0, 10.0, 0.0, VertexFreedom::FIXED);
    mesh->addVertex(0.0, 0.0, -20.0, VertexFreedom::Z_ONLY);
    mesh->addVertex(10.0, 0.0, -20.0, VertexFreedom::Z_ONLY);
    mesh->addVertex(10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);
    mesh->addVertex(0.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    mesh->addTriangle(0, 1, 2); mesh->addTriangle(0, 2, 3);
    mesh->addTriangle(4, 7, 6); mesh->addTriangle(4, 6, 5);
    mesh->addTriangle(3, 2, 6); mesh->addTriangle(3, 6, 7);
    mesh->addTriangle(0, 4, 5); mesh->addTriangle(0, 5, 1);
    mesh->addTriangle(1, 5, 6); mesh->addTriangle(1, 6, 2);
    mesh->addTriangle(0, 3, 7); mesh->addTriangle(0, 7, 4);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    // Capture positions of FIXED vertices before inversion
    std::vector<Vector3d> fixedBefore;
    for (uint32_t vi = 0; vi < std::min<uint32_t>(4, mesh->vertexCount()); ++vi) {
        fixedBefore.push_back(mesh->vertex(vi).position);
    }

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 3;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    runner.run();

    // FIXED vertices should not have moved
    for (uint32_t vi = 0; vi < std::min<uint32_t>(4, mesh->vertexCount()); ++vi) {
        const auto& pos = mesh->vertex(vi).position;
        REQUIRE(pos.x() == Approx(fixedBefore[vi].x()).margin(1e-8));
        REQUIRE(pos.y() == Approx(fixedBefore[vi].y()).margin(1e-8));
        REQUIRE(pos.z() == Approx(fixedBefore[vi].z()).margin(1e-8));
    }
}

TEST_CASE("Z_ONLY vertices change only Z coordinate", "[inversion][movement]") {
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    // Capture initial X/Y values
    std::vector<std::pair<double, double>> initialXY;
    for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
        initialXY.push_back({mesh->vertex(vi).position.x(), mesh->vertex(vi).position.y()});
    }

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 3;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    runner.run();

    // X and Y should remain unchanged for all vertices
    for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
        REQUIRE(mesh->vertex(vi).position.x() == Approx(initialXY[vi].first).margin(1e-8));
        REQUIRE(mesh->vertex(vi).position.y() == Approx(initialXY[vi].second).margin(1e-8));
    }
}

TEST_CASE("Vertices do not teleport — movement is bounded", "[inversion][movement]") {
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    auto before = capturePositions(*model);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 3;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    runner.run();

    auto after = capturePositions(*model);

    // No vertex should have moved more than 100m in Z (reasonable bound)
    for (size_t i = 0; i < before.size(); ++i) {
        double dz = std::abs(before[i].z() - after[i].z());
        REQUIRE(dz < 100.0);
    }
}

TEST_CASE("Shared vertices between groups remain synchronized after inversion", "[inversion][movement]") {
    auto model = std::make_shared<LithologyModel>();

    // Two adjacent boxes — shared face at z=-20
    auto mesh0 = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);
    auto mesh1 = makeBoxMesh(0.0, 0.0, -20.0, 10.0, 10.0, -50.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "upper"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "lower"; g1.density = 3.00;
    model->addGroup(g0);
    model->addGroup(g1);
    model->setGroupMesh(0, mesh0);
    model->setGroupMesh(1, mesh1);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 3;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    runner.run();

    // Shared face vertices: mesh0 bottom = {4,5,6,7}, mesh1 top = {0,1,2,3}
    // They should still be at the same positions
    struct Pair { uint32_t m0v, m1v; };
    Pair pairs[] = {{4, 0}, {5, 1}, {6, 2}, {7, 3}};

    for (const auto& [a, b] : pairs) {
        const auto& pa = mesh0->vertex(a).position;
        const auto& pb = mesh1->vertex(b).position;
        REQUIRE(pa.x() == Approx(pb.x()).margin(1e-6));
        REQUIRE(pa.y() == Approx(pb.y()).margin(1e-6));
        REQUIRE(pa.z() == Approx(pb.z()).margin(1e-6));
    }
}

} // namespace
} // namespace litho_invert
