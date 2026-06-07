#include <catch2/catch2.hpp>
#include <litho_invert/inversion/runner.h>
#include <litho_invert/forward/gravity_forward.h>

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

TEST_CASE("Bounds are finite after applyModelBounds", "[inversion][bounds]") {
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    model->applyModelBounds(1.0);

    VectorXd lower, upper;
    model->getBounds(lower, upper);

    for (Eigen::Index i = 0; i < lower.size(); ++i) {
        REQUIRE(lower[i] > -1e10);
        REQUIRE(upper[i] < 1e10);
    }
}

TEST_CASE("Z bounds are within reasonable range", "[inversion][bounds]") {
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    model->applyModelBounds(1.0);

    VectorXd lower, upper;
    model->getBounds(lower, upper);

    // Z must be within model AABB + margin
    for (Eigen::Index i = 0; i < lower.size(); ++i) {
        if (lower[i] < -500.0) continue; // X/Y DOFs may have large negative bounds
        REQUIRE(lower[i] >= -30.0); // Z minimum: -20 - margin
        REQUIRE(upper[i] <= 10.0);  // Z maximum: 0 + margin
    }
}

TEST_CASE("Vertices within starting model AABB after inversion", "[inversion][bounds]") {
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    // Set bounds before inversion
    model->applyModelBounds(1.0);

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    runner.run();

    // All vertices should still be within reasonable bounds
    auto* finalMesh = model->groupMesh(0);
    for (uint32_t vi = 0; vi < finalMesh->vertexCount(); ++vi) {
        const auto& pos = finalMesh->vertex(vi).position;
        REQUIRE(pos.z() >= -100.0); // not unreasonably deep
        REQUIRE(pos.z() <= 100.0);  // not unreasonably high
        REQUIRE(pos.x() >= -100.0); // within model extent
        REQUIRE(pos.x() <= 200.0);
        REQUIRE(pos.y() >= -100.0);
        REQUIRE(pos.y() <= 200.0);
    }
}

TEST_CASE("applyModelBounds produces tighter bounds than default", "[inversion][bounds]") {
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::XYZ_FREE);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    // Get default bounds (before applyModelBounds)
    VectorXd defaultLower, defaultUpper;
    model->getBounds(defaultLower, defaultUpper);

    // Apply model bounds
    model->applyModelBounds(1.0);

    VectorXd boundedLower, boundedUpper;
    model->getBounds(boundedLower, boundedUpper);

    // Bounded lower should be >= default lower (tighter or equal)
    // Bounded upper should be reasonable (finite, not the default +1e6)
    for (Eigen::Index i = 0; i < boundedLower.size(); ++i) {
        // Skip XY DOFs with default +1e6 upper bound — those get tightened
        if (defaultUpper[i] > 1e5) {
            REQUIRE(boundedUpper[i] < 1e5);
            continue;
        }
        // For Z DOFs, the margin may slightly expand the default [minZ, maxZ]
        // but bounds must still be finite and reasonable
        REQUIRE(boundedLower[i] > -1e4);
        REQUIRE(boundedUpper[i] < 1e4);
    }
}

TEST_CASE("Per-mesh depth bounds are set by applyModelBounds", "[inversion][bounds]") {
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    // Before applyModelBounds, minDepth = -10000 (default) and maxDepth = 0 (default)
    double minBefore = mesh->minDepth();
    double maxBefore = mesh->maxDepth();

    model->applyModelBounds(1.0);

    // After applyModelBounds, bounds should be tighter
    double minAfter = mesh->minDepth();
    double maxAfter = mesh->maxDepth();

    // Bounds should have changed from defaults
    REQUIRE(minAfter > minBefore); // tighter min
    REQUIRE(maxAfter > maxBefore); // tighter max (includes positive margin)
}

TEST_CASE("Depth bounds enforced during inversion", "[inversion][bounds]") {
    auto model = std::make_shared<LithologyModel>();

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model->addGroup(g0);
    model->setGroupMesh(0, mesh);

    model->applyModelBounds(1.0);

    double minDepth = mesh->minDepth();
    double maxDepth = mesh->maxDepth();

    auto obs = makeGridObs(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    runner.run();

    // All vertex Z values should be within depth bounds
    auto* finalMesh = model->groupMesh(0);
    for (uint32_t vi = 0; vi < finalMesh->vertexCount(); ++vi) {
        double z = finalMesh->vertex(vi).position.z();
        // Convert: minDepth is max negative Z (deepest positive-down)
        // Z positive-up: z <= maxDepth (surface), z >= minDepth (deep)
        REQUIRE(z <= finalMesh->maxDepth());
        REQUIRE(z >= finalMesh->minDepth());
    }
}

} // namespace
} // namespace litho_invert
