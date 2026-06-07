#include <catch2/catch2.hpp>
#include <litho_invert/litho/lithology_model.h>

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

    mesh->addTriangle(0, 1, 2); mesh->addTriangle(0, 2, 3); // +Z top
    mesh->addTriangle(4, 7, 6); mesh->addTriangle(4, 6, 5); // -Z bottom
    mesh->addTriangle(3, 2, 6); mesh->addTriangle(3, 6, 7); // +Y front
    mesh->addTriangle(0, 4, 5); mesh->addTriangle(0, 5, 1); // -Y back
    mesh->addTriangle(1, 5, 6); mesh->addTriangle(1, 6, 2); // +X right
    mesh->addTriangle(0, 3, 7); mesh->addTriangle(0, 7, 4); // -X left

    return mesh;
}

TEST_CASE("applyModelBounds sets per-mesh depth bounds", "[model][apply_bounds]") {
    LithologyModel model;

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::XYZ_FREE);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model.addGroup(g0);
    model.setGroupMesh(0, mesh);

    // Z bounds initially set to default (-10000, 0) by SurfaceMesh
    model.applyModelBounds(1.0);

    // After applyModelBounds, mesh Z bounds should match the model AABB + margin
    // AABB: z ∈ [-20, 0]
    // Margin ≈ max(10, 10, 20) * 0.01 * 1.0 = 0.2
    REQUIRE(mesh->minDepth() < -19.0);
    REQUIRE(mesh->maxDepth() > 0.0);
}

TEST_CASE("applyModelBounds sets X/Y bounds in global DOFs", "[model][apply_bounds]") {
    LithologyModel model;

    auto mesh0 = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::XYZ_FREE);
    auto mesh1 = makeBoxMesh(0.0, 0.0, -20.0, 10.0, 10.0, -50.0, VertexFreedom::XYZ_FREE);

    LithoGroup g0; g0.id = 0; g0.name = "upper"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "lower"; g1.density = 3.00;

    model.addGroup(g0);
    model.addGroup(g1);
    model.setGroupMesh(0, mesh0);
    model.setGroupMesh(1, mesh1);

    model.applyModelBounds(1.0);

    VectorXd lower, upper;
    model.getBounds(lower, upper);

    // X/Y bounds should be finite (not ±1e12) after applyModelBounds
    for (Eigen::Index i = 0; i < lower.size(); ++i) {
        REQUIRE(lower[i] > -1e10);
        REQUIRE(upper[i] < 1e10);
    }

    // Bounds should be within a reasonable range
    // margin = max(dimX=10, dimY=10, dimZ=50) * 0.01 * 1.0 = 0.5
    // X bounds: [0-0.5, 10+0.5] = [-0.5, 10.5]
    // Z bounds for depth axes: [-50.5, 0.5]
    for (Eigen::Index i = 0; i < lower.size(); ++i) {
        REQUIRE(lower[i] >= -60.0); // Z minimum ~ -50 - margin, with slack
        REQUIRE(upper[i] <= 11.0);  // X/Y upper bound ~10.5, Z upper bound ~0.5
    }
}

TEST_CASE("applyModelBounds with margin multiplier", "[model][apply_bounds]") {
    LithologyModel model;

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 100.0, 100.0, -200.0, VertexFreedom::XYZ_FREE);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model.addGroup(g0);
    model.setGroupMesh(0, mesh);

    // Large margin multiplier
    model.applyModelBounds(5.0);

    // With marginMultiplier=5, margin should be ~ 200 * 0.01 * 5 = 10m
    // Z bounds should extend by roughly 10m beyond the model bbox
    REQUIRE(mesh->minDepth() < -205.0); // -200 - margin
    REQUIRE(mesh->maxDepth() > 5.0);    // 0 + margin
}

TEST_CASE("applyModelBounds is safe on empty model", "[model][apply_bounds]") {
    LithologyModel model;
    REQUIRE_NOTHROW(model.applyModelBounds(1.0));
}

TEST_CASE("applyModelBounds handles model with zero extent", "[model][apply_bounds]") {
    LithologyModel model;

    // Degenerate: all vertices at same position
    auto mesh = std::make_shared<SurfaceMesh>();
    mesh->addVertex(0.0, 0.0, 0.0, VertexFreedom::XYZ_FREE);
    mesh->addVertex(0.0, 0.0, 0.0, VertexFreedom::XYZ_FREE);
    mesh->addVertex(0.0, 0.0, 0.0, VertexFreedom::XYZ_FREE);
    mesh->addTriangle(0, 1, 2);

    LithoGroup g0; g0.id = 0; g0.name = "point"; g0.density = 2.67;
    model.addGroup(g0);
    model.setGroupMesh(0, mesh);

    REQUIRE_NOTHROW(model.applyModelBounds(1.0));
    // Margin should be at least 0.01 (1.0 * 0.01 * 1.0)
}

} // namespace
} // namespace litho_invert
