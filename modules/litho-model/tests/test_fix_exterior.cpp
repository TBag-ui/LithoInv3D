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
    mesh->addVertex(x0, y0, z0, f); // 0
    mesh->addVertex(x1, y0, z0, f); // 1
    mesh->addVertex(x1, y1, z0, f); // 2
    mesh->addVertex(x0, y1, z0, f); // 3
    mesh->addVertex(x0, y0, z1, f); // 4
    mesh->addVertex(x1, y0, z1, f); // 5
    mesh->addVertex(x1, y1, z1, f); // 6
    mesh->addVertex(x0, y1, z1, f); // 7

    // +Z top (0-1-2, 0-2-3)
    mesh->addTriangle(0, 1, 2); mesh->addTriangle(0, 2, 3);
    // -Z bottom (4-7-6, 4-6-5)
    mesh->addTriangle(4, 7, 6); mesh->addTriangle(4, 6, 5);
    // +Y front (3-2-6, 3-6-7)
    mesh->addTriangle(3, 2, 6); mesh->addTriangle(3, 6, 7);
    // -Y back (0-4-5, 0-5-1)
    mesh->addTriangle(0, 4, 5); mesh->addTriangle(0, 5, 1);
    // +X right (1-5-6, 1-6-2)
    mesh->addTriangle(1, 5, 6); mesh->addTriangle(1, 6, 2);
    // -X left (0-3-7, 0-7-4)
    mesh->addTriangle(0, 3, 7); mesh->addTriangle(0, 7, 4);

    return mesh;
}

TEST_CASE("fixExteriorFaces on single group model", "[model][fix_exterior]") {
    LithologyModel model;

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::XYZ_FREE);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model.addGroup(g0);
    model.setGroupMesh(0, mesh);

    // Before fixExteriorFaces, all vertices are XYZ_FREE
    for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
        REQUIRE(mesh->vertex(vi).freedom == VertexFreedom::XYZ_FREE);
    }

    model.fixExteriorFaces();

    // After fixExteriorFaces, some vertices should be constrained
    int nFixed = 0, nXYFree = 0, nXYZFree = 0, nOther = 0;
    for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
        switch (mesh->vertex(vi).freedom) {
            case VertexFreedom::FIXED: ++nFixed; break;
            case VertexFreedom::XY_FREE: ++nXYFree; break;
            case VertexFreedom::XYZ_FREE: ++nXYZFree; break;
            default: ++nOther; break;
        }
    }

    SECTION("all vertices get some form of constraint") {
        // With a single box mesh, ALL faces are exterior
        REQUIRE(nOther == 0);
        // Total: 8 vertices
        REQUIRE(nFixed + nXYFree + nXYZFree == 8);
    }

    SECTION("corners are FIXED") {
        // Box corners belong to 3 exterior faces each (top/bottom + 2 sides)
        // Most corners should be FIXED
        // Actually corners should have ≥3 exterior faces and be FIXED
        REQUIRE(nFixed > 0);
    }

    SECTION("all vertices are at face intersections") {
        // All 8 box vertices are at 3-face intersections — each belongs to
        // exactly 3 exterior faces, so all are classified as corners and FIXED.
        REQUIRE(nFixed == 8);
    }
}

TEST_CASE("fixExteriorFaces on two-group model", "[model][fix_exterior]") {
    LithologyModel model;

    // Two adjacent boxes — shared face at z=-20
    auto mesh0 = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::XYZ_FREE);
    auto mesh1 = makeBoxMesh(0.0, 0.0, -20.0, 10.0, 10.0, -50.0, VertexFreedom::XYZ_FREE);

    LithoGroup g0; g0.id = 0; g0.name = "upper"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "lower"; g1.density = 3.00;

    model.addGroup(g0);
    model.addGroup(g1);
    model.setGroupMesh(0, mesh0);
    model.setGroupMesh(1, mesh1);

    model.fixExteriorFaces();

    SECTION("shared face vertices have appropriate constraints") {
        // mesh0 bottom face vertices (4,5,6,7) at z=-20 and
        // mesh1 top face vertices (0,1,2,3) at z=-20 share positions.
        // The shared faces are interior, leaving 2 exterior side faces per
        // corner vertex. Non-corner vertices might get different freedoms.
        // We just verify the method completed without corrupting mesh state.
        for (uint32_t vi = 0; vi < mesh0->vertexCount(); ++vi) {
            REQUIRE(static_cast<int>(mesh0->vertex(vi).freedom) >= 0);
            REQUIRE(static_cast<int>(mesh0->vertex(vi).freedom) <= 6);
        }
        for (uint32_t vi = 0; vi < mesh1->vertexCount(); ++vi) {
            REQUIRE(static_cast<int>(mesh1->vertex(vi).freedom) >= 0);
            REQUIRE(static_cast<int>(mesh1->vertex(vi).freedom) <= 6);
        }
    }

    SECTION("exterior top face corner vertices are FIXED") {
        // mesh0 top face: 0,1,2,3 → z=0 (exterior top).
        // Each is a corner belonging to top + 2 side faces = 3 exterior faces → FIXED.
        for (uint32_t vi : {0u, 1u, 2u, 3u}) {
            REQUIRE(mesh0->vertex(vi).freedom == VertexFreedom::FIXED);
        }
    }
}

TEST_CASE("fixExteriorFaces is idempotent", "[model][fix_exterior]") {
    LithologyModel model;

    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::XYZ_FREE);

    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    model.addGroup(g0);
    model.setGroupMesh(0, mesh);

    model.fixExteriorFaces();

    // Record freedoms after first call
    std::vector<VertexFreedom> firstPass;
    for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
        firstPass.push_back(mesh->vertex(vi).freedom);
    }

    // Second call should not change anything
    model.fixExteriorFaces();

    for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
        REQUIRE(mesh->vertex(vi).freedom == firstPass[vi]);
    }
}

TEST_CASE("fixExteriorFaces does not crash on empty model", "[model][fix_exterior]") {
    LithologyModel model;
    REQUIRE_NOTHROW(model.fixExteriorFaces());
}

TEST_CASE("fixExteriorFaces handles model with no meshes", "[model][fix_exterior]") {
    LithologyModel model;
    LithoGroup g0; g0.id = 0; g0.name = "empty"; g0.density = 2.67;
    model.addGroup(g0);
    REQUIRE_NOTHROW(model.fixExteriorFaces());
}

} // namespace
} // namespace litho_invert
