#include <catch2/catch2.hpp>
#include <litho_invert/surface/surface_mesh.h>

namespace litho_invert {
namespace {

TEST_CASE("classifySide on simple surface", "[surface][normals]") {
    SurfaceMesh mesh;

    // Create a simple 2x2 grid (4 vertices, 2 triangles) at z=0
    // v0---v1
    // | \   |
    // |  \  |
    // |   \ |
    // v2---v3
    mesh.addVertex(0.0, 0.0, 0.0);
    mesh.addVertex(1.0, 0.0, 0.0);
    mesh.addVertex(0.0, 1.0, 0.0);
    mesh.addVertex(1.0, 1.0, 0.0);

    // CCW winding when viewed from above (+Z)
    mesh.addTriangle(0, 1, 2); // CCW looking from +Z: v0→v1→v2
    mesh.addTriangle(1, 3, 2); // CCW looking from +Z: v1→v3→v2

    SECTION("point above surface returns +1") {
        int side = mesh.classifySide(Vector3d(0.5, 0.5, 10.0));
        REQUIRE(side == 1);
    }

    SECTION("point below surface returns -1") {
        int side = mesh.classifySide(Vector3d(0.5, 0.5, -10.0));
        REQUIRE(side == -1);
    }

    SECTION("point on surface returns 0") {
        int side = mesh.classifySide(Vector3d(0.5, 0.5, 0.0));
        REQUIRE(side == 0);
    }
}

TEST_CASE("classifySide outside XY extent", "[surface][normals]") {
    SurfaceMesh mesh;
    mesh.addVertex(0.0, 0.0, 0.0);
    mesh.addVertex(1.0, 0.0, 0.0);
    mesh.addVertex(0.0, 1.0, 0.0);
    mesh.addVertex(1.0, 1.0, 0.0);
    mesh.addTriangle(0, 1, 2);
    mesh.addTriangle(1, 3, 2);

    // Point outside the XY extent of the mesh
    SECTION("outside point above surface") {
        int side = mesh.classifySide(Vector3d(5.0, 5.0, 100.0));
        REQUIRE(side == 1);
    }

    SECTION("outside point below surface") {
        int side = mesh.classifySide(Vector3d(5.0, 5.0, -100.0));
        REQUIRE(side == -1);
    }
}

TEST_CASE("Triangle winding produces correct normal orientation", "[surface][normals]") {
    // Two triangles with opposite winding should produce opposite side
    // classification results for the same point

    SurfaceMesh mesh1;
    mesh1.addVertex(0.0, 0.0, 0.0);
    mesh1.addVertex(1.0, 0.0, 0.0);
    mesh1.addVertex(0.0, 1.0, 0.0);
    mesh1.addTriangle(0, 1, 2); // CCW when viewed from +Z

    SurfaceMesh mesh2;
    mesh2.addVertex(0.0, 0.0, 0.0);
    mesh2.addVertex(1.0, 0.0, 0.0);
    mesh2.addVertex(0.0, 1.0, 0.0);
    mesh2.addTriangle(0, 2, 1); // CW when viewed from +Z (opposite winding)

    // A point well above
    Vector3d above(0.3, 0.3, 10.0);

    int side1 = mesh1.classifySide(above);
    int side2 = mesh2.classifySide(above);

    // Opposite winding → opposite side classification
    REQUIRE(side1 == -side2);
}

TEST_CASE("isValid mesh validation", "[surface][validation]") {
    SECTION("empty mesh is invalid") {
        SurfaceMesh mesh;
        REQUIRE_FALSE(mesh.isValid());
    }

    SECTION("mesh with vertices but no triangles is invalid") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0);
        mesh.addVertex(1.0, 0.0, 0.0);
        mesh.addVertex(0.0, 1.0, 0.0);
        REQUIRE_FALSE(mesh.isValid());
    }

    SECTION("mesh with out-of-range triangle is invalid") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0);
        mesh.addVertex(1.0, 0.0, 0.0);
        mesh.addVertex(0.0, 1.0, 0.0);
        // v2=999 which doesn't exist
        mesh.addTriangle(0, 1, 999);
        REQUIRE_FALSE(mesh.isValid());
    }

    SECTION("valid mesh passes validation") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0);
        mesh.addVertex(1.0, 0.0, 0.0);
        mesh.addVertex(0.0, 1.0, 0.0);
        mesh.addTriangle(0, 1, 2);
        REQUIRE(mesh.isValid());
    }
}

TEST_CASE("Vertex count and naming", "[surface]") {
    SurfaceMesh mesh;
    mesh.addVertex(0.0, 0.0, 0.0);
    mesh.addVertex(1.0, 0.0, 0.0);

    REQUIRE(mesh.vertexCount() == 2);

    mesh.setName("test_surface");
    REQUIRE(mesh.name() == "test_surface");
}

TEST_CASE("neighbors build correctly on simple mesh", "[surface][neighbors]") {
    SurfaceMesh mesh;
    // 2x2 grid with 2 triangles
    mesh.addVertex(0.0, 0.0, 0.0); // 0
    mesh.addVertex(1.0, 0.0, 0.0); // 1
    mesh.addVertex(0.0, 1.0, 0.0); // 2
    mesh.addVertex(1.0, 1.0, 0.0); // 3
    mesh.addTriangle(0, 1, 2);
    mesh.addTriangle(1, 3, 2);

    mesh.buildNeighbors();

    // Vertex 0 should be neighbors with 1 and 2
    const auto& n0 = mesh.neighborVertices(0);
    REQUIRE(n0.size() == 2);
    REQUIRE((n0[0] == 1 || n0[1] == 1));
    REQUIRE((n0[0] == 2 || n0[1] == 2));

    // Vertex 3 should be neighbors with 1 and 2
    const auto& n3 = mesh.neighborVertices(3);
    REQUIRE(n3.size() == 2);
}
} // namespace
} // namespace litho_invert
