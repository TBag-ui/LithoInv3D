#include <catch2/catch2.hpp>
#include <litho_invert/surface/surface_mesh.h>

namespace litho_invert {
namespace {

TEST_CASE("VertexFreedom DOF counting", "[surface][dof]") {
    SECTION("FIXED vertex contributes 0 DOFs") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0, VertexFreedom::FIXED);
        mesh.addTriangle(0, 0, 0); // dummy, DOF system doesn't check triangles
        REQUIRE(mesh.dofCount() == 0);
    }

    SECTION("Z_ONLY vertex contributes 1 DOF") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0, VertexFreedom::Z_ONLY);
        mesh.addTriangle(0, 0, 0);
        REQUIRE(mesh.dofCount() == 1);
    }

    SECTION("X_ONLY vertex contributes 1 DOF") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0, VertexFreedom::X_ONLY);
        mesh.addTriangle(0, 0, 0);
        REQUIRE(mesh.dofCount() == 1);
    }

    SECTION("Y_ONLY vertex contributes 1 DOF") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0, VertexFreedom::Y_ONLY);
        mesh.addTriangle(0, 0, 0);
        REQUIRE(mesh.dofCount() == 1);
    }

    SECTION("XY_FREE vertex contributes 2 DOFs") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0, VertexFreedom::XY_FREE);
        mesh.addTriangle(0, 0, 0);
        REQUIRE(mesh.dofCount() == 2);
    }

    SECTION("XYZ_FREE vertex contributes 3 DOFs") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0, VertexFreedom::XYZ_FREE);
        mesh.addTriangle(0, 0, 0);
        REQUIRE(mesh.dofCount() == 3);
    }

    SECTION("ALONG_VECTOR vertex contributes 1 DOF") {
        SurfaceMesh mesh;
        mesh.addVertex(Vector3d(0.0, 0.0, 0.0), VertexFreedom::ALONG_VECTOR, Vector3d(0.0, 0.0, 1.0));
        mesh.addTriangle(0, 0, 0);
        REQUIRE(mesh.dofCount() == 1);
    }

    SECTION("Mixed freedoms produce correct total") {
        SurfaceMesh mesh;
        mesh.addVertex(0.0, 0.0, 0.0, VertexFreedom::Z_ONLY);   // 1
        mesh.addVertex(1.0, 0.0, 0.0, VertexFreedom::XY_FREE);   // 2
        mesh.addVertex(2.0, 0.0, 0.0, VertexFreedom::XYZ_FREE);  // 3
        mesh.addVertex(3.0, 0.0, 0.0, VertexFreedom::FIXED);     // 0
        mesh.addVertex(4.0, 0.0, 0.0, VertexFreedom::Z_ONLY);    // 1
        mesh.addTriangle(0, 1, 2);
        REQUIRE(mesh.dofCount() == 7);
    }
}

TEST_CASE("VertexFreedom DOF mappings correctness", "[surface][dof]") {
    SECTION("Z_ONLY DOF maps to axis 2 (Z)") {
        SurfaceMesh mesh;
        mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::Z_ONLY);
        mesh.addTriangle(0, 0, 0);
        REQUIRE(mesh.dofCount() == 1);
        const auto& mappings = mesh.dofMappings();
        REQUIRE(mappings[0].vertexIndex == 0);
        REQUIRE(mappings[0].axis == 2);
    }

    SECTION("X_ONLY DOF maps to axis 0 (X)") {
        SurfaceMesh mesh;
        mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::X_ONLY);
        mesh.addTriangle(0, 0, 0);
        const auto& mappings = mesh.dofMappings();
        REQUIRE(mappings[0].axis == 0);
    }

    SECTION("Y_ONLY DOF maps to axis 1 (Y)") {
        SurfaceMesh mesh;
        mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::Y_ONLY);
        mesh.addTriangle(0, 0, 0);
        const auto& mappings = mesh.dofMappings();
        REQUIRE(mappings[0].axis == 1);
    }

    SECTION("XY_FREE creates separate mappings for X and Y") {
        SurfaceMesh mesh;
        mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::XY_FREE);
        mesh.addTriangle(0, 0, 0);
        const auto& mappings = mesh.dofMappings();
        REQUIRE(mappings.size() == 2);
        int axes = (1 << mappings[0].axis) | (1 << mappings[1].axis);
        REQUIRE(axes == ((1 << 0) | (1 << 1))); // axes 0 and 1, no axis 2
    }

    SECTION("XYZ_FREE creates mappings for all three axes") {
        SurfaceMesh mesh;
        mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::XYZ_FREE);
        mesh.addTriangle(0, 0, 0);
        const auto& mappings = mesh.dofMappings();
        REQUIRE(mappings.size() == 3);
        int axes = (1 << mappings[0].axis) | (1 << mappings[1].axis) | (1 << mappings[2].axis);
        REQUIRE(axes == 7); // all 3 axes
    }
}

TEST_CASE("VertexFreedom changes rebuild DOF mappings", "[surface][dof]") {
    SurfaceMesh mesh;
    mesh.addVertex(0.0, 0.0, 0.0, VertexFreedom::Z_ONLY);
    mesh.addVertex(1.0, 0.0, 0.0, VertexFreedom::Z_ONLY);
    mesh.addTriangle(0, 0, 0);

    REQUIRE(mesh.dofCount() == 2);

    SECTION("setVertexFreedom triggers rebuild") {
        mesh.setVertexFreedom(0, VertexFreedom::XYZ_FREE);
        // Now: vertex 0 = 3 DOFs, vertex 1 = 1 DOF
        REQUIRE(mesh.dofCount() == 4);
    }

    SECTION("setAllFreedom triggers rebuild") {
        mesh.setAllFreedom(VertexFreedom::XY_FREE);
        // Both vertices = 2 DOFs each
        REQUIRE(mesh.dofCount() == 4);
    }
}

TEST_CASE("Padding vertices excluded from DOF count", "[surface][dof][padding]") {
    SurfaceMesh mesh;
    // Create a 3x3 grid (9 vertices) with 1 padding ring around a 1x1 interior
    // interiorDim=1, paddingRings=1 → fullDim=3 → 9 vertices
    // Interior: only vertex 4 (center) = index 1*3+1 = 4
    // Padding: all others (0,1,2,3,5,6,7,8)

    // Build 3x3 grid
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            mesh.addVertex(static_cast<double>(c), static_cast<double>(r), 0.0, VertexFreedom::Z_ONLY);
        }
    }

    mesh.setInteriorGrid(1, 1);  // 1x1 interior, 1 padding ring

    // Without padding enabled, all vertices get DOFs: 9 Z_ONLY = 9 DOFs
    // With padding, only interior vertex 4 gets a DOF: 1 DOF
    // Interior vertex: r=1, c=1 → index 1*3+1 = 4
    REQUIRE(mesh.dofCount() == 1);
    REQUIRE(mesh.isInteriorVertex(4));
    REQUIRE_FALSE(mesh.isPaddingVertex(4));

    // Verify a corner padding vertex
    REQUIRE(mesh.isPaddingVertex(0));
    REQUIRE_FALSE(mesh.isInteriorVertex(0));
}

} // namespace
} // namespace litho_invert
