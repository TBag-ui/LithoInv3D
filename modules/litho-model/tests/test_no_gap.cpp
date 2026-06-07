#include <catch2/catch2.hpp>
#include <litho_invert/litho/lithology_model.h>

using Catch::Approx;

namespace litho_invert {
namespace {

// Helper: build a simple closed box mesh (8 vertices, 12 triangles for 6 faces)
static std::shared_ptr<SurfaceMesh> makeBoxMesh(
    double x0, double y0, double z0,
    double x1, double y1, double z1,
    VertexFreedom freedom = VertexFreedom::XYZ_FREE)
{
    auto mesh = std::make_shared<SurfaceMesh>();
    // 8 vertices
    mesh->addVertex(x0, y0, z0, freedom); // 0
    mesh->addVertex(x1, y0, z0, freedom); // 1
    mesh->addVertex(x1, y1, z0, freedom); // 2
    mesh->addVertex(x0, y1, z0, freedom); // 3
    mesh->addVertex(x0, y0, z1, freedom); // 4
    mesh->addVertex(x1, y0, z1, freedom); // 5
    mesh->addVertex(x1, y1, z1, freedom); // 6
    mesh->addVertex(x0, y1, z1, freedom); // 7

    // 6 faces, 2 triangles each, CCW facing outward
    // +Z top:   0-1-2, 0-2-3
    mesh->addTriangle(0, 1, 2); mesh->addTriangle(0, 2, 3);
    // -Z bottom: 4-7-6, 4-6-5 (CCW looking from -Z)
    mesh->addTriangle(4, 7, 6); mesh->addTriangle(4, 6, 5);
    // +Y front: 3-2-6, 3-6-7
    mesh->addTriangle(3, 2, 6); mesh->addTriangle(3, 6, 7);
    // -Y back:  0-4-5, 0-5-1
    mesh->addTriangle(0, 4, 5); mesh->addTriangle(0, 5, 1);
    // +X right: 1-5-6, 1-6-2
    mesh->addTriangle(1, 5, 6); mesh->addTriangle(1, 6, 2);
    // -X left:  0-3-7, 0-7-4
    mesh->addTriangle(0, 3, 7); mesh->addTriangle(0, 7, 4);

    return mesh;
}

TEST_CASE("Shared vertices between adjacent meshes", "[model][no_gap]") {
    LithologyModel model;

    // Group 0: box from (0,0,0) to (1,1,-2) — background
    auto mesh0 = makeBoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, -2.0, VertexFreedom::XYZ_FREE);

    // Group 1: box from (0,0,-2) to (1,1,-5) — anomaly
    // Shares the z=-2 face with mesh0
    auto mesh1 = makeBoxMesh(0.0, 0.0, -2.0, 1.0, 1.0, -5.0, VertexFreedom::XYZ_FREE);

    LithoGroup g0; g0.id = 0; g0.name = "background"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "anomaly"; g1.density = 3.00;

    model.addGroup(g0);
    model.addGroup(g1);
    model.setGroupMesh(0, mesh0);
    model.setGroupMesh(1, mesh1);

    SECTION("shared position vertices map to single global DOF") {
        // The top face of mesh1 (= z=-2) shares 4 vertices with
        // the bottom face of mesh0. After rebuild, shared DOFs should
        // be deduplicated.
        uint32_t nDofs = model.totalDofCount();

        // Without sharing: 8 verts * 3 (XYZ_FREE) = 24 per mesh = 48 total
        // With sharing: 4 shared vertices × 3 axes = 12 shared DOFs
        //    + 4 unique vertices × 3 axes × 2 meshes = 24 unique DOFs
        //    Total: 24 unique + 12 shared = 36 DOFs
        REQUIRE(nDofs > 0);
        REQUIRE(nDofs < 48); // Must be deduplicated
    }

    SECTION("shared vertices have identical positions after applyParameterVector") {
        model.fixExteriorFaces();

        VectorXd params = model.assembleParameterVector();
        model.applyParameterVector(params);

        // Check a shared vertex: mesh0 vertex 4 (bottom corner) should equal
        // mesh1 vertex 0 (top corner) — both at (0, 0, -2)
        // mesh0 vertices: 4=(0,0,-2), 5=(1,0,-2), 6=(1,1,-2), 7=(0,1,-2)
        // mesh1 vertices: 0=(0,0,-2), 1=(1,0,-2), 2=(1,1,-2), 3=(0,1,-2)
        const auto& m0v4 = mesh0->vertex(4).position;
        const auto& m1v0 = mesh1->vertex(0).position;
        REQUIRE(m0v4.x() == Approx(m1v0.x()));
        REQUIRE(m0v4.y() == Approx(m1v0.y()));
        REQUIRE(m0v4.z() == Approx(m1v0.z()));
    }

    SECTION("shared position is maintained after parameter perturbation") {
        // Perturb the initial assembled parameters (no fixExteriorFaces).
        // Since shared vertices at identical positions share global DOFs,
        // they must stay synchronized after any parameter change.
        VectorXd params = model.assembleParameterVector();

        for (Eigen::Index i = 0; i < params.size(); ++i) {
            params[i] *= 1.1;
        }

        VectorXd lower, upper;
        model.getBounds(lower, upper);
        for (Eigen::Index i = 0; i < params.size(); ++i) {
            params[i] = std::max(lower[i], std::min(upper[i], params[i]));
        }

        model.applyParameterVector(params);

        const auto& m0v4 = mesh0->vertex(4).position;
        const auto& m1v0 = mesh1->vertex(0).position;
        REQUIRE(m0v4.x() == Approx(m1v0.x()).margin(1e-6));
        REQUIRE(m0v4.y() == Approx(m1v0.y()).margin(1e-6));
        REQUIRE(m0v4.z() == Approx(m1v0.z()).margin(1e-6));
    }

    SECTION("all four shared corners stay synchronized") {
        model.fixExteriorFaces();
        // mesh0 bottom: 4=(0,0,-2), 5=(1,0,-2), 6=(1,1,-2), 7=(0,1,-2)
        // mesh1 top:    0=(0,0,-2), 1=(1,0,-2), 2=(1,1,-2), 3=(0,1,-2)

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
}

TEST_CASE("Three meshes sharing a vertex", "[model][no_gap]") {
    LithologyModel model;

    // Three adjacent boxes sharing a face
    auto mesh0 = makeBoxMesh(0.0, 0.0, 0.0, 1.0, 1.0, -2.0, VertexFreedom::XYZ_FREE);
    auto mesh1 = makeBoxMesh(0.0, 0.0, -2.0, 1.0, 1.0, -4.0, VertexFreedom::XYZ_FREE);
    auto mesh2 = makeBoxMesh(0.0, 0.0, -4.0, 1.0, 1.0, -6.0, VertexFreedom::XYZ_FREE);

    LithoGroup g0; g0.id = 0; g0.name = "g0"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "g1"; g1.density = 3.00;
    LithoGroup g2; g2.id = 2; g2.name = "g2"; g2.density = 3.30;

    model.addGroup(g0);
    model.addGroup(g1);
    model.addGroup(g2);
    model.setGroupMesh(0, mesh0);
    model.setGroupMesh(1, mesh1);
    model.setGroupMesh(2, mesh2);

    SECTION("total DOF count is less than sum of individual mesh DOFs") {
        uint32_t total = model.totalDofCount();
        uint32_t sumIndividual = 8 * 3 * 3; // 8 vertices × 3 axes × 3 meshes = 72
        REQUIRE(total > 0);
        REQUIRE(total < sumIndividual);
    }

    SECTION("shared corners are identical across all three meshes") {
        model.fixExteriorFaces();
        VectorXd params = model.assembleParameterVector();
        model.applyParameterVector(params);

        // mesh0 bottom = (0,0,-2), mesh1 top = (0,0,-2), mesh1 bottom = (0,0,-4), mesh2 top = (0,0,-4)
        const auto& p_m0b = mesh0->vertex(4).position; // (0,0,-2)
        const auto& p_m1t = mesh1->vertex(0).position; // (0,0,-2)
        const auto& p_m1b = mesh1->vertex(4).position; // (0,0,-4)
        const auto& p_m2t = mesh2->vertex(0).position; // (0,0,-4)

        REQUIRE(p_m0b.x() == Approx(p_m1t.x()).margin(1e-6));
        REQUIRE(p_m0b.y() == Approx(p_m1t.y()).margin(1e-6));
        REQUIRE(p_m0b.z() == Approx(p_m1t.z()).margin(1e-6));

        REQUIRE(p_m1b.x() == Approx(p_m2t.x()).margin(1e-6));
        REQUIRE(p_m1b.y() == Approx(p_m2t.y()).margin(1e-6));
        REQUIRE(p_m1b.z() == Approx(p_m2t.z()).margin(1e-6));
    }
}

} // namespace
} // namespace litho_invert
