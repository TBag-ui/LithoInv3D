#include <catch2/catch2.hpp>
#include <litho_invert/litho/lithology_model.h>

using Catch::Approx;

namespace litho_invert {
namespace {

static std::shared_ptr<SurfaceMesh> makeSimpleMesh(VertexFreedom f) {
    auto mesh = std::make_shared<SurfaceMesh>();
    // 4 vertices forming a simple surface at z=-100
    mesh->addVertex(0.0, 0.0, -100.0, f);
    mesh->addVertex(1.0, 0.0, -100.0, f);
    mesh->addVertex(1.0, 1.0, -100.0, f);
    mesh->addVertex(0.0, 1.0, -100.0, f);
    mesh->addTriangle(0, 1, 2);
    mesh->addTriangle(0, 2, 3);
    return mesh;
}

TEST_CASE("Global DOF deduplication", "[model][global_dof]") {
    LithologyModel model;

    // Two meshes with identical vertex positions → shared DOFs
    auto mesh0 = makeSimpleMesh(VertexFreedom::Z_ONLY);
    auto mesh1 = makeSimpleMesh(VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "g0"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "g1"; g1.density = 3.00;

    model.addGroup(g0);
    model.addGroup(g1);
    model.setGroupMesh(0, mesh0);
    model.setGroupMesh(1, mesh1);

    SECTION("same position and axis share a single global DOF") {
        // Each mesh has 4 Z_ONLY vertices = 4 DOFs
        // Since positions are identical, global DOFs = 4 (not 8)
        uint32_t total = model.totalDofCount();
        REQUIRE(total == 4);
    }

    SECTION("globalDofIndex returns valid index for each local DOF") {
        for (int mi = 0; mi < 2; ++mi) {
            for (uint32_t li = 0; li < 4; ++li) {
                int gIdx = model.globalDofIndex(mi, static_cast<int>(li));
                REQUIRE(gIdx >= 0);
                REQUIRE(gIdx < 4);
            }
        }
    }

    SECTION("shared local DOFs map to same global DOF") {
        int g0_idx0 = model.globalDofIndex(0, 0); // mesh 0, local DOF 0
        int g1_idx0 = model.globalDofIndex(1, 0); // mesh 1, local DOF 0
        REQUIRE(g0_idx0 == g1_idx0);
    }
}

TEST_CASE("Total DOF count correctness", "[model][global_dof]") {
    LithologyModel model;

    // Mesh 0: 4 XYZ_FREE vertices at unique positions
    // Mesh 1: 4 Z_ONLY vertices at DIFFERENT positions
    auto mesh0 = std::make_shared<SurfaceMesh>();
    mesh0->addVertex(0.0, 0.0, -50.0, VertexFreedom::XYZ_FREE);
    mesh0->addVertex(1.0, 0.0, -50.0, VertexFreedom::XYZ_FREE);
    mesh0->addVertex(1.0, 1.0, -50.0, VertexFreedom::XYZ_FREE);
    mesh0->addVertex(0.0, 1.0, -50.0, VertexFreedom::XYZ_FREE);
    mesh0->addTriangle(0, 1, 2); mesh0->addTriangle(0, 2, 3);

    auto mesh1 = std::make_shared<SurfaceMesh>();
    mesh1->addVertex(0.0, 0.0, -150.0, VertexFreedom::Z_ONLY);
    mesh1->addVertex(1.0, 0.0, -150.0, VertexFreedom::Z_ONLY);
    mesh1->addVertex(1.0, 1.0, -150.0, VertexFreedom::Z_ONLY);
    mesh1->addVertex(0.0, 1.0, -150.0, VertexFreedom::Z_ONLY);
    mesh1->addTriangle(0, 1, 2); mesh1->addTriangle(0, 2, 3);

    LithoGroup g0; g0.id = 0; g0.name = "g0"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "g1"; g1.density = 3.00;

    model.addGroup(g0);
    model.addGroup(g1);
    model.setGroupMesh(0, mesh0);
    model.setGroupMesh(1, mesh1);

    uint32_t total = model.totalDofCount();
    // Mesh 0: 4 × 3 (XYZ) = 12 DOFs at unique positions
    // Mesh 1: 4 × 1 (Z_only) = 4 DOFs at different positions
    // No sharing (different positions) → 16 total
    REQUIRE(total == 16);
}

TEST_CASE("assembleParameterVector roundtrip", "[model][global_dof]") {
    LithologyModel model;

    auto mesh0 = makeSimpleMesh(VertexFreedom::Z_ONLY);
    auto mesh1 = makeSimpleMesh(VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "g0"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "g1"; g1.density = 3.00;

    model.addGroup(g0);
    model.addGroup(g1);
    model.setGroupMesh(0, mesh0);
    model.setGroupMesh(1, mesh1);

    // Assemble initial params (all Z = -100)
    VectorXd params = model.assembleParameterVector();
    for (Eigen::Index i = 0; i < params.size(); ++i) {
        params[i] = -200.0; // change all Z to -200
    }

    model.applyParameterVector(params);

    // Re-assemble — should get -200 back
    VectorXd params2 = model.assembleParameterVector();
    for (Eigen::Index i = 0; i < params2.size(); ++i) {
        REQUIRE(params2[i] == Approx(-200.0));
    }
}

} // namespace
} // namespace litho_invert
