#include <catch2/catch2.hpp>
#include <litho_invert/surface/surface_mesh.h>
#include <cmath>

using Catch::Approx;

namespace litho_invert {
namespace {

TEST_CASE("applyParams updates vertex positions", "[surface][dof]") {
    SurfaceMesh mesh;
    mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::XYZ_FREE);
    mesh.addTriangle(0, 0, 0);

    VectorXd params(3);
    params << 100.0, 200.0, 300.0;
    mesh.applyParams(params, 0);

    const auto& v = mesh.vertex(0);
    REQUIRE(v.position.x() == Approx(100.0));
    REQUIRE(v.position.y() == Approx(200.0));
    REQUIRE(v.position.z() == Approx(300.0));
}

TEST_CASE("applyParams respects vertex freedom", "[surface][dof]") {
    SECTION("Z_ONLY vertex only updates Z") {
        SurfaceMesh mesh;
        mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::Z_ONLY);
        mesh.addTriangle(0, 0, 0);

        VectorXd params(1);
        params << 500.0;
        mesh.applyParams(params, 0);

        const auto& v = mesh.vertex(0);
        REQUIRE(v.position.x() == Approx(10.0));
        REQUIRE(v.position.y() == Approx(20.0));
        REQUIRE(v.position.z() == Approx(500.0));
    }

    SECTION("XY_FREE vertex updates X and Y but not Z") {
        SurfaceMesh mesh;
        mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::XY_FREE);
        mesh.addTriangle(0, 0, 0);

        VectorXd params(2);
        params << 100.0, 200.0;
        mesh.applyParams(params, 0);

        const auto& v = mesh.vertex(0);
        REQUIRE(v.position.x() == Approx(100.0));
        REQUIRE(v.position.y() == Approx(200.0));
        REQUIRE(v.position.z() == Approx(30.0)); // unchanged
    }

    SECTION("X_ONLY vertex updates only X") {
        SurfaceMesh mesh;
        mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::X_ONLY);
        mesh.addTriangle(0, 0, 0);

        VectorXd params(1);
        params << 50.0;
        mesh.applyParams(params, 0);

        const auto& v = mesh.vertex(0);
        REQUIRE(v.position.x() == Approx(50.0));
        REQUIRE(v.position.y() == Approx(20.0));
        REQUIRE(v.position.z() == Approx(30.0));
    }

    SECTION("FIXED vertex does not change") {
        SurfaceMesh mesh;
        mesh.addVertex(10.0, 20.0, 30.0, VertexFreedom::FIXED);
        mesh.addTriangle(0, 0, 0);

        VectorXd params(1); // empty params for this mesh
        mesh.applyParams(params, 0);

        const auto& v = mesh.vertex(0);
        REQUIRE(v.position.x() == Approx(10.0));
        REQUIRE(v.position.y() == Approx(20.0));
        REQUIRE(v.position.z() == Approx(30.0));
    }
}

TEST_CASE("extractParams roundtrip", "[surface][dof]") {
    SurfaceMesh mesh;
    mesh.addVertex(100.0, 200.0, 300.0, VertexFreedom::XYZ_FREE);
    mesh.addVertex(400.0, 500.0, 600.0, VertexFreedom::Z_ONLY);
    mesh.addTriangle(0, 1, 0);

    // 3 DOFs (XYZ) + 1 DOF (Z only) = 4 total
    REQUIRE(mesh.dofCount() == 4);

    VectorXd extracted(4);
    mesh.extractParams(extracted, 0);

    // First dof: X of vertex 0
    REQUIRE(extracted[0] == Approx(100.0));
    // Second dof: Y of vertex 0
    REQUIRE(extracted[1] == Approx(200.0));
    // Third dof: Z of vertex 0
    REQUIRE(extracted[2] == Approx(300.0));
    // Fourth dof: Z of vertex 1
    REQUIRE(extracted[3] == Approx(600.0));

    // Now apply different params and re-extract
    SECTION("apply-extract round trip") {
        VectorXd new_params(4);
        new_params << 150.0, 250.0, 350.0, 650.0;
        mesh.applyParams(new_params, 0);

        VectorXd re_extracted(4);
        mesh.extractParams(re_extracted, 0);

        for (int i = 0; i < 4; ++i) {
            REQUIRE(re_extracted[i] == Approx(new_params[i]));
        }
    }
}

TEST_CASE("getBounds respects depth limits", "[surface][dof][bounds]") {
    SurfaceMesh mesh;
    mesh.addVertex(0.0, 0.0, -100.0, VertexFreedom::Z_ONLY);
    mesh.addVertex(0.0, 0.0, -200.0, VertexFreedom::XYZ_FREE);
    mesh.addTriangle(0, 1, 0);

    mesh.setBounds(-500.0, 50.0); // depth range: -500 to 50

    VectorXd lower, upper;
    uint32_t ndof = mesh.dofCount();
    lower.resize(ndof);
    upper.resize(ndof);
    mesh.getBounds(lower, upper, 0);

    // 1 Z_ONLY (Z) + 3 XYZ_FREE (X,Y,Z) = 4 DOFs
    REQUIRE(lower.size() == 4);
    REQUIRE(upper.size() == 4);

    // Z DOFs should respect depth bounds
    // DOF 0: Z_ONLY → Z axis → lower = -500, upper = 50
    REQUIRE(lower[0] == Approx(-500.0));
    REQUIRE(upper[0] == Approx(50.0));

    // DOF 1: XYZ_FREE → X axis → large bounds
    REQUIRE(lower[1] == Approx(-1e6));
    // DOF 2: XYZ_FREE → Y axis → large bounds
    REQUIRE(lower[2] == Approx(-1e6));
    // DOF 3: XYZ_FREE → Z axis → depth bounds
    REQUIRE(lower[3] == Approx(-500.0));
    REQUIRE(upper[3] == Approx(50.0));
}

TEST_CASE("getBounds with axis scaling", "[surface][dof][bounds]") {
    SurfaceMesh mesh;
    mesh.addVertex(0.0, 0.0, -100.0, VertexFreedom::Z_ONLY);
    mesh.addTriangle(0, 0, 0);

    mesh.setBounds(-1000.0, 100.0);
    mesh.setAxisScale(2, 10.0); // Z axis scale = 10

    VectorXd lower, upper;
    uint32_t ndof2 = mesh.dofCount();
    lower.resize(ndof2);
    upper.resize(ndof2);
    mesh.getBounds(lower, upper, 0);

    // Z bound scaled: lower = minZ / scale = -1000 / 10 = -100
    REQUIRE(lower[0] == Approx(-100.0));
    REQUIRE(upper[0] == Approx(10.0)); // 100 / 10
}

TEST_CASE("applyParams with axis scaling", "[surface][dof]") {
    SurfaceMesh mesh;
    mesh.addVertex(0.0, 0.0, -100.0, VertexFreedom::Z_ONLY);
    mesh.addTriangle(0, 0, 0);

    mesh.setAxisScale(2, 10.0); // Z axis scale = 10

    VectorXd params(1);
    params << 0.0; // param value 0 → position = 0 * 10 = 0
    mesh.applyParams(params, 0);

    REQUIRE(mesh.vertex(0).position.z() == Approx(0.0));

    params << 5.0; // param value 5 → position = 5 * 10 = 50
    mesh.applyParams(params, 0);
    REQUIRE(mesh.vertex(0).position.z() == Approx(50.0));
}

TEST_CASE("setBounds stores depth range", "[surface][bounds]") {
    SurfaceMesh mesh;
    mesh.addVertex(0.0, 0.0, 0.0);
    mesh.addTriangle(0, 0, 0);

    mesh.setBounds(-1500.0, 200.0);
    REQUIRE(mesh.minDepth() == Approx(-1500.0));
    REQUIRE(mesh.maxDepth() == Approx(200.0));
}

} // namespace
} // namespace litho_invert
