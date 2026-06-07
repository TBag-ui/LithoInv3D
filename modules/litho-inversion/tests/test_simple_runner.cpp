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

static LithologyModel makeSingleGroupModel() {
    LithologyModel model;

    // Box from (0, 0, 0) to (10, 10, -20)
    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);

    LithoGroup g0;
    g0.id = 0;
    g0.name = "body";
    g0.density = 2.67;
    model.addGroup(g0);
    model.setGroupMesh(0, mesh);

    return model;
}

static GravityData makeGridObservations(
    double x0, double x1, int nx,
    double y0, double y1, int ny,
    double zObs)
{
    GravityData data;
    double dx = (nx > 1) ? (x1 - x0) / (nx - 1) : 0.0;
    double dy = (ny > 1) ? (y1 - y0) / (ny - 1) : 0.0;

    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            double x = x0 + ix * dx;
            double y = y0 + iy * dy;
            data.push_back(GravityPoint(Vector3d(x, y, zObs), 0.0, 0.01));
        }
    }
    return data;
}

TEST_CASE("InversionRunner can be constructed with valid config", "[inversion][runner]") {
    auto model = std::make_shared<LithologyModel>(makeSingleGroupModel());
    auto obs = makeGridObservations(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 5;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";

    REQUIRE_NOTHROW(InversionRunner(config));
}

TEST_CASE("InversionRunner runs gravity-only inversion", "[inversion][runner]") {
    auto model = std::make_shared<LithologyModel>(makeSingleGroupModel());
    auto obs = makeGridObservations(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    auto result = runner.run();

    // Runner should complete without crashing
    REQUIRE(result.finalModel != nullptr);
    REQUIRE(result.totalIterations > 0);
}

TEST_CASE("InversionRunner records iteration history", "[inversion][runner]") {
    auto model = std::make_shared<LithologyModel>(makeSingleGroupModel());
    auto obs = makeGridObservations(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    auto result = runner.run();

    // History should have at least initial iteration (iter 0) + completed iterations
    REQUIRE(result.history.size() >= 1);
    // first entry is iteration 0
    REQUIRE(result.history[0].iteration == 0);
}

TEST_CASE("InversionRunner invokes iteration callback", "[inversion][runner]") {
    auto model = std::make_shared<LithologyModel>(makeSingleGroupModel());
    auto obs = makeGridObservations(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);

    int callbackCount = 0;
    runner.setIterationCallback([&callbackCount](const InversionIteration&) {
        ++callbackCount;
    });

    runner.run();

    // Callback should fire at least once (for initial iteration 0)
    REQUIRE(callbackCount >= 1);
}

TEST_CASE("InversionRunner final model is valid", "[inversion][runner]") {
    auto model = std::make_shared<LithologyModel>(makeSingleGroupModel());
    auto obs = makeGridObservations(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    auto result = runner.run();

    REQUIRE(result.finalModel != nullptr);
    REQUIRE(result.finalModel->groupMeshCount() > 0);
    auto* mesh = result.finalModel->groupMesh(0);
    REQUIRE(mesh != nullptr);
    REQUIRE(mesh->vertexCount() > 0);
    REQUIRE(mesh->isValid());
}

TEST_CASE("InversionRunner with XYZ_FREE vertices", "[inversion][runner]") {
    auto model = std::make_shared<LithologyModel>();

    // Override: create model with XYZ_FREE
    auto mesh = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::XYZ_FREE);

    LithologyModel xyzModel;
    LithoGroup g0; g0.id = 0; g0.name = "body"; g0.density = 2.67;
    xyzModel.addGroup(g0);
    xyzModel.setGroupMesh(0, mesh);

    auto modelPtr = std::make_shared<LithologyModel>(xyzModel);
    auto obs = makeGridObservations(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = modelPtr;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::XYZ_FREE;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    auto result = runner.run();

    REQUIRE(result.finalModel != nullptr);
    REQUIRE(result.totalIterations > 0);
}

TEST_CASE("InversionRunner preserves DOF count during inversion", "[inversion][runner]") {
    auto model = std::make_shared<LithologyModel>(makeSingleGroupModel());
    auto obs = makeGridObservations(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    uint32_t dofsBefore = model->totalDofCount();

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

    uint32_t dofsAfter = model->totalDofCount();
    REQUIRE(dofsBefore == dofsAfter);
}

TEST_CASE("InversionRunner handles two-group model", "[inversion][runner]") {
    auto model = std::make_shared<LithologyModel>();

    // Two stacked boxes sharing a face at z=-20
    auto mesh0 = makeBoxMesh(0.0, 0.0, 0.0, 10.0, 10.0, -20.0, VertexFreedom::Z_ONLY);
    auto mesh1 = makeBoxMesh(0.0, 0.0, -20.0, 10.0, 10.0, -50.0, VertexFreedom::Z_ONLY);

    LithoGroup g0; g0.id = 0; g0.name = "upper"; g0.density = 2.67;
    LithoGroup g1; g1.id = 1; g1.name = "lower"; g1.density = 3.00;

    model->addGroup(g0);
    model->addGroup(g1);
    model->setGroupMesh(0, mesh0);
    model->setGroupMesh(1, mesh1);

    auto obs = makeGridObservations(2.5, 7.5, 2, 2.5, 7.5, 2, 0.0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.maxIterations = 2;
    config.vertexFreedom = VertexFreedom::Z_ONLY;
    config.solver = "lbfgsb";
    config.lambda = 0.01;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    auto result = runner.run();

    REQUIRE(result.finalModel != nullptr);
    REQUIRE(result.finalModel->groupMeshCount() == 2);
}

} // namespace
} // namespace litho_invert
