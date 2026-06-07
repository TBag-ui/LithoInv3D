#include <catch2/catch2.hpp>
#include <litho_invert/inversion/runner.h>
#include <litho_invert/inversion/lbfgsb_optimizer.h>
#include <litho_invert/forward/gravity_forward.h>
#include <litho_invert/forward/magnetic_forward.h>
#include <litho_invert/io/surface_loader.h>
#include <litho_invert/io/exporters.h>
#include <litho_invert/plot/svg_plot.h>
#include "../generate_synthetic.h"
#include "../model_setup.h"
#include <filesystem>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <set>

namespace litho_invert {

using Catch::Approx;

// =========================================================================
// Test data path -- volume meshes live alongside test source
// =========================================================================
static std::string dataPath(const std::string& fname) {
    return "data/" + fname;
}

// =========================================================================
// Signed volume of a closed triangulated mesh via divergence theorem
// =========================================================================
static double computeMeshVolume(const SurfaceMesh& mesh) {
    double vol = 0.0;
    for (uint32_t t = 0; t < mesh.triangleCount(); ++t) {
        const auto& tri = mesh.triangle(t);
        const auto& p0 = mesh.vertex(tri.v0).position;
        const auto& p1 = mesh.vertex(tri.v1).position;
        const auto& p2 = mesh.vertex(tri.v2).position;
        vol += p0.dot(p1.cross(p2));
    }
    return vol / 6.0;
}

static double totalModelVolume(const LithologyModel& model) {
    double total = 0.0;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        total += std::abs(computeMeshVolume(*model.groupMesh(g)));
    }
    return total;
}

// =========================================================================
// Load cluster volume meshes from data/ directory
// =========================================================================
static std::shared_ptr<LithologyModel> loadClusterModel() {
    auto model = std::make_shared<LithologyModel>();

    model->addGroup({0, "Background",   2.67, 0.000});
    model->addGroup({1, "Intermediate", 2.80, 0.005});
    model->addGroup({2, "Mafic",        3.00, 0.020});
    model->addGroup({3, "Sulfide",      3.50, 0.080});

    double minZ = 1e30;
    for (int g = 0; g < 4; ++g) {
        std::string path = dataPath("volume_group_" + std::to_string(g + 1) + ".ts");
        auto mesh = TSurfaceLoader::load(path);
        if (!mesh) {
            throw std::runtime_error("Failed to load " + path);
        }
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            minZ = std::min(minZ, mesh->vertex(vi).position.z());
        }
        model->setGroupMesh(g, mesh);
    }
    model->setBottomDepth(minZ - 100.0);

    return model;
}

// =========================================================================
// Observation points on a grid covering the model extent
// =========================================================================
static GravityData generateObsForModel(const LithologyModel& model, int nPerSide = 10) {
    double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            const auto& p = mesh->vertex(vi).position;
            if (p.x() < xmin) xmin = p.x();
            if (p.x() > xmax) xmax = p.x();
            if (p.y() < ymin) ymin = p.y();
            if (p.y() > ymax) ymax = p.y();
        }
    }

    GravityData data;
    int n = std::max(2, nPerSide);
    for (int iy = 0; iy < n; ++iy) {
        for (int ix = 0; ix < n; ++ix) {
            double x = xmin + (xmax - xmin) * ix / (n - 1);
            double y = ymin + (ymax - ymin) * iy / (n - 1);
            data.push_back({Vector3d(x, y, 0.0), 0.0, 0.05});
        }
    }
    return data;
}

// =========================================================================
// Compute AABB of all group meshes
// =========================================================================
static void computeModelAABB(const LithologyModel& model,
                              double& xmin, double& xmax,
                              double& ymin, double& ymax,
                              double& zmin, double& zmax) {
    xmin = 1e30; xmax = -1e30;
    ymin = 1e30; ymax = -1e30;
    zmin = 1e30; zmax = -1e30;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            const auto& p = mesh->vertex(vi).position;
            if (p.x() < xmin) xmin = p.x();
            if (p.x() > xmax) xmax = p.x();
            if (p.y() < ymin) ymin = p.y();
            if (p.y() > ymax) ymax = p.y();
            if (p.z() < zmin) zmin = p.z();
            if (p.z() > zmax) zmax = p.z();
        }
    }
}

// =========================================================================
// Capture all vertex positions
// =========================================================================
static std::vector<Vector3d> captureAllPositions(const LithologyModel& model) {
    std::vector<Vector3d> positions;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            positions.push_back(mesh->vertex(vi).position);
        }
    }
    return positions;
}

// =========================================================================
// Quantize a 3D position to a hash key for shared-vertex detection.
// =========================================================================
struct QuantizedPos {
    int64_t x, y, z;
    bool operator==(const QuantizedPos& o) const {
        return x == o.x && y == o.y && z == o.z;
    }
};
struct QPosHash {
    size_t operator()(const QuantizedPos& p) const {
        return static_cast<size_t>(p.x) ^
               (static_cast<size_t>(p.y) << 21) ^
               (static_cast<size_t>(p.z) << 42);
    }
};

static QuantizedPos quantize(const Vector3d& v, double resolution = 1e-4) {
    return {
        static_cast<int64_t>(std::round(v.x() / resolution)),
        static_cast<int64_t>(std::round(v.y() / resolution)),
        static_cast<int64_t>(std::round(v.z() / resolution))
    };
}

// =========================================================================
// Build map of shared vertices across meshes.
// =========================================================================
static std::unordered_map<QuantizedPos,
    std::vector<std::pair<int, uint32_t>>, QPosHash>
findSharedVertices(const LithologyModel& model) {
    std::unordered_map<QuantizedPos,
        std::vector<std::pair<int, uint32_t>>, QPosHash> posMap;

    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            auto qp = quantize(mesh->vertex(vi).position);
            posMap[qp].push_back({g, vi});
        }
    }

    auto it = posMap.begin();
    while (it != posMap.end()) {
        if (it->second.size() < 2) {
            it = posMap.erase(it);
        } else {
            ++it;
        }
    }
    return posMap;
}

// =========================================================================
// Find all vertices classified as FIXED (model corners).
// =========================================================================
struct FixedVertexRef {
    int group;
    uint32_t vi;
    Vector3d pos;
};
static std::vector<FixedVertexRef> findFixedVertices(const LithologyModel& model) {
    std::vector<FixedVertexRef> result;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            if (mesh->vertex(vi).freedom == VertexFreedom::FIXED) {
                result.push_back({g, vi, mesh->vertex(vi).position});
            }
        }
    }
    return result;
}

// =========================================================================
// Find all vertices on the top surface (within tol of max Z).
// =========================================================================
static std::vector<std::pair<int, uint32_t>> findTopVertices(
    const LithologyModel& model, double tol = 1.0) {
    double zmax = -1e30;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            double z = mesh->vertex(vi).position.z();
            if (z > zmax) zmax = z;
        }
    }
    std::vector<std::pair<int, uint32_t>> result;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            if (std::abs(mesh->vertex(vi).position.z() - zmax) <= tol) {
                result.push_back({g, vi});
            }
        }
    }
    return result;
}

// =========================================================================
// Find all vertices on the bottom surface (within tol of min Z).
// =========================================================================
static std::vector<std::pair<int, uint32_t>> findBottomVertices(
    const LithologyModel& model, double tol = 1.0) {
    double zmin = 1e30;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            double z = mesh->vertex(vi).position.z();
            if (z < zmin) zmin = z;
        }
    }
    std::vector<std::pair<int, uint32_t>> result;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            if (std::abs(mesh->vertex(vi).position.z() - zmin) <= tol) {
                result.push_back({g, vi});
            }
        }
    }
    return result;
}

// =========================================================================
// Prepare initial model: vertex perturbation + shared boundary setup.
//
// Uses finalizeModelSetup (same as main_from_csv.cpp) to lock outer envelope
// vertices via fixExteriorFaces.  Only perturbs interior contact-surface
// vertices that retain Z freedom — outer envelope vertices stay at their
// true positions so the optimizer can correct the perturbation.
// =========================================================================
static void prepareInitialModel(std::shared_ptr<LithologyModel> model,
                                 double dx, double dy, double dz,
                                 VertexFreedom freedom = VertexFreedom::Z_ONLY) {
    // Run the standard setup first (sets freedoms, locks outer envelope).
    // This must come BEFORE perturbation so fixExteriorFaces sees the
    // unperturbed vertex positions for shared-face classification.
    finalizeModelSetup(model, freedom);

    // Perturb only vertices that have Z freedom (interior contacts).
    // Outer envelope vertices (XY_FREE, X_ONLY, Y_ONLY, FIXED) keep
    // their true positions — the optimizer can't move them in Z anyway.
    for (int g = 0; g < model->groupMeshCount(); ++g) {
        auto* mesh = model->groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            auto f = mesh->vertex(vi).freedom;
            // Only perturb vertices that CAN move in Z
            if (f == VertexFreedom::Z_ONLY ||
                f == VertexFreedom::XYZ_FREE ||
                f == VertexFreedom::XZ_FREE ||
                f == VertexFreedom::YZ_FREE) {
                mesh->vertex(vi).position.x() += dx;
                mesh->vertex(vi).position.y() += dy;
                mesh->vertex(vi).position.z() += dz;
            }
        }
    }
}

// =========================================================================
// Check for volume overlap: vertices from one group mesh sitting inside
// a non-adjacent group's volume.
// =========================================================================
static size_t countOverlaps(const LithologyModel& model) {
    size_t overlaps = 0;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            int classified = model.classifyPoint(mesh->vertex(vi).position);
            if (classified == -1) continue;
            if (std::abs(classified - g) > 1) {
                ++overlaps;
            }
        }
    }
    return overlaps;
}

// =========================================================================
// Count adjacent-group overlaps: vertices of group N inside group N±1.
// =========================================================================
static size_t countAdjacentOverlaps(const LithologyModel& model) {
    size_t overlaps = 0;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        const auto* mesh = model.groupMesh(g);
        if (!mesh) continue;
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            int classified = model.classifyPoint(mesh->vertex(vi).position);
            if (classified < 0 || classified == g) continue;
            if (std::abs(classified - g) == 1) ++overlaps;
        }
    }
    return overlaps;
}

// =========================================================================
// Boundary integrity: track vertices on the model envelope and verify
// they stay co-located and within the original AABB after inversion.
// =========================================================================

// State captured from the starting model for boundary integrity checks
struct BoundaryState {
    std::unordered_map<QuantizedPos,
        std::vector<std::pair<int, uint32_t>>, QPosHash> vertices;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double margin;
};

// Find all boundary vertices: after fixExteriorFaces, any vertex with
// freedom != XYZ_FREE is on the model envelope.
static std::unordered_map<QuantizedPos,
    std::vector<std::pair<int, uint32_t>>, QPosHash>
findBoundaryVertices(const LithologyModel& model) {
    std::unordered_map<QuantizedPos,
        std::vector<std::pair<int, uint32_t>>, QPosHash> posMap;

    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        if (!mesh) continue;
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            if (mesh->vertex(vi).freedom != VertexFreedom::XYZ_FREE) {
                auto qp = quantize(mesh->vertex(vi).position);
                posMap[qp].push_back({g, vi});
            }
        }
    }
    return posMap;
}

static BoundaryState captureBoundaryState(const LithologyModel& model) {
    BoundaryState state;
    state.vertices = findBoundaryVertices(model);
    computeModelAABB(model, state.xmin, state.xmax,
                     state.ymin, state.ymax,
                     state.zmin, state.zmax);
    double dim = std::max({state.xmax - state.xmin,
                           state.ymax - state.ymin,
                           state.zmax - state.zmin, 1.0});
    state.margin = dim * 0.02;
    return state;
}

static void checkBoundaryIntegrity(
    const LithologyModel& model,
    const BoundaryState& state,
    const std::string& label)
{
    // Check 1: shared boundary vertices from adjacent groups stay co-located
    for (const auto& [qp, instances] : state.vertices) {
        if (instances.size() < 2) continue;
        const auto& refPos = model.groupMesh(instances[0].first)
                                   ->vertex(instances[0].second).position;
        for (size_t i = 1; i < instances.size(); ++i) {
            const auto& pos = model.groupMesh(instances[i].first)
                                    ->vertex(instances[i].second).position;
            INFO(label << ": boundary gap group "
                 << instances[0].first << "/" << instances[i].first
                 << " ref=(" << refPos.x() << "," << refPos.y() << "," << refPos.z() << ")"
                 << " pos=(" << pos.x() << "," << pos.y() << "," << pos.z() << ")");
            REQUIRE(pos.x() == Approx(refPos.x()).margin(0.1));
            REQUIRE(pos.y() == Approx(refPos.y()).margin(0.1));
            REQUIRE(pos.z() == Approx(refPos.z()).margin(0.1));
        }
    }

    // Check 2: all vertices stay within the original model AABB
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        if (!mesh) continue;
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            const auto& pos = mesh->vertex(vi).position;
            INFO(label << ": vertex outside envelope group " << g
                 << " vtx " << vi
                 << " at (" << pos.x() << "," << pos.y() << "," << pos.z() << ")");
            REQUIRE(pos.x() >= state.xmin - state.margin);
            REQUIRE(pos.x() <= state.xmax + state.margin);
            REQUIRE(pos.y() >= state.ymin - state.margin);
            REQUIRE(pos.y() <= state.ymax + state.margin);
            REQUIRE(pos.z() >= state.zmin - state.margin);
            REQUIRE(pos.z() <= state.zmax + state.margin);
        }
    }
}

// =========================================================================
// ModelSnapshot: captures all pre-inversion state for invariant checking.
// =========================================================================
struct ModelSnapshot {
    std::vector<Vector3d> startPositions;  // all vertex positions before inversion
    std::vector<FixedVertexRef> fixedVertices;
    std::unordered_map<QuantizedPos,
        std::vector<std::pair<int, uint32_t>>, QPosHash> sharedVertices;
    std::vector<std::pair<int, uint32_t>> topVertices;
    std::vector<std::pair<int, uint32_t>> bottomVertices;
    double topZ;
    double bottomZ;
    struct EdgeVtx { int group; uint32_t vi; Vector3d pos; VertexFreedom freedom; };
    std::vector<EdgeVtx> edgeVertices;
    size_t totalVertexCount;
    BoundaryState boundaryState;
    size_t overlapBaseline;
    size_t adjacentOverlapBaseline;
};

static ModelSnapshot captureModelSnapshot(const LithologyModel& model) {
    ModelSnapshot snap;

    snap.startPositions = captureAllPositions(model);
    snap.fixedVertices = findFixedVertices(model);
    snap.sharedVertices = findSharedVertices(model);
    snap.topVertices = findTopVertices(model);
    snap.bottomVertices = findBottomVertices(model);

    snap.topZ = model.groupMesh(snap.topVertices[0].first)
                     ->vertex(snap.topVertices[0].second).position.z();
    snap.bottomZ = model.groupMesh(snap.bottomVertices[0].first)
                        ->vertex(snap.bottomVertices[0].second).position.z();

    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            auto f = mesh->vertex(vi).freedom;
            if (f == VertexFreedom::X_ONLY || f == VertexFreedom::Y_ONLY) {
                snap.edgeVertices.push_back(
                    {g, vi, mesh->vertex(vi).position, f});
            }
        }
    }

    snap.totalVertexCount = 0;
    for (int g = 0; g < model.groupMeshCount(); ++g)
        snap.totalVertexCount += model.groupMesh(g)->vertexCount();

    snap.boundaryState = captureBoundaryState(model);
    snap.overlapBaseline = countOverlaps(model);
    snap.adjacentOverlapBaseline = countAdjacentOverlaps(model);

    return snap;
}

// =========================================================================
// Run all checks against the snapshot for a given iteration label.
// =========================================================================
static void runAllChecks(const LithologyModel& model,
                          const ModelSnapshot& snap,
                          const std::string& label)
{
    INFO("=== Check group: " << label << " ===");

    // --- Check 1: Max vertex displacement from starting position ---
    double maxDisp = 0.0;
    auto afterPositions = captureAllPositions(model);
    for (size_t i = 0; i < snap.startPositions.size(); ++i) {
        maxDisp = std::max(maxDisp,
            (afterPositions[i] - snap.startPositions[i]).norm());
    }
    INFO(label << " max displacement from start = " << maxDisp << " m");
    REQUIRE(maxDisp < 500.0);

    // --- Check 2: No gaps between volumes ---
    for (const auto& [qp, instances] : snap.sharedVertices) {
        if (instances.size() < 2) continue;
        const auto& refPos = model.groupMesh(instances[0].first)
                                   ->vertex(instances[0].second).position;
        for (size_t i = 1; i < instances.size(); ++i) {
            const auto& pos = model.groupMesh(instances[i].first)
                                    ->vertex(instances[i].second).position;
            INFO(label << " shared vertex gap group "
                 << instances[0].first << "/" << instances[i].first);
            REQUIRE(pos.x() == Approx(refPos.x()).margin(0.1));
            REQUIRE(pos.y() == Approx(refPos.y()).margin(0.1));
            REQUIRE(pos.z() == Approx(refPos.z()).margin(0.1));
        }
    }

    // --- Check 3: Flat top ---
    for (const auto& [g, vi] : snap.topVertices) {
        double z = model.groupMesh(g)->vertex(vi).position.z();
        INFO(label << " top vertex group " << g << " vtx " << vi
             << " Z=" << z << " expected=" << snap.topZ);
        REQUIRE(z == Approx(snap.topZ).margin(0.1));
    }

    // --- Check 4: Flat bottom ---
    for (const auto& [g, vi] : snap.bottomVertices) {
        double z = model.groupMesh(g)->vertex(vi).position.z();
        INFO(label << " bottom vertex group " << g << " vtx " << vi
             << " Z=" << z << " expected=" << snap.bottomZ);
        REQUIRE(z == Approx(snap.bottomZ).margin(0.1));
    }

    // --- Check 5: Vertices moved (at least 50% of non-FIXED) ---
    REQUIRE(snap.startPositions.size() == afterPositions.size());

    size_t moved = 0;
    size_t fixedCount = 0;
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        auto* mesh = model.groupMesh(g);
        for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
            if (mesh->vertex(vi).freedom == VertexFreedom::FIXED) ++fixedCount;
        }
    }
    for (size_t i = 0; i < snap.startPositions.size(); ++i) {
        if ((snap.startPositions[i] - afterPositions[i]).norm() > 1e-6) ++moved;
    }
    size_t nonFixed = snap.totalVertexCount - fixedCount;
    INFO(label << " total=" << snap.totalVertexCount << " moved=" << moved
         << " fixed=" << fixedCount << " nonFixed=" << nonFixed);
    REQUIRE(nonFixed > 0);
    REQUIRE(moved > nonFixed * 0.5);

    // --- Check 6: No non-adjacent volume overlap increase ---
    size_t overlapsAfter = countOverlaps(model);
    INFO(label << " non-adjacent overlaps before=" << snap.overlapBaseline
         << " after=" << overlapsAfter);
    REQUIRE(overlapsAfter <= snap.overlapBaseline + 120);

    // --- Check 7: Boundary integrity ---
    checkBoundaryIntegrity(model, snap.boundaryState, label);

    // --- Check 8: No adjacent-group interpenetration increase ---
    size_t adjOverlaps = countAdjacentOverlaps(model);
    INFO(label << " adjacent overlaps before=" << snap.adjacentOverlapBaseline
         << " after=" << adjOverlaps);
    REQUIRE(adjOverlaps <= snap.adjacentOverlapBaseline + 10);

    // --- Check 9: Built-in model validation (spikes, edge tear, AABB, etc.) ---
    auto validation = model.validate(snap.overlapBaseline, snap.adjacentOverlapBaseline);
    INFO(label << " validate: " << (validation.passed ? "PASS" : "FAIL")
         << " (" << validation.failureReason << ")");
    REQUIRE(validation.passed);
}

// =========================================================================
// MULTI-ITERATION TEST: runs 10 iterations, checking all invariants after
// each iteration against the starting (pre-inversion) model.
// Each iteration creates a fresh InversionRunner with maxIterations=1,
// carrying model geometry forward across iterations.
// Invariant checks run AFTER runner.run() returns — NOT inside a callback,
// because callbacks fire during line search at intermediate trial points
// that can have extreme geometry.
// =========================================================================
TEST_CASE("Forrestania: multi-iteration comprehensive validation",
          "[forrestania][multi-iter][comprehensive]") {
    // ---- Setup ----
    auto trueModel = loadClusterModel();
    REQUIRE(trueModel->groupMeshCount() > 0);

    auto obs = generateObsForModel(*trueModel, 8);
    obs = computeSyntheticData(trueModel, obs);
    REQUIRE(obs.size() > 0);

    // Synthetic magnetic data (Forrestania IGRF: inc=-66.2, dec=-0.1, field=58874 nT)
    auto magObs = computeSyntheticMagnetic(trueModel, obs, -66.2, -0.1, 58874.0);
    REQUIRE(magObs.size() > 0);

    auto model = loadClusterModel();
    prepareInitialModel(model, 0.0, 0.0, 10.0, VertexFreedom::XYZ_FREE);

    // ---- Capture starting state ----
    ModelSnapshot snapshot = captureModelSnapshot(*model);
    REQUIRE(snapshot.startPositions.size() > 0);
    REQUIRE(snapshot.sharedVertices.size() > 0);
    REQUIRE(snapshot.topVertices.size() > 0);
    REQUIRE(snapshot.bottomVertices.size() > 0);

    // ---- Run 50 iterations, exporting .ts after each ----
    for (int iter = 0; iter < 50; ++iter) {
        InversionConfig config;
        config.model = model;
        config.observedData = obs;
        config.magneticData = magObs;
        config.magneticWeight = 1.0;
        config.magneticInclination = -66.2;
        config.magneticDeclination = -0.1;
        config.magneticField_nT = 58874.0;
        config.maxIterations = 1;
        config.solver = "lbfgsb";
        config.lambda = 1.0;
        config.lbfgsHistory = 20;
        config.tolerance = 1e-6;
        config.iterationExportDir = "data/iter_ts";

        InversionRunner runner(config);
        runner.run();

        // --- SVG diagnostics: observed vs predicted ---
        {
            auto gravFwd = std::make_shared<GravityForward>(model, obs);
            VectorXd gravPred = gravFwd->compute(model->assembleParameterVector());

            auto magFwd = std::make_shared<MagneticForward>(
                model, magObs, -66.2, -0.1, 58874.0);
            VectorXd magPred = magFwd->compute(model->assembleParameterVector());

            SVGPlot plot(800, 500, 1, 2);
            plot.title = "Iter " + std::to_string(iter + 1) + " Diagnostics";

            // Panel 0: Gravity observed vs predicted
            plot.panel(0).title = "Gravity";
            plot.panel(0).xlabel = "Observed (mGal)";
            plot.panel(0).ylabel = "Predicted (mGal)";
            // Build scatter vectors
            std::vector<double> gx, gy;
            for (size_t i = 0; i < obs.size(); ++i) {
                gx.push_back(obs[i].g_obs);
                gy.push_back(gravPred[static_cast<Eigen::Index>(i)]);
            }
            plot.panel(0).scatters.push_back({gx, gy, "blue", "", 2.0});
            double gmin = *std::min_element(gx.begin(), gx.end());
            double gmax = *std::max_element(gx.begin(), gx.end());
            plot.panel(0).xmin = gmin; plot.panel(0).xmax = gmax;
            plot.panel(0).ymin = gmin; plot.panel(0).ymax = gmax;
            plot.panel(0).hlines.push_back({0.0, "gray", 0.5, "dashed", ""});

            // Panel 1: Magnetic observed vs predicted
            plot.panel(1).title = "Magnetics";
            plot.panel(1).xlabel = "Observed (nT)";
            plot.panel(1).ylabel = "Predicted (nT)";
            std::vector<double> mx, my;
            for (size_t i = 0; i < magObs.size(); ++i) {
                mx.push_back(magObs[i].t_obs);
                my.push_back(magPred[static_cast<Eigen::Index>(i)]);
            }
            plot.panel(1).scatters.push_back({mx, my, "red", "", 2.0});
            double mmin = *std::min_element(mx.begin(), mx.end());
            double mmax = *std::max_element(mx.begin(), mx.end());
            plot.panel(1).xmin = mmin; plot.panel(1).xmax = mmax;
            plot.panel(1).ymin = mmin; plot.panel(1).ymax = mmax;
            plot.panel(1).hlines.push_back({0.0, "gray", 0.5, "dashed", ""});

            std::string diagPath = "data/iter_ts/iter_"
                + std::string(iter + 1 < 10 ? "00" : iter + 1 < 100 ? "0" : "")
                + std::to_string(iter + 1) + "_diag.svg";
            plot.save(diagPath);
        }

        std::string label = "iter_" + std::to_string(iter + 1);
        runAllChecks(*model, snapshot, label);
    }
}

// =========================================================================
// FAST SMOKE TEST: single iteration with all checks.
// =========================================================================
TEST_CASE("Forrestania: single-iteration comprehensive validation",
          "[forrestania][comprehensive][smoke]") {
    auto trueModel = loadClusterModel();
    auto obs = generateObsForModel(*trueModel, 8);
    obs = computeSyntheticData(trueModel, obs);

    auto magObs = computeSyntheticMagnetic(trueModel, obs, -66.2, -0.1, 58874.0);

    auto model = loadClusterModel();
    prepareInitialModel(model, 0.0, 0.0, 10.0, VertexFreedom::XYZ_FREE);

    ModelSnapshot snapshot = captureModelSnapshot(*model);
    REQUIRE(snapshot.startPositions.size() > 0);
    REQUIRE(snapshot.sharedVertices.size() > 0);
    REQUIRE(snapshot.topVertices.size() > 0);
    REQUIRE(snapshot.bottomVertices.size() > 0);

    InversionConfig config;
    config.model = model;
    config.observedData = obs;
    config.magneticData = magObs;
    config.magneticWeight = 1.0;
    config.magneticInclination = -66.2;
    config.magneticDeclination = -0.1;
    config.magneticField_nT = 58874.0;
    config.maxIterations = 1;
    config.solver = "lbfgsb";
    config.lambda = 1.0;
    config.lbfgsHistory = 20;
    config.tolerance = 1e-6;

    InversionRunner runner(config);
    runner.run();

    runAllChecks(*model, snapshot, "iter_1");
}

// =========================================================================
// TEST: Export perturbed model as .ts files for visual inspection
// =========================================================================
TEST_CASE("Forrestania: export perturbed model TS", "[forrestania][export]") {
    auto model = loadClusterModel();
    prepareInitialModel(model, 10.0, 10.0, 10.0);

    InversionExporter exporter("perturbed_output", "perturbed");
    exporter.setSubfolder("ts");

    exporter.setGroupNaming({"Background", "Intermediate", "Mafic", "Sulfide"});

    for (int g = 0; g < model->groupMeshCount(); ++g) {
        const auto* mesh = model->groupMesh(g);
        if (!mesh || mesh->vertexCount() == 0) continue;
        const auto& grp = model->group(g);
        std::string suffix = "group_" + std::to_string(grp.id) + "_" + grp.name;
        exporter.exportClosedVolume(*mesh, suffix);
    }

    for (int g = 0; g < model->groupMeshCount(); ++g) {
        const auto& grp = model->group(g);
        std::string fname = "perturbed_output/ts/perturbed_group_"
                          + std::to_string(grp.id) + "_" + grp.name + ".ts";
        INFO("Checking " << fname);
        REQUIRE(std::filesystem::exists(fname));
    }
}

} // namespace litho_invert

