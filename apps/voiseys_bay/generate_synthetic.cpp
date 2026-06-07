#include "generate_synthetic.h"
#include <litho_invert/forward/gravity_forward.h>
#include <litho_invert/forward/magnetic_forward.h>
#include <litho_invert/em/em_active_forward.h>
#include <litho_invert/em/em_mt_forward.h>
#include <cmath>
#include <iostream>

namespace litho_invert {

namespace {

constexpr double GRID_SIZE = 500.0;    // total grid extent in x and y (m)
constexpr double CELL_SIZE = 50.0;     // grid cell size (m)
constexpr int N_CELLS = 10;            // cells per dimension
constexpr int N_VERTS = 11;            // vertices per dimension

// Ellipsoid parameters for Surface 0 (anorthosite/troctolite boundary)
constexpr double LENS_CX = 0.0;
constexpr double LENS_CY = 0.0;
constexpr double LENS_CZ = -300.0;
constexpr double LENS_A = 150.0;   // x half-axis
constexpr double LENS_B = 75.0;    // y half-axis
constexpr double LENS_C = 50.0;    // z half-axis (vertical)

// Ellipsoid parameters for Surface 1 (troctolite/sulfide boundary)
constexpr double POD_CX = 0.0;
constexpr double POD_CY = 0.0;
constexpr double POD_CZ = -350.0;
constexpr double POD_A = 40.0;
constexpr double POD_B = 20.0;
constexpr double POD_C = 15.0;

constexpr double REGIONAL_S0 = -600.0;  // regional depth for surface 0
constexpr double REGIONAL_S1 = -600.0;  // regional depth for surface 1

// Returns true if (x,y) is inside the ellipsoid footprint
bool insideEllipsoid(double x, double y,
                     double cx, double cy,
                     double a, double b) {
    double dx = x - cx;
    double dy = y - cy;
    double val = (dx / a) * (dx / a) + (dy / b) * (dy / b);
    return val <= 1.0;
}

// Compute z on the ellipsoid surface for given (x,y).
// Returns the TOP of the ellipsoid (z = cz + c * sqrt(...)).
double ellipsoidTopZ(double x, double y,
                     double cx, double cy, double cz,
                     double a, double b, double c) {
    double dx = x - cx;
    double dy = y - cy;
    double val = (dx / a) * (dx / a) + (dy / b) * (dy / b);
    if (val > 1.0) return 0.0;  // outside footprint (should not be called)
    return cz + c * std::sqrt(1.0 - val);
}

// Create the grid vertex positions for a 500m x 500m grid with 50m spacing.
// Returns the x and y coordinates as vectors.
void makeGridCoords(std::vector<double>& xs, std::vector<double>& ys) {
    xs.resize(N_VERTS);
    ys.resize(N_VERTS);
    double half = GRID_SIZE / 2.0;
    for (int i = 0; i < N_VERTS; ++i) {
        xs[i] = -half + i * CELL_SIZE;
        ys[i] = -half + i * CELL_SIZE;
    }
}

} // anonymous namespace

// -----------------------------------------------------------------------
// Generate true synthetic model
// -----------------------------------------------------------------------
std::shared_ptr<LithologyModel> generateTrueModel() {
    auto model = std::make_shared<LithologyModel>();

    // Add lithology groups in order from surface to depth
    model->addGroup({0, "Anorthosite", 2.70, 0.001});
    model->addGroup({1, "Troctolite", 2.95, 0.02});
    model->addGroup({2, "Massive Sulfide", 4.00, 0.3});
    model->setBottomDepth(-5000.0);

    // Generate grid coordinates
    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    // --- Surface 0: Anorthosite/Troctolite boundary ---
    // Ellipsoidal lens; inside: bulges UP toward -250, outside: -600
    auto surf0 = std::make_shared<SurfaceMesh>();
    surf0->setName("anorthosite_troctolite");
    surf0->setBounds(-1000.0, 0.0);

    for (int iy = 0; iy < N_VERTS; ++iy) {
        for (int ix = 0; ix < N_VERTS; ++ix) {
            double x = xs[ix];
            double y = ys[iy];
            double z;
            if (insideEllipsoid(x, y, LENS_CX, LENS_CY, LENS_A, LENS_B)) {
                z = ellipsoidTopZ(x, y, LENS_CX, LENS_CY, LENS_CZ,
                                  LENS_A, LENS_B, LENS_C);
            } else {
                z = REGIONAL_S0;
            }
            surf0->addVertex(x, y, z, VertexFreedom::Z_ONLY);
        }
    }

    // Triangulate: two triangles per cell
    for (int iy = 0; iy < N_CELLS; ++iy) {
        for (int ix = 0; ix < N_CELLS; ++ix) {
            uint32_t i0 = static_cast<uint32_t>(iy * N_VERTS + ix);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = static_cast<uint32_t>((iy + 1) * N_VERTS + ix);
            uint32_t i3 = i2 + 1;
            surf0->addTriangle(i0, i1, i2);
            surf0->addTriangle(i1, i3, i2);
        }
    }

    // --- Surface 1: Troctolite/Sulfide boundary ---
    // Ellipsoidal pod; inside: bulges UP toward -335, outside: -600
    auto surf1 = std::make_shared<SurfaceMesh>();
    surf1->setName("troctolite_sulfide");
    surf1->setBounds(-1000.0, 0.0);

    for (int iy = 0; iy < N_VERTS; ++iy) {
        for (int ix = 0; ix < N_VERTS; ++ix) {
            double x = xs[ix];
            double y = ys[iy];
            double z;
            if (insideEllipsoid(x, y, POD_CX, POD_CY, POD_A, POD_B)) {
                z = ellipsoidTopZ(x, y, POD_CX, POD_CY, POD_CZ,
                                  POD_A, POD_B, POD_C);
            } else {
                z = REGIONAL_S1;
            }
            surf1->addVertex(x, y, z, VertexFreedom::Z_ONLY);
        }
    }

    // Same triangulation as surface 0
    for (int iy = 0; iy < N_CELLS; ++iy) {
        for (int ix = 0; ix < N_CELLS; ++ix) {
            uint32_t i0 = static_cast<uint32_t>(iy * N_VERTS + ix);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = static_cast<uint32_t>((iy + 1) * N_VERTS + ix);
            uint32_t i3 = i2 + 1;
            surf1->addTriangle(i0, i1, i2);
            surf1->addTriangle(i1, i3, i2);
        }
    }

    model->addSurface(surf0);
    model->addSurface(surf1);

    std::cout << "True model: " << model->groupCount() << " groups, "
              << model->surfaceCount() << " surfaces" << std::endl;
    return model;
}

// -----------------------------------------------------------------------
// Generate initial flat surfaces for inversion
// -----------------------------------------------------------------------
std::shared_ptr<LithologyModel> generateInitialModel() {
    auto model = std::make_shared<LithologyModel>();

    model->addGroup({0, "Anorthosite", 2.70, 0.001});
    model->addGroup({1, "Troctolite", 2.95, 0.02});
    model->addGroup({2, "Massive Sulfide", 4.00, 0.3});
    model->setBottomDepth(-5000.0);

    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    // Surface 0: flat at -400m (above surface 1, preserving correct layer order)
    auto surf0 = std::make_shared<SurfaceMesh>();
    surf0->setName("anorthosite_troctolite");
    surf0->setBounds(-1000.0, 0.0);

    for (int iy = 0; iy < N_VERTS; ++iy) {
        for (int ix = 0; ix < N_VERTS; ++ix) {
            surf0->addVertex(xs[ix], ys[iy], -400.0, VertexFreedom::Z_ONLY);
        }
    }

    for (int iy = 0; iy < N_CELLS; ++iy) {
        for (int ix = 0; ix < N_CELLS; ++ix) {
            uint32_t i0 = static_cast<uint32_t>(iy * N_VERTS + ix);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = static_cast<uint32_t>((iy + 1) * N_VERTS + ix);
            uint32_t i3 = i2 + 1;
            surf0->addTriangle(i0, i1, i2);
            surf0->addTriangle(i1, i3, i2);
        }
    }

    // Surface 1: flat at -500m (below surface 0)
    auto surf1 = std::make_shared<SurfaceMesh>();
    surf1->setName("troctolite_sulfide");
    surf1->setBounds(-1000.0, 0.0);

    for (int iy = 0; iy < N_VERTS; ++iy) {
        for (int ix = 0; ix < N_VERTS; ++ix) {
            surf1->addVertex(xs[ix], ys[iy], -500.0, VertexFreedom::Z_ONLY);
        }
    }

    for (int iy = 0; iy < N_CELLS; ++iy) {
        for (int ix = 0; ix < N_CELLS; ++ix) {
            uint32_t i0 = static_cast<uint32_t>(iy * N_VERTS + ix);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = static_cast<uint32_t>((iy + 1) * N_VERTS + ix);
            uint32_t i3 = i2 + 1;
            surf1->addTriangle(i0, i1, i2);
            surf1->addTriangle(i1, i3, i2);
        }
    }

    model->addSurface(surf0);
    model->addSurface(surf1);

    return model;
}

// -----------------------------------------------------------------------
// Generate observation points
// -----------------------------------------------------------------------
GravityData generateObservationPoints() {
    GravityData data;
    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    // Observations at vertex positions, z=0
    for (int iy = 0; iy < N_VERTS; ++iy) {
        for (int ix = 0; ix < N_VERTS; ++ix) {
            data.push_back({Vector3d(xs[ix], ys[iy], 0.0), 0.0, 0.0});
        }
    }

    return data;
}

// -----------------------------------------------------------------------
// Compute synthetic data
// -----------------------------------------------------------------------
GravityData computeSyntheticData(std::shared_ptr<LithologyModel> trueModel,
                                  const GravityData& observationPoints,
                                  double paddingDensity) {
    GravityForward forward(trueModel, observationPoints);

    // Include padding contribution so synthetic data is consistent with
    // the inversion forward model when padding is enabled.
    if (paddingDensity > 0.0) {
        forward.enablePadding(true, -100000.0);
        forward.setPaddingDensity(paddingDensity);
    }

    VectorXd params = trueModel->assembleParameterVector();
    VectorXd gravity = forward.compute(params);

    GravityData result = observationPoints;
    for (size_t i = 0; i < result.size(); ++i) {
        result[i].g_obs = gravity[static_cast<Index>(i)];
    }

    return result;
}

// -----------------------------------------------------------------------
// Compute synthetic magnetic data
// -----------------------------------------------------------------------
MagneticData computeSyntheticMagnetic(std::shared_ptr<LithologyModel> trueModel,
                                       const GravityData& observationPoints,
                                       double inc_deg, double dec_deg,
                                       double field_nT,
                                       double paddingSusceptibility) {
    MagneticData magObs;
    magObs.reserve(observationPoints.size());
    for (const auto& gp : observationPoints) {
        magObs.push_back({gp.position, 0.0, 50.0}); // 50 nT uncertainty per point
    }

    MagneticForward forward(trueModel, magObs, inc_deg, dec_deg, field_nT);

    if (paddingSusceptibility != 0.0) {
        forward.enablePadding(true, -100000.0);
        forward.setPaddingSusceptibility(paddingSusceptibility);
    }

    VectorXd params = trueModel->assembleParameterVector();
    VectorXd magResponse = forward.compute(params);

    MagneticData result = magObs;
    for (size_t i = 0; i < result.size(); ++i) {
        result[i].t_obs = magResponse[static_cast<Index>(i)];
    }

    return result;
}

// -----------------------------------------------------------------------
// Build airborne TEM geometry for testing
// -----------------------------------------------------------------------
void buildAirborneTEMGeometry(std::vector<EMSource>& sources,
                               std::vector<EMReceiver>& receivers,
                               EMConfig& config) {
    // Single coaxial source-receiver pair centered over the survey area.
    // Vertical magnetic dipole at 30m height (typical helicopter TEM bird).
    EMSource src;
    src.position = {0.0, 0.0, 30.0};
    src.type = EMSurveyType::AirborneTEM;
    src.waveform = "step";
    src.loopRadius_m = 0.0;       // treat as vertical magnetic dipole
    src.loopHeight_m = 30.0;
    sources.push_back(src);

    // Coaxial receiver (z-component, same position as source)
    EMReceiver recv;
    recv.position = {0.0, 0.0, 30.0};
    recv.component = "z";
    receivers.push_back(recv);

    // Configure time gates: early to late (seconds after turn-off)
    // Typical VTEM gates covering ~50 µs to ~10 ms
    config.solverMethod = "ie";
    config.timeGates_s = {5e-5, 1e-4, 2e-4, 4e-4, 8e-4, 1.6e-3, 3.2e-3, 6.4e-3};
    config.autoComputeSkinDepth = true;
}

// -----------------------------------------------------------------------
// Compute synthetic active EM data
// -----------------------------------------------------------------------
ActiveEMData computeSyntheticActiveEM(std::shared_ptr<LithologyModel> trueModel,
                                       const GravityData& observationPoints,
                                       const std::vector<EMSource>& sources,
                                       const std::vector<EMReceiver>& receivers,
                                       const EMConfig& config,
                                       double paddingConductivity) {
    ActiveEMData data;

    // Build one ActiveEMPoint per observation position per time gate.
    // Each uses source=0, receiver=0 (the single airborne bird).
    int nGates = static_cast<int>(config.timeGates_s.size());
    for (const auto& gp : observationPoints) {
        for (int g = 0; g < nGates; ++g) {
            ActiveEMPoint pt;
            pt.position = gp.position;
            pt.sourceIndex = 0;
            pt.receiverIndex = 0;
            pt.gateIndex = g;
            pt.timeGate_s = config.timeGates_s[g];
            pt.isFrequencyDomain = false;
            pt.d_obs = 0.0;       // filled below
            pt.d_std = 1e-12;     // small uncertainty (nT/s scale)
            pt.residualNorm = EMResidualNorm::Log10;  // wide dynamic range
            data.push_back(pt);
        }
    }

    // Build the forward model and compute synthetic responses
    EMActiveForward forward(trueModel, sources, receivers, data, config);
    if (paddingConductivity > 0.0) {
        forward.enablePadding(true, -100000.0);
        forward.setPaddingConductivity(paddingConductivity);
    }

    VectorXd params = trueModel->assembleParameterVector();
    VectorXd emResponse = forward.compute(params);

    // Write responses back into the data array
    for (size_t i = 0; i < data.size(); ++i) {
        data[i].d_obs = emResponse[static_cast<Index>(i)];
    }

    return data;
}

// -----------------------------------------------------------------------
// Build MT stations for testing
// -----------------------------------------------------------------------
void buildMTStations(std::vector<MTStation>& stations,
                      std::vector<double>& frequencies_Hz) {
    // Four stations at accessible locations (simulating ground survey).
    // Placed at survey corners — ridges/dry ground, not marshes.
    frequencies_Hz = {0.001, 0.01, 0.1, 1.0, 10.0, 100.0};  // 6 decades

    MTStation s1;
    s1.position = {-200.0, -200.0, 0.0};
    s1.name = "MT01_SW";
    s1.frequencies_Hz = frequencies_Hz;
    stations.push_back(s1);

    MTStation s2;
    s2.position = {200.0, -200.0, 0.0};
    s2.name = "MT02_SE";
    s2.frequencies_Hz = frequencies_Hz;
    stations.push_back(s2);

    MTStation s3;
    s3.position = {200.0, 200.0, 0.0};
    s3.name = "MT03_NE";
    s3.frequencies_Hz = frequencies_Hz;
    stations.push_back(s3);

    MTStation s4;
    s4.position = {-200.0, 200.0, 0.0};
    s4.name = "MT04_NW";
    s4.frequencies_Hz = frequencies_Hz;
    stations.push_back(s4);
}

// -----------------------------------------------------------------------
// Compute synthetic MT data
// -----------------------------------------------------------------------
MTData computeSyntheticMT(std::shared_ptr<LithologyModel> trueModel,
                           const std::vector<MTStation>& stations,
                           const EMConfig& config,
                           double paddingConductivity) {
    MTData data;

    // Build one MTImpedanceElement per station per frequency per tensor element.
    // We record Zxy (real and imaginary) as the primary off-diagonal element.
    // Zxx and Zyy are near-zero for a 1D Earth; Zyx = -Zxy.
    for (const auto& stn : stations) {
        for (size_t fi = 0; fi < stn.frequencies_Hz.size(); ++fi) {
            // Zxy — real
            MTImpedanceElement zxyReal;
            zxyReal.stationIndex = static_cast<int>(&stn - &stations[0]);
            zxyReal.frequencyIndex = static_cast<int>(fi);
            zxyReal.iComp = 0;  // Ex
            zxyReal.jComp = 1;  // Hy
            zxyReal.z_std = 0.05; // 5% uncertainty
            data.push_back(zxyReal);
        }
    }

    EMMTForward forward(trueModel, stations, data, config);
    if (paddingConductivity > 0.0) {
        forward.enablePadding(true, -100000.0);
        forward.setPaddingConductivity(paddingConductivity);
    }

    VectorXd params = trueModel->assembleParameterVector();
    VectorXd mtResponse = forward.compute(params);

    // Write Zxy real responses back
    for (size_t i = 0; i < data.size(); ++i) {
        data[i].zReal_obs = mtResponse[static_cast<Index>(i)];
    }

    return data;
}

} // namespace litho_invert
