#include "generate_synthetic.h"
#include <litho_invert/forward/gravity_forward.h>
#include <litho_invert/forward/magnetic_forward.h>
#include <map>
#include <set>
#include <functional>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

namespace litho_invert {

namespace {

// ---- Grid ----
constexpr double GRID_HALF = 1000.0;
constexpr double CELL_SIZE = 100.0;
constexpr int N_CELLS = 20;
constexpr int N_VERTS = 21;

// ---- Observation grid ----
constexpr double OBS_HALF = 950.0;
constexpr double OBS_CELL = 100.0;
constexpr int OBS_N = 20;

// ---- Lithology properties from Lithosplitter_out/ clustering ----
constexpr double RHO_GRANITEGNEISS = 2.68;
constexpr double CHI_GRANITEGNEISS = 0.000;
constexpr double RHO_MAFIC         = 2.93;
constexpr double CHI_MAFIC         = 0.016;
constexpr double RHO_SULFIDE       = 4.17;
constexpr double CHI_SULFIDE       = 0.080;

// ---- Surface 0: MaficComplex / GraniteGneiss ----
//
// Bowl/lopolith shape.  MaficComplex thick at Centre (0,0) and South (-115,0),
// pinching out toward edges where GraniteGneiss reaches surface.
//
// Borehole control (from Lithosplitter classified intervals):
//   BH-01_Centre (0,0):      Mafic 0–153 m  → S0 at z ≈ -153 m
//   BH-03_South  (-115,0):   Mafic 0–15 m   → S0 at z ≈ -15 m (then Sulfide)
//   BH-02_East   (300,0):    Mafic 0–27 m   → S0 at z ≈ -27 m (then GraniteGneiss)
//   BH-04_FarWest (-500,0):  GraniteGneiss at surface → S0 near z=0
//   BH-05_North  (0,200):    GraniteGneiss at surface → S0 near z=0
// ----
constexpr double S0_CX = -57.5;
constexpr double S0_CY = 0.0;
constexpr double S0_INNER_DEPTH = 150.0;  // max Mafic thickness at centre
constexpr double S0_A  = 350.0;
constexpr double S0_B  = 180.0;
constexpr double S0_REGIONAL = -5.0;       // near-surface GraniteGneiss outside bowl

// ---- Surface 1: GraniteGneiss / MassiveSulfide ----
//
// Elongated ellipsoidal pod between Centre and South marking the sulfide body.
//
// Borehole control:
//   BH-01_Centre (0,0):    Sulfide 162–180 m → S1 top near -162 m
//   BH-03_South  (-115,0): Sulfide 15–54 m and 99–105 m → complex
//   Other holes:            no sulfide → S1 below -600 m
// ----
constexpr double S1_CX = -20.0;
constexpr double S1_CY = 0.0;
constexpr double S1_CZ = -175.0;
constexpr double S1_A  = 150.0;
constexpr double S1_B  = 35.0;
constexpr double S1_C  = 15.0;
constexpr double S1_REGIONAL = -600.0;

// ---- Ellipsoid helpers ----
bool insideEllipsoid(double x, double y,
                     double cx, double cy, double a, double b) {
    double dx = x - cx;
    double dy = y - cy;
    return (dx*dx)/(a*a) + (dy*dy)/(b*b) <= 1.0;
}

double ellipsoidTopZ(double x, double y,
                     double cx, double cy, double cz,
                     double a, double b, double c) {
    double dx = x - cx;
    double dy = y - cy;
    double val = (dx*dx)/(a*a) + (dy*dy)/(b*b);
    if (val > 1.0) return 0.0;
    return cz + c * std::sqrt(1.0 - val);
}

double bowlZ(double x, double y,
            double cx, double cy, double innerDepth,
            double a, double b, double regionalZ) {
    double dx = x - cx;
    double dy = y - cy;
    double val = (dx*dx)/(a*a) + (dy*dy)/(b*b);
    if (val > 1.0) return regionalZ;
    return regionalZ - innerDepth * std::sqrt(1.0 - val);
}

void makeGridCoords(std::vector<double>& xs, std::vector<double>& ys) {
    xs.resize(N_VERTS);
    ys.resize(N_VERTS);
    for (int i = 0; i < N_VERTS; ++i) {
        xs[i] = -GRID_HALF + i * CELL_SIZE;
        ys[i] = -GRID_HALF + i * CELL_SIZE;
    }
}

std::shared_ptr<SurfaceMesh> makeTriangulatedGrid(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    const std::string& name,
    double minZ, double maxZ,
    std::function<double(double,double)> zFunc,
    VertexFreedom vf = VertexFreedom::Z_ONLY) {
    auto surf = std::make_shared<SurfaceMesh>();
    surf->setName(name);
    surf->setBounds(minZ, maxZ);

    for (int iy = 0; iy < N_VERTS; ++iy) {
        for (int ix = 0; ix < N_VERTS; ++ix) {
            double x = xs[ix];
            double y = ys[iy];
            double z = zFunc(x, y);
            surf->addVertex(x, y, z, vf);
        }
    }

    for (int iy = 0; iy < N_CELLS; ++iy) {
        for (int ix = 0; ix < N_CELLS; ++ix) {
            uint32_t i0 = static_cast<uint32_t>(iy * N_VERTS + ix);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = static_cast<uint32_t>((iy + 1) * N_VERTS + ix);
            uint32_t i3 = i2 + 1;
            surf->addTriangle(i0, i1, i2);
            surf->addTriangle(i1, i3, i2);
        }
    }
    return surf;
}

} // anonymous namespace

// ---- buildClosedMesh: convert top/bottom surface pair to a single closed mesh ----
std::shared_ptr<SurfaceMesh> buildClosedMesh(
    const SurfaceMesh& top, const SurfaceMesh& bottom)
{
    auto closed = std::make_shared<SurfaceMesh>();
    closed->setName(top.name() + "_closed");
    closed->setBounds(std::min(top.minDepth(), bottom.minDepth()),
                      std::max(top.maxDepth(), bottom.maxDepth()));

    uint32_t nv = top.vertexCount();

    // Copy vertices: top first, then bottom
    for (uint32_t v = 0; v < nv; ++v)
        closed->addVertex(top.vertex(v).position, top.vertex(v).freedom);
    for (uint32_t v = 0; v < nv; ++v)
        closed->addVertex(bottom.vertex(v).position, bottom.vertex(v).freedom);

    // Top faces (outward normal points up)
    for (uint32_t t = 0; t < top.triangleCount(); ++t) {
        const auto& tri = top.triangle(t);
        closed->addTriangle(tri.v0, tri.v1, tri.v2);
    }

    // Bottom faces (reversed — outward normal points down)
    for (uint32_t t = 0; t < bottom.triangleCount(); ++t) {
        const auto& tri = bottom.triangle(t);
        closed->addTriangle(nv + tri.v0, nv + tri.v2, nv + tri.v1);
    }

    // Side walls: find boundary edges of the top surface and connect to bottom
    std::map<std::pair<uint32_t,uint32_t>, int> orientedCount;
    for (uint32_t t = 0; t < top.triangleCount(); ++t) {
        const auto& tri = top.triangle(t);
        orientedCount[{tri.v0, tri.v1}]++;
        orientedCount[{tri.v1, tri.v2}]++;
        orientedCount[{tri.v2, tri.v0}]++;
    }

    std::set<std::pair<uint32_t,uint32_t>> processed;
    for (const auto& kv : orientedCount) {
        uint32_t a = kv.first.first, b = kv.first.second;
        int count = kv.second;
        auto undirected = std::make_pair(std::min(a, b), std::max(a, b));
        if (processed.count(undirected)) continue;
        processed.insert(undirected);

        auto revIt = orientedCount.find({b, a});
        int revCount = (revIt != orientedCount.end()) ? revIt->second : 0;

        if (count + revCount == 1) {
            // Boundary edge. Use the oriented version that exists.
            uint32_t at, bt;
            if (count == 1) { at = a; bt = b; }
            else            { at = b; bt = a; }
            uint32_t ab = nv + at, bb = nv + bt;
            closed->addTriangle(at, bt, bb);
            closed->addTriangle(at, bb, ab);
        }
    }

    closed->buildNeighbors();
    return closed;
}

// ---- buildFlatMesh: create a flat surface matching the reference topology ----
std::shared_ptr<SurfaceMesh> buildFlatMesh(
    const SurfaceMesh& ref, double z, const std::string& name)
{
    auto flat = std::make_shared<SurfaceMesh>();
    flat->setName(name);
    flat->setBounds(ref.minDepth(), ref.maxDepth());
    for (uint32_t v = 0; v < ref.vertexCount(); ++v) {
        const auto& pos = ref.vertex(v).position;
        flat->addVertex(pos.x(), pos.y(), z, ref.vertex(v).freedom);
    }
    for (uint32_t t = 0; t < ref.triangleCount(); ++t) {
        const auto& tri = ref.triangle(t);
        flat->addTriangle(tri.v0, tri.v1, tri.v2);
    }
    return flat;
}

// =======================================================================
// True model
// =======================================================================
std::shared_ptr<LithologyModel> generateTrueModel() {
    auto model = std::make_shared<LithologyModel>();

    // 4 groups: upper Mafic → Granite → lower Mafic → Sulfide
    // Lower Mafic has same properties as upper Mafic (same lithology, two layers)
    model->addGroup({0, "MaficComplex_1", RHO_MAFIC,         CHI_MAFIC});
    model->addGroup({1, "GraniteGneiss",  RHO_GRANITEGNEISS, CHI_GRANITEGNEISS});
    model->addGroup({2, "MaficComplex_2", RHO_MAFIC,         CHI_MAFIC});
    model->addGroup({3, "MassiveSulfide", RHO_SULFIDE,       CHI_SULFIDE});
    model->setBottomDepth(-5000.0);

    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    // Surface 0: MaficComplex_1 / GraniteGneiss (bowl — Mafic thick at centre)
    auto surf0 = makeTriangulatedGrid(xs, ys,
        "mafic1_granitegneiss", -10000.0, 100.0,
        [](double x, double y) -> double {
            return bowlZ(x, y, S0_CX, S0_CY, S0_INNER_DEPTH,
                        S0_A, S0_B, S0_REGIONAL);
        });

    // Surface 1: GraniteGneiss / MaficComplex_2 (bottom of Granite)
    auto surf1 = makeTriangulatedGrid(xs, ys,
        "granitegneiss_mafic2", -10000.0, 100.0,
        [](double x, double y) -> double {
            if (insideEllipsoid(x, y, S1_CX, S1_CY, S1_A, S1_B)) {
                return ellipsoidTopZ(x, y, S1_CX, S1_CY, S1_CZ,
                                     S1_A, S1_B, S1_C);
            }
            return S1_REGIONAL;
        });

    // Surface 2: MaficComplex_2 / MassiveSulfide
    // In the true model, lower Mafic is pinched out (S2 == S1 everywhere)
    auto surf2 = makeTriangulatedGrid(xs, ys,
        "mafic2_sulfide", -10000.0, 100.0,
        [&](double x, double y) -> double {
            // Copy S1 — lower Mafic has zero thickness in true model
            if (insideEllipsoid(x, y, S1_CX, S1_CY, S1_A, S1_B)) {
                return ellipsoidTopZ(x, y, S1_CX, S1_CY, S1_CZ,
                                     S1_A, S1_B, S1_C);
            }
            return S1_REGIONAL;
        });

    auto flatTop = buildFlatMesh(*surf0, 0.0, "flat_top");
    auto flatBottom = buildFlatMesh(*surf0, model->bottomDepth(), "flat_bottom");
    model->setGroupMesh(0, buildClosedMesh(*flatTop, *surf0));
    model->setGroupMesh(1, buildClosedMesh(*surf0, *surf1));
    model->setGroupMesh(2, buildClosedMesh(*surf1, *surf2));
    model->setGroupMesh(3, buildClosedMesh(*surf2, *flatBottom));

    // Diagnostic
    std::cout << "True model: " << model->groupCount() << " groups, "
              << model->groupMeshCount() << " group meshes" << std::endl;
    std::cout << "Group classification at borehole positions:" << std::endl;
    struct BhPos { const char* name; double x, y; };
    for (auto& bh : {BhPos{"BH-01_Centre", 0.0, 0.0},
                     BhPos{"BH-02_East", 300.0, 0.0},
                     BhPos{"BH-03_South", -115.0, 0.0},
                     BhPos{"BH-04_FarWest", -500.0, 0.0},
                     BhPos{"BH-05_North", 0.0, 200.0}}) {
        int gSurf = model->classifyPoint(Vector3d(bh.x, bh.y, -5.0));
        int gMid  = model->classifyPoint(Vector3d(bh.x, bh.y, -80.0));
        int gDeep = model->classifyPoint(Vector3d(bh.x, bh.y, -170.0));
        std::cout << "  " << bh.name << " (" << bh.x << "," << bh.y
                  << "): z=-5m group=" << gSurf
                  << " z=-80m group=" << gMid
                  << " z=-170m group=" << gDeep << std::endl;
    }

    return model;
}

// =======================================================================
// True model with remanent magnetization
// =======================================================================
std::shared_ptr<LithologyModel> generateTrueModelWithRemanence() {
    auto model = std::make_shared<LithologyModel>();

    LithoGroup g0(0, "MaficComplex_1", RHO_MAFIC, CHI_MAFIC);
    g0.remanence_magnitude = 0.5;    // A/m, Q ≈ 1.5
    g0.remanence_inclination = 65.0; // slightly shallower than induced
    g0.remanence_declination = -30.0; // more westerly
    model->addGroup(g0);

    LithoGroup g1(1, "GraniteGneiss",  RHO_GRANITEGNEISS, CHI_GRANITEGNEISS);
    g1.remanence_magnitude = 0.0;  // non-magnetic, no remanence
    g1.remanence_inclination = 75.0;
    g1.remanence_declination = -20.0;
    model->addGroup(g1);

    // Lower Mafic: same remanence as upper Mafic
    LithoGroup g2(2, "MaficComplex_2", RHO_MAFIC, CHI_MAFIC);
    g2.remanence_magnitude = 0.5;
    g2.remanence_inclination = 65.0;
    g2.remanence_declination = -30.0;
    model->addGroup(g2);

    LithoGroup g3(3, "MassiveSulfide", RHO_SULFIDE, CHI_SULFIDE);
    g3.remanence_magnitude = 4.0;    // A/m, Q ≈ 4.5
    g3.remanence_inclination = 45.0; // distinctly shallower (different paleomag pole)
    g3.remanence_declination = -60.0; // strongly westerly
    model->addGroup(g3);

    model->setBottomDepth(-5000.0);

    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    // Surface 0: MaficComplex_1 / GraniteGneiss (bowl — Mafic thick at centre)
    auto surf0 = makeTriangulatedGrid(xs, ys,
        "mafic1_granitegneiss", -10000.0, 100.0,
        [](double x, double y) -> double {
            return bowlZ(x, y, S0_CX, S0_CY, S0_INNER_DEPTH,
                        S0_A, S0_B, S0_REGIONAL);
        });

    // Surface 1: GraniteGneiss / MaficComplex_2 (bottom of Granite)
    auto surf1 = makeTriangulatedGrid(xs, ys,
        "granitegneiss_mafic2", -10000.0, 100.0,
        [](double x, double y) -> double {
            if (insideEllipsoid(x, y, S1_CX, S1_CY, S1_A, S1_B)) {
                return ellipsoidTopZ(x, y, S1_CX, S1_CY, S1_CZ,
                                     S1_A, S1_B, S1_C);
            }
            return S1_REGIONAL;
        });

    // Surface 2: MaficComplex_2 / MassiveSulfide (pinched in true model)
    auto surf2 = makeTriangulatedGrid(xs, ys,
        "mafic2_sulfide", -10000.0, 100.0,
        [](double x, double y) -> double {
            if (insideEllipsoid(x, y, S1_CX, S1_CY, S1_A, S1_B)) {
                return ellipsoidTopZ(x, y, S1_CX, S1_CY, S1_CZ,
                                     S1_A, S1_B, S1_C);
            }
            return S1_REGIONAL;
        });

    auto flatTop = buildFlatMesh(*surf0, 0.0, "flat_top");
    auto flatBottom = buildFlatMesh(*surf0, model->bottomDepth(), "flat_bottom");
    model->setGroupMesh(0, buildClosedMesh(*flatTop, *surf0));
    model->setGroupMesh(1, buildClosedMesh(*surf0, *surf1));
    model->setGroupMesh(2, buildClosedMesh(*surf1, *surf2));
    model->setGroupMesh(3, buildClosedMesh(*surf2, *flatBottom));

    std::cout << "True model with remanence: " << model->groupCount() << " groups, "
              << model->groupMeshCount() << " group meshes" << std::endl;
    for (int g = 0; g < model->groupCount(); ++g) {
        const auto& grp = model->group(g);
        std::cout << "  [" << grp.id << "] " << grp.name
                  << ": chi=" << std::scientific << std::setprecision(3) << grp.susceptibility
                  << "  M_rem=" << std::fixed << std::setprecision(1) << grp.remanence_magnitude
                  << " A/m  I_rem=" << grp.remanence_inclination
                  << "  D_rem=" << grp.remanence_declination << std::endl;
    }

    return model;
}

// =======================================================================
// Initial flat model
// =======================================================================
std::shared_ptr<LithologyModel> generateInitialModel() {
    auto model = std::make_shared<LithologyModel>();

    model->addGroup({0, "MaficComplex_1", RHO_MAFIC,         CHI_MAFIC});
    model->addGroup({1, "GraniteGneiss",  RHO_GRANITEGNEISS, CHI_GRANITEGNEISS});
    model->addGroup({2, "MaficComplex_2", RHO_MAFIC,         CHI_MAFIC});
    model->addGroup({3, "MassiveSulfide", RHO_SULFIDE,       CHI_SULFIDE});
    model->setBottomDepth(-5000.0);

    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    // Surface 0: flat at -100 m (Mafic_1 → GraniteGneiss)
    auto surf0 = makeTriangulatedGrid(xs, ys,
        "mafic1_granitegneiss", -10000.0, 100.0,
        [](double, double) { return -100.0; });

    // Surface 1: flat at -200 m (GraniteGneiss → Mafic_2)
    auto surf1 = makeTriangulatedGrid(xs, ys,
        "granitegneiss_mafic2", -10000.0, 100.0,
        [](double, double) { return -200.0; });

    // Surface 2: flat at -300 m (Mafic_2 → Sulfide)
    auto surf2 = makeTriangulatedGrid(xs, ys,
        "mafic2_sulfide", -10000.0, 100.0,
        [](double, double) { return -300.0; });

    auto flatTop = buildFlatMesh(*surf0, 0.0, "flat_top");
    auto flatBottom = buildFlatMesh(*surf0, model->bottomDepth(), "flat_bottom");
    model->setGroupMesh(0, buildClosedMesh(*flatTop, *surf0));
    model->setGroupMesh(1, buildClosedMesh(*surf0, *surf1));
    model->setGroupMesh(2, buildClosedMesh(*surf1, *surf2));
    model->setGroupMesh(3, buildClosedMesh(*surf2, *flatBottom));

    return model;
}

// =======================================================================
// Interpolated initial model from borehole contact depths
// =======================================================================
std::shared_ptr<LithologyModel> generateInterpolatedInitialModel(
    int interiorDim, double cellSize, double gridHalf)
{
    auto model = std::make_shared<LithologyModel>();

    model->addGroup({0, "MaficComplex_1", RHO_MAFIC,         CHI_MAFIC});
    model->addGroup({1, "GraniteGneiss",  RHO_GRANITEGNEISS, CHI_GRANITEGNEISS});
    model->addGroup({2, "MaficComplex_2", RHO_MAFIC,         CHI_MAFIC});
    model->addGroup({3, "MassiveSulfide", RHO_SULFIDE,       CHI_SULFIDE});
    model->setBottomDepth(-5000.0);

    // Build vertex positions
    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    struct ContactPt { double x, y, z; };

    // ---- Surface 0: MaficComplex_1 / GraniteGneiss (top of Granite) ----
    // Same as before — depth where GraniteGneiss first appears below upper Mafic.
    std::vector<ContactPt> s0Contacts = {
        {   0.0,   0.0, -153.0},  // BH-01_Centre: GraniteGneiss at 153–162 m
        { 300.0,   0.0,  -27.0},  // BH-02_East:   GraniteGneiss at 27–81 m
        {-115.0,   0.0,  -15.0},  // BH-03_South:  Mafic ends at 15 m
        {-500.0,   0.0,   -5.0},  // BH-04_FarWest: GraniteGneiss at surface
        {   0.0, 200.0,   -5.0},  // BH-05_North:  GraniteGneiss at surface
        {-1000.0, -1000.0, -5.0},
        { 1000.0, -1000.0, -5.0},
        {-1000.0,  1000.0, -5.0},
        { 1000.0,  1000.0, -5.0},
    };

    // ---- Surface 1: GraniteGneiss / MaficComplex_2 (bottom of Granite) ----
    // Where lower Mafic exists (BH-02), this is the top of lower Mafic.
    // Where lower Mafic is absent, S1 meets S2 at the Sulfide top.
    std::vector<ContactPt> s1Contacts = {
        {   0.0,   0.0, -162.0},  // BH-01: Granite bottom = Sulfide top (no lower Mafic)
        { 300.0,   0.0,  -81.0},  // BH-02: Granite bottom = lower Mafic top
        {-115.0,   0.0,  -57.0},  // BH-03: first Granite ends at 57 m
        {-500.0,   0.0, -500.0},  // BH-04: GraniteGneiss only, deep
        {   0.0, 200.0, -500.0},  // BH-05: GraniteGneiss only, deep
        {-1000.0, -1000.0, -500.0},
        { 1000.0, -1000.0, -500.0},
        {-1000.0,  1000.0, -500.0},
        { 1000.0,  1000.0, -500.0},
    };

    // ---- Surface 2: MaficComplex_2 / MassiveSulfide (top of Sulfide) ----
    // Where lower Mafic is absent, S2 = S1 (pinched, Sulfide directly below Granite).
    // Where lower Mafic exists (BH-02), S2 is the base of lower Mafic.
    std::vector<ContactPt> s2Contacts = {
        {   0.0,   0.0, -162.0},  // BH-01: Sulfide top (S2 = S1, lower Mafic absent)
        { 300.0,   0.0, -500.0},  // BH-02: borehole ends in Mafic at 180 m, Sulfide deeper
        {-115.0,   0.0, -105.0},  // BH-03: second Sulfide at 99–105 m
        {-500.0,   0.0, -500.0},  // BH-04: no Sulfide
        {   0.0, 200.0, -500.0},  // BH-05: no Sulfide
        {-1000.0, -1000.0, -500.0},
        { 1000.0, -1000.0, -500.0},
        {-1000.0,  1000.0, -500.0},
        { 1000.0,  1000.0, -500.0},
    };

    // IDW interpolation with power=2
    auto idw = [](double x, double y, const std::vector<ContactPt>& contacts) -> double {
        double sumW = 0.0, sumWZ = 0.0;
        constexpr double minDist = 1.0;
        for (const auto& c : contacts) {
            double dx = x - c.x;
            double dy = y - c.y;
            double d2 = dx*dx + dy*dy;
            if (d2 < minDist*minDist) d2 = minDist*minDist;
            double w = 1.0 / d2;
            sumW += w;
            sumWZ += w * c.z;
        }
        return sumWZ / sumW;
    };

    // Build surfaces
    auto surf0 = makeTriangulatedGrid(xs, ys,
        "mafic1_granitegneiss", -10000.0, 100.0,
        [&](double x, double y) { return idw(x, y, s0Contacts); });

    auto surf1 = makeTriangulatedGrid(xs, ys,
        "granitegneiss_mafic2", -10000.0, 100.0,
        [&](double x, double y) { return idw(x, y, s1Contacts); });

    auto surf2 = makeTriangulatedGrid(xs, ys,
        "mafic2_sulfide", -10000.0, 100.0,
        [&](double x, double y) { return idw(x, y, s2Contacts); });

    auto flatTop = buildFlatMesh(*surf0, 0.0, "flat_top");
    auto flatBottom = buildFlatMesh(*surf0, model->bottomDepth(), "flat_bottom");
    model->setGroupMesh(0, buildClosedMesh(*flatTop, *surf0));
    model->setGroupMesh(1, buildClosedMesh(*surf0, *surf1));
    model->setGroupMesh(2, buildClosedMesh(*surf1, *surf2));
    model->setGroupMesh(3, buildClosedMesh(*surf2, *flatBottom));

    std::cout << "Interpolated initial model: " << model->groupCount() << " groups, "
              << model->groupMeshCount() << " group meshes" << std::endl;
    std::cout << "  Surface 0 (mafic1_granitegneiss) from "
              << s0Contacts.size() << " control points" << std::endl;
    std::cout << "  Surface 1 (granitegneiss_mafic2) from "
              << s1Contacts.size() << " control points" << std::endl;
    std::cout << "  Surface 2 (mafic2_sulfide) from "
              << s2Contacts.size() << " control points" << std::endl;

    // Diagnostic: surface depths at boreholes
    std::cout << "Surface depths at borehole positions:" << std::endl;
    struct Bh { double x, y; const char* name; };
    for (auto& bh : {Bh{0.0, 0.0, "BH-01_Centre"},
                     Bh{300.0, 0.0, "BH-02_East"},
                     Bh{-115.0, 0.0, "BH-03_South"},
                     Bh{-500.0, 0.0, "BH-04_FarWest"},
                     Bh{0.0, 200.0, "BH-05_North"}}) {
        double z0 = idw(bh.x, bh.y, s0Contacts);
        double z1 = idw(bh.x, bh.y, s1Contacts);
        double z2 = idw(bh.x, bh.y, s2Contacts);
        std::cout << "  " << bh.name << " (" << bh.x << "," << bh.y
                  << "): S0=" << std::fixed << std::setprecision(1) << z0
                  << "m  S1=" << z1 << "m  S2=" << z2 << "m" << std::endl;
    }

    return model;
}

// =======================================================================
// Observation points
// =======================================================================
GravityData generateObservationPoints() {
    GravityData data;
    for (int iy = 0; iy < OBS_N; ++iy) {
        for (int ix = 0; ix < OBS_N; ++ix) {
            double x = -OBS_HALF + ix * OBS_CELL;
            double y = -OBS_HALF + iy * OBS_CELL;
            data.push_back({Vector3d(x, y, 0.0), 0.0, 0.01});
        }
    }
    return data;
}

// =======================================================================
// Synthetic gravity
// =======================================================================
GravityData computeSyntheticData(std::shared_ptr<LithologyModel> trueModel,
                                  const GravityData& observationPoints) {
    GravityForward forward(trueModel, observationPoints);
    VectorXd params = trueModel->assembleParameterVector();
    VectorXd gravity = forward.compute(params);

    GravityData result = observationPoints;
    for (size_t i = 0; i < result.size(); ++i) {
        result[i].g_obs = gravity[static_cast<Index>(i)];
    }

    return result;
}

// =======================================================================
// Synthetic magnetics
// =======================================================================
MagneticData computeSyntheticMagnetic(std::shared_ptr<LithologyModel> trueModel,
                                       const GravityData& observationPoints,
                                       double inc_deg, double dec_deg,
                                       double field_nT,
                                       RemanentMagnetizationMode mode) {
    MagneticData magObs;
    magObs.reserve(observationPoints.size());
    for (const auto& gp : observationPoints) {
        magObs.push_back({gp.position, 0.0, 50.0});
    }

    MagneticForward forward(trueModel, magObs, inc_deg, dec_deg, field_nT);
    forward.setRemanenceMode(mode);

    VectorXd params = trueModel->assembleParameterVector();
    VectorXd magResponse = forward.compute(params);

    MagneticData result = magObs;
    for (size_t i = 0; i < result.size(); ++i) {
        result[i].t_obs = magResponse[static_cast<Index>(i)];
    }

    return result;
}

// =======================================================================
// Lithosplitter constraints — derived from real cluster intervals
// =======================================================================
//
// cluster_id → litho_group_id mapping (4 groups with MaficLower repeated):
//   2 (Intermediate,Norite+Gabbro) → 0 MaficComplex_1 (top)
//   3 (Intermediate,Troctolite)    → 0 MaficComplex_1
//   1 (Felsic Granite)    → 1 GraniteGneiss (middle)
//   2 (Intermediate,Norite+Gabbro) → 2 MaficComplex_2 (lower, same props)
//   3 (Intermediate,Troctolite)    → 2 MaficComplex_2
//   0 (MassiveSulfide)    → 3 MassiveSulfide (bottom)
//   4 (MassiveSulfide)    → 3 MassiveSulfide
//
// Constraints extracted from Lithosplitter_out/ BH-*_classified.csv
// by merging consecutive intervals with the same mapped group.

std::vector<Constraint> generateLithosplitterConstraints() {
    std::vector<Constraint> constraints;

    auto addConstraint = [&](double x, double y,
                              double z_top, double z_bottom,
                              int litho_group_id) {
        Constraint c;
        c.position = Vector3d(x, y, 0.0);
        c.z_top = z_top;
        c.z_bottom = z_bottom;
        c.litho_group_id = litho_group_id;
        constraints.push_back(c);
    };

    // ===================================================================
    // BH-01_Centre  (0, 0)  — max depth 180 m
    //
    //  From BH-1_classified.csv:
    //    0–153 m: C2+C3 → MaficComplex_1 (0) [Norite+Troctolite]
    //    153–162 m: C1 → GraniteGneiss (1) [Felsic Granite]
    //    162–180 m: C0+C4 → MassiveSulfide (3)
    // ===================================================================
    addConstraint(0.0, 0.0, 0.0, 153.0, 0);    // MaficComplex_1
    addConstraint(0.0, 0.0, 153.0, 162.0, 1);  // GraniteGneiss
    addConstraint(0.0, 0.0, 162.0, 180.0, 3);  // MassiveSulfide

    // ===================================================================
    // BH-02_East  (300, 0)  — max depth 180 m
    //
    //  From BH-2_classified.csv:
    //    0–27 m: C2 → MaficComplex_1 (0) [Gabbro+Norite]
    //    27–81 m: C1 → GraniteGneiss (1) [Felsic Granite]
    //    81–180 m: C2+C3 → MaficComplex_2 (2) [Norite+Troctolite]
    // ===================================================================
    addConstraint(300.0, 0.0, 0.0, 27.0, 0);     // MaficComplex_1
    addConstraint(300.0, 0.0, 27.0, 81.0, 1);    // GraniteGneiss
    addConstraint(300.0, 0.0, 81.0, 180.0, 2);   // MaficComplex_2

    // ===================================================================
    // BH-03_South  (-115, 0)  — max depth 180 m
    //
    //  From BH-3_classified.csv (complex alternating lithologies):
    //    0–15 m: C2 → MaficComplex_1 (0) [Norite]
    //    15–54 m: C0+C4 → MassiveSulfide (3)
    //    54–57 m: C1 → GraniteGneiss (1)
    //    57–78 m: C3 → MaficComplex_2 (2) [Troctolite, below Granite]
    //    78–99 m: C1 → GraniteGneiss (1)
    //    99–105 m: C0+C4 → MassiveSulfide (3)
    //    105–114 m: C2 → MaficComplex_2 (2)
    //    114–150 m: C1 → GraniteGneiss (1)
    //    150–180 m: C2+C3 → MaficComplex_2 (2)
    // ===================================================================
    addConstraint(-115.0, 0.0, 0.0, 15.0, 0);     // MaficComplex_1
    addConstraint(-115.0, 0.0, 15.0, 54.0, 3);    // MassiveSulfide
    addConstraint(-115.0, 0.0, 54.0, 57.0, 1);    // GraniteGneiss
    addConstraint(-115.0, 0.0, 57.0, 78.0, 2);    // MaficComplex_2 (below Granite)
    addConstraint(-115.0, 0.0, 78.0, 99.0, 1);    // GraniteGneiss
    addConstraint(-115.0, 0.0, 99.0, 105.0, 3);   // MassiveSulfide
    addConstraint(-115.0, 0.0, 105.0, 114.0, 2);  // MaficComplex_2
    addConstraint(-115.0, 0.0, 114.0, 150.0, 1);  // GraniteGneiss
    addConstraint(-115.0, 0.0, 150.0, 180.0, 2);  // MaficComplex_2

    // ===================================================================
    // BH-04_FarWest  (-500, 0)  — max depth 180 m
    //
    //  No classified data available; inferred from regional mapping.
    //  0–180 m: Felsic Granite → GraniteGneiss (1)
    // ===================================================================
    addConstraint(-500.0, 0.0, 0.0, 180.0, 1);  // GraniteGneiss

    // ===================================================================
    // BH-05_North  (0, 200)  — max depth 180 m
    //
    //  No classified data available; inferred from regional mapping.
    //  0–180 m: Felsic Granite → GraniteGneiss (1)
    // ===================================================================
    addConstraint(0.0, 200.0, 0.0, 180.0, 1);  // GraniteGneiss

    // Report
    std::cout << "Lithosplitter constraints: " << constraints.size()
              << " intervals from 5 boreholes" << std::endl;
    std::cout << "  BH-01_Centre: 3 intervals (Mafic+Felsic+Sulfide)"
              << std::endl;
    std::cout << "  BH-02_East:   3 intervals (Mafic+Felsic+Mafic)"
              << std::endl;
    std::cout << "  BH-03_South:  9 intervals (alternating lithologies)"
              << std::endl;
    std::cout << "  BH-04_FarWest: 1 interval (GraniteGneiss)"
              << std::endl;
    std::cout << "  BH-05_North:  1 interval  (GraniteGneiss)"
              << std::endl;

    return constraints;
}

// =========================================================================
// makePaddedGrid — regular N×N grid with optional padding rings
// =========================================================================

std::shared_ptr<SurfaceMesh> makePaddedGrid(
    const std::string& name,
    double flatZ,
    int interiorDim,
    int paddingRings,
    double cellSize,
    double gridHalf,
    VertexFreedom vf)
{
    auto surf = std::make_shared<SurfaceMesh>();
    surf->setName(name);
    surf->setBounds(-10000.0, 100.0);

    int fullDim = interiorDim + 2 * paddingRings;
    double fullHalf = gridHalf + paddingRings * cellSize;
    int nVerts = fullDim;
    int nCells = fullDim - 1;

    // Build grid positions
    std::vector<double> xs(nVerts), ys(nVerts);
    for (int i = 0; i < nVerts; ++i) {
        xs[i] = -fullHalf + i * cellSize;
        ys[i] = -fullHalf + i * cellSize;
    }

    // Create vertices: interior = Z_ONLY, padding = FIXED
    for (int iy = 0; iy < nVerts; ++iy) {
        for (int ix = 0; ix < nVerts; ++ix) {
            bool isPad = (ix < paddingRings || ix >= paddingRings + interiorDim
                       || iy < paddingRings || iy >= paddingRings + interiorDim);
            VertexFreedom f = isPad ? VertexFreedom::FIXED : vf;
            surf->addVertex(xs[ix], ys[iy], flatZ, f);
        }
    }

    // Create triangles (same pattern for padded and non-padded)
    for (int iy = 0; iy < nCells; ++iy) {
        for (int ix = 0; ix < nCells; ++ix) {
            uint32_t i0 = static_cast<uint32_t>(iy * nVerts + ix);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = static_cast<uint32_t>((iy + 1) * nVerts + ix);
            uint32_t i3 = i2 + 1;
            surf->addTriangle(i0, i1, i2);
            surf->addTriangle(i1, i3, i2);
        }
    }

    // Register padding info
    surf->setInteriorGrid(interiorDim, paddingRings);

    return surf;
}

// =======================================================================
// Cluster ID constraints — read from BH-*_classified.csv, remap to
// geological group indices so litho_group_id matches the model.
//
// BH (Lithosplitter) cluster_id → geological group index:
//   BH 0 (Troctolite, rho~2.90)     → geo 0 (inv c3, top mafic)
//   BH 1 (GraniteGneiss, rho~2.65)  → geo 2 (inv c1, middle)
//   BH 2 (MassiveSulfide, rho~4.1)  → geo 3 (inv c0, sulfide hi-chi)
//   BH 3 (Gabbro, rho~2.89)         → geo 1 (inv c2, lower mafic)
//   BH 4 (GraniteGneiss, rho~2.68)  → geo 2 (inv c1, middle)
// =======================================================================
std::vector<Constraint> generateClusterIdConstraints(
    const std::vector<ClusterProperties>& invClusters,
    const std::vector<int>& ordering)
{
    std::vector<Constraint> constraints;

    // Geological ordering (top→bottom): from boreholes, or cluster_id ascending
    std::vector<int> kGeoToInv = ordering;
    if (kGeoToInv.empty()) {
        std::set<int> uniqueIds;
        for (const auto& c : invClusters) uniqueIds.insert(c.cluster_id);
        kGeoToInv.assign(uniqueIds.begin(), uniqueIds.end());
    }

    // Step 1: collect per-cluster density/susceptibility from BH CSVs
    struct BhAccum { double rho_sum = 0, chi_sum = 0; int count = 0; };
    BhAccum bhAccum[5];

    const char* bhPaths[] = {
        "../../Lithosplitter_out/BH-01_Centre_classified.csv",
        "../../Lithosplitter_out/BH-02_East_classified.csv",
        "../../Lithosplitter_out/BH-03_South_classified.csv",
    };

    for (const char* fpath : bhPaths) {
        std::ifstream file(fpath);
        if (!file.is_open()) continue;
        std::string line;
        int lineno = 0;
        while (std::getline(file, line)) {
            ++lineno;
            if (lineno == 1 || line.empty() || line[0] == '#') continue;
            // Lithosplitter_out columns: 0=hole_id,1=from_m,2=to_m,
            // 3=mid_x,4=mid_y,5=mid_z,6=lith,7=cluster_id,...
            // 42=density_gcc, 43=susceptibility_si
            int col = 0, cid = -1;
            double rho = -1, chi = -1;
            std::stringstream ss(line);
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                switch (col) {
                case 7:  cid = (tok.find('.') == std::string::npos) ? std::stoi(tok) : -1; break;
                case 42: rho = std::stod(tok); break;
                case 43: chi = std::stod(tok); break;
                }
                ++col;
            }
            if (cid >= 0 && cid < 5 && rho > 0 && chi >= 0) {
                bhAccum[cid].rho_sum += rho;
                bhAccum[cid].chi_sum += chi;
                bhAccum[cid].count++;
            }
        }
    }

    // Step 2: build inv cluster lookup
    double invRho[5] = {}, invChi[5] = {};
    for (const auto& c : invClusters) {
        if (c.cluster_id >= 0 && c.cluster_id < 5) {
            invRho[c.cluster_id] = c.density_median;
            invChi[c.cluster_id] = c.susceptibility_median;
        }
    }

    // Step 3: match each BH cluster to closest inv cluster
    //         distance = |Δrho| + 10·|Δchi|  (normalises χ to ρ scale)
    int bhToInv[5] = {-1, -1, -1, -1, -1};
    for (int bh = 0; bh < 5; ++bh) {
        if (bhAccum[bh].count == 0) continue;
        double bhRho = bhAccum[bh].rho_sum / bhAccum[bh].count;
        double bhChi = bhAccum[bh].chi_sum / bhAccum[bh].count;
        int best = -1;
        double bestDist = 1e18;
        for (int inv = 0; inv < 5; ++inv) {
            if (invRho[inv] == 0.0) continue;
            double dist = std::abs(bhRho - invRho[inv])
                        + 10.0 * std::abs(bhChi - invChi[inv]);
            if (dist < bestDist) { bestDist = dist; best = inv; }
        }
        bhToInv[bh] = best;
        std::cout << "  BH c" << bh << " (rho=" << std::fixed << std::setprecision(3)
                  << bhRho << " chi=" << bhChi << ") → inv c"
                  << best << " (rho=" << invRho[best] << " chi=" << invChi[best]
                  << ") dist=" << bestDist << std::endl;
    }

    // Step 4: convert inv cluster → geo group
    int maxInvCid = 0;
    for (int cid : kGeoToInv) if (cid > maxInvCid) maxInvCid = cid;
    std::vector<int> invToGeo(maxInvCid + 1, 0);
    for (int gi = 0; gi < static_cast<int>(kGeoToInv.size()); ++gi)
        invToGeo[kGeoToInv[gi]] = gi;

    int nGeo = static_cast<int>(kGeoToInv.size());
    int defaultGeo = nGeo > 2 ? nGeo - 3 : 0; // middle-ish group
    std::vector<int> kBhToGeo(5, defaultGeo);
    for (int bh = 0; bh < 5; ++bh) {
        if (bhToInv[bh] >= 0 && bhToInv[bh] < static_cast<int>(invToGeo.size()))
            kBhToGeo[bh] = invToGeo[bhToInv[bh]];
    }

    struct BoreholeInfo {
        double x, y;
        std::string filepath;
    };

    // Try relative paths from build/release
    std::vector<BoreholeInfo> boreholes = {
        {   0.0,   0.0, "../../Lithosplitter_out/BH-01_Centre_classified.csv"},
        { 300.0,   0.0, "../../Lithosplitter_out/BH-02_East_classified.csv"},
        {-115.0,   0.0, "../../Lithosplitter_out/BH-03_South_classified.csv"},
        {-500.0,   0.0, "../../Lithosplitter_out/BH-04_FarWest_classified.csv"},
        {   0.0, 200.0, "../../Lithosplitter_out/BH-05_North_classified.csv"},
    };

    for (const auto& bh : boreholes) {
        std::ifstream file(bh.filepath);
        if (!file.is_open()) {
            // BH-04 and BH-05 have no classified data — use regional Felsic Granite
            if (bh.x == -500.0 || bh.y == 200.0) {
                Constraint c;
                c.position = Vector3d(bh.x, bh.y, 0.0);
                c.z_top = 0.0;
                c.z_bottom = 180.0;
                c.litho_group_id = kBhToGeo[1]; // BH c1 → geo group 2 (inv c1, Felsic Granite)
                constraints.push_back(c);
            } else {
                std::cerr << "generateClusterIdConstraints: Could not open "
                          << bh.filepath << std::endl;
            }
            continue;
        }

        std::string line;
        int lineNum = 0;
        struct Interval { double from_m, to_m; int cluster_id; };
        std::vector<Interval> intervals;

        while (std::getline(file, line)) {
            ++lineNum;
            if (lineNum == 1 || line.empty() || line[0] == '#') continue;

            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, ',')) {
                tokens.push_back(token);
            }

            // Columns: 0=hole_id, 1=from_m, 2=to_m, ..., 7=cluster_id
            if (tokens.size() < 8) continue;

            try {
                double from_m = std::stod(tokens[1]);
                double to_m   = std::stod(tokens[2]);
                int cluster_id = std::stoi(tokens[7]);
                intervals.push_back({from_m, to_m, cluster_id});
            } catch (...) {
                continue;
            }
        }

        if (intervals.empty()) continue;

        // Merge consecutive intervals with the same cluster_id
        Interval merged = intervals[0];
        for (size_t i = 1; i < intervals.size(); ++i) {
            if (intervals[i].cluster_id == merged.cluster_id) {
                merged.to_m = intervals[i].to_m;
            } else {
                Constraint c;
                c.position = Vector3d(bh.x, bh.y, 0.0);
                c.z_top = merged.from_m;
                c.z_bottom = merged.to_m;
                c.litho_group_id = kBhToGeo[merged.cluster_id];
                constraints.push_back(c);
                merged = intervals[i];
            }
        }
        // Last interval
        {
            Constraint c;
            c.position = Vector3d(bh.x, bh.y, 0.0);
            c.z_top = merged.from_m;
            c.z_bottom = merged.to_m;
            c.litho_group_id = kBhToGeo[merged.cluster_id];
            constraints.push_back(c);
        }
    }

    std::cout << "Cluster ID constraints: " << constraints.size()
              << " intervals from " << boreholes.size() << " boreholes" << std::endl;
    return constraints;
}

// =======================================================================
// Cluster ID true model — groups in cluster_id ascending order,
// bowl+ellipsoid geometry (Voisey's Bay primitives).
// =======================================================================
std::shared_ptr<LithologyModel> generateClusterIdTrueModel(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<int>& ordering)
{
    std::vector<int> kGeoToInv = ordering;
    if (kGeoToInv.empty()) {
        std::set<int> uniqueIds;
        for (const auto& c : clusters) uniqueIds.insert(c.cluster_id);
        kGeoToInv.assign(uniqueIds.begin(), uniqueIds.end());
    }

    auto findCluster = [&](int cid) -> const ClusterProperties* {
        for (const auto& c : clusters) {
            if (c.cluster_id == cid) return &c;
        }
        return nullptr;
    };

    auto model = std::make_shared<LithologyModel>();

    // Use CSV remanence values where available; fall back to sensible defaults.
    for (int gi = 0; gi < static_cast<int>(kGeoToInv.size()); ++gi) {
        int invCid = kGeoToInv[gi];
        const auto* cp = findCluster(invCid);

        double rho = cp ? cp->density_median : 2.67;
        double chi = cp ? cp->susceptibility_median : 0.0;
        std::string name = "Cluster_ID_" + std::to_string(invCid);

        LithoGroup g(gi, name, rho, chi);

        if (cp) {
            g.remanence_magnitude    = cp->remanence_magnitude;
            g.remanence_inclination  = cp->remanence_inclination;
            g.remanence_declination  = cp->remanence_declination;
        }

        model->addGroup(g);
    }

    model->setBottomDepth(-5000.0);

    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    // Build surface names from the actual cluster ordering
    std::string s0name = "contact_cluster_id_" + std::to_string(kGeoToInv[0])
                       + "_cluster_id_" + std::to_string(kGeoToInv[1]);
    std::string s1name = "contact_cluster_id_" + std::to_string(kGeoToInv[1])
                       + "_cluster_id_" + std::to_string(kGeoToInv[2]);
    std::string s2name = "contact_cluster_id_" + std::to_string(kGeoToInv[2])
                       + "_cluster_id_" + std::to_string(kGeoToInv[3]);
    std::string s3name = "contact_cluster_id_" + std::to_string(kGeoToInv[3])
                       + "_cluster_id_" + std::to_string(kGeoToInv[4]);

    // ---- Surface 0: kGeoToInv[0] / kGeoToInv[1] ----
    auto surf0 = makeTriangulatedGrid(xs, ys,
        s0name, -10000.0, 100.0,
        [](double, double) { return -30.0; });

    // ---- Surface 1: kGeoToInv[1] / kGeoToInv[2] ----
    auto surf1 = makeTriangulatedGrid(xs, ys,
        s1name, -10000.0, 100.0,
        [](double x, double y) -> double {
            return bowlZ(x, y, S0_CX, S0_CY, S0_INNER_DEPTH,
                        S0_A, S0_B, S0_REGIONAL);
        });

    // ---- Surface 2: kGeoToInv[2] / kGeoToInv[3] ----
    auto surf2 = makeTriangulatedGrid(xs, ys,
        s2name, -10000.0, 100.0,
        [](double x, double y) -> double {
            if (insideEllipsoid(x, y, S1_CX, S1_CY, S1_A, S1_B)) {
                return ellipsoidTopZ(x, y, S1_CX, S1_CY, S1_CZ,
                                     S1_A, S1_B, S1_C);
            }
            return S1_REGIONAL;
        });

    // ---- Surface 3: kGeoToInv[3] / kGeoToInv[4] ----
    auto surf3 = makeTriangulatedGrid(xs, ys,
        s3name, -10000.0, 100.0,
        [](double x, double y) -> double {
            if (insideEllipsoid(x, y, S1_CX, S1_CY, S1_A, S1_B)) {
                return ellipsoidTopZ(x, y, S1_CX, S1_CY, S1_CZ,
                                     S1_A, S1_B, S1_C) - 20.0;
            }
            return S1_REGIONAL;
        });

    auto flatTop = buildFlatMesh(*surf0, 0.0, "flat_top");
    auto flatBottom = buildFlatMesh(*surf0, model->bottomDepth(), "flat_bottom");
    model->setGroupMesh(0, buildClosedMesh(*flatTop, *surf0));
    model->setGroupMesh(1, buildClosedMesh(*surf0, *surf1));
    model->setGroupMesh(2, buildClosedMesh(*surf1, *surf2));
    model->setGroupMesh(3, buildClosedMesh(*surf2, *surf3));
    model->setGroupMesh(4, buildClosedMesh(*surf3, *flatBottom));

    // Diagnostic
    std::cout << "Cluster ID true model: " << model->groupCount() << " groups, "
              << model->groupMeshCount() << " group meshes" << std::endl;
    for (int g = 0; g < model->groupCount(); ++g) {
        const auto& grp = model->group(g);
        std::cout << "  [" << grp.id << "] " << grp.name
                  << ": rho=" << std::fixed << std::setprecision(3) << grp.density
                  << "  chi=" << std::scientific << std::setprecision(4) << grp.susceptibility
                  << "  M_rem=" << std::fixed << std::setprecision(1)
                  << grp.remanence_magnitude << " A/m"
                  << "  I_rem=" << grp.remanence_inclination
                  << "  D_rem=" << grp.remanence_declination << std::endl;
    }
    std::cout << "Group classification at borehole positions:" << std::endl;
    struct BhPos { const char* name; double x, y; };
    for (auto& bh : {BhPos{"BH-01_Centre", 0.0, 0.0},
                     BhPos{"BH-02_East", 300.0, 0.0},
                     BhPos{"BH-03_South", -115.0, 0.0},
                     BhPos{"BH-04_FarWest", -500.0, 0.0},
                     BhPos{"BH-05_North", 0.0, 200.0}}) {
        int gSurf = model->classifyPoint(Vector3d(bh.x, bh.y, -5.0));
        int gMid  = model->classifyPoint(Vector3d(bh.x, bh.y, -80.0));
        int gDeep = model->classifyPoint(Vector3d(bh.x, bh.y, -170.0));
        int gBot  = model->classifyPoint(Vector3d(bh.x, bh.y, -200.0));
        std::cout << "  " << bh.name << " (" << bh.x << "," << bh.y
                  << "): z=-5m g=" << gSurf
                  << "  z=-80m g=" << gMid
                  << "  z=-170m g=" << gDeep
                  << "  z=-200m g=" << gBot << std::endl;
    }

    return model;
}

// =======================================================================
// Derive stratigraphic ordering from borehole classified intervals
// =======================================================================
std::vector<int> deriveOrderingFromBoreholes(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<std::string>& classifiedCsvPaths)
{
    int nClusters = static_cast<int>(clusters.size());

    // Find the shallowest depth each cluster_id first appears in any borehole
    std::vector<double> minDepth(nClusters, 1e18);
    std::vector<int> allCids;
    for (const auto& c : clusters) {
        allCids.push_back(c.cluster_id);
    }

    int filesRead = 0;
    for (const auto& fpath : classifiedCsvPaths) {
        std::ifstream file(fpath);
        if (!file.is_open()) continue;
        ++filesRead;

        std::string line;
        int lineno = 0;
        while (std::getline(file, line)) {
            ++lineno;
            if (lineno == 1 || line.empty() || line[0] == '#') continue;

            std::stringstream ss(line);
            std::string tok;
            int col = 0, cid = -1;
            double from_m = -1;
            while (std::getline(ss, tok, ',')) {
                if (col == 1) from_m = std::stod(tok);
                else if (col == 7) {
                    try { cid = std::stoi(tok); } catch (...) { cid = -1; }
                }
                ++col;
            }
            if (cid >= 0 && from_m >= 0) {
                // Map BH cluster_id to its position in allCids
                for (int i = 0; i < nClusters; ++i) {
                    if (allCids[i] == cid && from_m < minDepth[i]) {
                        minDepth[i] = from_m;
                    }
                }
            }
        }
    }

    // Sort cluster_ids by min depth (shallowest = top = group 0)
    std::vector<std::pair<int, double>> cidDepth;
    for (int i = 0; i < nClusters; ++i) {
        cidDepth.emplace_back(allCids[i], minDepth[i]);
    }
    std::sort(cidDepth.begin(), cidDepth.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });

    std::vector<int> ordering;
    for (const auto& [cid, depth] : cidDepth) {
        ordering.push_back(cid);
    }

    if (filesRead > 0) {
        std::cout << "deriveOrderingFromBoreholes: " << filesRead << " boreholes, "
                  << ordering.size() << " clusters (top→bottom):";
    } else {
        std::cout << "deriveOrderingFromBoreholes: no borehole data, "
                  << ordering.size() << " clusters (cluster_id ascending):";
    }
    for (int cid : ordering) std::cout << " c" << cid;
    std::cout << std::endl;

    return ordering;
}

// =======================================================================
// Layered true model — flat horizontal surfaces
// =======================================================================
std::shared_ptr<LithologyModel> generateLayeredTrueModel(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<double>& surfaceDepths,
    const std::vector<int>& ordering)
{
    std::vector<int> kGeoToInv = ordering;
    if (kGeoToInv.empty()) {
        std::set<int> uniqueIds;
        for (const auto& c : clusters) uniqueIds.insert(c.cluster_id);
        kGeoToInv.assign(uniqueIds.begin(), uniqueIds.end());
    }
    int nGroups = static_cast<int>(kGeoToInv.size());
    int nSurfs = nGroups - 1;

    if (static_cast<int>(surfaceDepths.size()) != nSurfs) {
        std::cerr << "generateLayeredTrueModel: expected " << nSurfs
                  << " surface depths for " << nGroups << " groups, got "
                  << surfaceDepths.size() << std::endl;
        return nullptr;
    }

    auto findCluster = [&](int cid) -> const ClusterProperties* {
        for (const auto& c : clusters) {
            if (c.cluster_id == cid) return &c;
        }
        return nullptr;
    };

    auto model = std::make_shared<LithologyModel>();
    model->setBottomDepth(-5000.0);

    for (int gi = 0; gi < nGroups; ++gi) {
        int invCid = kGeoToInv[gi];
        const auto* cp = findCluster(invCid);

        double rho = cp ? cp->density_median : 2.67;
        double chi = cp ? cp->susceptibility_median : 0.0;
        std::string name = "Cluster_ID_" + std::to_string(invCid);

        LithoGroup g(gi, name, rho, chi);
        if (cp) {
            g.remanence_magnitude    = cp->remanence_magnitude;
            g.remanence_inclination  = cp->remanence_inclination;
            g.remanence_declination  = cp->remanence_declination;
        }
        model->addGroup(g);
    }

    // Build grid coordinates
    constexpr int INTERIOR_DIM = 21;
    constexpr double GRID_HALF = 1000.0;
    constexpr double CELL_SIZE = 100.0;
    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    // Build flat horizontal contact surfaces
    std::vector<std::shared_ptr<SurfaceMesh>> surfaces;
    for (int si = 0; si < nSurfs; ++si) {
        std::string sname = "contact_c" + std::to_string(kGeoToInv[si])
                          + "_c" + std::to_string(kGeoToInv[si + 1]);
        auto surf = makeTriangulatedGrid(xs, ys, sname, -10000.0, 100.0,
            [z = surfaceDepths[si]](double, double) { return z; });
        surfaces.push_back(surf);
    }

    auto flatTop = buildFlatMesh(*surfaces[0], 0.0, "flat_top");
    auto flatBottom = buildFlatMesh(*surfaces[0], model->bottomDepth(), "flat_bottom");
    for (int g = 0; g < nGroups; ++g) {
        const SurfaceMesh* top = (g == 0) ? flatTop.get() : surfaces[g - 1].get();
        const SurfaceMesh* bot = (g == nGroups - 1) ? flatBottom.get() : surfaces[g].get();
        model->setGroupMesh(g, buildClosedMesh(*top, *bot));
    }

    // Diagnostic
    std::cout << "Layered true model: " << nGroups << " groups, "
              << model->groupMeshCount() << " group meshes" << std::endl;
    for (int g = 0; g < nGroups; ++g) {
        const auto& grp = model->group(g);
        std::cout << "  [" << grp.id << "] " << grp.name
                  << ": rho=" << std::fixed << std::setprecision(3) << grp.density
                  << "  chi=" << std::scientific << std::setprecision(4) << grp.susceptibility
                  << "  M_rem=" << std::fixed << std::setprecision(1)
                  << grp.remanence_magnitude << " A/m"
                  << "  I_rem=" << grp.remanence_inclination
                  << "  D_rem=" << grp.remanence_declination << std::endl;
    }
    for (int si = 0; si < nSurfs; ++si) {
        std::cout << "  Surface " << si << ": z=" << surfaceDepths[si] << "m"
                  << std::endl;
    }

    return model;
}

// =======================================================================
// Dipping layered true model — planar surfaces with consistent dip
// =======================================================================
std::shared_ptr<LithologyModel> generateDippingLayeredTrueModel(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<double>& surfaceDepths,
    double dipAngleDeg,
    double dipDirectionDeg,
    const std::vector<int>& ordering)
{
    std::vector<int> kGeoToInv = ordering;
    if (kGeoToInv.empty()) {
        std::set<int> uniqueIds;
        for (const auto& c : clusters) uniqueIds.insert(c.cluster_id);
        kGeoToInv.assign(uniqueIds.begin(), uniqueIds.end());
    }
    int nGroups = static_cast<int>(kGeoToInv.size());
    int nSurfs = nGroups - 1;

    if (static_cast<int>(surfaceDepths.size()) != nSurfs) {
        std::cerr << "generateDippingLayeredTrueModel: expected " << nSurfs
                  << " surface depths for " << nGroups << " groups, got "
                  << surfaceDepths.size() << std::endl;
        return nullptr;
    }

    auto findCluster = [&](int cid) -> const ClusterProperties* {
        for (const auto& c : clusters) {
            if (c.cluster_id == cid) return &c;
        }
        return nullptr;
    };

    auto model = std::make_shared<LithologyModel>();
    model->setBottomDepth(-5000.0);

    for (int gi = 0; gi < nGroups; ++gi) {
        int invCid = kGeoToInv[gi];
        const auto* cp = findCluster(invCid);

        double rho = cp ? cp->density_median : 2.67;
        double chi = cp ? cp->susceptibility_median : 0.0;
        std::string name = "Cluster_ID_" + std::to_string(invCid);

        LithoGroup g(gi, name, rho, chi);
        if (cp) {
            g.remanence_magnitude    = cp->remanence_magnitude;
            g.remanence_inclination  = cp->remanence_inclination;
            g.remanence_declination  = cp->remanence_declination;
        }
        model->addGroup(g);
    }

    // Precompute dip vector components
    const double pi = std::acos(-1.0);
    double dipRad = dipAngleDeg * pi / 180.0;
    double dirRad = dipDirectionDeg * pi / 180.0;
    double tanDip = std::tan(dipRad);
    double dxDip = std::sin(dirRad) * tanDip;  // east component
    double dyDip = std::cos(dirRad) * tanDip;  // north component

    // Build grid coordinates
    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);

    // Build dipping planar contact surfaces
    std::vector<std::shared_ptr<SurfaceMesh>> surfaces;
    for (int si = 0; si < nSurfs; ++si) {
        double z0 = surfaceDepths[si];
        std::string sname = "contact_c" + std::to_string(kGeoToInv[si])
                          + "_c" + std::to_string(kGeoToInv[si + 1]);
        auto surf = makeTriangulatedGrid(xs, ys, sname, -10000.0, 100.0,
            [z0, dxDip, dyDip](double x, double y) {
                return z0 + x * dxDip + y * dyDip;
            });
        surfaces.push_back(surf);
    }

    auto flatTop = buildFlatMesh(*surfaces[0], 0.0, "flat_top");
    auto flatBottom = buildFlatMesh(*surfaces[0], model->bottomDepth(), "flat_bottom");
    for (int g = 0; g < nGroups; ++g) {
        const SurfaceMesh* top = (g == 0) ? flatTop.get() : surfaces[g - 1].get();
        const SurfaceMesh* bot = (g == nGroups - 1) ? flatBottom.get() : surfaces[g].get();
        model->setGroupMesh(g, buildClosedMesh(*top, *bot));
    }

    // Diagnostic
    std::cout << "Dipping layered true model: " << nGroups << " groups, "
              << model->groupMeshCount() << " group meshes"
              << "  dip=" << dipAngleDeg << "°  direction=" << dipDirectionDeg << "°"
              << std::endl;
    for (int g = 0; g < nGroups; ++g) {
        const auto& grp = model->group(g);
        std::cout << "  [" << grp.id << "] " << grp.name
                  << ": rho=" << std::fixed << std::setprecision(3) << grp.density
                  << "  chi=" << std::scientific << std::setprecision(4) << grp.susceptibility
                  << "  M_rem=" << std::fixed << std::setprecision(1)
                  << grp.remanence_magnitude << " A/m"
                  << "  I_rem=" << grp.remanence_inclination
                  << "  D_rem=" << grp.remanence_declination << std::endl;
    }
    for (int si = 0; si < nSurfs; ++si) {
        std::cout << "  Surface " << si << ": z@(0,0)=" << surfaceDepths[si] << "m"
                  << std::endl;
    }

    return model;
}

} // namespace litho_invert

