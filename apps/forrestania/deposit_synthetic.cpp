#include "deposit_synthetic.h"
#include "generate_synthetic.h"
#include <litho_invert/forward/gravity_forward.h>
#include <litho_invert/forward/magnetic_forward.h>
#include <litho_invert/surface/surface_mesh.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <set>
#include <filesystem>

namespace litho_invert {

namespace {

constexpr double GRID_HALF = 1000.0;
constexpr double CELL_SIZE = 100.0;
constexpr int N_VERTS = 21;
constexpr int N_CELLS = 20;

constexpr double OBS_HALF = 950.0;
constexpr double OBS_CELL = 100.0;
constexpr int OBS_N = 20;

// ---- CSV helpers ----
std::vector<std::string> splitCSV(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool inQuotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '"') { inQuotes = !inQuotes; continue; }
        if (c == ',' && !inQuotes) {
            fields.push_back(field);
            field.clear();
        } else {
            field += c;
        }
    }
    fields.push_back(field);
    return fields;
}

std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, b - a + 1);
}

double parseDouble(const std::string& s) {
    try { return std::stod(trim(s)); }
    catch (...) { return 0.0; }
}

int parseInt(const std::string& s) {
    try { return std::stoi(trim(s)); }
    catch (...) { return -1; }
}

// ---- Ellipsoid helpers ----
bool insideEllipsoid(double x, double y,
                     double cx, double cy, double a, double b) {
    double dx = x - cx, dy = y - cy;
    return (dx*dx)/(a*a) + (dy*dy)/(b*b) <= 1.0;
}

double ellipsoidTopZ(double x, double y,
                     double cx, double cy, double cz,
                     double a, double b, double c) {
    double dx = x - cx, dy = y - cy;
    double val = (dx*dx)/(a*a) + (dy*dy)/(b*b);
    if (val > 1.0) return cz; // outside ellipsoid: regional depth
    return cz + c * std::sqrt(1.0 - val);
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
    std::function<double(double,double)> zFunc)
{
    auto surf = std::make_shared<SurfaceMesh>();
    surf->setName(name);
    surf->setBounds(minZ, maxZ);

    for (int iy = 0; iy < N_VERTS; ++iy) {
        for (int ix = 0; ix < N_VERTS; ++ix) {
            double x = xs[ix], y = ys[iy];
            double z = zFunc(x, y);
            surf->addVertex(x, y, z, VertexFreedom::Z_ONLY);
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

// =========================================================================
// Load deposit data
// =========================================================================
DepositData loadDepositData(const std::string& depositPath) {
    DepositData data;

    // Assign default borehole positions (triangular pattern)
    data.holePositions["BH-1"] = {0.0, 0.0};
    data.holePositions["BH-2"] = {300.0, 0.0};
    data.holePositions["BH-3"] = {-150.0, 200.0};

    // Read synthetic_dataset.csv
    std::string datasetPath = depositPath + "/Dataset/synthetic_dataset.csv";
    std::ifstream dsFile(datasetPath);
    if (!dsFile.is_open()) {
        std::cerr << "WARNING: Cannot open " << datasetPath << std::endl;
        return data;
    }

    std::string line;
    std::getline(dsFile, line); // skip header

    // Find column indices
    auto header = splitCSV(line);
    int colHole = -1, colFrom = -1, colTo = -1, colLith = -1;
    int colDens = -1, colSusc = -1;
    for (size_t i = 0; i < header.size(); ++i) {
        std::string h = trim(header[i]);
        if (h == "hole_id") colHole = static_cast<int>(i);
        else if (h == "from_m") colFrom = static_cast<int>(i);
        else if (h == "to_m") colTo = static_cast<int>(i);
        else if (h == "ground_truth_lithology") colLith = static_cast<int>(i);
        else if (h == "density_gcc") colDens = static_cast<int>(i);
        else if (h == "susceptibility_si") colSusc = static_cast<int>(i);
    }

    std::set<std::string> holeSet;
    while (std::getline(dsFile, line)) {
        auto fields = splitCSV(line);
        if (static_cast<int>(fields.size()) <= std::max({colHole, colFrom, colTo, colDens}))
            continue;

        DepositData::Interval iv;
        iv.hole_id = trim(fields[colHole]);
        iv.from_m = parseDouble(fields[colFrom]);
        iv.to_m = parseDouble(fields[colTo]);
        iv.lithology = (colLith >= 0) ? trim(fields[colLith]) : "";
        iv.density = parseDouble(fields[colDens]);
        iv.susceptibility = (colSusc >= 0) ? parseDouble(fields[colSusc]) : 0.0;
        data.intervals.push_back(iv);
        holeSet.insert(iv.hole_id);
    }

    data.holeIds.assign(holeSet.begin(), holeSet.end());

    // Read inversion_start.csv
    std::string clusterPath = depositPath + "/Results/inversion_start.csv";
    data.clusters = loadClusterProperties(clusterPath);

    std::cout << "loadDepositData: " << data.intervals.size() << " intervals from "
              << data.holeIds.size() << " boreholes, "
              << data.clusters.size() << " clusters" << std::endl;

    return data;
}

// =========================================================================
// Derive cluster-to-group mapping
// =========================================================================
std::vector<ClusterGroupMapping> deriveClusterMapping(const DepositData& data, int nGroups) {
    std::vector<ClusterGroupMapping> mapping;

    if (data.clusters.empty()) return mapping;

    // Sort clusters by density
    std::vector<int> idx(data.clusters.size());
    for (size_t i = 0; i < idx.size(); ++i) idx[i] = static_cast<int>(i);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        return data.clusters[a].density_median < data.clusters[b].density_median;
    });

    // Quantile-based group assignment: ensures every group gets at least
    // one cluster, avoiding gaps that cause "no clusters map to group N".
    int n = static_cast<int>(data.clusters.size());
    for (size_t pos = 0; pos < data.clusters.size(); ++pos) {
        int groupId;
        if (nGroups == 2) {
            groupId = (pos >= n / 2) ? 1 : 0;
        } else {
            // 3 groups: split sorted clusters into thirds (by index, not density)
            if (n >= 3) {
                double frac = static_cast<double>(pos) / (n - 1);
                if (frac >= 0.67)      groupId = 2;
                else if (frac >= 0.33) groupId = 1;
                else                   groupId = 0;
            } else if (n == 2) {
                groupId = (pos == 0) ? 0 : 2;
            } else {
                groupId = 1; // single cluster → middle group
            }
        }
        mapping.push_back({data.clusters[idx[pos]].cluster_id, groupId});
    }

    return mapping;
}

// =========================================================================
// Generate deposit true model
// =========================================================================
std::shared_ptr<LithologyModel> generateDepositTrueModel(
    const DepositData& data,
    const std::vector<ClusterGroupMapping>& clusterMap)
{
    auto model = std::make_shared<LithologyModel>();

    // Determine number of groups and their properties
    int nGroups = 0;
    for (const auto& m : clusterMap)
        nGroups = std::max(nGroups, m.litho_group_id + 1);

    // Accumulate properties per group (weighted by sample count)
    std::vector<double> sumDens(nGroups, 0.0), sumSusc(nGroups, 0.0);
    std::vector<int> sumN(nGroups, 0);
    std::vector<std::string> groupNames(nGroups);
    std::vector<std::string> allNames;

    for (const auto& cluster : data.clusters) {
        int gid = -1;
        for (const auto& m : clusterMap) {
            if (m.cluster_id == cluster.cluster_id) { gid = m.litho_group_id; break; }
        }
        if (gid < 0) continue;
        double w = static_cast<double>(cluster.sample_count);
        sumDens[gid] += cluster.density_median * w;
        sumSusc[gid] += cluster.susceptibility_median * w;
        sumN[gid] += cluster.sample_count;
        allNames.push_back(cluster.working_name);
    }

    const char* defaultNames[] = {"Background", "HostRock", "OreBody"};
    for (int g = 0; g < nGroups; ++g) {
        if (sumN[g] > 0) {
            sumDens[g] /= static_cast<double>(sumN[g]);
            sumSusc[g] /= static_cast<double>(sumN[g]);
        }
        groupNames[g] = (g < 3) ? defaultNames[g] : ("Group" + std::to_string(g));
    }

    // Add groups
    for (int g = 0; g < nGroups; ++g) {
        LithoGroup grp(g, groupNames[g], sumDens[g], sumSusc[g]);
        model->addGroup(grp);
    }
    model->setBottomDepth(-5000.0);

    // ---- Auto-generate surface geometry ----
    //
    // Map borehole intervals to groups and find contact depths.
    // Contact between group A (above) and group B (below) → surface i depth.
    // Surface i separates group i from group i+1.
    //
    // For each borehole, we find the shallowest depth where group transitions
    // from group i to i+1.  We also find where group i+1 appears.
    //
    // Strategy:
    //   1. For each borehole, walk intervals top-to-bottom, assign groups
    //   2. Record first occurrence depth of each group
    //   3. Surface i = contact between group i and i+1 = first occurrence
    //      depth of group i+1
    //
    // The regional depth for surface i is the median contact depth across
    // all boreholes.  Boreholes where group i+1 is absent get the regional depth.
    //
    // An ellipsoid perturbation is placed at any borehole where group i+1
    // (ore/host) occurs significantly shallower than the regional.

    // Assign groups to each interval
    struct BhInterval {
        std::string holeId;
        double from, to;
        double density;
        int groupId;
    };
    std::vector<BhInterval> bhIntervals;
    for (const auto& iv : data.intervals) {
        // Find nearest cluster by density
        int bestCluster = -1;
        double bestDiff = 1e9;
        for (const auto& cl : data.clusters) {
            double diff = std::abs(iv.density - cl.density_median);
            if (diff < bestDiff) { bestDiff = diff; bestCluster = cl.cluster_id; }
        }
        int gid = -1;
        for (const auto& m : clusterMap) {
            if (m.cluster_id == bestCluster) { gid = m.litho_group_id; break; }
        }
        if (gid < 0) gid = 0;
        bhIntervals.push_back({iv.hole_id, iv.from_m, iv.to_m, iv.density, gid});
    }

    // For each borehole, find first occurrence depth of each group
    // (depth positive-down in borehole convention)
    std::map<std::string, std::vector<double>> groupFirstDepth; // per hole, per group
    for (const auto& holeId : data.holeIds) {
        groupFirstDepth[holeId].resize(nGroups, 9999.0);
        for (const auto& bi : bhIntervals) {
            if (bi.holeId != holeId) continue;
            int g = bi.groupId;
            if (bi.from < groupFirstDepth[holeId][g])
                groupFirstDepth[holeId][g] = bi.from;
        }
    }

    // Compute regional depth for each surface as median of contacts
    std::vector<double> regionalDepth(nGroups - 1, -500.0);
    for (int si = 0; si < nGroups - 1; ++si) {
        std::vector<double> depths;
        for (const auto& holeId : data.holeIds) {
            double d = groupFirstDepth[holeId][si + 1];
            if (d < 9000.0) depths.push_back(-d); // convert to z-up
        }
        if (!depths.empty()) {
            std::sort(depths.begin(), depths.end());
            regionalDepth[si] = depths[depths.size() / 2];
        }
    }

    // Build contact surfaces
    std::vector<double> xs, ys;
    makeGridCoords(xs, ys);
    std::vector<std::shared_ptr<SurfaceMesh>> surfaces;

    for (int si = 0; si < nGroups - 1; ++si) {
        // Find "anomaly" boreholes: where contact is significantly above regional
        // (group i+1 appears shallower than regional → ellipsoid lifts surface up)
        struct AnomalyInfo {
            double cx, cy, depth; // z-up depth (negative)
        };
        std::vector<AnomalyInfo> anomalies;
        double regional = regionalDepth[si];

        for (const auto& holeId : data.holeIds) {
            double d = groupFirstDepth[holeId][si + 1];
            if (d < 9000.0 && -d > regional + 50.0) { // at least 50m shallower
                auto pos = data.holePositions.find(holeId);
                if (pos != data.holePositions.end()) {
                    anomalies.push_back({pos->second.first, pos->second.second, -d});
                }
            }
        }

        // If no strong anomaly, find ANY borehole with group i+1 present
        if (anomalies.empty()) {
            for (const auto& holeId : data.holeIds) {
                double d = groupFirstDepth[holeId][si + 1];
                if (d < 9000.0 && -d > regional + 10.0) {
                    auto pos = data.holePositions.find(holeId);
                    if (pos != data.holePositions.end()) {
                        anomalies.push_back({pos->second.first, pos->second.second, -d});
                    }
                }
            }
        }

        // Build ellipsoid parameters
        struct EllipsoidParams {
            double cx, cy, cz, a, b, c;
        };
        std::vector<EllipsoidParams> ellipsoids;

        for (const auto& anom : anomalies) {
            // Ellipsoid size: proportional to thickness of anomalous unit
            // (how much above regional)
            double thickness = regional - anom.depth; // positive = above regional
            double a = std::max(100.0, thickness * 3.0);
            double b = std::max(50.0, thickness * 1.5);
            double c = std::max(20.0, thickness * 0.8);
            // Center the ellipsoid so its top touches the anomaly depth
            double cz = anom.depth - c;
            ellipsoids.push_back({anom.cx, anom.cy, cz, a, b, c});
        }

        if (ellipsoids.empty()) {
            // Create a small default perturbation at first borehole
            double cx = 0, cy = 0;
            if (!data.holeIds.empty()) {
                auto it = data.holePositions.find(data.holeIds[0]);
                if (it != data.holePositions.end()) {
                    cx = it->second.first;
                    cy = it->second.second;
                }
            }
            ellipsoids.push_back({cx, cy, regional - 50.0, 200.0, 100.0, 50.0});
        }

        // Build surface z-function
        auto zFunc = [regional, ellipsoids](double x, double y) -> double {
            double bestZ = regional; // default to regional
            for (const auto& ep : ellipsoids) {
                if (insideEllipsoid(x, y, ep.cx, ep.cy, ep.a, ep.b)) {
                    double z = ellipsoidTopZ(x, y, ep.cx, ep.cy, ep.cz,
                                             ep.a, ep.b, ep.c);
                    bestZ = std::max(bestZ, z); // highest (shallowest) wins
                }
            }
            return bestZ;
        };

        std::string sname = "surface_" + std::to_string(si);
        auto surf = makeTriangulatedGrid(xs, ys, sname, -10000.0, 100.0, zFunc);
        surfaces.push_back(surf);
    }

    // Build closed group meshes from contact surfaces
    auto flatTop = buildFlatMesh(*surfaces[0], 0.0, "flat_top");
    auto flatBottom = buildFlatMesh(*surfaces[0], model->bottomDepth(), "flat_bottom");
    for (int g = 0; g < nGroups; ++g) {
        const SurfaceMesh* top = (g == 0) ? flatTop.get() : surfaces[g - 1].get();
        const SurfaceMesh* bot = (g == nGroups - 1) ? flatBottom.get() : surfaces[g].get();
        model->setGroupMesh(g, buildClosedMesh(*top, *bot));
    }

    // Diagnostic
    std::cout << "True model: " << model->groupCount() << " groups, "
              << model->groupMeshCount() << " group meshes" << std::endl;
    for (int g = 0; g < model->groupCount(); ++g) {
        const auto& grp = model->group(g);
        std::cout << "  [" << grp.id << "] " << grp.name
                  << ": rho=" << std::fixed << std::setprecision(2) << grp.density
                  << "  chi=" << std::scientific << std::setprecision(3) << grp.susceptibility
                  << std::endl;
    }

    return model;
}

// =========================================================================
// Observation grid
// =========================================================================
GravityData generateDepositObservationPoints(const DepositData& data) {
    GravityData pts;
    for (int iy = 0; iy < OBS_N; ++iy) {
        for (int ix = 0; ix < OBS_N; ++ix) {
            double x = -OBS_HALF + ix * OBS_CELL;
            double y = -OBS_HALF + iy * OBS_CELL;
            pts.push_back({Vector3d(x, y, 0.0), 0.0, 0.01});
        }
    }
    return pts;
}

// =========================================================================
// Constraints from borehole intervals
// =========================================================================
std::vector<Constraint> generateDepositConstraints(
    const DepositData& data,
    const std::vector<ClusterGroupMapping>& clusterMap,
    const std::vector<ClusterProperties>& clusters)
{
    std::vector<Constraint> constraints;

    // Map cluster_id → group_id
    std::map<int, int> c2g;
    for (const auto& m : clusterMap)
        c2g[m.cluster_id] = m.litho_group_id;

    // For each borehole, merge consecutive intervals with same group
    for (const auto& holeId : data.holeIds) {
        double x = 0, y = 0;
        auto it = data.holePositions.find(holeId);
        if (it != data.holePositions.end()) {
            x = it->second.first;
            y = it->second.second;
        }

        // Collect intervals for this hole sorted by depth
        struct HoleIv { double from, to; int gid; };
        std::vector<HoleIv> holeIvs;
        for (const auto& iv : data.intervals) {
            if (iv.hole_id != holeId) continue;
            // Find nearest cluster
            int bestCl = -1;
            double bestDiff = 1e9;
            for (const auto& cl : clusters) {
                double diff = std::abs(iv.density - cl.density_median);
                if (diff < bestDiff) { bestDiff = diff; bestCl = cl.cluster_id; }
            }
            int gid = (c2g.count(bestCl)) ? c2g[bestCl] : 0;
            holeIvs.push_back({iv.from_m, iv.to_m, gid});
        }

        // Sort by depth
        std::sort(holeIvs.begin(), holeIvs.end(),
                  [](const HoleIv& a, const HoleIv& b) { return a.from < b.from; });

        // Merge consecutive intervals with same group
        for (size_t i = 0; i < holeIvs.size(); ) {
            double from = holeIvs[i].from;
            double to = holeIvs[i].to;
            int gid = holeIvs[i].gid;
            size_t j = i + 1;
            while (j < holeIvs.size() && holeIvs[j].gid == gid
                   && std::abs(holeIvs[j].from - to) < 0.1) {
                to = holeIvs[j].to;
                ++j;
            }
            Constraint c;
            c.position = Vector3d(x, y, 0.0);
            c.z_top = from;
            c.z_bottom = to;
            c.litho_group_id = gid;
            constraints.push_back(c);
            i = j;
        }
    }

    std::cout << "Constraints: " << constraints.size() << " intervals from "
              << data.holeIds.size() << " boreholes" << std::endl;

    return constraints;
}

// =========================================================================
// Generate INI config
// =========================================================================
std::string generateDepositIni(const DepositData& data,
                                const std::string& clusterCsvPath,
                                const std::string& outputDir,
                                int stride,
                                int maxIterations,
                                double propertyMismatch)
{
    // property_inversion_interval: fire property inversion every N geometry
    // iterations.  Set to 10 so it triggers multiple times within a 50-iteration
    // run but doesn't dominate the geometry optimisation.
    int propInterval = std::min(10, maxIterations / 3);
    if (propInterval < 5) propInterval = 5;

    std::ostringstream oss;
    oss << "# Auto-generated test inversion config\n";
    oss << "# Deposit: " << outputDir << "\n";
    oss << "# property_mismatch = " << propertyMismatch << "\n";
    oss << "[inversion]\n";
    oss << "max_iterations = " << maxIterations << "\n";
    oss << "tolerance = 1e-5\n";
    oss << "lambda = 0.0\n";
    oss << "omega = 10.0\n";
    oss << "control_point_stride = " << stride << "\n";
    oss << "lbfgs_history = 10\n";
    oss << "enable_property_inversion = true\n";
    oss << "property_inversion_interval = " << propInterval << "\n";
    oss << "property_inversion_max_iter = 10\n";
    oss << "\n[gravity]\n";
    oss << "density_min = 1.5\n";
    oss << "density_max = 6.0\n";
    oss << "\n[magnetic]\n";
    oss << "magnetic_weight = 0.1\n";
    oss << "inclination = 75.0\n";
    oss << "declination = -20.0\n";
    oss << "field_nT = 55000.0\n";
    oss << "susceptibility_min = 0.0\n";
    oss << "susceptibility_max = 1.0\n";
    oss << "remanence_min = 0.0\n";
    oss << "remanence_max = 10.0\n";
    oss << "\n[topography]\n";
    oss << "mode = none\n";
    oss << "dem_file =\n";
    oss << "datum_elevation = 0.0\n";
    oss << "bouguer_density = 2.67\n";
    oss << "padding_rings = 0\n";
    oss << "padding_cell_size = 0.0\n";
    oss << "invert_halfspace_properties = false\n";
    oss << "\n[padding]\n";
    oss << "enable_padding_group = true\n";
    oss << "padding_density_initial = 2.68\n";
    oss << "padding_density_lower = 1.5\n";
    oss << "padding_density_upper = 6.0\n";
    oss << "padding_depth = -100000.0\n";
    oss << "\n[data]\n";
    oss << "cluster_csv = " << clusterCsvPath << "\n";
    oss << "\n[output]\n";
    oss << "output_dir = " << outputDir << "\n";

    return oss.str();
}

void perturbModelProperties(LithologyModel& model, double factor) {
    // Deterministic seed per group so results are reproducible.
    // Hash the group name to get a consistent seed.
    for (int g = 0; g < model.groupCount(); ++g) {
        auto& grp = model.group(g);
        // Simple hash of group name for reproducibility
        uint32_t hash = 0;
        for (char c : grp.name) hash = hash * 31 + static_cast<uint32_t>(c);
        // Convert hash to uniform [-1, 1]
        double u = static_cast<double>(hash % 10000) / 10000.0;
        double perturbation = factor * (2.0 * u - 1.0);

        double origRho = grp.density;
        double origChi = grp.susceptibility;
        grp.density = origRho * (1.0 + perturbation);
        grp.susceptibility = origChi * (1.0 + perturbation);

        std::cout << "  Property mismatch [" << grp.name << "]: "
                  << "rho " << std::fixed << std::setprecision(3) << origRho
                  << " -> " << grp.density
                  << " (" << std::showpos << perturbation * 100.0 << "%"
                  << std::noshowpos << "), "
                  << "chi " << std::scientific << std::setprecision(3) << origChi
                  << " -> " << grp.susceptibility
                  << std::endl;
    }
}

} // namespace litho_invert
