#include <litho_invert/inversion/runner.h>
#include <litho_invert/io/io.h>
#include <litho_invert/io/ini_config.h>
// dem_loader.h removed — topography/padding no longer supported in volumetric architecture
#include <litho_invert/forward/gravity_forward.h>
#include <litho_invert/forward/magnetic_forward.h>
#include <litho_invert/plot/svg_plot.h>
#include "generate_synthetic.h"
#include "deposit_synthetic.h"
#include "cluster_loader.h"
#include "mesh_components.h"
#include "model_setup.h"
#include <litho_invert/io/gravity_loader.h>
#include <litho_invert/io/magnetic_loader.h>
#include <litho_invert/io/surface_loader.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <memory>
#include <filesystem>
#include <cmath>
#include <set>
#include <map>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <algorithm>

using namespace litho_invert;

int main(int argc, char* argv[]) {
    std::string outputDir = ".";
    std::string csvPath = "Lithosplitter_out/inversion_start.csv";
    std::string iniPath;
    std::string depositPath; // --deposit <path> mode

    // CLI overrides (applied after INI config is loaded)
    int cliMaxIter = 0;
    int cliStride = 0;
    double cliPropertyMismatch = -1.0; // negative = use default
    std::string cliTrueModelMode;
    std::string cliTrueModelDepths;
    double cliTrueModelDip = -1.0;
    double cliTrueModelDipDir = -1.0;

    // Mode detection: --deposit <path> | *.ini | output_dir [cluster_csv]
    bool useIni = false;
    bool useDeposit = false;
    if (argc > 1) {
        std::string arg1 = argv[1];
        if (arg1 == "--deposit" && argc > 2) {
            depositPath = argv[2];
            useDeposit = true;
        } else if (arg1.size() >= 4 && arg1.compare(arg1.size() - 4, 4, ".ini") == 0) {
            iniPath = arg1;
            useIni = true;
        } else {
            outputDir = argv[1];
        }
    }
    if (argc > 2 && !useDeposit && !useIni) csvPath = argv[2];

    // Parse optional --key=value (or --key value) overrides
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto tryVal = [&](const std::string& key) -> std::string {
            // --key=value format
            if (arg.rfind("--" + key + "=", 0) == 0)
                return arg.substr(key.size() + 3);
            // --key value format
            if (arg == "--" + key && i + 1 < argc)
                return argv[++i];
            return "";
        };
        std::string val;
        if (!(val = tryVal("max-iter")).empty())
            cliMaxIter = std::stoi(val);
        else if (!(val = tryVal("stride")).empty())
            cliStride = std::stoi(val);
        else if (!(val = tryVal("property-mismatch")).empty())
            cliPropertyMismatch = std::stod(val);
        else if (!(val = tryVal("true-model")).empty())
            cliTrueModelMode = val;
        else if (!(val = tryVal("depths")).empty())
            cliTrueModelDepths = val;
        else if (!(val = tryVal("dip")).empty())
            cliTrueModelDip = std::stod(val);
        else if (!(val = tryVal("dip-dir")).empty())
            cliTrueModelDipDir = std::stod(val);
    }

    std::filesystem::create_directories(outputDir);

    // ---- Load INI config if provided ----
    InversionConfig config;
    int paddingRings = 0;
    std::string groupColumn = "litho_group_id";
    std::vector<std::vector<int>> groupClusters; // per-group list of composing cluster IDs
    std::vector<ClusterGroupMapping> clusterMap;
    std::vector<ClusterProperties> clusters; // loaded from cluster CSV

    // Borehole classified CSV paths for ordering derivation
    std::vector<std::string> boreholeCsvPaths;

    // Real data mode (populated from INI, defaults for non-INI paths)
    bool useRealData = false;
    std::string obsGravCsv;
    std::string obsMagCsv;
    double localOriginX = 0.0;
    double localOriginY = 0.0;
    std::vector<std::string> surfaceCsvPaths;
    std::vector<std::string> contactMeshPaths;
    std::vector<std::string> groupMeshPaths;     // direct closed-volume meshes (new volumetric)
    std::string labelGridPath;                   // 3D label grid binary (new volumetric)
    double zDatum = 0.0;

    // Magnetic field constants (INI or defaults)
    double magInc = 75.0, magDec = -20.0, magField = 55000.0;

    // True model configuration (populated from INI or CLI)
    std::string trueModelMode = "voiseys";
    std::vector<double> trueModelDepths;
    double trueModelDip = 62.0;
    double trueModelDipDir = 90.0;

    if (useDeposit) {
        std::cout << "=== Deposit Test — " << depositPath << " ===" << std::endl;
        csvPath = depositPath + "/Results/inversion_start.csv";
    } else if (useIni) {
        std::cout << "=== Voisey's Bay — INI-Driven Joint Inversion ===" << std::endl;
        std::cout << "Config: " << iniPath << std::endl;
        std::cout << "Output: " << outputDir << std::endl;

        IniConfig ini;
        if (!ini.load(iniPath)) {
            std::cerr << "FATAL: Could not load INI config from " << iniPath << std::endl;
            return 1;
        }

        // Resolve relative paths against INI directory
        std::filesystem::path iniDir = std::filesystem::path(iniPath).parent_path();
        auto resolvePath = [&](const std::string& p) -> std::string {
            if (p.empty()) return p;
            std::filesystem::path fp(p);
            if (fp.is_absolute()) return p;
            // If the path already resolves from CWD, use as-is
            // (avoids double-resolving when INI already has resolved paths)
            if (std::filesystem::exists(fp)) return p;
            std::filesystem::path resolved = iniDir / fp;
            if (std::filesystem::exists(resolved)) return resolved.string();
            // File doesn't exist at either location — prefer INI-relative
            return resolved.string();
        };

        config = InversionConfig::fromIni(ini);
        outputDir = ini.getString("output", "output_dir", outputDir);
        outputDir = resolvePath(outputDir);
        csvPath = resolvePath(ini.getString("data", "cluster_csv", csvPath));
        paddingRings = config.topography.paddingRings;
        std::filesystem::create_directories(outputDir);

        // Read group configuration
        groupColumn = ini.getString("data", "group_column", "litho_group_id");

        // Parse borehole classified CSV paths for ordering derivation
        std::string bhCsvsStr = ini.getString("data", "borehole_classified_csvs", "");
        if (!bhCsvsStr.empty()) {
            boreholeCsvPaths.clear();
            std::stringstream bss(bhCsvsStr);
            std::string btok;
            while (std::getline(bss, btok, ',')) {
                btok.erase(0, btok.find_first_not_of(" \t"));
                btok.erase(btok.find_last_not_of(" \t") + 1);
                if (!btok.empty()) boreholeCsvPaths.push_back(resolvePath(btok));
            }
        }

        // Parse real data mode keys
        useRealData = ini.getBool("data", "real_data", false);
        obsGravCsv = resolvePath(ini.getString("data", "observed_gravity_csv", ""));
        obsMagCsv = resolvePath(ini.getString("data", "observed_magnetic_csv", ""));
        localOriginX = ini.getDouble("data", "local_origin_x", 0.0);
        localOriginY = ini.getDouble("data", "local_origin_y", 0.0);
        zDatum = ini.getDouble("data", "z_datum", 0.0);

        // Parse surface CSV paths for external starting geometry
        std::string surfCsvsStr = ini.getString("data", "surface_csvs", "");
        if (!surfCsvsStr.empty()) {
            surfaceCsvPaths.clear();
            std::stringstream scs(surfCsvsStr);
            std::string stok;
            while (std::getline(scs, stok, ',')) {
                stok.erase(0, stok.find_first_not_of(" \t"));
                stok.erase(stok.find_last_not_of(" \t") + 1);
                if (!stok.empty()) surfaceCsvPaths.push_back(resolvePath(stok));
            }
        }

        // Parse contact mesh paths (TS/OBJ) for 3D mesh pipeline (layered)
        std::string contactMeshesStr = ini.getString("data", "contact_meshes", "");
        if (!contactMeshesStr.empty()) {
            std::stringstream cms(contactMeshesStr);
            std::string ctok;
            while (std::getline(cms, ctok, ',')) {
                ctok.erase(0, ctok.find_first_not_of(" \t"));
                ctok.erase(ctok.find_last_not_of(" \t") + 1);
                if (!ctok.empty()) contactMeshPaths.push_back(resolvePath(ctok));
            }
        }

        // Parse group mesh paths (TS/OBJ) for volumetric pipeline — one
        // closed boundary mesh per lithology group (no flat top/bottom needed).
        std::string groupMeshesStr = ini.getString("data", "group_meshes", "");
        if (!groupMeshesStr.empty()) {
            std::stringstream gms(groupMeshesStr);
            std::string gtok;
            while (std::getline(gms, gtok, ',')) {
                gtok.erase(0, gtok.find_first_not_of(" \t"));
                gtok.erase(gtok.find_last_not_of(" \t") + 1);
                if (!gtok.empty()) groupMeshPaths.push_back(resolvePath(gtok));
            }
        }

        // Parse label grid path (binary 3D label array from pipeline)
        labelGridPath = ini.getString("data", "label_grid", "");
        if (!labelGridPath.empty()) {
            labelGridPath = resolvePath(labelGridPath);
        }

        // Parse [group_mapping] section when using litho_group_id mode
        if (groupColumn == "litho_group_id") {
            int numGroups = ini.getInt("group_mapping", "num_groups", 0);
            if (numGroups > 0) {
                for (int g = 0; g < numGroups; ++g) {
                    std::string key = "group_" + std::to_string(g) + "_clusters";
                    std::string csv = ini.getString("group_mapping", key, "");
                    std::vector<int> cids;
                    std::stringstream ss(csv);
                    std::string tok;
                    while (std::getline(ss, tok, ',')) {
                        tok.erase(0, tok.find_first_not_of(" \t"));
                        tok.erase(tok.find_last_not_of(" \t") + 1);
                        if (!tok.empty()) {
                            int cid = std::stoi(tok);
                            cids.push_back(cid);
                            clusterMap.push_back({cid, g});
                        }
                    }
                    groupClusters.push_back(cids);
                }
                std::cout << "Group mapping: " << numGroups << " groups from INI [group_mapping]"
                          << std::endl;
            } else {
                std::cerr << "FATAL: group_column=litho_group_id but no [group_mapping] section"
                          << std::endl;
                return 1;
            }
        }

        // Parse [true_model] section
        trueModelMode = ini.getString("true_model", "mode", "voiseys");
        trueModelDepths.clear();
        std::string depthsStr = ini.getString("true_model", "surface_depths", "");
        if (!depthsStr.empty()) {
            std::stringstream dss(depthsStr);
            std::string dtok;
            while (std::getline(dss, dtok, ',')) {
                dtok.erase(0, dtok.find_first_not_of(" \t"));
                dtok.erase(dtok.find_last_not_of(" \t") + 1);
                if (!dtok.empty()) trueModelDepths.push_back(std::stod(dtok));
            }
        }
        trueModelDip = ini.getDouble("true_model", "dip_angle", 0.0);
        trueModelDipDir = ini.getDouble("true_model", "dip_direction", 90.0);

        // Save resolved config for reproducibility
        IniConfig resolved = config.toIni();
        resolved.setString("data", "group_column", groupColumn);
        resolved.setBool("data", "real_data", useRealData);
        resolved.setString("data", "observed_gravity_csv", obsGravCsv);
        resolved.setString("data", "observed_magnetic_csv", obsMagCsv);
        resolved.setDouble("data", "local_origin_x", localOriginX);
        resolved.setDouble("data", "local_origin_y", localOriginY);
        resolved.setDouble("data", "z_datum", zDatum);
        if (!csvPath.empty()) {
            resolved.setString("data", "cluster_csv", csvPath);
        }
        if (!groupMeshPaths.empty()) {
            std::stringstream gmss;
            for (size_t gi = 0; gi < groupMeshPaths.size(); ++gi) {
                if (gi > 0) gmss << ",";
                gmss << groupMeshPaths[gi];
            }
            resolved.setString("data", "group_meshes", gmss.str());
        }
        if (!contactMeshPaths.empty()) {
            std::stringstream cmss;
            for (size_t ci = 0; ci < contactMeshPaths.size(); ++ci) {
                if (ci > 0) cmss << ",";
                cmss << contactMeshPaths[ci];
            }
            resolved.setString("data", "contact_meshes", cmss.str());
        }
        if (!labelGridPath.empty()) {
            resolved.setString("data", "label_grid", labelGridPath);
        }
        if (!surfaceCsvPaths.empty()) {
            std::stringstream scss;
            for (size_t si = 0; si < surfaceCsvPaths.size(); ++si) {
                if (si > 0) scss << ",";
                scss << surfaceCsvPaths[si];
            }
            resolved.setString("data", "surface_csvs", scss.str());
        }
        resolved.setString("true_model", "mode", trueModelMode);
        if (!trueModelDepths.empty()) {
            std::stringstream dss;
            for (size_t i = 0; i < trueModelDepths.size(); ++i) {
                if (i > 0) dss << ",";
                dss << trueModelDepths[i];
            }
            resolved.setString("true_model", "surface_depths", dss.str());
        }
        resolved.setDouble("true_model", "dip_angle", trueModelDip);
        resolved.setDouble("true_model", "dip_direction", trueModelDipDir);
        resolved.save(outputDir + "/resolved_config.ini");
    } else {
        std::cout << "=== Voisey's Bay — CSV-Driven Joint Inversion ===" << std::endl;
    }
    if (!csvPath.empty()) {
        std::cout << "Cluster properties: " << csvPath << std::endl;
    } else {
        std::cout << "Properties: original hardcoded (no cluster CSV)" << std::endl;
    }
    std::cout << "Output: " << outputDir << std::endl;
    if (paddingRings > 0) {
        std::cout << "Padding rings: " << paddingRings << std::endl;
    }

    // Apply CLI overrides for true model (applies to all modes)
    if (!cliTrueModelMode.empty()) trueModelMode = cliTrueModelMode;
    if (!cliTrueModelDepths.empty()) {
        trueModelDepths.clear();
        std::stringstream dss(cliTrueModelDepths);
        std::string dtok;
        while (std::getline(dss, dtok, ',')) {
            dtok.erase(0, dtok.find_first_not_of(" \t"));
            dtok.erase(dtok.find_last_not_of(" \t") + 1);
            if (!dtok.empty()) trueModelDepths.push_back(std::stod(dtok));
        }
    }
    if (cliTrueModelDip >= 0) trueModelDip = cliTrueModelDip;
    if (cliTrueModelDipDir >= 0) trueModelDipDir = cliTrueModelDipDir;

    std::cout << "True model mode: " << trueModelMode;
    if (!trueModelDepths.empty()) {
        std::cout << " depths=[";
        for (size_t i = 0; i < trueModelDepths.size(); ++i) {
            if (i > 0) std::cout << ",";
            std::cout << trueModelDepths[i];
        }
        std::cout << "]";
    }
    if (trueModelMode == "dipping") {
        std::cout << " dip=" << trueModelDip << "° dir=" << trueModelDipDir << "°";
    }
    std::cout << std::endl;

    // =====================================================================
    // 1. Determine mode: cluster CSV vs original hardcoded properties
    // =====================================================================
    bool useClusterCsv = !csvPath.empty();
    if (useClusterCsv && !std::filesystem::exists(csvPath)) {
        std::cerr << "WARNING: cluster CSV not found at " << csvPath
                  << " — falling back to original hardcoded properties." << std::endl;
        useClusterCsv = false;
    }

    // Load cluster properties early — needed for both true model and initial model
    if (useClusterCsv) {
        clusters = loadClusterProperties(csvPath);
        if (clusters.empty()) {
            std::cerr << "FATAL: No cluster properties loaded from " << csvPath << std::endl;
            return 1;
        }
    }

    // Compute stratigraphic ordering from borehole data (or cluster_id ascending fallback)
    std::vector<int> ordering;
    if (useClusterCsv) {
        ordering = deriveOrderingFromBoreholes(clusters, boreholeCsvPaths);
    }

    // =====================================================================
    // 2. Generate true model and synthetic data
    // =====================================================================
    DepositData depositData;
    if (useDeposit) {
        depositData = loadDepositData(depositPath);
        clusterMap = deriveClusterMapping(depositData, 3);
        outputDir = depositPath + "/Inversion";
        std::filesystem::create_directories(outputDir);

        // Resolve parameters: CLI overrides > defaults
        int stride = (cliStride > 0) ? cliStride : 10;
        int maxIter = (cliMaxIter > 0) ? cliMaxIter : 50;
        double propMismatch = (cliPropertyMismatch >= 0.0) ? cliPropertyMismatch : 0.1;

        // Generate and load INI for reproducible deposit config
        std::string iniStr = generateDepositIni(depositData, csvPath, outputDir,
                                                stride, maxIter, propMismatch);
        std::string iniWritePath = outputDir + "/test_config.ini";
        {
            std::ofstream ofs(iniWritePath);
            ofs << iniStr;
        }
        IniConfig ini;
        ini.load(iniWritePath);
        config = InversionConfig::fromIni(ini);
        iniPath = iniWritePath;  // needed for remanence mode fixup below
        useIni = true; // route through INI config path below

        std::cout << "Config: max_iter=" << maxIter << " stride=" << stride
                  << " property_mismatch=" << propMismatch << std::endl;
    }

    std::shared_ptr<LithologyModel> trueModel;
    GravityData obsPoints;
    GravityData syntheticData;
    MagneticData syntheticMagData;
    std::vector<Constraint> constraints;

    if (useRealData) {
        // =====================================================================
        // Real observed data mode — skip true model, load from prepped CSV
        // =====================================================================
        std::cout << "\n=== Real Observed Data Mode ===" << std::endl;

        syntheticData = CSVGravityLoader::load(obsGravCsv);
        std::cout << "Loaded " << syntheticData.size() << " gravity points" << std::endl;

        if (!obsMagCsv.empty()) {
            syntheticMagData = CSVMagneticLoader::load(obsMagCsv);
            std::cout << "Loaded " << syntheticMagData.size() << " magnetic points" << std::endl;
        }

        // Apply coordinate transform (MGA50 → local)
        for (auto& pt : syntheticData) {
            pt.position.x() -= localOriginX;
            pt.position.y() -= localOriginY;
            pt.position.z() -= zDatum;
        }
        for (auto& pt : syntheticMagData) {
            pt.position.x() -= localOriginX;
            pt.position.y() -= localOriginY;
            pt.position.z() -= zDatum;
        }

        // No borehole constraints with real data
        constraints.clear();

    } else {

    // =====================================================================
    // 2. Generate true model and synthetic data (synthetic mode)
    // =====================================================================

    std::cout << "\n--- True model ---" << std::endl;
    if (useDeposit) {
        trueModel = generateDepositTrueModel(depositData, clusterMap);
        obsPoints = generateDepositObservationPoints(depositData);

        // Apply property mismatch for robustness testing.
        // The true model's density and susceptibility are perturbed so the
        // synthetic data reflects a different earth than the initial model
        // assumes — a realistic scenario where physical properties have error.
        double propMismatch = (cliPropertyMismatch >= 0.0) ? cliPropertyMismatch : 0.1;
        if (propMismatch > 0.0) {
            std::cout << "\n--- Applying property mismatch (±"
                      << (propMismatch * 100.0) << "%) ---" << std::endl;
            perturbModelProperties(*trueModel, propMismatch);
        }
    } else if (useClusterCsv) {
        if (groupColumn == "cluster_id") {
            if (trueModelMode == "layered") {
                if (trueModelDepths.empty()) {
                    // Auto-generate evenly-spaced depths
                    int nGroups = static_cast<int>(clusters.size());
                    int nSurfs = nGroups - 1;
                    double spacing = 5000.0 / nGroups;
                    for (int i = 0; i < nSurfs; ++i)
                        trueModelDepths.push_back(-(i + 1) * spacing);
                    std::cout << "Auto-generated " << nSurfs
                              << " evenly-spaced surface depths" << std::endl;
                }
                trueModel = generateLayeredTrueModel(clusters, trueModelDepths, ordering);
            } else if (trueModelMode == "dipping") {
                if (trueModelDepths.empty()) {
                    int nGroups = static_cast<int>(clusters.size());
                    int nSurfs = nGroups - 1;
                    double spacing = 5000.0 / nGroups;
                    for (int i = 0; i < nSurfs; ++i)
                        trueModelDepths.push_back(-(i + 1) * spacing);
                    std::cout << "Auto-generated " << nSurfs
                              << " evenly-spaced surface depths" << std::endl;
                }
                trueModel = generateDippingLayeredTrueModel(
                    clusters, trueModelDepths, trueModelDip, trueModelDipDir, ordering);
            } else {
                trueModel = generateClusterIdTrueModel(clusters, ordering);
            }
        } else {
            trueModel = generateTrueModelWithRemanence();
        }
        obsPoints = generateObservationPoints();
    } else {
        trueModel = generateTrueModel();
        obsPoints = generateObservationPoints();
    }

    // =====================================================================
    // 3. Synthetic data from true model
    // =====================================================================

    {
        syntheticData = computeSyntheticData(trueModel, obsPoints);

        RemanentMagnetizationMode synthMagMode = useClusterCsv
            ? RemanentMagnetizationMode::FixedVectorPerGroup
            : RemanentMagnetizationMode::EffectiveSusceptibility;
        syntheticMagData = computeSyntheticMagnetic(
            trueModel, obsPoints, magInc, magDec, magField, synthMagMode);
    }

    std::cout << "Gravity points: " << syntheticData.size() << std::endl;
    std::cout << "Magnetic points: " << syntheticMagData.size() << std::endl;

    if (useDeposit) {
        constraints = generateDepositConstraints(depositData, clusterMap, depositData.clusters);
    } else if (groupColumn == "cluster_id") {
        constraints = generateClusterIdConstraints(clusters, ordering);
    } else {
        constraints = generateLithosplitterConstraints();
    }

    } // end else (synthetic mode)

    // =====================================================================
    // 4. Build initial model — either from cluster CSV or original hardcoded
    // =====================================================================
    std::shared_ptr<LithologyModel> initialModel;
    RemanentMagnetizationMode remanenceMode;
    bool hasUnknownRemanence = false;
    std::vector<std::string> groupExportNames; // per-group name for export files

    constexpr int INTERIOR_DIM = 21;
    constexpr double GRID_HALF = 1000.0;
    constexpr double CELL_SIZE = 100.0;

    if (useClusterCsv) {
        // ---- Cluster CSV mode (clusters already loaded above) ----

        if (groupColumn == "cluster_id") {
            // ---- cluster_id mode: each cluster = one group (geological order) ----
            std::cout << "\n--- Initial model (cluster_id mode, " << clusters.size()
                      << " clusters) ---" << std::endl;
            initialModel = buildModelFromClustersDirect(clusters, -5000.0, ordering);

            // Build groupExportNames from ordering
            groupExportNames.clear();
            for (int cid : ordering) {
                groupExportNames.push_back("cluster_id_" + std::to_string(cid));
            }

            int nGroups = initialModel->groupCount();
            int nSurfs = nGroups - 1;

            // Build vertex grid positions
            std::vector<double> vxs, vys;
            for (int i = 0; i < INTERIOR_DIM; ++i) {
                vxs.push_back(-GRID_HALF + i * CELL_SIZE);
                vys.push_back(-GRID_HALF + i * CELL_SIZE);
            }

            auto buildFlatSurf = [&](const std::string& name, double z) {
                auto s = std::make_shared<SurfaceMesh>();
                s->setName(name);
                s->setBounds(-10000.0, 100.0);
                for (int iy = 0; iy < INTERIOR_DIM; ++iy)
                    for (int ix = 0; ix < INTERIOR_DIM; ++ix)
                        s->addVertex(vxs[ix], vys[iy], z, config.vertexFreedom);
                int nc = INTERIOR_DIM - 1;
                for (int iy = 0; iy < nc; ++iy)
                    for (int ix = 0; ix < nc; ++ix) {
                        uint32_t i0 = iy * INTERIOR_DIM + ix;
                        uint32_t i1 = i0 + 1;
                        uint32_t i2 = (iy + 1) * INTERIOR_DIM + ix;
                        uint32_t i3 = i2 + 1;
                        s->addTriangle(i0, i1, i2);
                        s->addTriangle(i1, i3, i2);
                    }
                return s;
            };

            auto buildDippingSurf = [&](const std::string& name, double z0,
                                         double dipAngle, double dipDir) {
                auto s = std::make_shared<SurfaceMesh>();
                s->setName(name);
                s->setBounds(-10000.0, 100.0);
                double tanDip = std::tan(dipAngle * M_PI / 180.0);
                double dipRad = dipDir * M_PI / 180.0;
                double sx = std::sin(dipRad);
                double sy = std::cos(dipRad);
                for (int iy = 0; iy < INTERIOR_DIM; ++iy)
                    for (int ix = 0; ix < INTERIOR_DIM; ++ix) {
                        double z = z0 - tanDip * (vxs[ix] * sx + vys[iy] * sy);
                        s->addVertex(vxs[ix], vys[iy], z, config.vertexFreedom);
                    }
                int nc = INTERIOR_DIM - 1;
                for (int iy = 0; iy < nc; ++iy)
                    for (int ix = 0; ix < nc; ++ix) {
                        uint32_t i0 = iy * INTERIOR_DIM + ix;
                        uint32_t i1 = i0 + 1;
                        uint32_t i2 = (iy + 1) * INTERIOR_DIM + ix;
                        uint32_t i3 = i2 + 1;
                        s->addTriangle(i0, i1, i2);
                        s->addTriangle(i1, i3, i2);
                    }
                return s;
            };

            // Build surface from CSV point cloud using IDW interpolation
            auto buildSurfFromCSV = [&](const std::string& name, const std::string& csvPath) {
                std::vector<Vector3d> pts;
                std::ifstream file(csvPath);
                if (!file.is_open()) {
                    std::cerr << "WARNING: Cannot open surface CSV: " << csvPath << std::endl;
                    return std::shared_ptr<SurfaceMesh>(nullptr);
                }
                std::string line;
                int lineno = 0;
                while (std::getline(file, line)) {
                    ++lineno;
                    if (lineno == 1 || line.empty() || line[0] == '#') continue;
                    std::stringstream ss(line);
                    std::string tok;
                    double x = 0, y = 0, z = 0;
                    int col = 0;
                    while (std::getline(ss, tok, ',') && col < 3) {
                        try {
                            double v = std::stod(tok);
                            if (col == 0) x = v;
                            else if (col == 1) y = v;
                            else z = v;
                        } catch (...) {}
                        ++col;
                    }
                    pts.push_back(Vector3d(x, y, z));
                }
                if (pts.empty()) {
                    std::cerr << "WARNING: No points in surface CSV: " << csvPath << std::endl;
                    return std::shared_ptr<SurfaceMesh>(nullptr);
                }

                // IDW interpolation of z at each grid point
                auto idwZ = [&](double gx, double gy, const std::vector<Vector3d>& p) {
                    double sw = 0.0, swz = 0.0;
                    for (auto& pt : p) {
                        double d2 = (gx - pt.x())*(gx - pt.x()) + (gy - pt.y())*(gy - pt.y());
                        if (d2 < 1.0) d2 = 1.0;
                        double w = 1.0 / d2;
                        sw += w; swz += w * pt.z();
                    }
                    return swz / sw;
                };

                auto s = std::make_shared<SurfaceMesh>();
                s->setName(name);
                s->setBounds(-10000.0, 100.0);
                for (int iy = 0; iy < INTERIOR_DIM; ++iy)
                    for (int ix = 0; ix < INTERIOR_DIM; ++ix)
                        s->addVertex(vxs[ix], vys[iy],
                                     idwZ(vxs[ix], vys[iy], pts),
                                     config.vertexFreedom);
                int nc = INTERIOR_DIM - 1;
                for (int iy = 0; iy < nc; ++iy)
                    for (int ix = 0; ix < nc; ++ix) {
                        uint32_t i0 = iy * INTERIOR_DIM + ix;
                        uint32_t i1 = i0 + 1;
                        uint32_t i2 = (iy + 1) * INTERIOR_DIM + ix;
                        uint32_t i3 = i2 + 1;
                        s->addTriangle(i0, i1, i2);
                        s->addTriangle(i1, i3, i2);
                    }
                std::cout << "  Loaded " << pts.size() << " points from " << csvPath << std::endl;
                return s;
            };

            // Determine initial surface depths
            std::vector<double> initSurfZs(nSurfs, 0.0);
            if (!trueModelDepths.empty() && static_cast<int>(trueModelDepths.size()) == nSurfs) {
                // Use true model surface depths (offset slightly so optimizer has work to do)
                for (int i = 0; i < nSurfs; ++i)
                    initSurfZs[i] = trueModelDepths[i] - 20.0;
            } else if (!boreholeCsvPaths.empty()) {
                // Derive surface depths from borehole first-appearance data
                std::vector<double> minDepth(nGroups, 1e18);
                std::vector<int> allCids = ordering;
                for (const auto& fpath : boreholeCsvPaths) {
                    std::ifstream file(fpath);
                    if (!file.is_open()) continue;
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
                            for (int gi = 0; gi < nGroups; ++gi) {
                                if (allCids[gi] == cid && from_m < minDepth[gi])
                                    minDepth[gi] = from_m;
                            }
                        }
                    }
                }
                // Surface i separates group i (above) from group i+1 (below).
                // Use the shallowest depth of the lower group as the contact.
                for (int i = 0; i < nSurfs; ++i) {
                    if (minDepth[i + 1] < 1e17)
                        initSurfZs[i] = -(minDepth[i + 1] + 5.0); // slightly below
                    else
                        initSurfZs[i] = -(50.0 + i * 100.0); // fallback
                }
            } else {
                // Shallow evenly-spaced fallback
                for (int i = 0; i < nSurfs; ++i)
                    initSurfZs[i] = -(50.0 + i * 100.0);
            }

            // Build N-1 contact surfaces.
            // Map contact mesh paths by group pair (1-based from Python filenames
            // like contact_1_2.ts → groups (0,1) after 0-based conversion).
            std::map<std::pair<int,int>, std::string> contactMeshMap;
            for (const auto& p : contactMeshPaths) {
                std::filesystem::path fp(p);
                std::string stem = fp.stem().string(); // e.g. "contact_1_2"
                // Parse "contact_GA_GB"
                int ga = -1, gb = -1;
                size_t u1 = stem.rfind('_');
                if (u1 != std::string::npos && u1 > 0) {
                    size_t u0 = stem.rfind('_', u1 - 1);
                    if (u0 != std::string::npos) {
                        try {
                            ga = std::stoi(stem.substr(u0 + 1, u1 - u0 - 1)) - 1; // 0-based
                            gb = std::stoi(stem.substr(u1 + 1)) - 1;
                        } catch (...) { ga = -1; gb = -1; }
                    }
                }
                if (ga >= 0 && gb >= 0)
                    contactMeshMap[{ga, gb}] = p;
            }

            bool useGroupMeshes = !groupMeshPaths.empty()
                                && static_cast<int>(groupMeshPaths.size()) == nGroups;
            bool useContactMeshes = !contactMeshMap.empty();
            bool useSurfaceCsvs = (!useGroupMeshes && !useContactMeshes && !surfaceCsvPaths.empty()
                                   && static_cast<int>(surfaceCsvPaths.size()) == nSurfs);

            if (useGroupMeshes) {
                // ---- Volumetric path: load pre-built closed group volumes ----
                std::cout << "  Loading " << nGroups << " group volume meshes"
                          << " (volumetric pipeline)" << std::endl;

                // Helper: look up cluster properties by cluster_id
                auto findCluster = [&](int cid) -> const ClusterProperties* {
                    for (const auto& c : clusters)
                        if (c.cluster_id == cid) return &c;
                    return nullptr;
                };

                int origNumGroups = nGroups;
                for (int g = 0; g < origNumGroups; ++g) {
                    auto mesh = loadSurfaceMesh(groupMeshPaths[g]);
                    if (!mesh) {
                        std::cerr << "ERROR: Failed to load group mesh: "
                                  << groupMeshPaths[g] << std::endl;
                        return 1;
                    }

                    // Coordinate transform: world → local
                    for (uint32_t vi = 0; vi < mesh->vertexCount(); ++vi) {
                        auto& v = mesh->vertex(vi);
                        v.position.x() -= localOriginX;
                        v.position.y() -= localOriginY;
                        v.position.z() -= zDatum;
                    }

                    mesh->buildNeighbors();

                    // Detect disconnected components — split into separate
                    // groups so property inversion can assign independent
                    // density/susceptibility to each disconnected volume.
                    auto components = splitDisconnectedComponents(*mesh);
                    int nComp = static_cast<int>(components.size());

                    // Assign first component to the original group
                    components[0]->setName("group_" + groupExportNames[g]);
                    initialModel->setGroupMesh(g, components[0]);
                    std::cout << "    Group " << g << " (" << groupExportNames[g] << "): "
                              << components[0]->vertexCount() << " vertices, "
                              << components[0]->triangleCount() << " triangles, "
                              << components[0]->dofCount() << " DOFs";
                    if (nComp > 1) {
                        std::cout << " [" << nComp << " disconnected components]";
                    }
                    std::cout << std::endl;

                    // Create additional groups for extra disconnected components
                    for (int c = 1; c < nComp; ++c) {
                        int newG = initialModel->groupCount();
                        int origClusterId = ordering[g];
                        const ClusterProperties* cp = findCluster(origClusterId);

                        std::string newName = "Cluster_ID_" + std::to_string(origClusterId)
                                            + "_vol" + std::to_string(c);
                        LithoGroup newGroup(newG, newName,
                                          cp ? cp->density_median : 2.67,
                                          cp ? cp->susceptibility_median : 0.0);
                        if (cp) {
                            newGroup.remanence_magnitude = cp->remanence_magnitude;
                            newGroup.remanence_inclination = cp->remanence_inclination;
                            newGroup.remanence_declination = cp->remanence_declination;
                        }
                        initialModel->addGroup(newGroup);

                        std::string exportName = groupExportNames[g] + "_vol" + std::to_string(c);
                        groupExportNames.push_back(exportName);

                        components[c]->setName("group_" + exportName);
                        initialModel->setGroupMesh(newG, components[c]);

                        ordering.push_back(origClusterId);

                        std::cout << "    Group " << newG << " (" << exportName
                                  << ") [split from cluster_id=" << origClusterId << "]: "
                                  << components[c]->vertexCount() << " vertices, "
                                  << components[c]->triangleCount() << " triangles, "
                                  << components[c]->dofCount() << " DOFs" << std::endl;
                    }
                }

                // Apply shared model setup: setBounds, setAllFreedom,
                // fixExteriorFaces, applyModelBounds, setBottomDepth
                finalizeModelSetup(initialModel, config.vertexFreedom);

                std::cout << "  Bottom depth: " << initialModel->bottomDepth()
                          << " m" << std::endl;

                // Load label grid if provided
                if (!labelGridPath.empty()) {
                    std::ifstream lgf(labelGridPath, std::ios::binary);
                    if (lgf.is_open()) {
                        int lnx, lny, lnz;
                        double lx0, ly0, lz0, ldx, ldy, ldz;
                        lgf >> lnx >> lny >> lnz >> lx0 >> ly0 >> lz0 >> ldx >> ldy >> ldz;
                        lgf.get(); // consume newline after header
                        std::vector<int> labelData(lnx * lny * lnz);
                        lgf.read(reinterpret_cast<char*>(labelData.data()),
                                 labelData.size() * sizeof(int));
                        if (lgf.gcount() == static_cast<std::streamsize>(
                                labelData.size() * sizeof(int))) {
                            initialModel->setLabelGrid(labelData, lnx, lny, lnz,
                                                       lx0, ly0, lz0, ldx, ldy, ldz);
                            std::cout << "  Label grid: " << lnx << "x" << lny << "x" << lnz
                                      << " cells loaded" << std::endl;
                        } else {
                            std::cerr << "WARNING: label grid file truncated (got "
                                      << lgf.gcount() << " bytes, expected "
                                      << labelData.size() * sizeof(int) << ")"
                                      << std::endl;
                        }
                    } else {
                        std::cerr << "WARNING: could not open label grid: "
                                  << labelGridPath << std::endl;
                    }
                }

            } else if (useContactMeshes) {
                std::cout << "  Loading " << contactMeshMap.size()
                          << " contact meshes (3D pipeline, " << nSurfs
                          << " surfaces needed)" << std::endl;

                int mshLoaded = 0, mshFallback = 0;
                std::vector<std::shared_ptr<SurfaceMesh>> contactSurfs(nSurfs);
                for (int i = 0; i < nSurfs; ++i) {
                    int ga = i;       // 0-based group above
                    int gb = i + 1;   // 0-based group below

                    std::shared_ptr<SurfaceMesh> surf;
                    auto it = contactMeshMap.find({ga, gb});
                    if (it != contactMeshMap.end()) {
                        surf = loadSurfaceMesh(it->second);
                        if (!surf) {
                            std::cerr << "ERROR: Failed to load contact mesh: "
                                      << it->second << std::endl;
                            return 1;
                        }
                        // Coordinate transform: world (MGA50) → local
                        for (uint32_t vi = 0; vi < surf->vertexCount(); ++vi) {
                            auto& v = surf->vertex(vi);
                            v.position.x() -= localOriginX;
                            v.position.y() -= localOriginY;
                            v.position.z() -= zDatum;
                        }
                        ++mshLoaded;
                    } else {
                        // Fallback: generate flat surface at estimated depth
                        double zFallback = initSurfZs[i];
                        surf = buildFlatSurf("fallback", zFallback);
                        ++mshFallback;
                    }

                    std::string sname = "contact_" + groupExportNames[i] + "_"
                                      + groupExportNames[i + 1];
                    surf->setName(sname);
                    surf->setBounds(-10000.0, 100.0);
                    surf->setAllFreedom(config.vertexFreedom);

                    if (config.meshBoundaryMode == "fixed") {
                        std::map<std::pair<uint32_t,uint32_t>, int> edgeCount;
                        for (uint32_t ti = 0; ti < surf->triangleCount(); ++ti) {
                            const auto& tri = surf->triangle(ti);
                            auto mkEdge = [](uint32_t a, uint32_t b) {
                                return std::make_pair(std::min(a,b), std::max(a,b));
                            };
                            edgeCount[mkEdge(tri.v0, tri.v1)]++;
                            edgeCount[mkEdge(tri.v1, tri.v2)]++;
                            edgeCount[mkEdge(tri.v2, tri.v0)]++;
                        }
                        int nFixed = 0;
                        for (const auto& [e, cnt] : edgeCount) {
                            if (cnt == 1) {
                                surf->setVertexFreedom(e.first, VertexFreedom::FIXED);
                                surf->setVertexFreedom(e.second, VertexFreedom::FIXED);
                                ++nFixed;
                            }
                        }
                        std::cout << "    Fixed boundary vertices on " << nFixed
                                  << " boundary edges" << std::endl;
                    }

                    surf->buildNeighbors();

                    std::cout << "    Surface " << i << " (" << sname << "): "
                              << surf->vertexCount() << " vertices, "
                              << surf->triangleCount() << " triangles, "
                              << surf->dofCount() << " DOFs" << std::endl;

                    contactSurfs[i] = surf;
                }

                std::cout << "  Loaded: " << mshLoaded
                          << " from mesh files, " << mshFallback << " fallback flat"
                          << std::endl;

                // Build closed group meshes from contact surfaces
                auto flatTop = buildFlatMesh(*contactSurfs[0], 0.0, "flat_top");
                auto flatBottom = buildFlatMesh(*contactSurfs[0], initialModel->bottomDepth(), "flat_bottom");
                for (int g = 0; g < nGroups; ++g) {
                    const SurfaceMesh* top = (g == 0) ? flatTop.get() : contactSurfs[g - 1].get();
                    const SurfaceMesh* bot = (g == nGroups - 1) ? flatBottom.get() : contactSurfs[g].get();
                    initialModel->setGroupMesh(g, buildClosedMesh(*top, *bot));
                }
            } else if (useSurfaceCsvs) {
                std::cout << "  Using surface CSVs for " << nSurfs << " contacts" << std::endl;
                std::vector<std::shared_ptr<SurfaceMesh>> csvSurfs(nSurfs);
                for (int i = 0; i < nSurfs; ++i) {
                    std::string sname = "contact_" + groupExportNames[i] + "_"
                                      + groupExportNames[i + 1];
                    auto s = buildSurfFromCSV(sname, surfaceCsvPaths[i]);
                    if (s) csvSurfs[i] = s;
                }
                auto flatTop = buildFlatMesh(*csvSurfs[0], 0.0, "flat_top");
                auto flatBottom = buildFlatMesh(*csvSurfs[0], initialModel->bottomDepth(), "flat_bottom");
                for (int g = 0; g < nGroups; ++g) {
                    const SurfaceMesh* top = (g == 0) ? flatTop.get() : csvSurfs[g - 1].get();
                    const SurfaceMesh* bot = (g == nGroups - 1) ? flatBottom.get() : csvSurfs[g].get();
                    initialModel->setGroupMesh(g, buildClosedMesh(*top, *bot));
                }
            } else {
                bool useDipping = (trueModelDip != 0.0);
                std::vector<std::shared_ptr<SurfaceMesh>> initSurfs(nSurfs);
                for (int i = 0; i < nSurfs; ++i) {
                    double z = initSurfZs[i];
                    std::string sname = "contact_" + groupExportNames[i] + "_"
                                      + groupExportNames[i + 1];
                    if (useDipping) {
                        initSurfs[i] = buildDippingSurf(sname, z, trueModelDip, trueModelDipDir);
                        std::cout << "  Surface " << i << " (" << sname << ") at z=" << z
                                  << "m (center), dipping " << trueModelDip << "deg toward "
                                  << trueModelDipDir << "deg" << std::endl;
                    } else {
                        initSurfs[i] = buildFlatSurf(sname, z);
                        std::cout << "  Surface " << i << " (" << sname << ") at z=" << z << "m"
                                  << std::endl;
                    }
                }
                auto flatTop = buildFlatMesh(*initSurfs[0], 0.0, "flat_top");
                auto flatBottom = buildFlatMesh(*initSurfs[0], initialModel->bottomDepth(), "flat_bottom");
                for (int g = 0; g < nGroups; ++g) {
                    const SurfaceMesh* top = (g == 0) ? flatTop.get() : initSurfs[g - 1].get();
                    const SurfaceMesh* bot = (g == nGroups - 1) ? flatBottom.get() : initSurfs[g].get();
                    initialModel->setGroupMesh(g, buildClosedMesh(*top, *bot));
                }
            }

            hasUnknownRemanence = configureRemanenceFromClustersDirect(
                initialModel, clusters, remanenceMode, ordering);

        } else {
            // ---- litho_group_id mode: clusters→groups via mapping ----
            std::cout << "\n--- Initial model (litho_group_id mode) ---" << std::endl;
            initialModel = buildModelFromClusters(clusters, clusterMap, -5000.0);

            // Build groupExportNames from cluster composition
            for (int g = 0; g < initialModel->groupCount(); ++g) {
                std::string name = "clusters";
                // Collect cluster IDs for this group
                std::vector<int> cids;
                for (int gi = 0; gi < static_cast<int>(groupClusters.size()); ++gi) {
                    if (gi == g) {
                        cids = groupClusters[gi];
                        break;
                    }
                }
                std::sort(cids.begin(), cids.end());
                for (int cid : cids) {
                    name += "_" + std::to_string(cid);
                }
                groupExportNames.push_back(name);
            }

            // Build interpolated surfaces from borehole contact depths (IDW)
            {
                std::vector<double> vxs, vys;
                for (int i = 0; i < INTERIOR_DIM; ++i) {
                    vxs.push_back(-GRID_HALF + i * CELL_SIZE);
                    vys.push_back(-GRID_HALF + i * CELL_SIZE);
                }

                struct Cpt { double x, y, z; };
                std::vector<Cpt> s0Pts = {
                    {   0.0,   0.0, -153.0},
                    { 300.0,   0.0,  -27.0},
                    {-115.0,   0.0,  -15.0},
                    {-500.0,   0.0,   -5.0},
                    {   0.0, 200.0,   -5.0},
                    {-1000.0, -1000.0, -5.0},
                    { 1000.0, -1000.0, -5.0},
                    {-1000.0,  1000.0, -5.0},
                    { 1000.0,  1000.0, -5.0},
                };
                std::vector<Cpt> s1Pts = {
                    {   0.0,   0.0, -162.0},
                    { 300.0,   0.0,  -81.0},
                    {-115.0,   0.0,  -57.0},
                    {-500.0,   0.0, -500.0},
                    {   0.0, 200.0, -500.0},
                    {-1000.0, -1000.0, -500.0},
                    { 1000.0, -1000.0, -500.0},
                    {-1000.0,  1000.0, -500.0},
                    { 1000.0,  1000.0, -500.0},
                };
                std::vector<Cpt> s2Pts = {
                    {   0.0,   0.0, -162.0},
                    { 300.0,   0.0, -500.0},
                    {-115.0,   0.0, -105.0},
                    {-500.0,   0.0, -500.0},
                    {   0.0, 200.0, -500.0},
                    {-1000.0, -1000.0, -500.0},
                    { 1000.0, -1000.0, -500.0},
                    {-1000.0,  1000.0, -500.0},
                    { 1000.0,  1000.0, -500.0},
                };

                auto idwZ = [](double x, double y, const std::vector<Cpt>& pts) {
                    double sw = 0.0, swz = 0.0;
                    for (auto& p : pts) {
                        double d2 = (x-p.x)*(x-p.x) + (y-p.y)*(y-p.y);
                        if (d2 < 1.0) d2 = 1.0;
                        double w = 1.0 / d2;
                        sw += w; swz += w * p.z;
                    }
                    return swz / sw;
                };

                auto buildSurf = [&](const std::string& name, const std::vector<Cpt>& pts) {
                    auto s = std::make_shared<SurfaceMesh>();
                    s->setName(name);
                    s->setBounds(-10000.0, 100.0);
                    for (int iy = 0; iy < INTERIOR_DIM; ++iy)
                        for (int ix = 0; ix < INTERIOR_DIM; ++ix)
                            s->addVertex(vxs[ix], vys[iy],
                                         idwZ(vxs[ix], vys[iy], pts),
                                         config.vertexFreedom);
                    int nc = INTERIOR_DIM - 1;
                    for (int iy = 0; iy < nc; ++iy)
                        for (int ix = 0; ix < nc; ++ix) {
                            uint32_t i0 = iy * INTERIOR_DIM + ix;
                            uint32_t i1 = i0 + 1;
                            uint32_t i2 = (iy + 1) * INTERIOR_DIM + ix;
                            uint32_t i3 = i2 + 1;
                            s->addTriangle(i0, i1, i2);
                            s->addTriangle(i1, i3, i2);
                        }
                    return s;
                };

                int nGroups = initialModel->groupCount();
                int nSurfs = nGroups - 1;
                // Use generic naming if groupExportNames is populated
                std::string s0name = "contact_" + groupExportNames[0] + "_" + groupExportNames[1];
                std::string s1name = "contact_" + groupExportNames[1] + "_" + groupExportNames[2];
                std::string s2name = "contact_" + groupExportNames[2] + "_" + groupExportNames[3];
                auto surf0 = buildSurf(s0name, s0Pts);
                auto surf1 = buildSurf(s1name, s1Pts);
                auto surf2 = buildSurf(s2name, s2Pts);

                auto flatTop = buildFlatMesh(*surf0, 0.0, "flat_top");
                auto flatBottom = buildFlatMesh(*surf0, initialModel->bottomDepth(), "flat_bottom");
                initialModel->setGroupMesh(0, buildClosedMesh(*flatTop, *surf0));
                initialModel->setGroupMesh(1, buildClosedMesh(*surf0, *surf1));
                initialModel->setGroupMesh(2, buildClosedMesh(*surf1, *surf2));
                initialModel->setGroupMesh(3, buildClosedMesh(*surf2, *flatBottom));

                std::cout << "  Surface depths at boreholes:" << std::endl;
                for (auto& bh : {Cpt{0.0, 0.0, 0}, Cpt{300.0, 0.0, 0},
                                 Cpt{-115.0, 0.0, 0}, Cpt{-500.0, 0.0, 0},
                                 Cpt{0.0, 200.0, 0}}) {
                    double z0 = idwZ(bh.x, bh.y, s0Pts);
                    double z1 = idwZ(bh.x, bh.y, s1Pts);
                    double z2 = idwZ(bh.x, bh.y, s2Pts);
                    std::cout << "    (" << bh.x << "," << bh.y
                              << "): S0=" << std::fixed << std::setprecision(1) << z0
                              << "m  S1=" << z1 << "m  S2=" << z2 << "m" << std::endl;
                }
            }

            hasUnknownRemanence = configureRemanenceFromClusters(
                initialModel, clusters, clusterMap, remanenceMode);
        }
    } else {
        // ---- Original hardcoded properties mode ----
        std::cout << "\n--- Initial model (original hardcoded properties) ---" << std::endl;

        if (paddingRings > 0) {
            // Build manually with padded grids
            initialModel = std::make_shared<LithologyModel>();
            initialModel->addGroup({0, "MaficComplex_1", 2.93, 0.016});
            initialModel->addGroup({1, "GraniteGneiss",  2.68, 0.000});
            initialModel->addGroup({2, "MaficComplex_2", 2.93, 0.016});
            initialModel->addGroup({3, "MassiveSulfide", 4.17, 0.080});
            initialModel->setBottomDepth(-5000.0);

            auto s0 = makePaddedGrid("mafic1_granitegneiss", -100.0, INTERIOR_DIM, paddingRings, CELL_SIZE, GRID_HALF, config.vertexFreedom);
            auto s1 = makePaddedGrid("granitegneiss_mafic2", -200.0, INTERIOR_DIM, paddingRings, CELL_SIZE, GRID_HALF, config.vertexFreedom);
            auto s2 = makePaddedGrid("mafic2_sulfide", -300.0, INTERIOR_DIM, paddingRings, CELL_SIZE, GRID_HALF, config.vertexFreedom);

            auto flatTop = buildFlatMesh(*s0, 0.0, "flat_top");
            auto flatBottom = buildFlatMesh(*s0, initialModel->bottomDepth(), "flat_bottom");
            initialModel->setGroupMesh(0, buildClosedMesh(*flatTop, *s0));
            initialModel->setGroupMesh(1, buildClosedMesh(*s0, *s1));
            initialModel->setGroupMesh(2, buildClosedMesh(*s1, *s2));
            initialModel->setGroupMesh(3, buildClosedMesh(*s2, *flatBottom));
        } else {
            initialModel = generateInitialModel();
        }

        // Original mode: induced-only (EffectiveSusceptibility)
        remanenceMode = RemanentMagnetizationMode::EffectiveSusceptibility;
    }

    std::cout << "\nInitial model:" << std::endl;
    for (int g = 0; g < initialModel->groupCount(); ++g) {
        const auto& grp = initialModel->group(g);
        std::cout << "  [" << grp.id << "] " << grp.name
                  << ": rho=" << std::fixed << std::setprecision(3) << grp.density
                  << "  chi=" << std::scientific << std::setprecision(4) << grp.susceptibility
                  << "  M_rem=" << std::fixed << std::setprecision(2)
                  << grp.remanence_magnitude << " A/m"
                  << "  I_rem=" << grp.remanence_inclination
                  << "  D_rem=" << grp.remanence_declination << std::endl;
    }

    // =====================================================================
    // 6. Export starting model
    // =====================================================================
    std::cout << "\n--- Exporting starting model ---" << std::endl;
    {
        InversionExporter initExp(outputDir, "csv_start");
        initExp.setSubfolder("starting");
        if (!groupExportNames.empty()) {
            initExp.setGroupNaming(groupExportNames);
        }

        initExp.exportStartingModel(*initialModel);
        std::cout << "  Wrote starting model to " << outputDir << "/starting/" << std::endl;
    }

    // =====================================================================
    // 7. Configure inversion
    // =====================================================================
    // Common config (both modes)
    config.model = initialModel;
    config.observedData = syntheticData;
    config.magneticData = syntheticMagData;
    config.constraints = constraints;

    if (!useIni) {
        // Legacy hardcoded config
        config.lambda = 0.0;
        config.maxIterations = 30;
        config.tolerance = 1e-4;
        config.lbfgsHistory = 10;
        config.controlPointStride = 4;

        config.enablePropertyInversion = true;
        config.propertyInversionInterval = 15;
        config.propertyInversionMaxIter = 10;
        config.propertyDensityMin = 1.5;
        config.propertyDensityMax = 6.0;

        config.enablePaddingGroup = true;
        config.paddingDensityInitial = 2.68;
        config.paddingDensityLower = 1.5;
        config.paddingDensityUpper = 6.0;
        config.paddingDepth = -100000.0;

        config.magneticWeight = 0.1;
        config.magneticInclination = magInc;
        config.magneticDeclination = magDec;
        config.magneticField_nT = magField;
        config.remanenceMode = remanenceMode;

        config.propertySusceptibilityMin = 0.0;
        config.propertySusceptibilityMax = 1.0;

        if (hasUnknownRemanence) {
            config.propertyRemanenceMin = 0.0;
            config.propertyRemanenceMax = 10.0;
        }

        config.omega = 1e1;
    } else {
        // INI mode: config already populated from INI, fill in runtime fields
        // If remanence mode wasn't set in INI, use the mode derived from the data source
        if (!iniPath.empty()) {
            IniConfig ini;
            ini.load(iniPath);
            if (ini.getString("magnetic", "remanence_mode", "").empty()) {
                config.remanenceMode = remanenceMode;
                if (hasUnknownRemanence) {
                    config.propertyRemanenceMin = 0.0;
                    config.propertyRemanenceMax = 10.0;
                }
            }
        }
    }

    // =====================================================================
    // 8. Export true model (synthetic mode only)
    // =====================================================================
    if (!useRealData && trueModel) {
        std::cout << "\n--- Exporting true model ---" << std::endl;
        {
            InversionExporter trueExp(outputDir, "csv_true");
            trueExp.setSubfolder("true_model");
            if (!groupExportNames.empty()) {
                trueExp.setGroupNaming(groupExportNames);
            }
            trueExp.exportStartingModel(*trueModel);
            std::cout << "  Wrote true model to " << outputDir << "/true_model/" << std::endl;
        }
    }

    // =====================================================================
    // 9. Run inversion
    // =====================================================================
    std::cout << "\n--- Running inversion (mode: "
              << (config.remanenceMode == RemanentMagnetizationMode::FixedVectorPerGroup
                  ? "FixedVectorPerGroup" : "EffectiveSusceptibility")
              << ") ---" << std::endl;

    // Create inprogress root directory
    std::string inprogressDir = outputDir + "/inprogress";
    std::filesystem::create_directories(inprogressDir);

    config.groupExportNames = groupExportNames;
    InversionRunner runner(config);
    runner.setIterationCallback([initialModel, inprogressDir,
                                  &syntheticData, &syntheticMagData,
                                  magInc, magDec, magField](
                                    const InversionIteration& iter) {
        std::cout << "  Iter " << std::setw(3) << iter.iteration
                  << ": RMS=" << std::fixed << std::setprecision(4)
                  << iter.rmsError
                  << "  obj=" << std::scientific << std::setprecision(3)
                  << iter.totalObjective
                  << "  DW_g=(" << std::fixed << std::setprecision(2)
                  << iter.dw_gravity_x << "," << iter.dw_gravity_y << ")";
        if (iter.magneticMisfit > 0.0) {
            std::cout << "  mag=" << std::scientific << std::setprecision(4)
                      << iter.magneticMisfit
                      << " DW_m=(" << std::fixed << std::setprecision(2)
                      << iter.dw_magnetic_x << "," << iter.dw_magnetic_y << ")";
        }
        if (iter.constraintPenalty > 0.0) {
            std::cout << "  pen=" << std::scientific << std::setprecision(3)
                      << iter.constraintPenalty;
        }
        std::cout << std::flush;

        // --- SVG diagnostics: observed vs predicted ---
        std::string iterDir = inprogressDir + "/iter_"
            + (iter.iteration < 10 ? "00" : iter.iteration < 100 ? "0" : "")
            + std::to_string(iter.iteration);
        std::filesystem::create_directories(iterDir);
        if (!syntheticData.empty()) {
            VectorXd params = initialModel->assembleParameterVector();
            auto gravFwd = std::make_shared<GravityForward>(initialModel, syntheticData);
            VectorXd gravPred = gravFwd->compute(params);

            SVGPlot plot(800, 600, 1, syntheticMagData.empty() ? 1 : 2);
            plot.title = "Iter " + std::to_string(iter.iteration) + " Diagnostics";

            // Panel 0: Gravity
            std::vector<double> gx, gy;
            for (size_t i = 0; i < syntheticData.size(); ++i) {
                gx.push_back(syntheticData[i].g_obs);
                gy.push_back(gravPred[static_cast<Eigen::Index>(i)]);
            }
            plot.panel(0).title = "Gravity";
            plot.panel(0).xlabel = "Observed (mGal)";
            plot.panel(0).ylabel = "Predicted (mGal)";
            plot.panel(0).scatters.push_back({gx, gy, "blue", "", 2.0});
            double gmin = *std::min_element(gx.begin(), gx.end());
            double gmax = *std::max_element(gx.begin(), gx.end());
            plot.panel(0).xmin = gmin; plot.panel(0).xmax = gmax;
            plot.panel(0).ymin = gmin; plot.panel(0).ymax = gmax;
            plot.panel(0).hlines.push_back({0.0, "gray", 0.5, "dashed", ""});

            // Panel 1: Magnetics
            if (!syntheticMagData.empty()) {
                auto magFwd = std::make_shared<MagneticForward>(
                    initialModel, syntheticMagData, magInc, magDec, magField);
                VectorXd magPred = magFwd->compute(params);
                std::vector<double> mx, my;
                for (size_t i = 0; i < syntheticMagData.size(); ++i) {
                    mx.push_back(syntheticMagData[i].t_obs);
                    my.push_back(magPred[static_cast<Eigen::Index>(i)]);
                }
                plot.panel(1).title = "Magnetics";
                plot.panel(1).xlabel = "Observed (nT)";
                plot.panel(1).ylabel = "Predicted (nT)";
                plot.panel(1).scatters.push_back({mx, my, "red", "", 2.0});
                double mmin = *std::min_element(mx.begin(), mx.end());
                double mmax = *std::max_element(mx.begin(), mx.end());
                plot.panel(1).xmin = mmin; plot.panel(1).xmax = mmax;
                plot.panel(1).ymin = mmin; plot.panel(1).ymax = mmax;
                plot.panel(1).hlines.push_back({0.0, "gray", 0.5, "dashed", ""});
            }

            std::string diagPath = iterDir + "/diagnostics.svg";
            plot.save(diagPath);
        }

        std::cout << std::endl;
    });

    auto result = runner.run();

    // =====================================================================
    // 9. Report
    // =====================================================================
    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "Converged: " << (result.converged ? "yes" : "no")
              << "  Iterations: " << result.totalIterations << std::endl;
    std::cout << "Final RMS: " << std::fixed << std::setprecision(4)
              << result.finalRMS << " mGal" << std::endl;

    if (trueModel) {
        std::cout << "\nProperties (recovered vs true vs cluster-CSV):" << std::endl;
        for (int g = 0; g < initialModel->groupCount(); ++g) {
            const auto& grp = initialModel->group(g);
            const auto& tgrp = trueModel->group(g);
            std::cout << "  " << grp.name << ":"
                      << " rho=" << std::fixed << std::setprecision(4)
                      << result.finalDensities[g]
                      << " (true " << tgrp.density
                      << ", csv " << grp.density << ")"
                      << "  chi=" << std::scientific << std::setprecision(4)
                      << result.finalSusceptibilities[g]
                      << " (true " << tgrp.susceptibility << ")";
            if (g < static_cast<int>(result.finalRemanences.size())) {
                std::cout << "  M_rem=" << std::fixed << std::setprecision(2)
                          << result.finalRemanences[g] << " A/m"
                          << " (true " << tgrp.remanence_magnitude << ")";
            }
            std::cout << std::endl;
        }
    } else {
        std::cout << "\nProperties (recovered vs initial):" << std::endl;
        for (int g = 0; g < initialModel->groupCount(); ++g) {
            const auto& grp = initialModel->group(g);
            std::cout << "  " << grp.name << ":"
                      << " rho=" << std::fixed << std::setprecision(4)
                      << result.finalDensities[g]
                      << " (initial " << grp.density << ")"
                      << "  chi=" << std::scientific << std::setprecision(4)
                      << result.finalSusceptibilities[g]
                      << " (initial " << grp.susceptibility << ")";
            if (g < static_cast<int>(result.finalRemanences.size())) {
                std::cout << "  M_rem=" << std::fixed << std::setprecision(2)
                          << result.finalRemanences[g] << " A/m"
                          << " (initial " << grp.remanence_magnitude << ")";
            }
            std::cout << std::endl;
        }
    }

    // =====================================================================
    // 10. Export inversion results
    // =====================================================================

    InversionExporter exporter(outputDir, "csv_inv");
    exporter.setSubfolder("final");
    if (!groupExportNames.empty()) {
        exporter.setGroupNaming(groupExportNames);
    }

    // Compute export grid bounds from data + model surfaces
    double xmin = 1e30, xmax = -1e30;
    double ymin = 1e30, ymax = -1e30;
    double zmin = 1e30, zmax = -1e30;

    for (const auto& pt : syntheticData) {
        xmin = std::min(xmin, pt.position.x());
        xmax = std::max(xmax, pt.position.x());
        ymin = std::min(ymin, pt.position.y());
        ymax = std::max(ymax, pt.position.y());
        zmin = std::min(zmin, pt.position.z());
        zmax = std::max(zmax, pt.position.z());
    }
    for (int si = 0; si < result.finalModel->groupMeshCount(); ++si) {
        const auto* s = result.finalModel->groupMesh(si);
        for (uint32_t vi = 0; vi < s->vertexCount(); ++vi) {
            const auto& v = s->vertex(vi);
            xmin = std::min(xmin, v.position.x());
            xmax = std::max(xmax, v.position.x());
            ymin = std::min(ymin, v.position.y());
            ymax = std::max(ymax, v.position.y());
            zmin = std::min(zmin, v.position.z());
            zmax = std::max(zmax, v.position.z());
        }
    }

    double exportCellSize = 50.0;
    double pad = exportCellSize;
    xmin = std::floor((xmin - pad) / exportCellSize) * exportCellSize;
    xmax = std::ceil((xmax + pad) / exportCellSize) * exportCellSize;
    ymin = std::floor((ymin - pad) / exportCellSize) * exportCellSize;
    ymax = std::ceil((ymax + pad) / exportCellSize) * exportCellSize;
    zmin = std::floor((zmin - pad) / exportCellSize) * exportCellSize;
    zmax = std::ceil((zmax + pad) / exportCellSize) * exportCellSize;

    std::cout << "[Export] Grid bounds: X[" << xmin << ", " << xmax
              << "] Y[" << ymin << ", " << ymax
              << "] Z[" << zmin << ", " << zmax
              << "] cell=" << exportCellSize << "m" << std::endl;

    exporter.exportAll(result, syntheticData,
                       xmin, xmax, ymin, ymax,
                       zmin, zmax, exportCellSize);

    // =====================================================================
    // 11. SVG convergence plot
    // =====================================================================
    std::cout << "\n--- Generating convergence plot ---" << std::endl;
    {
        const auto& h = result.history;
        int n = static_cast<int>(h.size());

        std::vector<double> iters(n), rms(n), obj(n), dw_gx(n), dw_gy(n),
                            dw_mx(n), dw_my(n);
        for (int i = 0; i < n; ++i) {
            iters[i]  = static_cast<double>(h[i].iteration);
            rms[i]    = h[i].rmsError;
            obj[i]    = h[i].totalObjective;
            dw_gx[i]  = h[i].dw_gravity_x;
            dw_gy[i]  = h[i].dw_gravity_y;
            dw_mx[i]  = h[i].dw_magnetic_x;
            dw_my[i]  = h[i].dw_magnetic_y;
        }

        SVGPlot plot(900, 650, 2, 2);
        plot.title = "Voisey's Bay — CSV-Driven Joint Inversion";

        {
            auto& p = plot.panel(0, 0);
            p.title = "Gravity RMS Misfit";
            p.xlabel = "Iteration";
            p.ylabel = "RMS (mGal)";
            plot.addLine(0, iters.data(), rms.data(), n, "#2166ac", "Gravity RMS", 2.0);
            plot.autoBounds(0);
        }

        {
            auto& p = plot.panel(0, 1);
            p.title = "Objective Function (log)";
            p.xlabel = "Iteration";
            p.ylabel = "Total Objective";
            p.logY = true;
            plot.addLine(1, iters.data(), obj.data(), n, "#b2182b", "Total obj", 2.0);
            plot.autoBounds(1);
        }

        {
            auto& p = plot.panel(1, 0);
            p.title = "Gravity Durbin-Watson";
            p.xlabel = "Iteration";
            p.ylabel = "Durbin-Watson";
            plot.addLine(2, iters.data(), dw_gx.data(), n, "#4daf4a", "DW_g (x)", 1.5);
            plot.addLine(2, iters.data(), dw_gy.data(), n, "#377eb8", "DW_g (y)", 1.5);
            plot.addHLine(2, 2.0, "#888888", 1.2, "dashed", "Random (2.0)");
            plot.autoBounds(2);
        }

        {
            auto& p = plot.panel(1, 1);
            p.title = "Magnetic Durbin-Watson";
            p.xlabel = "Iteration";
            p.ylabel = "Durbin-Watson";
            plot.addLine(3, iters.data(), dw_mx.data(), n, "#a65628", "DW_m (x)", 1.5);
            plot.addLine(3, iters.data(), dw_my.data(), n, "#f781bf", "DW_m (y)", 1.5);
            plot.addHLine(3, 2.0, "#888888", 1.2, "dashed", "Random (2.0)");
            plot.autoBounds(3);
        }

        std::string svgPath = outputDir + "/convergence_plot_csv.svg";
        plot.save(svgPath);
        std::cout << "  Wrote " << svgPath << std::endl;
        std::cout << "  (Open in any browser to view; print to PDF/PNG as needed)"
                  << std::endl;
    }

    std::cout << "\n=== Done ===" << std::endl;
    return 0;
}

