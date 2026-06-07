#include <litho_invert/inversion/runner.h>
#include <litho_invert/io/io.h>
#include <litho_invert/em/em_active_forward.h>
#include <litho_invert/em/em_mt_forward.h>
#include "generate_synthetic.h"
#include <iostream>
#include <iomanip>
#include <memory>
#include <filesystem>

using namespace litho_invert;

int main(int argc, char* argv[]) {
    std::string outputDir = ".";
    if (argc > 1) {
        outputDir = argv[1];
    }

    // Create output directory if it doesn't exist
    std::filesystem::create_directories(outputDir);

    std::cout << "=== Voisey's Bay Synthetic Inversion Example ===" << std::endl;
    std::cout << "Output directory: " << outputDir << std::endl;

    // 1. Generate or load true model and synthetic data
    std::cout << "\n--- Generating synthetic model ---" << std::endl;
    auto trueModel = generateTrueModel();
    auto obsPoints = generateObservationPoints();

    // Set conductivities on the true model.
    // Anorthosite:      1e-4 S/m (resistive host rock)
    // Troctolite:       1e-2 S/m (moderately conductive mafic rock)
    // Massive Sulfide:  1.0  S/m (highly conductive — pyrrhotite/chalcopyrite)
    trueModel->group(0).conductivity = 1e-4;
    trueModel->group(1).conductivity = 1e-2;
    trueModel->group(2).conductivity = 1.0;
    std::cout << "True model conductivities (S/m):" << std::endl;
    for (int g = 0; g < trueModel->groupCount(); ++g) {
        std::cout << "  " << trueModel->group(g).name << ": "
                  << std::scientific << std::setprecision(2)
                  << trueModel->group(g).conductivity << std::endl;
    }
    // Use same padding density as initial model so observed/predicted are consistent.
    // The padding contribution will be ~2.17 mGal in both synthetic data and predictions.
    // Property inversion can adjust the padding density to test the mechanism.
    auto syntheticData = computeSyntheticData(trueModel, obsPoints, 2.67);

    // Generate synthetic magnetic data
    std::cout << "\n--- Generating synthetic magnetic data ---" << std::endl;
    const double magInc = 75.0;
    const double magDec = -20.0;
    const double magField = 55000.0;
    auto syntheticMagData = computeSyntheticMagnetic(trueModel, obsPoints,
                                                      magInc, magDec, magField, 0.0);
    {
        double tmin = 1e99, tmax = -1e99;
        for (const auto& pt : syntheticMagData) {
            if (pt.t_obs < tmin) tmin = pt.t_obs;
            if (pt.t_obs > tmax) tmax = pt.t_obs;
        }
        std::cout << "Magnetic points: " << syntheticMagData.size() << std::endl;
        std::cout << "Magnetic range: " << tmin << " to " << tmax << " nT" << std::endl;
    }

    // Report true model statistics
    {
        std::cout << "Observation points: " << syntheticData.size() << std::endl;
        double gmin = 1e99, gmax = -1e99;
        for (const auto& pt : syntheticData) {
            if (pt.g_obs < gmin) gmin = pt.g_obs;
            if (pt.g_obs > gmax) gmax = pt.g_obs;
        }
        std::cout << "Gravity range: " << gmin << " to " << gmax << " mGal" << std::endl;

        // True surface statistics
        for (int s = 0; s < trueModel->surfaceCount(); ++s) {
            const SurfaceMesh* surf = trueModel->surface(s);
            double zmin = 1e99, zmax = -1e99;
            for (uint32_t v = 0; v < surf->vertexCount(); ++v) {
                double z = surf->vertex(v).position.z();
                if (z < zmin) zmin = z;
                if (z > zmax) zmax = z;
            }
            std::cout << "True surface " << s << " z range: " << zmin << " to " << zmax << std::endl;
        }
    }

    // =====================================================================
    // EM forward-model test: run active EM and MT forward models standalone
    // against the TRUE model to verify they produce reasonable output.
    // =====================================================================
    std::cout << "\n--- EM forward-model diagnostic (true model) ---" << std::endl;

    // -- Active EM (airborne TEM) --
    std::vector<EMSource> emSources;
    std::vector<EMReceiver> emReceivers;
    EMConfig emConfig;
    buildAirborneTEMGeometry(emSources, emReceivers, emConfig);
    std::cout << "Active EM: " << emSources.size() << " source(s), "
              << emReceivers.size() << " receiver(s), "
              << emConfig.timeGates_s.size() << " time gates" << std::endl;

    auto syntheticEM = computeSyntheticActiveEM(
        trueModel, obsPoints, emSources, emReceivers, emConfig, 1e-4);

    {
        double emin = 1e99, emax = -1e99;
        for (const auto& pt : syntheticEM) {
            if (pt.d_obs < emin) emin = pt.d_obs;
            if (pt.d_obs > emax) emax = pt.d_obs;
        }
        std::cout << "Active EM data points: " << syntheticEM.size() << std::endl;
        std::cout << "Active EM response range: " << std::scientific << std::setprecision(4)
                  << emin << " to " << emax << std::endl;

        // Show a few representative responses at the center point
        std::cout << "Center-point response by gate:" << std::endl;
        int nGates = static_cast<int>(emConfig.timeGates_s.size());
        for (size_t i = 0; i < syntheticEM.size(); ++i) {
            const auto& pt = syntheticEM[i];
            // Center point is at (0,0,0) — grid index 60 (11x11 grid, (5,5))
            if (std::abs(pt.position.x()) < 1.0 && std::abs(pt.position.y()) < 1.0) {
                std::cout << "  gate " << pt.gateIndex
                          << " (t=" << std::scientific << std::setprecision(1)
                          << pt.timeGate_s << " s): "
                          << std::setprecision(4) << pt.d_obs << std::endl;
            }
        }
    }

    // -- MT --
    std::vector<MTStation> mtStations;
    std::vector<double> mtFrequencies;
    buildMTStations(mtStations, mtFrequencies);
    std::cout << "\nMT: " << mtStations.size() << " stations, "
              << mtFrequencies.size() << " frequencies" << std::endl;

    EMConfig mtConfig;
    mtConfig.solverMethod = "ie";
    auto syntheticMT = computeSyntheticMT(trueModel, mtStations, mtConfig, 1e-4);

    {
        double zmin = 1e99, zmax = -1e99;
        for (const auto& pt : syntheticMT) {
            if (pt.zReal_obs < zmin) zmin = pt.zReal_obs;
            if (pt.zReal_obs > zmax) zmax = pt.zReal_obs;
        }
        std::cout << "MT data points: " << syntheticMT.size() << std::endl;
        std::cout << "MT Zxy real range: " << std::fixed << std::setprecision(2)
                  << zmin << " to " << zmax << " ohm" << std::endl;

        // Show Zxy at station MT01_SW across frequencies
        std::cout << "Station MT01_SW Zxy vs frequency:" << std::endl;
        for (const auto& pt : syntheticMT) {
            if (pt.stationIndex == 0) {
                double f = mtStations[0].frequencies_Hz[pt.frequencyIndex];
                std::cout << "  " << std::scientific << std::setprecision(4) << f
                          << " Hz: Re(Zxy) = " << std::fixed << std::setprecision(3)
                          << pt.zReal_obs << " ohm" << std::endl;
            }
        }
    }

    // -- Sanity check: EM forward model against initial flat model --
    // This verifies the EM forward works with the initial model (which will
    // be used during inversion if EM is enabled).
    {
        auto initialEM = generateInitialModel();
        initialEM->group(0).conductivity = 1e-4;
        initialEM->group(1).conductivity = 1e-2;
        initialEM->group(2).conductivity = 1.0;
        for (int i = 0; i < initialEM->surfaceCount(); ++i) {
            initialEM->surface(i)->buildNeighbors();
        }

        auto initEMData = computeSyntheticActiveEM(
            initialEM, obsPoints, emSources, emReceivers, emConfig, 1e-4);

        double emin = 1e99, emax = -1e99;
        for (const auto& pt : initEMData) {
            if (pt.d_obs < emin) emin = pt.d_obs;
            if (pt.d_obs > emax) emax = pt.d_obs;
        }
        std::cout << "\nInitial flat-model Active EM range: "
                  << std::scientific << std::setprecision(4)
                  << emin << " to " << emax << std::endl;
    }

    // 2. Build initial model with flat surfaces
    std::cout << "\n--- Building initial model ---" << std::endl;
    auto initialModel = generateInitialModel();

    // Set conductivities on the initial model
    initialModel->group(0).conductivity = 1e-4;
    initialModel->group(1).conductivity = 1e-2;
    initialModel->group(2).conductivity = 1.0;

    // Build neighbors on the initial surfaces (for regularization)
    for (int i = 0; i < initialModel->surfaceCount(); ++i) {
        initialModel->surface(i)->buildNeighbors();
    }

    // Diagnostic: compute gravity misfit of initial model vs true model
    {
        GravityForward fwdTrue(trueModel, obsPoints);
        fwdTrue.enablePadding(true, -100000.0);
        fwdTrue.setPaddingDensity(2.67);
        VectorXd trueParams = trueModel->assembleParameterVector();
        VectorXd trueGrav = fwdTrue.compute(trueParams);

        GravityForward fwdInit(initialModel, obsPoints);
        fwdInit.enablePadding(true, -100000.0);
        fwdInit.setPaddingDensity(2.67);
        VectorXd initParams = initialModel->assembleParameterVector();
        VectorXd initGrav = fwdInit.compute(initParams);

        double trueRMS = 0.0, initRMS = 0.0;
        for (size_t i = 0; i < syntheticData.size(); ++i) {
            double rTrue = syntheticData[i].g_obs - trueGrav[static_cast<Index>(i)];
            double rInit = syntheticData[i].g_obs - initGrav[static_cast<Index>(i)];
            trueRMS += rTrue * rTrue;
            initRMS += rInit * rInit;
        }
        trueRMS = std::sqrt(trueRMS / syntheticData.size());
        initRMS = std::sqrt(initRMS / syntheticData.size());
        std::cout << "True model self-misfit RMS: " << trueRMS << " mGal"
                  << " (should be ~0)" << std::endl;
        std::cout << "Initial flat model misfit RMS: " << initRMS << " mGal"
                  << " (max possible improvement)" << std::endl;
        std::cout << "Initial model param count: " << initParams.size() << std::endl;
    }

    // 3. Configure and run inversion
    std::cout << "\n--- Running inversion ---" << std::endl;

    InversionConfig config;
    config.model = initialModel;
    config.observedData = syntheticData;
    config.lambda = 0.0;            // no regularization for maximum shape recovery
    config.maxIterations = 200;
    config.tolerance = 1e-5;
    config.lbfgsHistory = 10;

    // Downsample to 6x6=36 CPs per surface (72 total vs 121 data = overdetermined)
    config.controlPointStride = 2;

    // Property inversion: currently disabled — geometry-only to avoid
    // density adjustments short-circuiting surface shape recovery.
    config.enablePropertyInversion = false;
    config.propertyInversionInterval = 50;
    config.propertyInversionMaxIter = 15;
    config.propertyDensityMin = 1.5;
    config.propertyDensityMax = 6.0;

    // Padding group: deep half-space with free density to absorb regional trends
    config.enablePaddingGroup = true;
    config.paddingDensityInitial = 2.67;
    config.paddingDensityLower = 1.5;
    config.paddingDensityUpper = 6.0;
    config.paddingDepth = -100000.0;

    // Magnetic joint inversion
    config.magneticData = syntheticMagData;
    config.magneticWeight = 0.01;   // balanced for t_std=50 nT uncertainties
    config.magneticInclination = magInc;
    config.magneticDeclination = magDec;
    config.magneticField_nT = magField;
    config.propertySusceptibilityMin = 0.0;
    config.propertySusceptibilityMax = 1.0;

    // Active EM joint inversion — currently disabled for debugging crash.
    // The standalone EM forward-model diagnostic (above) confirms the
    // forward models work.  Re-enable once the joint gradient is fixed.
    bool enableEMInversion = true;  // re-enabled after parameterCount fix
    if (enableEMInversion) {
        config.activeEMData = syntheticEM;
        config.emSources = emSources;
        config.emReceivers = emReceivers;
        config.emConfig = emConfig;
        config.activeEMWeight = 1e-6;
        config.propertyConductivityMin = 1e-5;
        config.propertyConductivityMax = 1e2;

        config.mtData = syntheticMT;
        config.mtStations = mtStations;
        config.mtWeight = 1e-6;

        std::cout << "\nEM inversion: activeEM points=" << syntheticEM.size()
                  << " weight=" << std::scientific << config.activeEMWeight
                  << "  MT points=" << syntheticMT.size()
                  << " weight=" << config.mtWeight << std::endl;
    }

    InversionRunner runner(config);

    // Set iteration callback to print progress
    runner.setIterationCallback([](const InversionIteration& iter) {
        if (iter.iteration % 10 == 0 || iter.iteration == 1) {
            std::cout << "  Iter " << std::setw(3) << iter.iteration
                      << ": RMS=" << std::fixed << std::setprecision(4)
                      << iter.rmsError
                      << " DW_g=(" << std::fixed << std::setprecision(2)
                      << iter.dw_gravity_x << "," << iter.dw_gravity_y << ")";
            if (iter.magneticMisfit > 0.0) {
                std::cout << " mag=" << std::scientific << std::setprecision(4)
                          << iter.magneticMisfit;
            }
            if (iter.dw_magnetic_x != 2.0 || iter.dw_magnetic_y != 2.0) {
                std::cout << " DW_m=(" << std::fixed << std::setprecision(2)
                          << iter.dw_magnetic_x << "," << iter.dw_magnetic_y << ")";
            }
            if (iter.activeEMMisfit > 0.0) {
                std::cout << " aem=" << std::scientific << std::setprecision(4)
                          << iter.activeEMMisfit;
                std::cout << " DW_aem=(" << std::fixed << std::setprecision(2)
                          << iter.dw_activeEM_x << "," << iter.dw_activeEM_y << ")";
            }
            if (iter.mtMisfit > 0.0) {
                std::cout << " mt=" << std::scientific << std::setprecision(4)
                          << iter.mtMisfit;
                std::cout << " DW_mt=(" << std::fixed << std::setprecision(2)
                          << iter.dw_mt_x << "," << iter.dw_mt_y << ")";
            }
            std::cout << std::endl;
        }
    });

    auto result = runner.run();

    // 4. Export true model surfaces for comparison
    std::cout << "\n--- Exporting true surfaces ---" << std::endl;
    {
        InversionExporter trueExp(outputDir, "voiseys_bay_true");
        for (int i = 0; i < trueModel->surfaceCount(); ++i) {
            const SurfaceMesh* surf = trueModel->surface(i);
            std::string surfName = surf->name().empty()
                ? ("surface_" + std::to_string(i))
                : surf->name();
            trueExp.exportTS(*surf, surfName);
        }
    }

    // 5. Report inversion results
    std::cout << "\n--- Inversion Results ---" << std::endl;
    std::cout << "Converged: " << (result.converged ? "yes" : "no") << std::endl;
    std::cout << "Iterations: " << result.totalIterations << std::endl;
    std::cout << "Final objective: " << result.finalMisfit << std::endl;
    std::cout << "Final RMS (mGal): " << result.finalRMS << std::endl;

    // Recovered surface statistics
    std::cout << "\nRecovered surfaces:" << std::endl;
    for (int s = 0; s < initialModel->surfaceCount(); ++s) {
        const SurfaceMesh* surf = initialModel->surface(s);
        double zmin = 1e99, zmax = -1e99;
        for (uint32_t v = 0; v < surf->vertexCount(); ++v) {
            double z = surf->vertex(v).position.z();
            if (z < zmin) zmin = z;
            if (z > zmax) zmax = z;
        }
        std::cout << "  Surface " << s << " (" << surf->name()
                  << "): z range " << zmin << " to " << zmax
                  << ", vertices: " << surf->vertexCount()
                  << ", triangles: " << surf->triangleCount() << std::endl;
    }

    // Report recovered densities
    std::cout << "\nRecovered densities (g/cm^3):" << std::endl;
    for (int g = 0; g < initialModel->groupCount(); ++g) {
        std::cout << "  " << initialModel->group(g).name << ": "
                  << std::fixed << std::setprecision(4)
                  << result.finalDensities[g] << std::endl;
    }
    if (config.enablePaddingGroup) {
        std::cout << "  Padding (deep background): "
                  << std::fixed << std::setprecision(4)
                  << result.finalPaddingDensity << std::endl;
    }

    // Report recovered susceptibilities
    std::cout << "\nRecovered susceptibilities (SI):" << std::endl;
    for (int g = 0; g < initialModel->groupCount(); ++g) {
        std::cout << "  " << initialModel->group(g).name << ": "
                  << std::fixed << std::setprecision(4)
                  << result.finalSusceptibilities[g] << std::endl;
    }
    if (config.enablePaddingGroup) {
        std::cout << "  Padding (deep background): "
                  << std::fixed << std::setprecision(4)
                  << result.finalPaddingSusceptibility << std::endl;
    }

    // Report recovered conductivities
    std::cout << "\nRecovered conductivities (S/m):" << std::endl;
    for (int g = 0; g < initialModel->groupCount(); ++g) {
        std::cout << "  " << initialModel->group(g).name << ": "
                  << std::scientific << std::setprecision(4)
                  << result.finalConductivities[g] << std::endl;
    }
    if (config.enablePaddingGroup && result.finalPaddingConductivity > 0.0) {
        std::cout << "  Padding (deep background): "
                  << std::scientific << std::setprecision(4)
                  << result.finalPaddingConductivity << std::endl;
    }

    // 6. Export inversion results
    std::cout << "\n--- Exporting results ---" << std::endl;
    InversionExporter exporter(outputDir, "voiseys_bay");
    exporter.exportAll(result, syntheticData,
                       -250.0, 250.0,  // x range
                       -250.0, 250.0,  // y range
                       -600.0, 0.0,    // z range
                       50.0);          // cell size

    std::cout << "\n=== Done ===" << std::endl;
    return 0;
}
