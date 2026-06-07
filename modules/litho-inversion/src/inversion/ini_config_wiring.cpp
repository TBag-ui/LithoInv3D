#include <litho_invert/inversion/runner.h>
#include <litho_invert/io/ini_config.h>

namespace litho_invert {

InversionConfig InversionConfig::fromIni(const IniConfig& ini) {
    InversionConfig config;

    // --- [inversion] ---
    config.maxIterations          = ini.getInt("inversion", "max_iterations", 500);
    config.tolerance              = ini.getDouble("inversion", "tolerance", 1e-6);
    config.lambda                 = ini.getDouble("inversion", "lambda", 1.0);
    config.omega                  = ini.getDouble("inversion", "omega", 1e6);
    config.lbfgsHistory           = ini.getInt("inversion", "lbfgs_history", 10);
    config.fdStep                 = ini.getDouble("inversion", "fd_step", 1.0);
    config.armijoC1               = ini.getDouble("inversion", "armijo_c1", 1e-4);
    config.lineSearchMaxIter      = ini.getInt("inversion", "line_search_max_iter", 50);
    config.disableLineSearch      = ini.getBool("inversion", "disable_line_search", false);
    config.controlPointStride     = ini.getInt("inversion", "control_point_stride", 0);
    config.enablePropertyInversion = ini.getBool("inversion", "enable_property_inversion", false);
    config.enableGeometryInversion = ini.getBool("inversion", "enable_geometry_inversion", true);
    config.propertyInversionInterval = ini.getInt("inversion", "property_inversion_interval", 50);
    config.propertyInversionMaxIter  = ini.getInt("inversion", "property_inversion_max_iter", 20);
    config.propertyDamping          = ini.getDouble("inversion", "property_damping", 0.01);

    // Vertex freedom
    std::string vf = ini.getString("inversion", "vertex_freedom", "z_only");
    if (vf == "xyz_free") {
        config.vertexFreedom = VertexFreedom::XYZ_FREE;
    } else if (vf == "xy_free") {
        config.vertexFreedom = VertexFreedom::XY_FREE;
    } else {
        config.vertexFreedom = VertexFreedom::Z_ONLY;
    }

    // Mesh boundary mode
    config.meshBoundaryMode = ini.getString("inversion", "mesh_boundary_mode", "free");

    // Solver selection
    config.solver               = ini.getString("inversion", "solver", "lbfgsb");
    config.gncgCGMaxIter        = ini.getInt("inversion", "gncg_cg_max_iter", 50);
    config.gncgCGTolerance      = ini.getDouble("inversion", "gncg_cg_tolerance", 1e-6);
    config.enableEigenvalueScaling = ini.getBool("inversion", "enable_eigenvalue_scaling", true);

    // --- [gravity] ---
    config.propertyDensityMin = ini.getDouble("gravity", "density_min", 1.0);
    config.propertyDensityMax = ini.getDouble("gravity", "density_max", 6.0);
    config.gravityUncertainty = ini.getDouble("gravity", "gravity_uncertainty", 0.0);

    // --- [magnetic] ---
    config.magneticWeight          = ini.getDouble("magnetic", "magnetic_weight", 1.0);
    config.magneticInclination     = ini.getDouble("magnetic", "inclination", 75.0);
    config.magneticDeclination     = ini.getDouble("magnetic", "declination", -20.0);
    config.magneticField_nT        = ini.getDouble("magnetic", "field_nT", 55000.0);
    config.magneticUncertainty     = ini.getDouble("magnetic", "magnetic_uncertainty", 0.0);
    config.propertySusceptibilityMin = ini.getDouble("magnetic", "susceptibility_min", 0.0);
    config.propertySusceptibilityMax = ini.getDouble("magnetic", "susceptibility_max", 1.0);
    config.propertyRemanenceMin      = ini.getDouble("magnetic", "remanence_min", 0.0);
    config.propertyRemanenceMax      = ini.getDouble("magnetic", "remanence_max", 10.0);
    config.propertyRemanenceComponentMin = ini.getDouble("magnetic", "remanence_component_min", -10.0);
    config.propertyRemanenceComponentMax = ini.getDouble("magnetic", "remanence_component_max", 10.0);

    // Remanence mode
    std::string remMode = ini.getString("magnetic", "remanence_mode", "effective_susceptibility");
    if (remMode == "fixed_vector") {
        config.remanenceMode = RemanentMagnetizationMode::FixedVectorPerGroup;
    } else if (remMode == "vector") {
        config.remanenceMode = RemanentMagnetizationMode::VectorPerGroup;
    } else {
        config.remanenceMode = RemanentMagnetizationMode::EffectiveSusceptibility;
    }

    // --- [padding] ---
    config.enablePaddingGroup      = ini.getBool("padding", "enable_padding_group", false);
    config.paddingDensityInitial   = ini.getDouble("padding", "padding_density_initial", 2.67);
    config.paddingDensityLower     = ini.getDouble("padding", "padding_density_lower", 1.0);
    config.paddingDensityUpper     = ini.getDouble("padding", "padding_density_upper", 6.0);
    config.paddingDepth            = ini.getDouble("padding", "padding_depth", -100000.0);
    config.paddingSusceptibilityInitial = ini.getDouble("padding", "padding_susceptibility_initial", 0.0);
    config.paddingSusceptibilityLower   = ini.getDouble("padding", "padding_susceptibility_lower", 0.0);
    config.paddingSusceptibilityUpper   = ini.getDouble("padding", "padding_susceptibility_upper", 1.0);
    config.paddingConductivityInitial   = ini.getDouble("padding", "padding_conductivity_initial", 1e-4);
    config.paddingConductivityLower     = ini.getDouble("padding", "padding_conductivity_lower", 1e-5);
    config.paddingConductivityUpper     = ini.getDouble("padding", "padding_conductivity_upper", 1e2);

    // --- [topography] ---
    std::string topoMode = ini.getString("topography", "mode", "none");
    if (topoMode == "raw") {
        config.topography.mode = TopographyMode::Raw;
    } else if (topoMode == "terrain_corrected") {
        config.topography.mode = TopographyMode::TerrainCorrected;
    } else {
        config.topography.mode = TopographyMode::None;
    }
    config.topography.demFile        = ini.getString("topography", "dem_file", "");
    config.topography.datumElevation = ini.getDouble("topography", "datum_elevation", 0.0);
    config.topography.bouguerDensity = ini.getDouble("topography", "bouguer_density", 2.67);
    config.topography.paddingRings   = ini.getInt("topography", "padding_rings", 0);
    config.topography.paddingCellSize = ini.getDouble("topography", "padding_cell_size", 0.0);
    config.topography.invertHalfspaceProperties =
        ini.getBool("topography", "invert_halfspace_properties", false);

    // projected_beyond_survey_area: controls how closure surfaces are built
    // when the model extends beyond the survey footprint via padding rings.
    //   independent_flat          — closure surfaces are independent horizontal grids (default)
    //   projected_beyond_survey   — closure surfaces copy x,y from reference surface topology
    std::string projMode = ini.getString("topography", "projected_beyond_survey_area", "independent_flat");
    config.topography.paddingProjection =
        (projMode == "projected_beyond_survey") ? PaddingProjection::ProjectedBeyondSurvey
                                                : PaddingProjection::IndependentFlat;

    // --- [regularization] ---
    config.enableReferenceModel = ini.getBool("regularization", "enable_reference_model", false);
    config.lambdaRef = ini.getDouble("regularization", "lambda_ref", 0.1);
    config.gradientSmoothingWeight = ini.getDouble("regularization", "gradient_smoothing_weight", 0.0);

    // --- [bounds] ---
    config.enableDepthBounds = ini.getBool("bounds", "enable_depth_bounds", false);
    config.depthBoundMargin = ini.getDouble("bounds", "depth_bound_margin", 500.0);

    // --- [data] ---
    config.failOnMissingMeshes = ini.getBool("data", "fail_on_missing_meshes", true);
    // Note: actual data files are loaded by the caller (main), not here.
    // The INI only stores paths for reference / reproducibility.

    // --- [output] ---
    config.iterationExportDir = ini.getString("output", "iteration_export_dir", "");
    config.lineSearchExportDir = ini.getString("output", "line_search_export_dir", "");

    return config;
}

IniConfig InversionConfig::toIni() const {
    IniConfig ini;

    // --- [inversion] ---
    ini.setInt("inversion", "max_iterations", maxIterations);
    ini.setDouble("inversion", "tolerance", tolerance);
    ini.setDouble("inversion", "lambda", lambda);
    ini.setDouble("inversion", "omega", omega);
    ini.setInt("inversion", "lbfgs_history", lbfgsHistory);
    ini.setDouble("inversion", "fd_step", fdStep);
    ini.setDouble("inversion", "armijo_c1", armijoC1);
    ini.setInt("inversion", "line_search_max_iter", lineSearchMaxIter);
    ini.setBool("inversion", "disable_line_search", disableLineSearch);
    ini.setInt("inversion", "control_point_stride", controlPointStride);
    ini.setBool("inversion", "enable_property_inversion", enablePropertyInversion);
    ini.setBool("inversion", "enable_geometry_inversion", enableGeometryInversion);
    ini.setInt("inversion", "property_inversion_interval", propertyInversionInterval);
    ini.setInt("inversion", "property_inversion_max_iter", propertyInversionMaxIter);
    ini.setDouble("inversion", "property_damping", propertyDamping);

    switch (vertexFreedom) {
    case VertexFreedom::XYZ_FREE:
        ini.setString("inversion", "vertex_freedom", "xyz_free");
        break;
    case VertexFreedom::XY_FREE:
        ini.setString("inversion", "vertex_freedom", "xy_free");
        break;
    default:
        ini.setString("inversion", "vertex_freedom", "z_only");
        break;
    }

    ini.setString("inversion", "solver", solver);
    ini.setString("inversion", "mesh_boundary_mode", meshBoundaryMode);
    ini.setInt("inversion", "gncg_cg_max_iter", gncgCGMaxIter);
    ini.setDouble("inversion", "gncg_cg_tolerance", gncgCGTolerance);
    ini.setBool("inversion", "enable_eigenvalue_scaling", enableEigenvalueScaling);

    // --- [gravity] ---
    ini.setDouble("gravity", "density_min", propertyDensityMin);
    ini.setDouble("gravity", "density_max", propertyDensityMax);
    ini.setDouble("gravity", "gravity_uncertainty", gravityUncertainty);

    // --- [magnetic] ---
    ini.setDouble("magnetic", "magnetic_weight", magneticWeight);
    ini.setDouble("magnetic", "inclination", magneticInclination);
    ini.setDouble("magnetic", "declination", magneticDeclination);
    ini.setDouble("magnetic", "field_nT", magneticField_nT);
    ini.setDouble("magnetic", "magnetic_uncertainty", magneticUncertainty);
    ini.setDouble("magnetic", "susceptibility_min", propertySusceptibilityMin);
    ini.setDouble("magnetic", "susceptibility_max", propertySusceptibilityMax);
    ini.setDouble("magnetic", "remanence_min", propertyRemanenceMin);
    ini.setDouble("magnetic", "remanence_max", propertyRemanenceMax);
    ini.setDouble("magnetic", "remanence_component_min", propertyRemanenceComponentMin);
    ini.setDouble("magnetic", "remanence_component_max", propertyRemanenceComponentMax);

    switch (remanenceMode) {
    case RemanentMagnetizationMode::FixedVectorPerGroup:
        ini.setString("magnetic", "remanence_mode", "fixed_vector");
        break;
    case RemanentMagnetizationMode::VectorPerGroup:
        ini.setString("magnetic", "remanence_mode", "vector");
        break;
    default:
        ini.setString("magnetic", "remanence_mode", "effective_susceptibility");
        break;
    }

    // --- [padding] ---
    ini.setBool("padding", "enable_padding_group", enablePaddingGroup);
    ini.setDouble("padding", "padding_density_initial", paddingDensityInitial);
    ini.setDouble("padding", "padding_density_lower", paddingDensityLower);
    ini.setDouble("padding", "padding_density_upper", paddingDensityUpper);
    ini.setDouble("padding", "padding_depth", paddingDepth);
    ini.setDouble("padding", "padding_susceptibility_initial", paddingSusceptibilityInitial);
    ini.setDouble("padding", "padding_susceptibility_lower", paddingSusceptibilityLower);
    ini.setDouble("padding", "padding_susceptibility_upper", paddingSusceptibilityUpper);
    ini.setDouble("padding", "padding_conductivity_initial", paddingConductivityInitial);
    ini.setDouble("padding", "padding_conductivity_lower", paddingConductivityLower);
    ini.setDouble("padding", "padding_conductivity_upper", paddingConductivityUpper);

    // --- [topography] ---
    switch (topography.mode) {
    case TopographyMode::Raw:
        ini.setString("topography", "mode", "raw");
        break;
    case TopographyMode::TerrainCorrected:
        ini.setString("topography", "mode", "terrain_corrected");
        break;
    default:
        ini.setString("topography", "mode", "none");
        break;
    }
    ini.setString("topography", "dem_file", topography.demFile);
    ini.setDouble("topography", "datum_elevation", topography.datumElevation);
    ini.setDouble("topography", "bouguer_density", topography.bouguerDensity);
    ini.setInt("topography", "padding_rings", topography.paddingRings);
    ini.setDouble("topography", "padding_cell_size", topography.paddingCellSize);
    ini.setBool("topography", "invert_halfspace_properties", topography.invertHalfspaceProperties);
    ini.setString("topography", "projected_beyond_survey_area",
                  topography.paddingProjection == PaddingProjection::ProjectedBeyondSurvey
                      ? "projected_beyond_survey" : "independent_flat");

    // --- [regularization] ---
    ini.setBool("regularization", "enable_reference_model", enableReferenceModel);
    ini.setDouble("regularization", "lambda_ref", lambdaRef);
    ini.setDouble("regularization", "gradient_smoothing_weight", gradientSmoothingWeight);

    // --- [bounds] ---
    ini.setBool("bounds", "enable_depth_bounds", enableDepthBounds);
    ini.setDouble("bounds", "depth_bound_margin", depthBoundMargin);

    // --- [data] ---
    ini.setBool("data", "fail_on_missing_meshes", failOnMissingMeshes);

    // --- [output] ---
    if (!iterationExportDir.empty()) {
        ini.setString("output", "iteration_export_dir", iterationExportDir);
    }
    if (!lineSearchExportDir.empty()) {
        ini.setString("output", "line_search_export_dir", lineSearchExportDir);
    }

    return ini;
}

} // namespace litho_invert

