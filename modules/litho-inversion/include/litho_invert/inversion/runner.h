#pragma once
#include <litho_invert/core/common.h>
#include <litho_invert/surface/surface_mesh.h>
#include <litho_invert/litho/lithology_model.h>
#include <litho_invert/forward/gravity_forward.h>
#include <litho_invert/forward/magnetic_forward.h>
#include <litho_invert/em/em_data.h>
#include <litho_invert/em/em_active_forward.h>
#include <litho_invert/em/em_mt_forward.h>
#include <litho_invert/inversion/optimizer.h>
#include <litho_invert/inversion/lbfgsb_optimizer.h>
#include <litho_invert/inversion/objective.h>
#include <litho_invert/inversion/joint_objective.h>
#include <litho_invert/regularization/smoothness.h>
#include <litho_invert/regularization/reference_model.h>
#include <litho_invert/inversion/constraint_handler.h>
#include <litho_invert/inversion/inversion_result.h>
#include <litho_invert/io/exporters.h>
#include <memory>
#include <vector>
#include <functional>

namespace litho_invert {

class IniConfig;  // forward declaration

struct InversionConfig {
    std::shared_ptr<LithologyModel> model;
    GravityData observedData;
    std::vector<Constraint> constraints;
    double lambda = 1.0;         // regularization weight
    double omega = 1e6;          // constraint penalty weight
    int maxIterations = 500;
    double tolerance = 1e-6;
    int lbfgsHistory = 10;

    // Finite-difference step size for Jacobian computation (metres).
    // Larger values help with noisy geophysical gradients but reduce accuracy.
    double fdStep = 1.0;

    // Armijo constant for line search (c1 in sufficient-decrease condition).
    // Default 1e-4. Increase to 1e-3 or 1e-2 for noisy objective functions.
    double armijoC1 = 1e-4;

    // Max backtracking iterations in line search. Default 50.
    int lineSearchMaxIter = 50;

    // Solver selection
    // "lbfgsb" — Limited-memory BFGS with bounds (default)
    // "gncg"   — Gauss-Newton Conjugate Gradient (SimPEG-style)
    std::string solver = "lbfgsb";

    // GNCG-specific settings
    int gncgCGMaxIter = 50;         // max inner CG iterations
    double gncgCGTolerance = 1e-6;   // CG convergence tolerance

    // Control-point downsampling: when > 0, only every N-th grid vertex gets
    // a free DOF. Non-control vertices are bilinearly interpolated. This
    // reduces the parameter count for faster Jacobian computation.
    int controlPointStride = 0;  // 0 = use all vertices (no downsampling)

    // Property inversion (alternating geometry/density optimization)
    bool enablePropertyInversion = false;
    bool enableGeometryInversion = true;   // false = property-only mode
    int  propertyInversionInterval = 50;    // geometry iters between property steps
    int  propertyInversionMaxIter = 20;     // max iterations per property phase
    double propertyDamping = 0.01;          // Tikhonov damping for property inversion
    double propertyDensityMin = 1.0;        // lower bound for group densities (g/cm^3)
    double propertyDensityMax = 6.0;        // upper bound for group densities (g/cm^3)

    // Padding group (deep half-space with free density)
    bool   enablePaddingGroup = false;
    double paddingDensityInitial = 2.67;    // starting density (g/cm^3)
    double paddingDensityLower = 1.0;
    double paddingDensityUpper = 6.0;
    double paddingDepth = -100000.0;        // deep bottom z (positive-up)

    // Magnetic data (empty = gravity-only, backward compatible)
    MagneticData magneticData;
    double magneticWeight = 1.0;        // α — weight of magnetic misfit vs gravity
    double magneticInclination = 75.0;  // deg from horizontal, positive down
    double magneticDeclination = -20.0; // deg from north, positive east
    double magneticField_nT = 55000.0;  // inducing field strength

    // Per-data-type default uncertainties (applied when data point sigma=0).
    // SimPEG uses 0.015 mGal for gravity and 60 nT for airborne TMI.
    double gravityUncertainty = 0.0;     // mGal, 0 = unweighted
    double magneticUncertainty = 0.0;    // nT, 0 = unweighted

    // Eigenvalue-based Hessian scaling toggle (SimPEG ScalingMultipleDataMisfits_ByEig)
    bool enableEigenvalueScaling = true;

    // Magnetic property inversion bounds
    double propertySusceptibilityMin = 0.0;
    double propertySusceptibilityMax = 1.0;

    // Remanent magnetization mode (applies to all groups)
    RemanentMagnetizationMode remanenceMode =
        RemanentMagnetizationMode::EffectiveSusceptibility;

    // Remanence bounds — FixedVectorPerGroup mode (A/m)
    double propertyRemanenceMin = 0.0;
    double propertyRemanenceMax = 10.0;   // 10 A/m ≈ substantial remnant

    // Remanence bounds — VectorPerGroup mode (A/m per component)
    double propertyRemanenceComponentMin = -10.0;
    double propertyRemanenceComponentMax = 10.0;

    // Magnetic padding
    double paddingSusceptibilityInitial = 0.0;
    double paddingSusceptibilityLower = 0.0;
    double paddingSusceptibilityUpper = 1.0;

    // Reference model regularization: penalizes deviation from starting model
    // Use with reduced lambda to prevent depth drift while preserving features
    bool enableReferenceModel = false;
    double lambdaRef = 0.1;
    VectorXd referenceParams;  // when non-empty, use instead of capturing from model

    // When true, L-BFGS-B always accepts α=1.0 without backtracking trial evaluations.
    bool disableLineSearch = false;

    // Hard depth bounds: constrains DOFs to ±margin from starting position
    bool enableDepthBounds = false;
    double depthBoundMargin = 500.0;

    // Spatial gradient smoothing weight (0 = disabled, 0.3 typical).
    // Applies Laplacian smoothing to the joint gradient so the L-BFGS-B
    // search direction moves neighboring vertices coherently.
    double gradientSmoothingWeight = 0.0;

    // Surface vertex freedom (Z_ONLY, XYZ_FREE, XY_FREE)
    VertexFreedom vertexFreedom = VertexFreedom::Z_ONLY;

    // Mesh boundary mode for loaded contact meshes:
    //   "free"  — all vertices get vertexFreedom (default)
    //   "fixed" — boundary vertices (hull edges) get FIXED freedom
    std::string meshBoundaryMode = "free";

    // When true, missing contact meshes cause a fatal error instead of
    // silently generating flat fallback surfaces. Default true: the
    // starting model must be explicitly provided via contact meshes.
    bool failOnMissingMeshes = true;

    // Per-iteration .ts export: when non-empty, each geometry iteration
    // exports all group meshes as GOCAD TSurf files into this directory,
    // overwriting the previous iteration's files. Empty = disabled.
    std::string iterationExportDir;

    // Line-search diagnostic export: when non-empty, each trial evaluation
    // during L-BFGS-B line search exports all group meshes as closed-volume
    // .ts files into subdirectories (line_search/iter_NNN/trial_MM/).
    // Separate from iterationExportDir so trial geometry is clearly
    // distinguishable from validated iteration geometry.
    std::string lineSearchExportDir;

    // Per-group export names for iteration .ts files (e.g. "cluster_id_0").
    // When empty, uses "group_0", "group_1", etc.
    std::vector<std::string> groupExportNames;

    // Topography + lateral padding
    TopographyConfig topography;

    // =====================================================================
    // Active EM (airborne + large-loop)
    // Leave ActiveEMData empty for gravity/mag-only inversions.
    // =====================================================================
    ActiveEMData activeEMData;          // observation points (tagged with src/recv/gate)
    std::vector<EMSource> emSources;    // transmitter geometries
    std::vector<EMReceiver> emReceivers; // receiver geometries
    double activeEMWeight = 1.0;        // α_aem — weight vs gravity misfit
    EMConfig emConfig;                  // solver, subsetting, physics config

    // Active EM property inversion bounds
    double propertyConductivityMin = 1e-5;   // S/m
    double propertyConductivityMax = 1e2;    // S/m
    double paddingConductivityInitial = 1e-4; // S/m
    double paddingConductivityLower = 1e-5;
    double paddingConductivityUpper = 1e2;

    // =====================================================================
    // Magnetotelluric (passive, plane-wave)
    // Leave MTData empty for inversions without MT.
    // =====================================================================
    MTData mtData;
    std::vector<MTStation> mtStations;  // station descriptions
    double mtWeight = 1.0;              // α_mt — weight vs gravity misfit

    // --- Factory methods for INI config ---
    static InversionConfig fromIni(const IniConfig& ini);
    IniConfig toIni() const;
};

class InversionRunner {
public:
    explicit InversionRunner(const InversionConfig& config);

    InversionResult run();

    using IterationCallback = std::function<void(const InversionIteration&)>;
    void setIterationCallback(IterationCallback cb);

    std::shared_ptr<Optimizer> geometryOptimizer() { return m_geometryOptimizer; }

private:
    InversionConfig m_config;
    IterationCallback m_callback;
    std::shared_ptr<GravityForward> m_forward;
    std::shared_ptr<ObjectiveFunction> m_objective;
    std::shared_ptr<Optimizer> m_geometryOptimizer;
    std::shared_ptr<SurfaceSmoothness> m_smoothness;
    std::shared_ptr<ReferenceModelRegularization> m_referenceReg;
    std::shared_ptr<ConstraintHandler> m_constraintHandler;
    std::shared_ptr<MagneticForward> m_magneticForward;
    std::shared_ptr<JointObjective> m_jointObjective;
    std::shared_ptr<EMActiveForward> m_activeEMForward;
    std::shared_ptr<EMMTForward> m_mtForward;

    std::unique_ptr<InversionExporter> m_iterationExporter;
    int m_iterationCounter = 0;
    int m_lineSearchTrialCounter = 0;

    void setup();
    void runPropertyInversion(InversionResult& result);

    // Apply spatial Laplacian smoothing to the gradient so that neighboring
    // vertices have correlated search directions, preventing surface tearing
    // during L-BFGS-B line search with large step sizes.
    VectorXd smoothGradient(const VectorXd& g) const;
};

} // namespace litho_invert

