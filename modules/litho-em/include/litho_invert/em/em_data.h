#pragma once
#include <litho_invert/core/common.h>
#include <vector>
#include <string>

namespace litho_invert {

// =========================================================================
// EM residual norm — how misfit is computed per datum.
// Linear:     residual = d_obs - d_pred
// Log10:      residual = log10(d_obs) - log10(d_pred)
//   Log-space handles the wide dynamic range of EM decays (10^5 → 10^-2 nT/s
//   for time-domain, or pT–µT for frequency-domain).  Use this when the
//   amplitude span exceeds 2–3 decades.
// =========================================================================
enum class EMResidualNorm {
    Linear,  // residual = d_obs - d_pred        (unit: same as data)
    Log10    // residual = log10(d_obs) - log10(d_pred)  (dimensionless)
};

// =========================================================================
// EM survey geometry types.
// Determines how the source field is computed and what subsetting strategy
// the trust-region algorithm applies.
// =========================================================================
enum class EMSurveyType {
    AirborneFEM,   // helicopter-borne frequency-domain EM (e.g. RESOLVE, DIGHEM)
    AirborneTEM,   // helicopter or fixed-wing time-domain EM (e.g. VTEM, SkyTEM)
    LargeLoopTEM,  // ground-based large-loop time-domain EM (e.g. UTEM, FLYTEM)
    SparseMT       // magnetotelluric (passive, plane-wave natural source)
};

// =========================================================================
// EM transmitter / source description.
// For MT surveys this is unused (plane-wave source has no position).
// =========================================================================
struct EMSource {
    Vector3d position;          // (x, y, z) — transmitter center position
    EMSurveyType type = EMSurveyType::AirborneTEM;

    // Waveform parameters (time-domain only, ignored for frequency-domain)
    std::string waveform = "step";  // "step", "half-sine", "trapezoid"
    double pulseWidth_s = 0.004;    // pulse width or half-sine duration (seconds)
    double rampTime_s = 1e-5;       // turn-off ramp time (seconds)

    // Loop geometry (airborne bird or ground loop)
    double loopRadius_m = 0.0;      // 0 = treat as vertical magnetic dipole
    double loopHeight_m = 30.0;     // height above ground (airborne) or 0 (ground)
    Vector3d loopNormal = {0, 0, 1}; // loop axis direction (usually vertical)

    // For large-loop: polygon vertices defining the loop path
    std::vector<Vector3d> loopPath;

    EMSource() = default;
};

// =========================================================================
// EM receiver description.
// Each receiver measures the secondary field (dB/dt or B) at a position
// relative to the source. For airborne systems receivers are in the bird;
// for large-loop they may be on the ground.
// =========================================================================
struct EMReceiver {
    Vector3d position;  // (x, y, z) — receiver position
    std::string component = "z";  // "x", "y", "z" — measured component

    EMReceiver() = default;
};

// =========================================================================
// Active EM observation point (airborne or large-loop).
// Tagged with source index, receiver index, and time gate or frequency index
// so the forward model knows which transmitter-receiver geometry to use.
// =========================================================================
struct ActiveEMPoint {
    Vector3d position;          // observation location (for bookkeeping / DW)
    int sourceIndex = 0;        // index into the EMSource array
    int receiverIndex = 0;      // index into the EMReceiver array for this source
    int gateIndex = 0;          // time-gate or frequency index

    double d_obs = 0.0;         // observed data value
    double d_std = 0.0;         // standard deviation (0 = unweighted)

    bool isFrequencyDomain = false;  // true = Hz domain, false = time domain
    double frequency_Hz = 0.0;       // frequency (if frequency-domain)
    double timeGate_s = 0.0;         // gate center time after turn-off (if time-domain)

    EMResidualNorm residualNorm = EMResidualNorm::Linear;

    ActiveEMPoint() = default;
};

using ActiveEMData = std::vector<ActiveEMPoint>;

// =========================================================================
// MT impedance tensor element.
// Each MT station records a 2×2 complex impedance tensor at each frequency:
//   [ Zxx Zxy ]
//   [ Zyx Zyy ]
// where Z_ij = E_i / H_j (units: V/A = Ω).
//
// The forward model predicts Z_ij from the subsurface conductivity
// distribution under plane-wave excitation.
// =========================================================================
struct MTImpedanceElement {
    int stationIndex = 0;       // which MT station
    int frequencyIndex = 0;     // which frequency band
    int iComp = 0;              // E-field component: 0=x, 1=y
    int jComp = 0;              // H-field component: 0=x, 1=y

    double zReal_obs = 0.0;     // observed Re(Z_ij)  [Ω]
    double zImag_obs = 0.0;     // observed Im(Z_ij)  [Ω]
    double z_std = 0.0;         // uncertainty (0 = unweighted)

    MTImpedanceElement() = default;
};

// =========================================================================
// MT station description.
// Sparse stations at accessible locations (ridges, clearings — not marshes).
// =========================================================================
struct MTStation {
    Vector3d position;          // (x, y, z) — station location
    std::string name;           // label (e.g. "MT01")
    std::vector<double> frequencies_Hz;  // frequencies measured at this station

    MTStation() = default;
};

using MTData = std::vector<MTImpedanceElement>;

// =========================================================================
// EM forward configuration.
// All knobs the user may want to tune for EM forward modeling.
// Passed to EMSolver::solve() and used by trust-region subsetting in
// EMActiveForward.
// =========================================================================
struct EMConfig {
    // --- Solver selection ---
    std::string solverMethod = "ie";  // "ie" (integral equation), "fdfd", "fdtd"

    // --- Frequency / time-domain ---
    double frequency_Hz = 0.0;        // 0 = time-domain, >0 = frequency-domain
    std::vector<double> timeGates_s;  // gate center times for time-domain (seconds)

    // --- Trust-region active-EM subsetting ---
    // Early iterations: wide margin to capture gross conductor location.
    // Late iterations: narrow margin, focused on conductive groups only.
    // Periodic full-mesh checks catch shielding / out-of-region conductors.
    //
    // EXTENSION POINT: To add a new subsetting strategy, implement it in
    //   EMActiveForward::buildActiveSet() and add configuration fields here.
    int    subsettingIterationsWide = 20;        // use wide margin for first N iters
    double subsettingSkinDepthMultiplierWide = 3.0;   // margin multiplier, early
    double subsettingSkinDepthMultiplierNarrow = 1.5; // margin multiplier, late
    double subsettingConductivityThreshold = 1e-3;    // S/m — filter groups below this
    int    subsettingFullCheckInterval = 10;     // validate full mesh every N iters
    double subsettingFullCheckTolerance = 0.05;  // relative misfit difference trigger

    // --- Skin depth reference ---
    // δ = 503 / sqrt(σ · f)  [m, for S/m and Hz]
    // For time-domain, f_equiv = 1 / (2 · t_gate_center)
    bool   autoComputeSkinDepth = true;

    // --- Physics ---
    bool   includeDisplacementCurrents = false; // quasi-static for most EM prospecting
    bool   useTotalField = false;               // false = secondary field only (default)
    double backgroundConductivity = 1e-4;       // S/m — half-space reference

    // --- Numerical ---
    int    integrationPointsPerTriangle = 1;    // quadrature points for facet integration
    double nearFieldThreshold_m = 1.0;          // switch to analytical formula below this distance
    int    maxSubsurfaceBounces = 0;            // >0 for layered-Earth IE (future extension)

    EMConfig() = default;
};

// =========================================================================
// EXTENSION POINT: To add a new EM survey type (e.g. borehole EM, CSEM):
//   1. Add a value to EMSurveyType enum above
//   2. Add data struct below (or in a new header)
//   3. Implement a new forward class extending EMForwardModel
//   4. Register it in JointObjective via addObjective()
// =========================================================================

} // namespace litho_invert

