#pragma once
#include <litho_invert/core/common.h>
#include <litho_invert/em/em_data.h>
#include <litho_invert/litho/lithology_model.h>
#include <memory>

namespace litho_invert {

// =========================================================================
// EMSolver — stateless compute engine for electromagnetic responses.
//
// The solver takes geometry (LithologyModel), physical properties
// (conductivity per group), and survey parameters (EMConfig) and returns
// predicted data at receiver positions.
//
// This interface exists so the forward model can swap solvers without
// changing the inversion framework:
//   - "ie"  — integral equation on polyhedra (default, same mesh as grav/mag)
//   - "fdfd" — finite-difference frequency domain (voxelized sub-grid)
//   - "fdtd" — finite-difference time domain (full waveform)
//
// EXTENSION POINT: To add a new solver:
//   1. Implement this interface
//   2. Register a factory string in EMActiveForward::createSolver()
//   3. Add any solver-specific config to EMConfig (this file)
//
// The solver is designed to be swappable at runtime via the solverMethod
// string in EMConfig, so the same inversion can use different solvers
// for different survey types (e.g. IE for airborne, FDFD for large-loop).
// =========================================================================
class EMSolver {
public:
    virtual ~EMSolver() = default;

    // --- Active-source EM (airborne, large-loop) ---
    // Compute secondary-field response at receiver positions for a given
    // source. Returns one value per receiver per time gate / frequency.
    //
    // srcPos: transmitter center position
    // recvPositions: receiver positions
    // model: lithology model with current geometry and conductivities
    // config: solver + physics configuration
    //
    // Returns VectorXd sized recvPositions.size() * config.timeGates_s.size()
    // (or * 1 for single-frequency).  Layout: all gates for recv 0, then
    // all gates for recv 1, etc.
    virtual VectorXd solveActive(const Vector3d& srcPos,
                                 const std::vector<Vector3d>& recvPositions,
                                 const LithologyModel& model,
                                 const EMConfig& config) = 0;

    // --- Passive-source EM (MT) ---
    // Compute the 2×2 complex impedance tensor at each station for each
    // frequency. Returns 8 real values per station per frequency:
    //   Re(Zxx), Im(Zxx), Re(Zxy), Im(Zxy), Re(Zyx), Im(Zyx), Re(Zyy), Im(Zyy)
    //
    // stations: MT station positions
    // frequencies_Hz: frequencies to compute (same for all stations, or
    //                 the solver can interpolate internally)
    // model: lithology model with current geometry and conductivities
    // config: solver + physics configuration
    virtual VectorXd solveMT(const std::vector<Vector3d>& stationPositions,
                             const std::vector<double>& frequencies_Hz,
                             const LithologyModel& model,
                             const EMConfig& config) = 0;

    // --- Query solver capabilities ---
    virtual std::string methodName() const = 0;
    virtual bool supportsTimeDomain() const { return true; }
    virtual bool supportsFrequencyDomain() const { return true; }
    virtual bool supportsMT() const { return false; }
};

// =========================================================================
// EXTENSION POINT: Solver factory.
// Add new solver creation logic here when adding a new EMSolver
// implementation.  Called by EMActiveForward / EMMTForward during setup.
// =========================================================================
std::shared_ptr<EMSolver> createEMSolver(const std::string& method);

} // namespace litho_invert

