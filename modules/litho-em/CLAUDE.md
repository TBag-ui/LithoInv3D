# litho-em — EM Forward Models (Active EM + Magnetotelluric)

## Purpose

Active-source and magnetotelluric electromagnetic forward modeling for conductivity
constraints in the joint inversion. **IE solver is stubbed** — both forward models
currently use fallback approximations. This module is separated from litho-forward
so the IE solver can be implemented independently without touching the working
gravity/magnetic code.

## API

```cpp
#include <litho_invert/em/em_data.h>
#include <litho_invert/em/em_solver.h>
#include <litho_invert/em/em_forward_model.h>
#include <litho_invert/em/em_active_forward.h>
#include <litho_invert/em/em_mt_forward.h>

// === Data Types ===
struct EMSource { Vector3d position; EMSurveyType type; string waveform; ... };
struct EMReceiver { Vector3d position; string component; };
struct ActiveEMPoint { int sourceIndex, receiverIndex, gateIndex; double d_obs, d_std; ... };
struct MTImpedanceElement { int stationIndex, frequencyIndex, iComp, jComp; double zReal_obs, zImag_obs; ... };
struct MTStation { Vector3d position; string name; vector<double> frequencies_Hz; };

struct EMConfig {
    string solverMethod = "ie";     // "ie", "fdfd", "fdtd"
    double frequency_Hz = 0.0;      // 0 = time-domain
    vector<double> timeGates_s;
    // Trust-region subsetting
    int subsettingIterationsWide = 20;
    double subsettingSkinDepthMultiplierWide = 3.0;
    double subsettingSkinDepthMultiplierNarrow = 1.5;
    double subsettingConductivityThreshold = 1e-3;
    // Physics
    bool useTotalField = false;
    double backgroundConductivity = 1e-4;
    // Numerical
    int integrationPointsPerTriangle = 1;
    double nearFieldThreshold_m = 1.0;
};

// === Solver Interface (stubbed) ===
class EMSolver {
    virtual VectorXd solveActive(srcPos, recvPositions, model, config) = 0;
    virtual VectorXd solveMT(stationPositions, frequencies, model, config) = 0;
    virtual string methodName() const = 0;
};
shared_ptr<EMSolver> createEMSolver(method);  // returns nullptr for all methods

// === Active EM Forward ===
EMActiveForward(shared_ptr<LithologyModel>, vector<EMSource>, vector<EMReceiver>, ActiveEMData, EMConfig);
VectorXd compute(params);
MatrixXd computeJacobian(params);
VectorXd computeGroupUnitResponse(groupIndex);
unordered_set<int> buildActiveSet(sourceIndex, activeVolume);  // trust-region
void enablePadding(bool, depth);

// === MT Forward ===
EMMTForward(shared_ptr<LithologyModel>, vector<MTStation>, MTData, EMConfig);
VectorXd compute(params);          // 8·nStations·nFreq values
MatrixXd computeJacobian(params);
```

## Current Status

- **Active EM fallback**: conductive-dipole approx: response ∼ σ·δ³/R³
- **MT fallback**: 1D half-space impedance: Z = √(iωμ₀/σ)
- **Skin depth**: δ = 503/√(σ·f) [m]; time-domain: δ ≈ 503·√(t_gate/σ)
- **Trust-region subsetting**: functional — restricts active groups to near-source volume

## Extension: Implementing the IE Solver

1. Create `src/em/ie_solver.cpp` with class `IESolver : public EMSolver`
2. Implement `solveActive()` — integral equation over polyhedral cells
3. Implement `solveMT()` — plane-wave IE over full mesh
4. Register: `if (method == "ie") return make_shared<IESolver>();`
5. Remove fallback code in `computeSourceReceiver()` / `compute()`

## Build

```powershell
cd modules/litho-em
qmake litho-em.pro
nmake release
```

## Dependencies

- litho-core, litho-surface, litho-model
