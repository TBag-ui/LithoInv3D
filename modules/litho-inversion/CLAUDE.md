# litho-inversion — Optimizer, Objective, and Runner

## Purpose

The inversion engine. Wires up forward models, regularization, and constraints into
an objective function, then minimizes it with L-BFGS-B. The `InversionRunner` orchestrates
the main loop: alternating geometry optimization and property inversion phases.

## API

```cpp
#include <litho_invert/inversion/optimizer.h>
#include <litho_invert/inversion/lbfgsb_optimizer.h>
#include <litho_invert/inversion/objective.h>
#include <litho_invert/inversion/joint_objective.h>
#include <litho_invert/inversion/constraint_handler.h>
#include <litho_invert/inversion/runner.h>

// === Optimizer ===
struct OptimizerResult {
    VectorXd params; double value; int iterations; bool converged;
};

class LBFGSBOptimizer : public Optimizer {
    void setHistorySize(int m);       // default 10
    void setClearOnMinimize(bool);    // false = warm-start across property phases
    OptimizerResult minimize(obj, grad, x0, lower, upper) override;
};

// === Objective ===
class ObjectiveFunction {
    ObjectiveFunction(shared_ptr<ForwardModel>, const GravityData&);
    void addRegularization(shared_ptr<Regularization>);
    void setConstraintHandler(shared_ptr<ConstraintHandler>);
    double evaluate(params);
    VectorXd gradient(params);
    Components evaluateComponents(params);
    // Components: { dataMisfit, regularization, constraintPenalty, total }
};

// === Joint Objective ===
class JointObjective {
    JointObjective(gravityObj, magForward, magData, alpha_mag);
    void addActiveEM(activeEMForward, data, alpha);
    void addMT(mtForward, data, alpha);
    double evaluate(params);
    VectorXd gradient(params);
    JointComponents evaluateComponents(params);
    // Plus: gravityMisfit(), magneticMisfit(), activeEMMisfit(), mtMisfit()
    // Plus: gravityResiduals(), magneticResiduals(), etc.
};

// === Constraints ===
class ConstraintHandler {
    ConstraintHandler(shared_ptr<LithologyModel>, const vector<Constraint>&);
    double evaluate(params);
    VectorXd gradient(params);
    void setWeight(omega);  // default 1e6
};

// === Runner ===
struct InversionConfig {
    shared_ptr<LithologyModel> model;
    GravityData observedData; vector<Constraint> constraints;
    double lambda, omega; int maxIterations, controlPointStride;
    bool enablePropertyInversion, enableReferenceModel, enableDepthBounds;
    // Magnetic (empty = skip)
    MagneticData magneticData; double magneticWeight, inc, dec, field_nT;
    // EM (empty = skip)
    ActiveEMData activeEMData; EMConfig emConfig;
    MTData mtData; vector<MTStation> mtStations;
    // Padding, topography, property bounds...
};

struct InversionIteration {
    int iteration; double dataMisfit, regularization, constraintPenalty;
    double totalObjective, rmsError, magneticMisfit;
    double dw_gravity_x, dw_gravity_y;  // Durbin-Watson all data types
};

struct InversionResult {
    shared_ptr<LithologyModel> finalModel; vector<InversionIteration> history;
    VectorXd predictedData; bool converged; int totalIterations;
    double finalMisfit, finalRMS;
    VectorXd finalDensities, finalSusceptibilities, finalConductivities;
    shared_ptr<SurfaceMesh> closureTop, closureBottom, closureDeepBottom;
};

class InversionRunner {
    InversionRunner(const InversionConfig&);
    InversionResult run();
    void setIterationCallback(function<void(const InversionIteration&)>);
};
```

## Main Loop Pseudocode

```
While not converged:
  Phase A: Geometry Optimization
    optimizer.minimize(geoObj, geoGrad, params, lower, upper)
      geoObj → JointObjective::evaluateComponents + DW + callback
      geoGrad → JointObjective::gradient
    apply params → extrapolate padding
  Phase B: Property Inversion (if enabled)
    compute unit-response matrices U
    bounded LS solve: min ||W·(d_obs - U·ρ - U_pad·ρ_pad)||²
    write properties back to model groups
```

## Design Notes

- **Debug prints**: All `std::cerr` output guarded by `#ifdef LITHO_INVERT_DEBUG`
- **L-BFGS-B convergence**: ||proj(x − g) − x||_∞ < tolerance
- **Warm-start**: `setClearOnMinimize(false)` preserves history across property phases
- **Gravity analytical Jacobian**: Z-axis analytical, XY via finite differences — handled by `GravityForward::computeJacobian()`
- **Magnetic FD gradient**: via `computeJacobian()` (analytical not yet implemented)

## Build

```powershell
cd modules/litho-inversion
qmake litho-inversion.pro
nmake release
```

## Dependencies

ALL other modules: litho-core, litho-surface, litho-model, litho-forward,
litho-em, litho-io, litho-regularization

## Tests

```powershell
cd modules/litho-inversion/tests
qmake tests.pro && nmake release && release\tests.exe
```

Tests: quadratic bowl convergence, bounded Rosenbrock, constraint penalty,
objective evaluation, joint objective weighting, simple 1-DOF recovery.
