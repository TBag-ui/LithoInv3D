# litho-regularization — Surface Smoothness and Reference Model

## Purpose

Regularization terms for the inversion objective. Penalizes rough surfaces
(Laplacian curvature) and/or deviation from a reference starting model.

## API

```cpp
#include <litho_invert/regularization/regularization.h>
#include <litho_invert/regularization/smoothness.h>
#include <litho_invert/regularization/reference_model.h>

// === Abstract base ===
class Regularization {
public:
    virtual double evaluate(const VectorXd& params) = 0;
    virtual VectorXd gradient(const VectorXd& params) = 0;
    void setWeight(double lambda);
    double weight() const;
};

// === Surface Smoothness (Laplacian curvature) ===
class SurfaceSmoothness : public Regularization {
    explicit SurfaceSmoothness(shared_ptr<LithologyModel> model);
    // Φ_r = ½ Σ_v ||pos(v) − mean(neighbors(v))||²
    // Flat surface → zero penalty; folds/creases → quadratic penalty
};

// === Reference Model ===
class ReferenceModelRegularization : public Regularization {
    explicit ReferenceModelRegularization(shared_ptr<LithologyModel> model);
    void captureReference(const VectorXd& params);  // snapshot starting model
    // Φ_ref = ½·λ_ref·||m − m_start||²
    // Penalizes deviation from the captured starting model
};
```

## Mathematical Details

### Laplacian Curvature

For each free vertex, compute the Laplacian: L(v) = v − mean(neighbors(v)).
Penalty = ½ Σ_v ||L(v)||². Gradient = L(v) × ∂v/∂params (chain rule through
the parameter mapping, respecting control-point downsampling).

Flat grid surfaces have L(v) = 0 (neighbors are co-planar). The penalty grows
quadratically as vertices deviate from their local best-fit plane.

### Reference Model

Simple L² distance in parameter space. Captures the starting parameter vector
after setup (after control-point stride, padding extrapolation, etc.).
Controlled by `[regularization]` INI section: `enable_reference_model = true`,
`lambda_ref = 0.1`.

## Build

```powershell
cd modules/litho-regularization
qmake litho-regularization.pro
nmake release
```

## Dependencies

- litho-core (Eigen types)
- litho-surface (SurfaceMesh neighbors)
- litho-model (LithologyModel DOF system)
