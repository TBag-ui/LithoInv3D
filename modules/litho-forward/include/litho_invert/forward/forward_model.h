#pragma once
#include <litho_invert/core/common.h>

namespace litho_invert {

class ForwardModel {
public:
    virtual ~ForwardModel() = default;

    // Compute predicted data for given parameter vector
    virtual VectorXd compute(const VectorXd& params) = 0;

    // Compute Jacobian (d_pred / d_param). Default: central finite differences.
    virtual MatrixXd computeJacobian(const VectorXd& params);

    // Compute both predicted data and Jacobian together
    virtual void computeBoth(const VectorXd& params,
                             VectorXd& predicted,
                             MatrixXd& jacobian);

    // Finite-difference step size for Jacobian computation.
    // Default 0.01: small enough to give stable property derivatives (density,
    // susceptibility, remanence), but still large enough to avoid numerical noise
    // in geometry derivatives at typical inversion depths (100-5000 m).
    void setFiniteDifferenceStep(double h) { m_fdStep = h; }
    double finiteDifferenceStep() const { return m_fdStep; }

    virtual size_t dataCount() const = 0;
    virtual size_t parameterCount() const = 0;

protected:
    double m_fdStep = 0.01;  // finite-difference step for Jacobian
};

} // namespace litho_invert
