#pragma once
#include <litho_invert/core/common.h>
#include <functional>

namespace litho_invert {

struct OptimizerResult {
    VectorXd params;
    double finalValue = 0.0;
    int iterations = 0;
    bool converged = false;
    std::string message;
};

class Optimizer {
public:
    virtual ~Optimizer() = default;

    virtual OptimizerResult minimize(
        std::function<double(const VectorXd&)> objective,
        std::function<VectorXd(const VectorXd&)> gradient,
        const VectorXd& initialParams,
        const VectorXd& lowerBounds,
        const VectorXd& upperBounds) = 0;

    virtual void restart(const VectorXd& newParams) = 0;
    virtual void setMaxIterations(int n) = 0;
    virtual int iterationCount() const = 0;
    virtual void setTolerance(double ftol) { m_ftol = ftol; }
    virtual double tolerance() const { return m_ftol; }

protected:
    double m_ftol = 1e-6;
    int m_maxIter = 500;
    int m_iterCount = 0;
};

} // namespace litho_invert

