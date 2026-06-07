#pragma once
#include <litho_invert/inversion/optimizer.h>

namespace litho_invert {

class GNCGOptimizer : public Optimizer {
public:
    GNCGOptimizer();

    OptimizerResult minimize(
        std::function<double(const VectorXd&)> objective,
        std::function<VectorXd(const VectorXd&)> gradient,
        const VectorXd& initialParams,
        const VectorXd& lowerBounds,
        const VectorXd& upperBounds) override;

    void restart(const VectorXd& newParams) override;
    void setMaxIterations(int n) override { m_maxIter = n; }
    int iterationCount() const override { return m_iterCount; }

    // Inner CG solver settings
    void setCGMaxIter(int n) { m_cgMaxIter = n; }
    int cgMaxIter() const { return m_cgMaxIter; }
    void setCGTolerance(double tol) { m_cgTol = tol; }
    double cgTolerance() const { return m_cgTol; }

    // Finite-difference step for Hessian-vector products
    void setFDStep(double h) { m_fdStep = h; }
    double fdStep() const { return m_fdStep; }

    // Armijo constant for line search
    void setArmijoC1(double c1) { m_armijoC1 = c1; }
    double armijoC1() const { return m_armijoC1; }

    // Max backtracking iterations in line search
    void setLineSearchMaxIter(int n) { m_lineSearchMaxIter = n; }
    int lineSearchMaxIter() const { return m_lineSearchMaxIter; }

private:
    int m_cgMaxIter = 50;
    double m_cgTol = 1e-6;
    double m_fdStep = 2.0;
    double m_armijoC1 = 1e-4;
    int m_lineSearchMaxIter = 30;

    VectorXd m_x;

    // Conjugate gradient: solve H*x = b, where Hv(v) computes H*v
    VectorXd cgSolve(
        const std::function<VectorXd(const VectorXd&)>& Hv,
        const VectorXd& b);

    // Armijo backtracking line search with bound projection
    double lineSearch(
        std::function<double(const VectorXd&)> objective,
        const VectorXd& d,
        const VectorXd& lb,
        const VectorXd& ub,
        VectorXd& xNew);

    static VectorXd project(const VectorXd& x, const VectorXd& lb, const VectorXd& ub);
    static double projectedGradientNorm(const VectorXd& x, const VectorXd& g,
                                         const VectorXd& lb, const VectorXd& ub);
};

} // namespace litho_invert

