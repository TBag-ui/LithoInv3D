#include <litho_invert/inversion/gncg_optimizer.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>

namespace litho_invert {

GNCGOptimizer::GNCGOptimizer()
{
}

VectorXd GNCGOptimizer::project(const VectorXd& x,
                                 const VectorXd& lb,
                                 const VectorXd& ub)
{
    VectorXd xp = x;
    for (Index i = 0; i < x.size(); ++i) {
        xp[i] = std::max(lb[i], std::min(ub[i], x[i]));
    }
    return xp;
}

double GNCGOptimizer::projectedGradientNorm(const VectorXd& x,
                                              const VectorXd& g,
                                              const VectorXd& lb,
                                              const VectorXd& ub)
{
    double sum = 0.0;
    for (Index i = 0; i < x.size(); ++i) {
        double pg = g[i];
        if (x[i] <= lb[i] + 1e-12 && g[i] > 0) pg = 0.0;
        if (x[i] >= ub[i] - 1e-12 && g[i] < 0) pg = 0.0;
        sum += pg * pg;
    }
    return std::sqrt(sum);
}

VectorXd GNCGOptimizer::cgSolve(
    const std::function<VectorXd(const VectorXd&)>& Hv,
    const VectorXd& b)
{
    // Solve H*x = b via conjugate gradient (H is symmetric, approx PSD)
    int n = static_cast<int>(b.size());
    VectorXd x = VectorXd::Zero(n);
    VectorXd r = b;          // r = b - H*x, with x=0: r = b
    VectorXd p = r;
    double rsOld = r.squaredNorm();

    if (rsOld < 1e-30) return x;

    for (int j = 0; j < m_cgMaxIter; ++j) {
        VectorXd Hp = Hv(p);
        double p2 = p.squaredNorm();
        if (p2 < 1e-30) break;

        double alpha = rsOld / p2;
        // Use curvature condition: only advance if direction is a descent direction
        double denom = p.dot(Hp);
        if (denom > 1e-16) {
            alpha = rsOld / denom;
        }
        // else: skip CG update, use steepest descent scaling

        x += alpha * p;
        r -= alpha * Hp;
        double rsNew = r.squaredNorm();

        if (rsNew < m_cgTol * m_cgTol) break;

        double beta = rsNew / rsOld;
        p = r + beta * p;
        rsOld = rsNew;
    }

    return x;
}

double GNCGOptimizer::lineSearch(
    std::function<double(const VectorXd&)> objective,
    const VectorXd& d,
    const VectorXd& lb,
    const VectorXd& ub,
    VectorXd& xNew)
{
    double f0 = objective(m_x);
    double alpha = 1.0;
    const double rho = 0.5;

    for (int i = 0; i < m_lineSearchMaxIter; ++i) {
        xNew = project(m_x + alpha * d, lb, ub);
        double fNew = objective(xNew);

        if (fNew <= f0) {
            return alpha;
        }
        alpha *= rho;
    }

    // No improvement found — return current point
    xNew = m_x;
    return 0.0;
}

OptimizerResult GNCGOptimizer::minimize(
    std::function<double(const VectorXd&)> objective,
    std::function<VectorXd(const VectorXd&)> gradient,
    const VectorXd& initialParams,
    const VectorXd& lowerBounds,
    const VectorXd& upperBounds)
{
    OptimizerResult result;
    m_x = project(initialParams, lowerBounds, upperBounds);

    double pgNorm0 = 0.0;
    int stallCount = 0;

    for (int k = 0; k < m_maxIter; ++k) {
        m_iterCount = k + 1;

        VectorXd g = gradient(m_x);
        double pgNorm = projectedGradientNorm(m_x, g, lowerBounds, upperBounds);
        if (k == 0) pgNorm0 = pgNorm;

        if (k > 0) {
            double convergenceThreshold = m_ftol * std::max(1.0, pgNorm0);
            if (pgNorm < convergenceThreshold) {
                result.converged = true;
                result.message = "Converged: projected gradient below tolerance";
                break;
            }
        }

        // Gauss-Newton step via CG on (H + reg)*Δx = -g
        // Compute Hessian-vector product via FD on gradient:
        //   H(x)*v ≈ (g(x+h*v) - g(x-h*v)) / (2h)
        // This captures J^T*J + reg_hessian, i.e. the full Gauss-Newton Hessian.
        auto Hv = [&](const VectorXd& v) -> VectorXd {
            double vn = v.norm();
            if (vn < 1e-15) return VectorXd::Zero(v.size());
            double h = m_fdStep / vn;
            VectorXd xp = project(m_x + h * v, lowerBounds, upperBounds);
            VectorXd xm = project(m_x - h * v, lowerBounds, upperBounds);
            VectorXd gp = gradient(xp);
            VectorXd gm = gradient(xm);
            return (gp - gm) / (2.0 * h);
        };

        VectorXd d = cgSolve(Hv, -g);

        // Check if CG produced a descent direction
        double dg = d.dot(g);
        if (dg >= 0.0) {
            // CG failed — fall back to steepest descent
            d = -g;
            dg = -g.squaredNorm();
        }

        // Line search
        VectorXd xNew;
        double alpha = lineSearch(objective, d, lowerBounds, upperBounds, xNew);

        if (alpha == 0.0) {
            result.converged = true;
            result.message = "Converged: no step taken";
            break;
        }

        // Stall detection
        double stepNorm = (xNew - m_x).norm();
        if (stepNorm < 1e-8) {
            ++stallCount;
            if (stallCount >= 3) {
                result.converged = true;
                result.message = "Converged: step stalled";
                break;
            }
        } else {
            stallCount = 0;
        }

        m_x = xNew;
    }

    result.params = m_x;
    result.finalValue = objective(m_x);
    result.iterations = m_iterCount;

    if (!result.converged && m_iterCount >= m_maxIter) {
        result.message = "Maximum iterations reached";
    }

    return result;
}

void GNCGOptimizer::restart(const VectorXd& newParams)
{
    m_x = newParams;
}

} // namespace litho_invert

