#include <litho_invert/inversion/lbfgsb_optimizer.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace litho_invert {

LBFGSBOptimizer::LBFGSBOptimizer()
    : m_historySize(10)
{
}

VectorXd LBFGSBOptimizer::project(const VectorXd& x,
                                   const VectorXd& lb,
                                   const VectorXd& ub)
{
    VectorXd xp = x;
    for (Index i = 0; i < x.size(); ++i) {
        xp[i] = std::max(lb[i], std::min(ub[i], x[i]));
    }
    return xp;
}

void LBFGSBOptimizer::projectInPlace(VectorXd& x,
                                      const VectorXd& lb,
                                      const VectorXd& ub)
{
    for (Index i = 0; i < x.size(); ++i) {
        x[i] = std::max(lb[i], std::min(ub[i], x[i]));
    }
}

VectorXd LBFGSBOptimizer::computeSearchDirection(const VectorXd& g)
{
    if (m_sHistory.empty()) {
        // No history: steepest descent
        return -g;
    }

    // Compute gamma (scaling for initial Hessian approximation H0 = gamma * I)
    // gamma = s_{k-1}^T y_{k-1} / y_{k-1}^T y_{k-1}
    const VectorXd& sLast = m_sHistory.back();
    const VectorXd& yLast = m_yHistory.back();
    double sy = sLast.dot(yLast);
    double yy = yLast.squaredNorm();
    double gamma = (yy > std::numeric_limits<double>::epsilon()) ? (sy / yy) : 1.0;

    int m = static_cast<int>(m_sHistory.size());

    // Store alpha values from first loop
    std::vector<double> alpha(m);

    // First loop: forward through history (most recent to oldest)
    // q starts as g
    VectorXd q = g;
    for (int i = m - 1; i >= 0; --i) {
        alpha[i] = m_rhoHistory[i] * m_sHistory[i].dot(q);
        q -= alpha[i] * m_yHistory[i];
    }

    // r = H0 * q = gamma * q
    VectorXd r = gamma * q;

    // Second loop: backward through history (oldest to most recent)
    for (int i = 0; i < m; ++i) {
        double beta = m_rhoHistory[i] * m_yHistory[i].dot(r);
        r += m_sHistory[i] * (alpha[i] - beta);
    }

    // r = H_k * g_k, so search direction = -r
    return -r;
}

double LBFGSBOptimizer::projectedGradientNorm(const VectorXd& x,
                                                const VectorXd& g,
                                                const VectorXd& lb,
                                                const VectorXd& ub)
{
    // Compute norm of projected gradient: ||x - project(x - g, lb, ub)||
    VectorXd xMinusG = x - g;
    VectorXd proj = project(xMinusG, lb, ub);
    return (x - proj).norm();
}

double LBFGSBOptimizer::lineSearch(
    std::function<double(const VectorXd&)> objective,
    std::function<VectorXd(const VectorXd&)> gradient,
    const VectorXd& d, VectorXd& xNew, VectorXd& gNew)
{
    // Short-circuit: always accept the full step without backtracking.
    if (m_disableLineSearch) {
        xNew = m_x + d;
        projectInPlace(xNew, m_lb, m_ub);
        gNew = gradient(xNew);
        return 1.0;
    }

    double alpha = 1.0;

    double f0 = objective(m_x);
    double gd = m_g.dot(d);

    if (gd >= 0.0) {
        alpha = 0.0;
        xNew = m_x;
        gNew = m_g;
        return alpha;
    }

    const double minAlpha = 1e-16;

    for (int iter = 0; iter < m_lineSearchMaxIter; ++iter) {
        if (alpha < minAlpha) {
            break;
        }

        VectorXd xTry = m_x + alpha * d;
        projectInPlace(xTry, m_lb, m_ub);
        double fTry = objective(xTry);

        // Projected Armijo: use actual displacement after bound-clamping
        // rather than alpha*gd, which overestimates decrease for clamped DOFs.
        double projectedDecrease = m_g.dot(xTry - m_x);
        if (projectedDecrease >= 0.0) {
            alpha *= 0.5;
            continue;
        }
        double threshold = f0 + m_armijoC1 * projectedDecrease;

        if (fTry <= threshold) {
            xNew = xTry;
            gNew = gradient(xNew);
            return alpha;
        }

        alpha *= 0.5;
    }

    alpha = minAlpha;
    xNew = m_x + alpha * d;
    projectInPlace(xNew, m_lb, m_ub);
    gNew = gradient(xNew);
    return alpha;
}

void LBFGSBOptimizer::updateHistory(const VectorXd& s, const VectorXd& y)
{
    double rho = 1.0 / y.dot(s);

    if (static_cast<int>(m_sHistory.size()) >= m_historySize) {
        // Remove oldest entry (circular buffer behavior)
        m_sHistory.erase(m_sHistory.begin());
        m_yHistory.erase(m_yHistory.begin());
        m_rhoHistory.erase(m_rhoHistory.begin());
    }

    m_sHistory.push_back(s);
    m_yHistory.push_back(y);
    m_rhoHistory.push_back(rho);
}

OptimizerResult LBFGSBOptimizer::minimize(
    std::function<double(const VectorXd&)> objective,
    std::function<VectorXd(const VectorXd&)> gradient,
    const VectorXd& initialParams,
    const VectorXd& lowerBounds,
    const VectorXd& upperBounds)
{
    OptimizerResult result;
    m_iterCount = 0;

    // Initialize state
    m_lb = lowerBounds;
    m_ub = upperBounds;
    m_x = project(initialParams, m_lb, m_ub);
    m_g = gradient(m_x);

    // Clear history (skip when warm-starting from a previous run)
    if (m_clearOnMinimize) {
        m_sHistory.clear();
        m_yHistory.clear();
        m_rhoHistory.clear();
    }

    // Use the reference objective value for scale-aware convergence.
    // fRef is set by setReferenceObjective() from the caller — typically the
    // objective value at the initial guess. This avoids an extra full forward
    // evaluation inside the convergence check.
    double fRef = std::abs(m_fRef);

    double pgNorm0 = 0.0;  // initial projected gradient norm for relative convergence
    int stallCount = 0;
    for (int k = 0; k < m_maxIter; ++k) {
        m_iterCount = k + 1;

        // Convergence: projected gradient norm relative to initial value.
        // The first iteration is never converged — we must try at least one step.
        double pgNorm = projectedGradientNorm(m_x, m_g, m_lb, m_ub);
        if (k == 0) pgNorm0 = pgNorm;

        if (k > 0) {
            double convergenceThreshold = m_ftol * std::max(1.0, pgNorm0);
            if (pgNorm < convergenceThreshold) {
                result.converged = true;
                result.message = "Converged: projected gradient below tolerance";
                break;
            }
        }

        // Compute search direction
        VectorXd d = computeSearchDirection(m_g);

        // Line search
        VectorXd xNew, gNew;
        double alpha = lineSearch(objective, gradient, d, xNew, gNew);

        if (alpha == 0.0) {
            result.converged = true;
            result.message = "Converged: no step taken";
            break;
        }

        // Stall detection: if line search returns a near-zero step for
        // several consecutive iterations, the Armijo condition cannot be
        // satisfied at any scale — further iterations won't help.
        if (alpha < 1e-8) {
            ++stallCount;
            if (stallCount >= 3) {
                result.converged = true;
                result.message =
                    "Converged: line search stalled (alpha < 1e-8 for "
                    + std::to_string(stallCount) + " consecutive iterations)";
                break;
            }
        } else {
            stallCount = 0;
        }

        // Compute s and y
        VectorXd s = xNew - m_x;
        VectorXd y = gNew - m_g;

        // Update history if curvature condition is satisfied
        double sy = s.dot(y);
        if (sy > 1e-16) {
            updateHistory(s, y);
        }

        // Advance to next iterate
        m_x = xNew;
        m_g = gNew;
    }

    result.params = m_x;
    result.finalValue = objective(m_x);
    result.iterations = m_iterCount;

    if (!result.converged && m_iterCount >= m_maxIter) {
        result.message = "Maximum iterations reached";
    }

    return result;
}

void LBFGSBOptimizer::restart(const VectorXd& newParams)
{
    m_sHistory.clear();
    m_yHistory.clear();
    m_rhoHistory.clear();
    m_x = project(newParams, m_lb, m_ub);
    m_iterCount = 0;
}

void LBFGSBOptimizer::setHistory(const std::vector<VectorXd>& s,
                                  const std::vector<VectorXd>& y,
                                  const std::vector<double>& rho)
{
    m_sHistory = s;
    m_yHistory = y;
    m_rhoHistory = rho;
}

void LBFGSBOptimizer::getHistory(std::vector<VectorXd>& s,
                                  std::vector<VectorXd>& y,
                                  std::vector<double>& rho) const
{
    s = m_sHistory;
    y = m_yHistory;
    rho = m_rhoHistory;
}

} // namespace litho_invert

