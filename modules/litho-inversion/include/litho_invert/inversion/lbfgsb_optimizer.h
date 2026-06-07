#pragma once
#include <litho_invert/inversion/optimizer.h>

namespace litho_invert {

class LBFGSBOptimizer : public Optimizer {
public:
    LBFGSBOptimizer();

    OptimizerResult minimize(
        std::function<double(const VectorXd&)> objective,
        std::function<VectorXd(const VectorXd&)> gradient,
        const VectorXd& initialParams,
        const VectorXd& lowerBounds,
        const VectorXd& upperBounds) override;

    void restart(const VectorXd& newParams) override;
    void setMaxIterations(int n) override { m_maxIter = n; }
    int iterationCount() const override { return m_iterCount; }

    void setHistorySize(int m) { m_historySize = m; }
    int historySize() const { return m_historySize; }

    // Warm-start: when false, s/y/rho history is preserved across minimize()
    // calls so the next run starts with a Hessian approximation instead of
    // steepest descent. Default true (cold-start, backward compatible).
    void setClearOnMinimize(bool v) { m_clearOnMinimize = v; }
    bool clearOnMinimize() const { return m_clearOnMinimize; }

    // Transfer L-BFGS history between optimizer instances for warm-start.
    void setHistory(const std::vector<VectorXd>& s,
                    const std::vector<VectorXd>& y,
                    const std::vector<double>& rho);
    void getHistory(std::vector<VectorXd>& s,
                    std::vector<VectorXd>& y,
                    std::vector<double>& rho) const;
    bool hasHistory() const { return !m_sHistory.empty(); }

    // Armijo constant for line search sufficient-decrease condition.
    // Default 1e-4 (standard). Increase to 1e-2 for noisy geophysical gradients.
    void setArmijoC1(double c1) { m_armijoC1 = c1; }
    double armijoC1() const { return m_armijoC1; }

    // Max backtracking iterations in line search. Default 50.
    void setLineSearchMaxIter(int n) { m_lineSearchMaxIter = n; }
    int lineSearchMaxIter() const { return m_lineSearchMaxIter; }

    // When true, skip backtracking entirely and always accept α=1.0.
    void setDisableLineSearch(bool v) { m_disableLineSearch = v; }
    bool disableLineSearch() const { return m_disableLineSearch; }

    // Reference objective value for scale-aware convergence test.
    // Set to f(x0) before calling minimize() so the convergence threshold
    // scales with the problem magnitude rather than parameter magnitude.
    void setReferenceObjective(double fRef) { m_fRef = fRef; }
    double referenceObjective() const { return m_fRef; }

private:
    int m_historySize = 10;
    bool m_clearOnMinimize = true;  // cold-start by default
    double m_armijoC1 = 1e-4;
    int m_lineSearchMaxIter = 50;
    bool m_disableLineSearch = false;
    double m_fRef = 1.0;  // reference objective for convergence scaling

    // L-BFGS history buffers (stored oldest-first)
    std::vector<VectorXd> m_sHistory;  // s_k = x_{k+1} - x_k
    std::vector<VectorXd> m_yHistory;  // y_k = g_{k+1} - g_k
    std::vector<double> m_rhoHistory;  // 1/(y_k^T * s_k)

    // Current state
    VectorXd m_x;   // current iterate
    VectorXd m_g;   // current gradient
    VectorXd m_lb;  // lower bounds
    VectorXd m_ub;  // upper bounds

    // Clamp x to bounds
    static VectorXd project(const VectorXd& x, const VectorXd& lb, const VectorXd& ub);
    static void projectInPlace(VectorXd& x, const VectorXd& lb, const VectorXd& ub);

    // Two-loop recursion: compute search direction d = -H * g
    VectorXd computeSearchDirection(const VectorXd& g);

    // Projected gradient for convergence test
    static double projectedGradientNorm(const VectorXd& x, const VectorXd& g,
                                         const VectorXd& lb, const VectorXd& ub);

    // Line search with sufficient decrease condition and bounds
    double lineSearch(std::function<double(const VectorXd&)> objective,
                      std::function<VectorXd(const VectorXd&)> gradient,
                      const VectorXd& d, VectorXd& xNew, VectorXd& gNew);

    // Update L-BFGS history with new (s, y) pair
    void updateHistory(const VectorXd& s, const VectorXd& y);
};

} // namespace litho_invert

