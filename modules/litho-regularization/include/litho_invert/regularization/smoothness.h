#pragma once
#include <litho_invert/regularization/regularization.h>
#include <litho_invert/litho/lithology_model.h>
#include <memory>

namespace litho_invert {

// Penalizes Laplacian curvature: each vertex should be near the average of its neighbors.
// phi_r = 0.5 * SUM_{v} ||pos_v - mean(neighbors)||^2
class SurfaceSmoothness : public Regularization {
public:
    explicit SurfaceSmoothness(std::shared_ptr<LithologyModel> model);
    double evaluate(const VectorXd& params) override;
    VectorXd gradient(const VectorXd& params) override;
private:
    std::shared_ptr<LithologyModel> m_model;
};

} // namespace litho_invert

