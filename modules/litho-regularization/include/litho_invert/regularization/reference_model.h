#pragma once
#include <litho_invert/regularization/regularization.h>
#include <litho_invert/litho/lithology_model.h>
#include <memory>

namespace litho_invert {

class ReferenceModelRegularization : public Regularization {
public:
    explicit ReferenceModelRegularization(std::shared_ptr<LithologyModel> model);

    double evaluate(const VectorXd& params) override;
    VectorXd gradient(const VectorXd& params) override;

    void captureReference(const VectorXd& params);

private:
    std::shared_ptr<LithologyModel> m_model;
    VectorXd m_refParams;
    bool m_captured = false;
};

} // namespace litho_invert
