#include <litho_invert/regularization/reference_model.h>

namespace litho_invert {

ReferenceModelRegularization::ReferenceModelRegularization(
    std::shared_ptr<LithologyModel> model)
    : m_model(std::move(model))
{
}

void ReferenceModelRegularization::captureReference(const VectorXd& params) {
    m_refParams = params;
    m_captured = true;
}

double ReferenceModelRegularization::evaluate(const VectorXd& params) {
    if (!m_captured) return 0.0;
    VectorXd diff = params - m_refParams;
    return 0.5 * m_weight * diff.squaredNorm();
}

VectorXd ReferenceModelRegularization::gradient(const VectorXd& params) {
    if (!m_captured) return VectorXd::Zero(params.size());
    return m_weight * (params - m_refParams);
}

} // namespace litho_invert
