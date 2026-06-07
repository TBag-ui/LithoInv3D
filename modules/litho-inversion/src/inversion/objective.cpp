#include <litho_invert/inversion/objective.h>
#include <cmath>

namespace litho_invert {

ObjectiveFunction::ObjectiveFunction(
    std::shared_ptr<ForwardModel> forward,
    const GravityData& data)
    : m_forward(std::move(forward))
    , m_data(data)
{
}

void ObjectiveFunction::addRegularization(std::shared_ptr<Regularization> reg) {
    m_regs.push_back(std::move(reg));
}

void ObjectiveFunction::setConstraintHandler(std::shared_ptr<ConstraintHandler> handler) {
    m_constraints = std::move(handler);
}

double ObjectiveFunction::dataMisfit(const VectorXd& params) {
    VectorXd predicted = m_forward->compute(params);
    double misfit = 0.0;
    for (size_t i = 0; i < m_data.size(); ++i) {
        double w = (m_data[i].g_std > 0) ? (1.0 / m_data[i].g_std) : 1.0;
        double res = w * (m_data[i].g_obs - predicted[static_cast<Index>(i)]);
        misfit += 0.5 * res * res;
    }
    return misfit;
}

VectorXd ObjectiveFunction::residuals(const VectorXd& params) {
    VectorXd predicted = m_forward->compute(params);
    VectorXd resid(static_cast<Index>(m_data.size()));
    for (size_t i = 0; i < m_data.size(); ++i) {
        resid[static_cast<Index>(i)] = m_data[i].g_obs - predicted[static_cast<Index>(i)];
    }
    return resid;
}

VectorXd ObjectiveFunction::dataMisfitGradient(const VectorXd& params) {
    MatrixXd J = m_forward->computeJacobian(params);
    VectorXd predicted = m_forward->compute(params);

    // Weighted residuals: Wr = W^T * W * r
    VectorXd Wr(static_cast<Index>(m_data.size()));
    for (size_t i = 0; i < m_data.size(); ++i) {
        double w = (m_data[i].g_std > 0) ? (1.0 / m_data[i].g_std) : 1.0;
        double r = m_data[i].g_obs - predicted[static_cast<Index>(i)];
        Wr[static_cast<Index>(i)] = w * w * r;
    }

    return -J.transpose() * Wr;
}

double ObjectiveFunction::evaluate(const VectorXd& params) {
    return evaluateComponents(params).total;
}

VectorXd ObjectiveFunction::gradient(const VectorXd& params) {
    VectorXd grad = dataMisfitGradient(params);

    for (auto& reg : m_regs) {
        grad += reg->gradient(params);
    }

    if (m_constraints) {
        grad += m_constraints->gradient(params);
    }

    return grad;
}

size_t ObjectiveFunction::parameterCount() const {
    return m_forward->parameterCount();
}

size_t ObjectiveFunction::dataCount() const {
    return m_forward->dataCount();
}

ObjectiveFunction::Components ObjectiveFunction::evaluateComponents(const VectorXd& params) {
    Components c;
    VectorXd predicted = m_forward->compute(params);

    for (size_t i = 0; i < m_data.size(); ++i) {
        double w = (m_data[i].g_std > 0) ? (1.0 / m_data[i].g_std) : 1.0;
        double res = w * (m_data[i].g_obs - predicted[static_cast<Index>(i)]);
        c.dataMisfit += 0.5 * res * res;
    }

    for (auto& r : m_regs) {
        c.regularization += r->evaluate(params);
    }

    if (m_constraints) {
        c.constraintPenalty = m_constraints->evaluate(params);
    }

    c.total = c.dataMisfit + c.regularization + c.constraintPenalty;
    return c;
}

} // namespace litho_invert
