#include <litho_invert/inversion/constraint_handler.h>
#include <cmath>
#include <cstdint>

namespace litho_invert {

ConstraintHandler::ConstraintHandler(
    std::shared_ptr<LithologyModel> model,
    const std::vector<Constraint>& constraints)
    : m_model(std::move(model))
    , m_constraints(constraints)
{
}

double ConstraintHandler::computeViolation(
    const Constraint& c, const LithologyModel& model) const
{
    // Convert constraint depth (positive down) to model z (positive up, negative down)
    double midZ = -(c.z_top + c.z_bottom) / 2.0;
    Vector3d midpoint(c.position.x(), c.position.y(), midZ);

    int actualGroup = model.classifyPoint(midpoint);
    if (actualGroup == c.litho_group_id) {
        return 0.0;
    }

    // Violation proportional to constraint interval thickness
    return std::abs(c.z_bottom - c.z_top);
}

double ConstraintHandler::evaluate(const VectorXd& params) {
    m_model->applyParameterVector(params);

    double penalty = 0.0;
    for (const auto& c : m_constraints) {
        double v = computeViolation(c, *m_model);
        penalty += v * v;
    }

    return m_weight * penalty;
}

VectorXd ConstraintHandler::gradient(const VectorXd& params) {
    const uint32_t nDof = m_model->totalDofCount();
    VectorXd grad = VectorXd::Zero(static_cast<Index>(nDof));

    // Use central finite differences to approximate the gradient.
    // The classification function is discrete (nearest-neighbor sampling),
    // so the analytical gradient is zero almost everywhere. Finite differences
    // at a practical step size effectively "blur" the gradient over the step
    // distance, providing useful directional information.
    const double h = 1.0; // 1 m FD step

    for (uint32_t j = 0; j < nDof; ++j) {
        VectorXd paramsPert = params;
        paramsPert(static_cast<Index>(j)) += h;
        double fp = evaluate(paramsPert);

        paramsPert(static_cast<Index>(j)) -= 2.0 * h;
        double fm = evaluate(paramsPert);

        grad(static_cast<Index>(j)) = (fp - fm) / (2.0 * h);
    }

    return m_weight * grad;
}

} // namespace litho_invert

