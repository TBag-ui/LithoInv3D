#pragma once
#include <litho_invert/core/common.h>
#include <litho_invert/litho/lithology_model.h>
#include <memory>
#include <vector>

namespace litho_invert {

class ConstraintHandler {
public:
    ConstraintHandler(std::shared_ptr<LithologyModel> model,
                      const std::vector<Constraint>& constraints);

    double evaluate(const VectorXd& params);
    VectorXd gradient(const VectorXd& params);
    void setWeight(double omega) { m_weight = omega; }
    double weight() const { return m_weight; }
    size_t constraintCount() const { return m_constraints.size(); }

private:
    std::shared_ptr<LithologyModel> m_model;
    std::vector<Constraint> m_constraints;
    double m_weight = 1e6;  // large default near-hard constraint

    double computeViolation(const Constraint& c, const LithologyModel& model) const;
};

} // namespace litho_invert

