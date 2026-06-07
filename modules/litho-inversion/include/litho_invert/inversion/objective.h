#pragma once
#include <litho_invert/core/common.h>
#include <litho_invert/forward/forward_model.h>
#include <litho_invert/regularization/regularization.h>
#include <litho_invert/inversion/constraint_handler.h>
#include <memory>
#include <vector>

namespace litho_invert {

class ObjectiveFunction {
public:
    ObjectiveFunction(std::shared_ptr<ForwardModel> forward,
                      const GravityData& data);

    void addRegularization(std::shared_ptr<Regularization> reg);
    void setConstraintHandler(std::shared_ptr<ConstraintHandler> handler);

    double evaluate(const VectorXd& params);
    VectorXd gradient(const VectorXd& params);
    double dataMisfit(const VectorXd& params);
    VectorXd dataMisfitGradient(const VectorXd& params);
    VectorXd residuals(const VectorXd& params);  // obs - predicted

    size_t parameterCount() const;
    size_t dataCount() const;
    const GravityData& data() const { return m_data; }
    ForwardModel& forwardModel() { return *m_forward; }
    const ForwardModel& forwardModel() const { return *m_forward; }

    struct Components {
        double dataMisfit = 0.0;
        double regularization = 0.0;
        double constraintPenalty = 0.0;
        double total = 0.0;
    };
    Components evaluateComponents(const VectorXd& params);

private:
    std::shared_ptr<ForwardModel> m_forward;
    std::vector<std::shared_ptr<Regularization>> m_regs;
    std::shared_ptr<ConstraintHandler> m_constraints;
    const GravityData& m_data;
};

} // namespace litho_invert
