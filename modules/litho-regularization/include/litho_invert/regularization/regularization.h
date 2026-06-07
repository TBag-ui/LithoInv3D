#pragma once
#include <litho_invert/core/common.h>

namespace litho_invert {

class Regularization {
public:
    virtual ~Regularization() = default;
    virtual double evaluate(const VectorXd& params) = 0;
    virtual VectorXd gradient(const VectorXd& params) = 0;
    void setWeight(double lambda) { m_weight = lambda; }
    double weight() const { return m_weight; }
protected:
    double m_weight = 1.0;
};

} // namespace litho_invert

