#include <litho_invert/em/em_forward_model.h>
#include <cmath>

namespace litho_invert {

EMForwardModel::EMForwardModel(std::shared_ptr<LithologyModel> model,
                               const EMConfig& config)
    : m_model(std::move(model))
    , m_config(config)
{
}

double EMForwardModel::skinDepth(double conductivity_Sm, double frequency_Hz) {
    // δ = 503 / sqrt(σ · f)  [m, S/m, Hz]
    // Derived from δ = 1 / sqrt(π · f · μ₀ · σ) with μ₀ = 4π×10⁻⁷
    if (conductivity_Sm <= 0.0 || frequency_Hz <= 0.0) {
        return 1e9;  // effectively infinite — no attenuation
    }
    return 503.0 / std::sqrt(conductivity_Sm * frequency_Hz);
}

double EMForwardModel::skinDepthTimeDomain(double conductivity_Sm, double timeGate_s) {
    // Approximate effective frequency: f_equiv ≈ 1 / (2 · t_gate)
    // then use the frequency-domain skin depth formula.
    if (conductivity_Sm <= 0.0 || timeGate_s <= 0.0) {
        return 1e9;
    }
    double fEquiv = 1.0 / (2.0 * timeGate_s);
    return skinDepth(conductivity_Sm, fEquiv);
}

} // namespace litho_invert

