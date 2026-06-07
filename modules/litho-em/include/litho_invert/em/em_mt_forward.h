#pragma once
#include <litho_invert/em/em_forward_model.h>
#include <vector>

namespace litho_invert {

// =========================================================================
// EMMTForward — magnetotelluric forward model.
//
// Computes the 2×2 complex impedance tensor Z at each station for each
// frequency. The MT method uses natural plane-wave sources (no controlled
// transmitter), so there is no source geometry — only station positions
// and frequencies.
//
// The forward model predicts Re(Zxx), Im(Zxx), Re(Zxy), Im(Zxy),
// Re(Zyx), Im(Zyx), Re(Zyy), Im(Zyy) per station per frequency.
// These are compared with observed MTImpedanceElement values.
//
// Because MT stations are sparse (typically 10–50 stations accessible
// on foot from camp), there is no subsetting — the whole mesh contributes
// to every station.  This makes MT a good "regional constraint" that
// complements the high-resolution but localized active EM.
//
// EXTENSION POINT: To add tipper (vertical transfer function) or phase
//   tensor output, add columns to the predicted vector and extend
//   MTImpedanceElement with corresponding observed fields.
// =========================================================================
class EMMTForward : public EMForwardModel {
public:
    // stations: MT station descriptions (positions + frequency bands)
    // data: observed impedance tensor elements
    // model: lithology model
    // config: EM configuration (solver, frequencies)
    EMMTForward(std::shared_ptr<LithologyModel> model,
                const std::vector<MTStation>& stations,
                const MTData& data,
                const EMConfig& config);

    // --- ForwardModel interface ---
    VectorXd compute(const VectorXd& params) override;
    size_t dataCount() const override;
    size_t parameterCount() const override;

    // --- Unit response for property inversion ---
    VectorXd computeGroupUnitResponse(int groupIndex) override;

    // --- Padding ---
    void enablePadding(bool enabled, double paddingDepth = -100000.0);
    bool paddingEnabled() const { return m_paddingEnabled; }
    VectorXd computePaddingUnitResponse() const;

    // --- Accessors ---
    const std::vector<MTStation>& stations() const { return m_stations; }
    const MTData& mtData() const { return m_data; }

private:
    std::vector<MTStation> m_stations;
    const MTData& m_data;  // non-owning reference

    bool m_paddingEnabled = false;
    double m_paddingDepth = -100000.0;
    std::shared_ptr<SurfaceMesh> m_flatDeepBottom;

    void buildFlatDeepBottom();

    // Predicted data layout:
    // For each station, for each frequency: 8 values
    //   [Re(Zxx), Im(Zxx), Re(Zxy), Im(Zxy),
    //    Re(Zyx), Im(Zyx), Re(Zyy), Im(Zyy)]
    static constexpr int IMPEDANCE_COMPONENTS = 8;
};

} // namespace litho_invert

