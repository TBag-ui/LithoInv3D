#pragma once
#include <litho_invert/core/common.h>
#include <string>

namespace litho_invert {

enum class RemanentMagnetizationMode {
    EffectiveSusceptibility = 0,  // χ absorbs induced + remnant (current behaviour)
    FixedVectorPerGroup     = 1,  // user supplies remnant direction, invert for intensity
    VectorPerGroup          = 2   // invert for M_x, M_y, M_z per group
};

struct LithoGroup {
    int id = -1;
    std::string name;
    double density = 0.0;           // g/cm³
    double susceptibility = 0.0;    // SI (dimensionless) — induced only
    double conductivity = 1e-4;     // S/m — electrical conductivity

    // Remnant magnetization (SI volume magnetization, A/m).
    // Used only when remanenceMode != EffectiveSusceptibility.
    double remanence_magnitude = 0.0;      // A/m — inverted in FixedVector mode
    double remanence_inclination = 0.0;    // deg from horizontal, positive-down
    double remanence_declination = 0.0;    // deg from north, positive-east
    Vector3d magnetization;                // A/m — total remnant vector (Vector mode)

    LithoGroup() = default;
    LithoGroup(int id_, const std::string& name_, double density_,
               double susceptibility_ = 0.0, double conductivity_ = 1e-4)
        : id(id_), name(name_), density(density_),
          susceptibility(susceptibility_), conductivity(conductivity_) {}

    bool isValid() const { return id >= 0 && !name.empty() && density > 0.0; }
};

} // namespace litho_invert
