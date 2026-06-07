#pragma once
#include <litho_invert/litho/lithology_model.h>
#include <litho_invert/core/common.h>
#include <litho_invert/em/em_data.h>
#include <memory>

namespace litho_invert {

// Generate true synthetic Voisey's Bay model with 3 litho groups.
// Returns the true model (used for computing synthetic observations).
std::shared_ptr<LithologyModel> generateTrueModel();

// Generate initial flat surfaces for the inversion starting model.
// The model will have 2 surfaces with given vertex freedoms.
std::shared_ptr<LithologyModel> generateInitialModel();

// Generate observation points on a regular grid above the model.
// Returns GravityData with g_obs and g_std set to 0 (caller fills g_obs).
GravityData generateObservationPoints();

// Compute synthetic gravity data from the true model.
GravityData computeSyntheticData(std::shared_ptr<LithologyModel> trueModel,
                                  const GravityData& observationPoints,
                                  double paddingDensity = 0.0);

// Compute synthetic magnetic data from the true model.
// Uses the susceptibility values set in the true model groups.
MagneticData computeSyntheticMagnetic(std::shared_ptr<LithologyModel> trueModel,
                                       const GravityData& observationPoints,
                                       double inc_deg, double dec_deg,
                                       double field_nT,
                                       double paddingSusceptibility = 0.0);

// Build a simple airborne TEM source/receiver geometry for testing.
// Returns a single vertical-dipole source at center with one z-component
// receiver at the same position (coaxial bird).
void buildAirborneTEMGeometry(std::vector<EMSource>& sources,
                               std::vector<EMReceiver>& receivers,
                               EMConfig& config);

// Compute synthetic active EM data from the true model.
// Uses the conductivity values set in the true model groups.
ActiveEMData computeSyntheticActiveEM(std::shared_ptr<LithologyModel> trueModel,
                                       const GravityData& observationPoints,
                                       const std::vector<EMSource>& sources,
                                       const std::vector<EMReceiver>& receivers,
                                       const EMConfig& config,
                                       double paddingConductivity = 1e-4);

// Build sparse MT station geometry for testing.
// Returns 4 stations at accessible locations around the survey area.
void buildMTStations(std::vector<MTStation>& stations,
                      std::vector<double>& frequencies_Hz);

// Compute synthetic MT data from the true model.
MTData computeSyntheticMT(std::shared_ptr<LithologyModel> trueModel,
                           const std::vector<MTStation>& stations,
                           const EMConfig& config,
                           double paddingConductivity = 1e-4);

} // namespace litho_invert
