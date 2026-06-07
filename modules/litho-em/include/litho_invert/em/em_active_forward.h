#pragma once
#include <litho_invert/em/em_forward_model.h>
#include <litho_invert/surface/surface_mesh.h>
#include <vector>
#include <set>
#include <unordered_set>

namespace litho_invert {

// =========================================================================
// EMActiveForward — active-source EM forward model.
//
// Handles airborne TEM/FEM and large-loop TEM surveys.  Each source has
// one or more receivers; the forward model computes the secondary field
// response for every source-receiver-gate combination in the ActiveEMData.
//
// TRUST-REGION SUBSETTING:
//   Because active EM is sensitive to conductive bodies but can be
//   "shielded" (a strong conductor obscures deeper conductors behind it),
//   we compute responses on a dynamically-sized subset of the mesh:
//
//   - Early iterations: wide bounding box (3× skin depth), all litho groups
//   - Late iterations:  narrow bounding box (1.5× skin depth), conductive
//     groups only (σ ≥ threshold)
//   - Periodic full-mesh validation guards against missed conductors
//
//   The subsetting algorithm is in buildActiveSet() and is driven by the
//   configuration fields in EMConfig (prefixed "subsetting").
//
// EXTENSION POINT: To add a new active-EM survey type:
//   1. Add value to EMSurveyType enum (em_data.h)
//   2. If it has a different source geometry, add fields to EMSource
//   3. The solver (EMSolver::solveActive) handles the physics — this
//      class only manages geometry, subsetting, and property scaling
// =========================================================================
class EMActiveForward : public EMForwardModel {
public:
    // sources: array of transmitters (one per airborne line or ground loop)
    // receivers: array of receivers (per source: one for airborne bird,
    //            multiple for large-loop ground receivers)
    // data: observation points (each tagged with source/receiver/gate indices)
    // model: lithology model
    // config: EM configuration (solver, subsetting, physics)
    EMActiveForward(std::shared_ptr<LithologyModel> model,
                    const std::vector<EMSource>& sources,
                    const std::vector<EMReceiver>& receivers,
                    const ActiveEMData& data,
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
    const std::vector<EMSource>& sources() const { return m_sources; }
    const std::vector<EMReceiver>& receivers() const { return m_receivers; }
    const ActiveEMData& emData() const { return m_data; }

    // --- Subsetting ---
    // Returns the set of litho-group indices that are "active" for the
    // given source at the current iteration.  The forward model uses this
    // internally; exposed for diagnostics.
    //
    // activeVolume: output — the bounding box (xmin, xmax, ymin, ymax, zmin, zmax)
    //   that defines the active region for this source.
    std::unordered_set<int> buildActiveSet(int sourceIndex,
                                           std::vector<double>& activeVolume) const;

private:
    std::vector<EMSource> m_sources;
    std::vector<EMReceiver> m_receivers;
    const ActiveEMData& m_data;  // non-owning reference

    bool m_paddingEnabled = false;
    double m_paddingDepth = -100000.0;
    std::shared_ptr<SurfaceMesh> m_flatDeepBottom;

    void buildFlatDeepBottom();

    // Compute response for a single source-receiver pair at all time gates
    // or frequencies, optionally restricted to active groups.
    VectorXd computeSourceReceiver(int sourceIndex, int receiverIndex,
                                   const std::vector<int>& gateIndices,
                                   const std::unordered_set<int>* activeGroups) const;

    // Get the dominant conductivity in the model (used for skin-depth calc)
    double dominantConductivity() const;
};

} // namespace litho_invert

