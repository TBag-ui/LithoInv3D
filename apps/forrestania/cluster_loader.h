#pragma once
#include <litho_invert/core/common.h>
#include <litho_invert/litho/litho_group.h>
#include <litho_invert/litho/lithology_model.h>
#include <string>
#include <vector>
#include <memory>

namespace litho_invert {

// Physical properties loaded from Lithosplitter clustering output CSV.
// "n/a" or empty numeric fields are stored as NaN.
struct ClusterProperties {
    int cluster_id = -1;
    std::string working_name;
    int sample_count = 0;

    double density_median = 0.0;
    double density_p10 = 0.0;
    double density_p90 = 0.0;

    double susceptibility_median = 0.0;
    double susceptibility_p10 = 0.0;
    double susceptibility_p90 = 0.0;

    bool has_measured_density = false;
    bool has_measured_susceptibility = false;

    double remanence_magnitude = std::numeric_limits<double>::quiet_NaN();
    double remanence_inclination = std::numeric_limits<double>::quiet_NaN();
    double remanence_declination = std::numeric_limits<double>::quiet_NaN();
    bool has_measured_remanence = false;
};

// Load cluster properties from Lithosplitter inversion_start.csv.
// Format: cluster_id,working_name,sample_count,density_median_gcc,...,
//         susceptibility_median_SI,...,has_measured_density,
//         has_measured_susceptibility,remanence_magnitude_Am,
//         remanence_inclination_deg,remanence_declination_deg,
//         has_measured_remanence
std::vector<ClusterProperties> loadClusterProperties(const std::string& filepath);

// Build LithologyModel from cluster properties using cluster→group mapping.
// clusterToGroup maps cluster_id → litho_group_id.
// Properties are averaged across all clusters that map to the same group.
// Returns model with groups created but no surfaces (caller adds surfaces).
struct ClusterGroupMapping {
    int cluster_id;
    int litho_group_id;
};

std::shared_ptr<LithologyModel> buildModelFromClusters(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<ClusterGroupMapping>& mapping,
    double bottomDepth = -5000.0);

// Summarise remanence mode per group and initial values.
// Returns true if any group has unknown remanence (needs inversion).
bool configureRemanenceFromClusters(
    std::shared_ptr<LithologyModel> model,
    const std::vector<ClusterProperties>& clusters,
    const std::vector<ClusterGroupMapping>& mapping,
    RemanentMagnetizationMode& outMode);

// Build LithologyModel directly from cluster properties — one group per unique inv cluster_id.
// ordering: kGeoToInv — geo_group_index → cluster_id. If empty, uses cluster_id ascending.
// The inversion_start.csv must list clusters in stratigraphic order (top→bottom).
std::shared_ptr<LithologyModel> buildModelFromClustersDirect(
    const std::vector<ClusterProperties>& clusters,
    double bottomDepth = -5000.0,
    const std::vector<int>& ordering = {});

// Remanence config for direct mode — uses the same ordering.
bool configureRemanenceFromClustersDirect(
    std::shared_ptr<LithologyModel> model,
    const std::vector<ClusterProperties>& clusters,
    RemanentMagnetizationMode& outMode,
    const std::vector<int>& ordering = {});

} // namespace litho_invert
