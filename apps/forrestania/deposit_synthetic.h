#pragma once
#include <litho_invert/litho/lithology_model.h>
#include <litho_invert/core/common.h>
#include "cluster_loader.h"
#include <memory>
#include <vector>
#include <string>
#include <map>

namespace litho_invert {

struct DepositData {
    struct Interval {
        std::string hole_id;
        double from_m, to_m;
        std::string lithology;
        double density, susceptibility;
    };
    std::vector<Interval> intervals;
    std::vector<std::string> holeIds;
    std::map<std::string, std::pair<double, double>> holePositions;
    std::vector<ClusterProperties> clusters;
};

// Load deposit data from directory (reads synthetic_dataset.csv + inversion_start.csv)
DepositData loadDepositData(const std::string& depositPath);

// Auto-derive cluster-to-group mapping by sorting clusters by density.
// nGroups = 2 or 3.  Clusters are sorted by density and split into groups.
std::vector<ClusterGroupMapping> deriveClusterMapping(const DepositData& data, int nGroups);

// Generate a synthetic true model for a deposit.
// Creates groups from cluster mapping, and surfaces with ellipsoidal
// perturbations derived from borehole lithology contacts.
std::shared_ptr<LithologyModel> generateDepositTrueModel(
    const DepositData& data,
    const std::vector<ClusterGroupMapping>& clusterMap);

// Generate 20x20 observation grid centred on the boreholes.
GravityData generateDepositObservationPoints(const DepositData& data);

// Generate constraints from borehole intervals.
std::vector<Constraint> generateDepositConstraints(
    const DepositData& data,
    const std::vector<ClusterGroupMapping>& clusterMap,
    const std::vector<ClusterProperties>& clusters);

// Generate a simple INI config string for a deposit test.
std::string generateDepositIni(const DepositData& data,
                                const std::string& clusterCsvPath,
                                const std::string& outputDir,
                                int stride = 10,
                                int maxIterations = 50,
                                double propertyMismatch = 0.1);

// Perturb physical properties of a lithology model by a random factor.
// Each group's density and susceptibility are multiplied by
// (1 + factor * U(-1,1)), creating a controlled mismatch between
// the true model (used for synthetic data) and the initial model
// (used for inversion).
// @param model  The true model whose properties will be perturbed in-place.
// @param factor Relative perturbation magnitude (e.g., 0.1 = ±10%).
void perturbModelProperties(LithologyModel& model, double factor = 0.1);

} // namespace litho_invert
