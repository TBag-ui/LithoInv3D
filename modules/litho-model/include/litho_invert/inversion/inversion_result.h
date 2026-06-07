#pragma once
#include <litho_invert/core/common.h>
#include <litho_invert/surface/surface_mesh.h>
#include <litho_invert/litho/lithology_model.h>
#include <vector>
#include <memory>

namespace litho_invert {

struct InversionIteration {
    int iteration = 0;
    double dataMisfit = 0.0;
    double regularization = 0.0;
    double constraintPenalty = 0.0;
    double totalObjective = 0.0;
    double rmsError = 0.0;
    double magneticMisfit = 0.0;
    double dw_gravity_x = 2.0;
    double dw_gravity_y = 2.0;
    double dw_magnetic_x = 2.0;
    double dw_magnetic_y = 2.0;
    double activeEMMisfit = 0.0;
    double mtMisfit = 0.0;
    double dw_activeEM_x = 2.0;
    double dw_activeEM_y = 2.0;
    double dw_mt_x = 2.0;
    double dw_mt_y = 2.0;
};

struct InversionResult {
    std::shared_ptr<LithologyModel> finalModel;
    std::vector<InversionIteration> history;
    VectorXd predictedData;
    bool converged = false;
    int totalIterations = 0;
    double finalMisfit = 0.0;
    double finalRMS = 0.0;
    VectorXd finalDensities;
    double finalPaddingDensity = 0.0;
    VectorXd finalSusceptibilities;
    double finalPaddingSusceptibility = 0.0;
    VectorXd finalRemanences;
    std::vector<Vector3d> finalRemanenceVectors;
    VectorXd finalConductivities;
    double finalPaddingConductivity = 1e-4;

    std::shared_ptr<SurfaceMesh> closureTop;
    std::shared_ptr<SurfaceMesh> closureBottom;
    std::shared_ptr<SurfaceMesh> closureDeepBottom;
};

} // namespace litho_invert
