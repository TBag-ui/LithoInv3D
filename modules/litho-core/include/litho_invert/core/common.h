#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <string>
#include <cstdint>

namespace litho_invert {

using Vector3d = Eigen::Vector3d;
using Vector3i = Eigen::Vector3i;
using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;
using Index = Eigen::Index;

struct GravityPoint {
    Vector3d position;  // (x, y, z) observation location
    double g_obs = 0.0;       // observed gravity (mGal)
    double g_std = 0.0;       // standard deviation (mGal), 0 = unweighted

    GravityPoint() = default;
    GravityPoint(const Vector3d& pos, double g_obs_ = 0.0, double g_std_ = 0.0)
        : position(pos), g_obs(g_obs_), g_std(g_std_) {}
};

using GravityData = std::vector<GravityPoint>;

struct MagneticPoint {
    Vector3d position;     // (x, y, z) observation location
    double t_obs = 0.0;    // observed total-field anomaly (nT)
    double t_std = 0.0;    // standard deviation (nT), 0 = unweighted

    MagneticPoint() = default;
    MagneticPoint(const Vector3d& pos, double val = 0.0, double std = 0.0)
        : position(pos), t_obs(val), t_std(std) {}
};

using MagneticData = std::vector<MagneticPoint>;

struct Constraint {
    Vector3d position;          // (x, y) center of constraint column
    double z_top = 0.0;         // top depth (positive down)
    double z_bottom = 0.0;      // bottom depth (positive down)
    int litho_group_id = -1;    // known lithology in this interval
};

enum class TopographyMode {
    None,              // flat z=0 top (current behavior)
    Raw,               // topography included, no correction applied
    TerrainCorrected   // topography included, Bouguer slab subtracted
};

// ---------------------------------------------------------------------------
// PaddingProjection — controls how closure surfaces (top, bottom, deep) are
// built when padding rings extend the mesh beyond the survey area.
//
// IndependentFlat (default):
//   Closure surfaces are built as an independent, perfectly horizontal regular
//   grid whose extent matches the bounding box of the reference surface (the
//   topmost lithological horizon).  The x,y positions of closure-surface
//   vertices are computed from the bounding-box geometry alone — they do NOT
//   follow any topography-induced displacement in the reference surface.
//   Use this when:
//     - Topography is flat or not modelled
//     - You want synthetic data and inversion to agree on a strict horizontal
//       grid for the closure surfaces
//     - You are matching legacy behaviour (pre-padding or early padding builds)
//
// ProjectedBeyondSurvey:
//   Closure surfaces copy the x,y vertex positions from the reference surface
//   (surface[0]), then set z to the appropriate closure depth (topography
//   elevation, bottom depth, or deep-padding depth).  The resulting closure
//   mesh inherits the same triangulation as the reference surface.
//   Use this when:
//     - The reference surface has a non-regular grid (e.g. from an external
//       mesh with irregular triangulation)
//     - You want closure surfaces to exactly follow the same x,y topology as
//       the reference surface
//     - Padding rings were added and the closure surfaces must span the full
//       padded grid
//
// Both modes produce closure surfaces with the same number of vertices and
// triangles as the reference surface, so polyhedron construction is identical
// in the forward models.  The difference is only in the x,y placement of
// closure-surface vertices in the horizontal plane.
// ---------------------------------------------------------------------------
enum class PaddingProjection {
    IndependentFlat,          // closure surfaces are independent horizontal grids
    ProjectedBeyondSurvey     // closure surfaces copy x,y from reference surface
};

struct TopographyConfig {
    TopographyMode mode = TopographyMode::None;
    std::string demFile;
    double datumElevation = 0.0;
    double bouguerDensity = 2.67;   // g/cm³

    // Lateral padding
    int paddingRings = 0;
    double paddingCellSize = 0.0;   // 0 = same as interior

    // Closure surface construction — controls how closure surfaces extend
    // beyond the survey area when padding rings are added.
    //   IndependentFlat       — independent horizontal grid from bounding box
    //   ProjectedBeyondSurvey — copy x,y topology from reference surface
    PaddingProjection paddingProjection = PaddingProjection::IndependentFlat;

    // Half-space property inversion
    bool invertHalfspaceProperties = false;
};

} // namespace litho_invert
