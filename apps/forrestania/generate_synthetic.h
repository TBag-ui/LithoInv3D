#pragma once
#include <litho_invert/litho/lithology_model.h>
#include <litho_invert/core/common.h>
#include "cluster_loader.h"
#include <memory>
#include <vector>
#include <string>

namespace litho_invert {

// Lithology properties derived from geochemical clustering of
// Lithosplitter_out/ drill core data (cluster statistics):
//
//   cluster 0 (MassiveSulfide, n=9):      rho=4.316, chi=0.097
//   cluster 1 (Felsic Granite, n=40):     rho=2.647, chi=0.0005
//   cluster 2 (Intermediate Norite+Gabbro, n=55): rho=2.943, chi=0.0057
//   cluster 3 (Intermediate Troctolite, n=64):    rho=2.952, chi=0.0194
//   cluster 4 (MassiveSulfide, n=12):     rho=4.296, chi=0.0464
//
// Mapped to 4 litho groups (Mafic-on-top with MaficLower repeated):
//   MaficComplex_1 (clusters 2,3):    rho≈2.95  chi≈0.013  (top)
//   GraniteGneiss  (cluster 1):       rho≈2.65  chi≈0.0005
//   MaficComplex_2 (clusters 2,3):    rho≈2.95  chi≈0.013  (lower, same props)
//   MassiveSulfide (clusters 0,4):    rho≈4.31  chi≈0.072  (bottom)

// Generate true synthetic model honoring real borehole constraints from
// Lithosplitter_out/ data for all 5 boreholes.
std::shared_ptr<LithologyModel> generateTrueModel();

// Generate initial flat surfaces for the inversion starting model.
std::shared_ptr<LithologyModel> generateInitialModel();

// Generate interpolated starting model from borehole contact depths.
// Uses inverse-distance-weighted (IDW) interpolation so surfaces already
// approximate the borehole constraints before inversion begins.
// Group ordering (4 groups, 3 surfaces):
//   MaficComplex_1=0 (top), GraniteGneiss=1, MaficComplex_2=2, MassiveSulfide=3 (bottom)
std::shared_ptr<LithologyModel> generateInterpolatedInitialModel(int interiorDim = 21,
                                                                  double cellSize = 100.0,
                                                                  double gridHalf = 1000.0);

// Observation points on regular grid matching survey geometry.
GravityData generateObservationPoints();

// Build a closed triangulated mesh from two surfaces sharing the same topology.
// The top surface faces are kept as-is; bottom faces are reversed; side walls
// connect boundary edges. Used to convert layered (N-1 contact surfaces) to
// volumetric (N closed group meshes) for synthetic tests.
std::shared_ptr<SurfaceMesh> buildClosedMesh(
    const SurfaceMesh& top, const SurfaceMesh& bottom);

// Build a flat triangulated mesh matching the topology of a reference surface.
std::shared_ptr<SurfaceMesh> buildFlatMesh(
    const SurfaceMesh& ref, double z, const std::string& name);

// Compute synthetic gravity data from the true model.
GravityData computeSyntheticData(std::shared_ptr<LithologyModel> trueModel,
                                  const GravityData& observationPoints);

// Compute synthetic magnetic data from the true model.
MagneticData computeSyntheticMagnetic(std::shared_ptr<LithologyModel> trueModel,
                                       const GravityData& observationPoints,
                                       double inc_deg, double dec_deg,
                                       double field_nT,
                                       RemanentMagnetizationMode mode =
                                           RemanentMagnetizationMode::EffectiveSusceptibility);

// Generate true model with realistic remanent magnetization.
// Same geometry as generateTrueModel() but with remnant fields set
// on the magnetic litho groups.  Uses FixedVector mode directions
// representative of oriented core measurements.
//
// MassiveSulfide: Q ≈ 4.5, I_rem = 45°, D_rem = -60° (distinct from induced)
// MaficComplex:    Q ≈ 1.5, I_rem = 65°, D_rem = -30°
// GraniteGneiss:   no remanence (non-magnetic)
std::shared_ptr<LithologyModel> generateTrueModelWithRemanence();

// Build constraint list from real Lithosplitter_out/ cluster intervals.
// Uses all 5 boreholes.  Each constraint encodes (x, y, z_top, z_bottom,
// litho_group_id) where depths are positive-down metres.
//
// Top-of-hole constraint: the first (shallowest) interval from each
// borehole is always included regardless of its length, because the
// surface lithology is considered known from regional mapping.
std::vector<Constraint> generateLithosplitterConstraints();

// Build a list of constraint intervals from BH-*_classified.csv files.
// When groupColumn == "litho_group_id": uses the hardcoded 4-group mapping
// (generateLithosplitterConstraints).
// When groupColumn == "cluster_id": reads cluster_id directly from the
// classified CSV and remaps BH (Lithosplitter) cluster numbering to
// inversion_start cluster numbering via density-matched mapping.
// Constraints where litho_group_id = geological group index.
// BH (Lithosplitter) cluster_id → geo group mapping is computed automatically
// by matching density/susceptibility medians against inversion cluster properties.
std::vector<Constraint> generateClusterIdConstraints(
    const std::vector<ClusterProperties>& invClusters,
    const std::vector<int>& ordering = {});

std::shared_ptr<LithologyModel> generateClusterIdTrueModel(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<int>& ordering = {});

std::vector<int> deriveOrderingFromBoreholes(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<std::string>& classifiedCsvPaths);

std::shared_ptr<LithologyModel> generateLayeredTrueModel(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<double>& surfaceDepths,
    const std::vector<int>& ordering = {});

std::shared_ptr<LithologyModel> generateDippingLayeredTrueModel(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<double>& surfaceDepths,
    double dipAngleDeg,
    double dipDirectionDeg,
    const std::vector<int>& ordering = {});

// Build a regular N×N padded grid SurfaceMesh.
// interiorDim = N (invertible interior), paddingRings = P (extra rings each side).
// Total grid = (N + 2P) × (N + 2P).
// Interior vertices are Z_ONLY, padding vertices are FIXED.
// flatZ = initial z value for all vertices.
// cellSize = spacing between grid vertices.
// gridHalf = half-extent of the INTERIOR region (N-1)*cellSize/2.
std::shared_ptr<SurfaceMesh> makePaddedGrid(
    const std::string& name,
    double flatZ,
    int interiorDim,
    int paddingRings,
    double cellSize,
    double gridHalf,
    VertexFreedom vf = VertexFreedom::Z_ONLY);

} // namespace litho_invert
