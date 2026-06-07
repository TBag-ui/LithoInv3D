#pragma once
#include <litho_invert/core/common.h>
#include <litho_invert/surface/surface_mesh.h>
#include <litho_invert/litho/lithology_model.h>
#include <litho_invert/inversion/inversion_result.h>
#include <string>

namespace litho_invert {

class InversionExporter {
public:
    InversionExporter(const std::string& outputDir, const std::string& baseName);

    // Raw text (written FIRST for safety)
    void exportRawVertices(const SurfaceMesh& mesh, const std::string& suffix);
    void exportRawTriangles(const SurfaceMesh& mesh, const std::string& suffix);

    // OBJ export
    void exportOBJ(const SurfaceMesh& mesh, const std::string& suffix);

    // GOCAD TSurf export (.ts) — triangulated surface readable by Geoscience Analyst
    void exportTS(const SurfaceMesh& mesh, const std::string& suffix);

    // Interior-only exports (truncated to survey area, no padding vertices)
    void exportInteriorTS(const SurfaceMesh& mesh, const std::string& suffix);
    void exportInteriorOBJ(const SurfaceMesh& mesh, const std::string& suffix);

    // UBC-GIF mesh export (3D block model with litho codes)
    void exportUBCGIF(const LithologyModel& model,
                      double xmin, double xmax,
                      double ymin, double ymax,
                      double zmin, double zmax,
                      double cellSize);

    // Litho CSV export (x,y,z,litho_id,name,density)
    void exportLithoCSV(const LithologyModel& model,
                        double xmin, double xmax,
                        double ymin, double ymax,
                        double zmin, double zmax,
                        double cellSize);

    // Predicted vs observed CSV
    void exportPredictedCSV(const GravityData& observed, const VectorXd& predicted);

    // Inversion log (simple text summary)
    void exportLog(const InversionResult& result);

    // Surface CSV export (x, y, z plain text)
    void exportSurfaceCSV(const SurfaceMesh& mesh, const std::string& suffix);

    // Export a closed litho-unit volume as a single GOCAD TSurf.
    // In the unified volumetric architecture, each group mesh is already
    // a closed triangulated boundary — no top/bottom stitching needed.
    void exportClosedVolume(const SurfaceMesh& closedMesh, const std::string& suffix);

    // Export starting model: each group's closed volume mesh as TS and OBJ
    void exportStartingModel(const LithologyModel& model);

    // Export everything
    void exportAll(const InversionResult& result, const GravityData& observed,
                   double xmin, double xmax,
                   double ymin, double ymax,
                   double zmin, double zmax,
                   double cellSize);

    // Set a subfolder prefix (e.g., "starting", "final")
    void setSubfolder(const std::string& sub);

    // Set per-group export names (e.g., "clusters_2_3", "cluster_id_0")
    // When set, contact surfaces use "contact_" + names[i] + "_" + names[i+1]
    // and closed volumes use names[g]. When empty, falls back to legacy naming.
    void setGroupNaming(const std::vector<std::string>& groupExportNames);

private:
    std::string m_outDir;
    std::string m_base;
    std::string m_subfolder;
    std::vector<std::string> m_groupExportNames;
    std::string path(const std::string& suffix, const std::string& ext) const;
};

} // namespace litho_invert
