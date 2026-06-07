#include <litho_invert/io/exporters.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <map>
#include <set>
#include <utility>
#include <filesystem>
#include <algorithm>

namespace litho_invert {

InversionExporter::InversionExporter(const std::string& outputDir, const std::string& baseName)
    : m_outDir(outputDir), m_base(baseName)
{
}

std::string InversionExporter::path(const std::string& suffix, const std::string& ext) const {
    std::string dir = m_outDir;
    if (!m_subfolder.empty()) {
        dir += "/" + m_subfolder;
    }
    std::string filename = m_base;
    if (!suffix.empty()) {
        filename += "_" + suffix;
    }
    filename += "." + ext;
    return dir + "/" + filename;
}

void InversionExporter::setSubfolder(const std::string& sub) {
    m_subfolder = sub;
}

void InversionExporter::setGroupNaming(const std::vector<std::string>& groupExportNames) {
    m_groupExportNames = groupExportNames;
}

// -----------------------------------------------------------------------
// Raw vertices
// -----------------------------------------------------------------------
void InversionExporter::exportRawVertices(const SurfaceMesh& mesh, const std::string& suffix) {
    std::string filepath = path(suffix, "txt");
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "exporter: Could not write " << filepath << std::endl;
        return;
    }
    // Write header
    out << "# vertex_index x y z\n";
    for (uint32_t i = 0; i < mesh.vertexCount(); ++i) {
        const Vector3d& p = mesh.vertex(i).position;
        out << i << " "
            << std::fixed << std::setprecision(3) << p.x() << " "
            << std::fixed << std::setprecision(3) << p.y() << " "
            << std::fixed << std::setprecision(3) << p.z() << "\n";
    }
    std::cout << "exporter: Wrote " << mesh.vertexCount() << " vertices to "
              << filepath << std::endl;
}

// -----------------------------------------------------------------------
// Raw triangles
// -----------------------------------------------------------------------
void InversionExporter::exportRawTriangles(const SurfaceMesh& mesh, const std::string& suffix) {
    std::string filepath = path(suffix, "txt");
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "exporter: Could not write " << filepath << std::endl;
        return;
    }
    // Write header
    out << "# triangle_index v0 v1 v2\n";
    for (uint32_t i = 0; i < mesh.triangleCount(); ++i) {
        const Triangle& t = mesh.triangle(i);
        out << i << " " << t.v0 << " " << t.v1 << " " << t.v2 << "\n";
    }
    std::cout << "exporter: Wrote " << mesh.triangleCount() << " triangles to "
              << filepath << std::endl;
}

// -----------------------------------------------------------------------
// OBJ export
// -----------------------------------------------------------------------
void InversionExporter::exportOBJ(const SurfaceMesh& mesh, const std::string& suffix) {
    std::string filepath = path(suffix, "obj");
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "exporter: Could not write " << filepath << std::endl;
        return;
    }
    out << "# LithoInvert3D OBJ export - " << m_base << "_" << suffix << "\n";
    out << "# Vertices: " << mesh.vertexCount()
        << ", Triangles: " << mesh.triangleCount() << "\n";

    for (uint32_t i = 0; i < mesh.vertexCount(); ++i) {
        const Vector3d& p = mesh.vertex(i).position;
        out << "v " << std::fixed << std::setprecision(3)
            << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    for (uint32_t i = 0; i < mesh.triangleCount(); ++i) {
        const Triangle& t = mesh.triangle(i);
        // OBJ uses 1-based indices
        out << "f " << (t.v0 + 1) << " " << (t.v1 + 1) << " " << (t.v2 + 1) << "\n";
    }

    std::cout << "exporter: Wrote OBJ to " << filepath << std::endl;
}

// -----------------------------------------------------------------------
// GOCAD TSurf export (.ts) for Geoscience Analyst
// -----------------------------------------------------------------------
void InversionExporter::exportTS(const SurfaceMesh& mesh, const std::string& suffix) {
    std::string filepath = path(suffix, "ts");
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "exporter: Could not write " << filepath << std::endl;
        return;
    }

    // GOCAD TSurf header
    out << "GOCAD TSurf 1\n";
    out << "HEADER {\n";
    out << "  name: " << m_base << "_" << suffix << "\n";
    out << "}\n";
    out << "TFACE\n";

    // Vertices (1-based indexing for GOCAD)
    for (uint32_t i = 0; i < mesh.vertexCount(); ++i) {
        const Vector3d& p = mesh.vertex(i).position;
        out << "VRTX " << (i + 1) << " "
            << std::scientific << std::setprecision(8)
            << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // Triangles referencing vertices by 1-based index
    for (uint32_t i = 0; i < mesh.triangleCount(); ++i) {
        const Triangle& t = mesh.triangle(i);
        out << "TRGL " << (t.v0 + 1) << " " << (t.v1 + 1) << " " << (t.v2 + 1) << "\n";
    }

    out << "END\n";

    std::cout << "exporter: Wrote GOCAD TSurf to " << filepath
              << " (" << mesh.vertexCount() << " vertices, "
              << mesh.triangleCount() << " triangles)" << std::endl;
}
// -----------------------------------------------------------------------
// UBC-GIF mesh export
// -----------------------------------------------------------------------
void InversionExporter::exportUBCGIF(const LithologyModel& model,
                                     double xmin, double xmax,
                                     double ymin, double ymax,
                                     double zmin, double zmax,
                                     double cellSize)
{
    std::string filepath = path("ubc", "msh");
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "exporter: Could not write " << filepath << std::endl;
        return;
    }

    int nx = static_cast<int>(std::ceil((xmax - xmin) / cellSize));
    int ny = static_cast<int>(std::ceil((ymax - ymin) / cellSize));
    int nz = static_cast<int>(std::ceil((zmax - zmin) / cellSize));

    out << nx << " " << ny << " " << nz << "\n";

    // X edges (easting, increasing)
    out << std::fixed << std::setprecision(1);
    for (int i = 0; i <= nx; ++i) {
        out << (xmin + i * cellSize);
        if (i < nx) out << " ";
    }
    out << "\n";

    // Y edges (northing, increasing)
    for (int i = 0; i <= ny; ++i) {
        out << (ymin + i * cellSize);
        if (i < ny) out << " ";
    }
    out << "\n";

    // Z edges (depth, positive DOWN!)
    // Convert from positive-up (z_up) to positive-down (z_down)
    // z_up = -z_down, so zmin_up = -zmax_down, zmax_up = -zmin_down
    for (int i = 0; i <= nz; ++i) {
        double z_up = zmax - i * cellSize;   // from surface downward
        double z_down = -z_up;                // positive down
        out << std::fixed << std::setprecision(1) << z_down;
        if (i < nz) out << " ";
    }
    out << "\n";

    // Cell litho codes: fast-x, then y, then z
    out << std::setprecision(0);
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                double cx = xmin + (i + 0.5) * cellSize;
                double cy = ymin + (j + 0.5) * cellSize;
                // UBC-GIF z is positive down; model uses positive up
                double cz_up = zmax - (k + 0.5) * cellSize;

                int lithoId = model.classifyPoint(Vector3d(cx, cy, cz_up));
                out << lithoId;
                if (i < nx - 1) out << " ";
            }
            out << "\n";
        }
    }

    std::cout << "exporter: Wrote UBC-GIF mesh to " << filepath
              << " (" << nx << "x" << ny << "x" << nz << " cells)" << std::endl;
}

// -----------------------------------------------------------------------
// Litho CSV export
// -----------------------------------------------------------------------
void InversionExporter::exportLithoCSV(const LithologyModel& model,
                                       double xmin, double xmax,
                                       double ymin, double ymax,
                                       double zmin, double zmax,
                                       double cellSize)
{
    std::string filepath = path("litho", "csv");
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "exporter: Could not write " << filepath << std::endl;
        return;
    }

    int nx = static_cast<int>(std::ceil((xmax - xmin) / cellSize));
    int ny = static_cast<int>(std::ceil((ymax - ymin) / cellSize));
    int nz = static_cast<int>(std::ceil((zmax - zmin) / cellSize));

    // CSV header
    out << "x,y,z,litho_id,name,density\n";

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                double cx = xmin + (i + 0.5) * cellSize;
                double cy = ymin + (j + 0.5) * cellSize;
                double cz_up = zmax - (k + 0.5) * cellSize;

                int lithoId = model.classifyPoint(Vector3d(cx, cy, cz_up));
                std::string name = "";
                double density = 0.0;

                // Find the litho group by id
                for (int g = 0; g < model.groupCount(); ++g) {
                    if (model.group(g).id == lithoId) {
                        name = model.group(g).name;
                        density = model.group(g).density;
                        break;
                    }
                }

                out << std::fixed << std::setprecision(1) << cx << ","
                    << std::fixed << std::setprecision(1) << cy << ","
                    << std::fixed << std::setprecision(1) << cz_up << ","
                    << lithoId << ","
                    << name << ","
                    << std::fixed << std::setprecision(2) << density << "\n";
            }
        }
    }

    std::cout << "exporter: Wrote litho CSV to " << filepath
              << " (" << (nx*ny*nz) << " cells)" << std::endl;
}

// -----------------------------------------------------------------------
// Predicted vs observed CSV
// -----------------------------------------------------------------------
void InversionExporter::exportPredictedCSV(const GravityData& observed,
                                            const VectorXd& predicted)
{
    std::string filepath = path("predicted", "csv");
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "exporter: Could not write " << filepath << std::endl;
        return;
    }

    out << "x,y,z,observed,predicted,residual\n";
    size_t n = std::min(observed.size(), static_cast<size_t>(predicted.size()));
    for (size_t i = 0; i < n; ++i) {
        double res = observed[i].g_obs - predicted[static_cast<Index>(i)];
        out << std::fixed << std::setprecision(2)
            << observed[i].position.x() << ","
            << observed[i].position.y() << ","
            << observed[i].position.z() << ","
            << std::setprecision(6) << observed[i].g_obs << ","
            << std::setprecision(6) << predicted[static_cast<Index>(i)] << ","
            << std::setprecision(6) << res << "\n";
    }

    std::cout << "exporter: Wrote " << n << " predictions to "
              << filepath << std::endl;
}

// -----------------------------------------------------------------------
// Interior-only exports (truncated to survey area)
// -----------------------------------------------------------------------

namespace {

// Build an interior-only submesh from a padded SurfaceMesh.
// Returns a SurfaceMesh with only interior vertices and triangles where
// all three vertices are interior. If no padding is active, returns a
// copy of the full mesh.
std::shared_ptr<SurfaceMesh> buildInteriorSubmesh(const SurfaceMesh& mesh) {
    auto result = std::make_shared<SurfaceMesh>();
    result->setName(mesh.name());

    if (mesh.paddingRings() <= 0) {
        // No padding — just copy all vertices and triangles
        for (uint32_t i = 0; i < mesh.vertexCount(); ++i) {
            const auto& v = mesh.vertex(i);
            result->addVertex(v.position, VertexFreedom::FIXED);
        }
        for (uint32_t t = 0; t < mesh.triangleCount(); ++t) {
            const Triangle& tri = mesh.triangle(t);
            result->addTriangle(tri.v0, tri.v1, tri.v2);
        }
        return result;
    }

    // Build vertex index remapping: old index → new index (or -1 if padding)
    std::vector<int> remap(mesh.vertexCount(), -1);
    for (uint32_t vi = 0; vi < mesh.vertexCount(); ++vi) {
        if (mesh.isInteriorVertex(vi)) {
            remap[vi] = static_cast<int>(result->vertexCount());
            const auto& v = mesh.vertex(vi);
            result->addVertex(v.position, VertexFreedom::FIXED);
        }
    }

    // Keep only triangles where all three vertices are interior
    for (uint32_t t = 0; t < mesh.triangleCount(); ++t) {
        const Triangle& tri = mesh.triangle(t);
        int nv0 = remap[tri.v0];
        int nv1 = remap[tri.v1];
        int nv2 = remap[tri.v2];
        if (nv0 >= 0 && nv1 >= 0 && nv2 >= 0) {
            result->addTriangle(static_cast<uint32_t>(nv0),
                               static_cast<uint32_t>(nv1),
                               static_cast<uint32_t>(nv2));
        }
    }

    return result;
}

} // namespace

void InversionExporter::exportInteriorTS(const SurfaceMesh& mesh, const std::string& suffix) {
    auto interior = buildInteriorSubmesh(mesh);
    exportTS(*interior, suffix);
}

void InversionExporter::exportInteriorOBJ(const SurfaceMesh& mesh, const std::string& suffix) {
    auto interior = buildInteriorSubmesh(mesh);
    exportOBJ(*interior, suffix);
}

// -----------------------------------------------------------------------
// Inversion log
// -----------------------------------------------------------------------
void InversionExporter::exportLog(const InversionResult& result) {
    std::string filepath = path("log", "txt");
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "exporter: Could not write " << filepath << std::endl;
        return;
    }

    out << "=== LithoInvert3D Inversion Log ===\n";
    out << "Converged: " << (result.converged ? "yes" : "no") << "\n";
    out << "Total iterations: " << result.totalIterations << "\n";
    out << "Final objective value: " << result.finalMisfit << "\n";
    out << "Final RMS error (mGal): " << result.finalRMS << "\n";
    out << "History entries: " << result.history.size() << "\n";
    out << "\n";

    // Iteration history
    out << "Iteration history:\n";
    out << "iter\tdata_misfit\tmag_misfit\tregul\tconstraint\ttotal\trms\t"
        << "dw_gx\tdw_gy\tdw_mx\tdw_my\n";
    for (const auto& h : result.history) {
        out << h.iteration << "\t"
            << h.dataMisfit << "\t"
            << h.magneticMisfit << "\t"
            << h.regularization << "\t"
            << h.constraintPenalty << "\t"
            << h.totalObjective << "\t"
            << h.rmsError << "\t"
            << h.dw_gravity_x << "\t"
            << h.dw_gravity_y << "\t"
            << h.dw_magnetic_x << "\t"
            << h.dw_magnetic_y << "\n";
    }

    out << "\nGroup meshes:\n";
    if (result.finalModel) {
        for (int i = 0; i < result.finalModel->groupMeshCount(); ++i) {
            const SurfaceMesh* surf = result.finalModel->groupMesh(i);
            out << "Group " << i << " (" << surf->name() << "): "
                << surf->vertexCount() << " vertices, "
                << surf->triangleCount() << " triangles\n";
        }
    }

    std::cout << "exporter: Wrote inversion log to " << filepath << std::endl;
}

// -----------------------------------------------------------------------
// Surface CSV export (x, y, z plain text)
// -----------------------------------------------------------------------

void InversionExporter::exportSurfaceCSV(const SurfaceMesh& mesh,
                                          const std::string& suffix) {
    std::string filepath = path(suffix, "csv");
    std::ofstream ofs(filepath);
    if (!ofs) {
        std::cerr << "exporter: Cannot open " << filepath << std::endl;
        return;
    }

    // Header
    ofs << "x,y,z\n";

    const uint32_t nVerts = mesh.vertexCount();
    for (uint32_t i = 0; i < nVerts; ++i) {
        const Vector3d& p = mesh.vertex(i).position;
        ofs << std::fixed << std::setprecision(3)
            << p.x() << ',' << p.y() << ',' << p.z() << '\n';
    }

    std::cout << "exporter: Wrote " << nVerts << " vertices to "
              << filepath << std::endl;
}

// -----------------------------------------------------------------------
// Closed volume export
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// Grid resampling — projects any irregular surface onto a regular XY grid
// so all closed-volume surfaces share identical topology.
// -----------------------------------------------------------------------

static bool pointInTriangleXY(
    double px, double py,
    double x0, double y0,
    double x1, double y1,
    double x2, double y2,
    double& u, double& v, double& w)
{
    double d00 = x1 - x0, d01 = x2 - x0;
    double d10 = y1 - y0, d11 = y2 - y0;
    double det = d00 * d11 - d01 * d10;
    if (std::abs(det) < 1e-30) return false;

    double dx = px - x0, dy = py - y0;
    v = (dx * d11 - dy * d01) / det;
    w = (dy * d00 - dx * d10) / det;
    u = 1.0 - v - w;

    return (u >= -1e-10 && v >= -1e-10 && w >= -1e-10);
}

static double interpolateZ(const SurfaceMesh& src, double px, double py, double fallbackZ)
{
    double bestDist = 1e30;
    double bestZ = fallbackZ;

    for (uint32_t t = 0; t < src.triangleCount(); ++t) {
        const auto& tri = src.triangle(t);
        const auto& p0 = src.vertex(tri.v0).position;
        const auto& p1 = src.vertex(tri.v1).position;
        const auto& p2 = src.vertex(tri.v2).position;

        double u, v, w;
        if (pointInTriangleXY(px, py,
                              p0.x(), p0.y(),
                              p1.x(), p1.y(),
                              p2.x(), p2.y(), u, v, w)) {
            return u * p0.z() + v * p1.z() + w * p2.z();
        }

        double d0 = (px - p0.x()) * (px - p0.x()) + (py - p0.y()) * (py - p0.y());
        double d1 = (px - p1.x()) * (px - p1.x()) + (py - p1.y()) * (py - p1.y());
        double d2 = (px - p2.x()) * (px - p2.x()) + (py - p2.y()) * (py - p2.y());
        if (d0 < bestDist) { bestDist = d0; bestZ = p0.z(); }
        if (d1 < bestDist) { bestDist = d1; bestZ = p1.z(); }
        if (d2 < bestDist) { bestDist = d2; bestZ = p2.z(); }
    }

    // Also check vertices directly (handles vertex-only surfaces)
    for (uint32_t i = 0; i < src.vertexCount(); ++i) {
        const auto& p = src.vertex(i).position;
        double d = (px - p.x()) * (px - p.x()) + (py - p.y()) * (py - p.y());
        if (d < bestDist) { bestDist = d; bestZ = p.z(); }
    }

    return bestZ;
}

static std::shared_ptr<SurfaceMesh>
resampleToGrid(const SurfaceMesh& src,
               double xmin, double xmax, int nx,
               double ymin, double ymax, int ny,
               double fallbackZ = 0.0)
{
    auto grid = std::make_shared<SurfaceMesh>();

    double dx = (nx > 1) ? (xmax - xmin) / (nx - 1) : 0.0;
    double dy = (ny > 1) ? (ymax - ymin) / (ny - 1) : 0.0;

    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            double x = xmin + ix * dx;
            double y = ymin + iy * dy;
            double z = interpolateZ(src, x, y, fallbackZ);
            grid->addVertex(x, y, z);
        }
    }

    for (int iy = 0; iy < ny - 1; ++iy) {
        for (int ix = 0; ix < nx - 1; ++ix) {
            uint32_t i0 = static_cast<uint32_t>(iy * nx + ix);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = static_cast<uint32_t>((iy + 1) * nx + ix);
            uint32_t i3 = i2 + 1;
            grid->addTriangle(i0, i1, i2);
            grid->addTriangle(i1, i3, i2);
        }
    }

    return grid;
}

// Build a flat regular-grid surface (no source surface to interpolate from)
static std::shared_ptr<SurfaceMesh>
buildFlatGrid(double xmin, double xmax, int nx,
              double ymin, double ymax, int ny,
              double zValue)
{
    auto grid = std::make_shared<SurfaceMesh>();

    double dx = (nx > 1) ? (xmax - xmin) / (nx - 1) : 0.0;
    double dy = (ny > 1) ? (ymax - ymin) / (ny - 1) : 0.0;

    for (int iy = 0; iy < ny; ++iy) {
        for (int ix = 0; ix < nx; ++ix) {
            grid->addVertex(xmin + ix * dx, ymin + iy * dy, zValue);
        }
    }

    for (int iy = 0; iy < ny - 1; ++iy) {
        for (int ix = 0; ix < nx - 1; ++ix) {
            uint32_t i0 = static_cast<uint32_t>(iy * nx + ix);
            uint32_t i1 = i0 + 1;
            uint32_t i2 = static_cast<uint32_t>((iy + 1) * nx + ix);
            uint32_t i3 = i2 + 1;
            grid->addTriangle(i0, i1, i2);
            grid->addTriangle(i1, i3, i2);
        }
    }

    return grid;
}

// -----------------------------------------------------------------------
// Boundary loop extraction
// -----------------------------------------------------------------------

using Edge = std::pair<uint32_t, uint32_t>;
static Edge makeEdge(uint32_t a, uint32_t b) {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

static std::vector<std::vector<uint32_t>>
extractBoundaryLoops(const SurfaceMesh& mesh) {
    // Count how many triangles share each edge. Boundary edges appear once.
    std::map<Edge, int> edgeCount;
    for (uint32_t t = 0; t < mesh.triangleCount(); ++t) {
        const auto& tri = mesh.triangle(t);
        edgeCount[makeEdge(tri.v0, tri.v1)]++;
        edgeCount[makeEdge(tri.v1, tri.v2)]++;
        edgeCount[makeEdge(tri.v2, tri.v0)]++;
    }

    // Build adjacency from boundary edges only
    std::map<uint32_t, std::vector<uint32_t>> adj;
    for (const auto& [e, count] : edgeCount) {
        if (count == 1) {
            adj[e.first].push_back(e.second);
            adj[e.second].push_back(e.first);
        }
    }

    // Walk connected components of the boundary graph into ordered loops
    std::set<uint32_t> visited;
    std::vector<std::vector<uint32_t>> loops;
    for (const auto& [v, nbs] : adj) {
        if (visited.count(v) || nbs.size() != 2) continue;
        std::vector<uint32_t> loop;
        uint32_t curr = v;
        uint32_t prev = UINT32_MAX;
        while (!visited.count(curr)) {
            loop.push_back(curr);
            visited.insert(curr);
            uint32_t next = UINT32_MAX;
            for (uint32_t nb : adj[curr]) {
                if (nb != prev) { next = nb; break; }
            }
            if (next == UINT32_MAX) break;
            prev = curr;
            curr = next;
        }
        if (loop.size() >= 3) loops.push_back(std::move(loop));
    }
    return loops;
}

static void stitchLoop(
    std::ofstream& out,
    const std::vector<uint32_t>& topLoop,
    const std::vector<uint32_t>& botLoop,
    uint32_t topOff, uint32_t botOff,
    const SurfaceMesh& topSurf,
    const SurfaceMesh& botSurf)
{
    size_t nT = topLoop.size(), nB = botLoop.size();
    if (nT < 2 || nB < 2) return;

    size_t i = 0, j = 0;
    size_t steps = 0, maxSteps = nT + nB + 2;
    while ((i < nT || j < nB) && steps < maxSteps) {
        ++steps;
        uint32_t t0 = topLoop[i];
        uint32_t t1 = topLoop[(i + 1) % nT];
        uint32_t b0 = botLoop[j];
        uint32_t b1 = botLoop[(j + 1) % nB];

        double dx_t1b0 = topSurf.vertex(t1).position.x() - botSurf.vertex(b0).position.x();
        double dy_t1b0 = topSurf.vertex(t1).position.y() - botSurf.vertex(b0).position.y();
        double d1 = dx_t1b0 * dx_t1b0 + dy_t1b0 * dy_t1b0;

        double dx_t0b1 = topSurf.vertex(t0).position.x() - botSurf.vertex(b1).position.x();
        double dy_t0b1 = topSurf.vertex(t0).position.y() - botSurf.vertex(b1).position.y();
        double d2 = dx_t0b1 * dx_t0b1 + dy_t0b1 * dy_t0b1;

        bool advanceTop;
        if (i >= nT - 1 && j >= nB - 1) break;
        else if (i >= nT - 1) advanceTop = false;
        else if (j >= nB - 1) advanceTop = true;
        else advanceTop = (d1 < d2);

        if (advanceTop) {
            // Triangle: top[i] → top[i+1] → bot[j]
            out << "TRGL " << (t0 + 1 + topOff) << " "
                << (t1 + 1 + topOff) << " "
                << (b0 + 1 + botOff) << "\n";
            ++i;
        } else {
            // Triangle: top[i] → bot[j] → bot[j+1]
            out << "TRGL " << (t0 + 1 + topOff) << " "
                << (b0 + 1 + botOff) << " "
                << (b1 + 1 + botOff) << "\n";
            ++j;
        }
    }
}

// -----------------------------------------------------------------------
// Closed volume export helpers
// -----------------------------------------------------------------------

static void writeClosedVolumeTS(
    std::ofstream& out,
    const std::string& name,
    const SurfaceMesh& topSurf,
    const SurfaceMesh& botSurf)
{
    const uint32_t nVertsTop = topSurf.vertexCount();
    const uint32_t nTrisTop  = topSurf.triangleCount();
    const uint32_t nVertsBot = botSurf.vertexCount();
    const uint32_t nTrisBot  = botSurf.triangleCount();

    out << "GOCAD TSurf 1\n";
    out << "HEADER {\n";
    out << "  name: " << name << "\n";
    out << "}\n";
    out << "TFACE\n";

    // Top surface vertices (indices 1..nVertsTop)
    for (uint32_t i = 0; i < nVertsTop; ++i) {
        const auto& p = topSurf.vertex(i).position;
        out << "VRTX " << (i + 1) << " "
            << std::scientific << std::setprecision(8)
            << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // Bottom surface vertices (indices nVertsTop+1..nVertsTop+nVertsBot)
    for (uint32_t i = 0; i < nVertsBot; ++i) {
        const auto& p = botSurf.vertex(i).position;
        out << "VRTX " << (nVertsTop + i + 1) << " "
            << std::scientific << std::setprecision(8)
            << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    // Top triangles (original winding — normals up)
    for (uint32_t i = 0; i < nTrisTop; ++i) {
        const auto& t = topSurf.triangle(i);
        out << "TRGL " << (t.v0 + 1) << " " << (t.v1 + 1) << " " << (t.v2 + 1) << "\n";
    }

    // Bottom triangles (reversed winding — normals down)
    for (uint32_t i = 0; i < nTrisBot; ++i) {
        const auto& t = botSurf.triangle(i);
        out << "TRGL " << (t.v0 + 1 + nVertsTop) << " "
            << (t.v2 + 1 + nVertsTop) << " "
            << (t.v1 + 1 + nVertsTop) << "\n";
    }

    // Wall triangles — stitch boundary loops from top and bottom
    auto topLoops = extractBoundaryLoops(topSurf);
    auto botLoops = extractBoundaryLoops(botSurf);

    if (!topLoops.empty() && !botLoops.empty()) {
        // Match loops by XY centroid proximity
        for (auto& tl : topLoops) {
            double cxT = 0, cyT = 0;
            for (uint32_t vi : tl) {
                cxT += topSurf.vertex(vi).position.x();
                cyT += topSurf.vertex(vi).position.y();
            }
            cxT /= tl.size(); cyT /= tl.size();

            double bestDist = 1e30;
            size_t best = 0;
            for (size_t k = 0; k < botLoops.size(); ++k) {
                double cxB = 0, cyB = 0;
                for (uint32_t vi : botLoops[k]) {
                    cxB += botSurf.vertex(vi).position.x();
                    cyB += botSurf.vertex(vi).position.y();
                }
                cxB /= botLoops[k].size(); cyB /= botLoops[k].size();
                double dx = cxT - cxB, dy = cyT - cyB;
                double dist = dx * dx + dy * dy;
                if (dist < bestDist) { bestDist = dist; best = k; }
            }
            stitchLoop(out, tl, botLoops[best], 0, nVertsTop, topSurf, botSurf);
        }
    }

    out << "END\n";
}

void InversionExporter::exportClosedVolume(const SurfaceMesh& closedMesh,
                                            const std::string& suffix) {
    if (closedMesh.vertexCount() == 0) return;

    std::string filepath = path(suffix, "ts");

    if (!m_subfolder.empty()) {
        std::string dir = m_outDir + "/" + m_subfolder;
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
    }

    // The mesh is already a closed triangulated boundary — write directly
    std::ofstream out(filepath);
    if (!out.is_open()) {
        std::cerr << "exporter: Could not write " << filepath << std::endl;
        return;
    }

    out << "GOCAD TSurf 1\n";
    out << "HEADER {\n";
    out << "  name: " << m_base << "_" << suffix << "\n";
    out << "}\n";
    out << "TFACE\n";

    for (uint32_t i = 0; i < closedMesh.vertexCount(); ++i) {
        const auto& p = closedMesh.vertex(i).position;
        out << "VRTX " << (i + 1) << " "
            << std::scientific << std::setprecision(8)
            << p.x() << " " << p.y() << " " << p.z() << "\n";
    }

    for (uint32_t i = 0; i < closedMesh.triangleCount(); ++i) {
        const auto& t = closedMesh.triangle(i);
        out << "TRGL " << (t.v0 + 1) << " " << (t.v1 + 1) << " " << (t.v2 + 1) << "\n";
    }

    out << "END\n";

    std::cout << "exporter: Wrote closed volume GOCAD TSurf to " << filepath
              << " (" << closedMesh.vertexCount() << " vertices, "
              << closedMesh.triangleCount() << " triangles)" << std::endl;
}

void InversionExporter::exportStartingModel(const LithologyModel& model) {
    // Ensure output subdirectory exists
    if (!m_subfolder.empty()) {
        std::string dir = m_outDir + "/" + m_subfolder;
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
    }

    int nGroups = model.groupCount();
    bool useNaming = !m_groupExportNames.empty()
                     && static_cast<int>(m_groupExportNames.size()) == nGroups;

    // Export each group's closed volume mesh directly
    for (int g = 0; g < model.groupMeshCount(); ++g) {
        const SurfaceMesh* mesh = model.groupMesh(g);
        if (!mesh || mesh->vertexCount() == 0) continue;

        std::string gname;
        if (useNaming) {
            gname = m_groupExportNames[g];
        } else {
            gname = "group_" + std::to_string(g) + "_" + model.group(g).name;
            for (auto& ch : gname) {
                if (ch == ' ' || ch == '-' || ch == '.' || ch == '(' || ch == ')')
                    ch = '_';
            }
        }
        exportTS(*mesh, gname);
        exportOBJ(*mesh, gname);
        exportClosedVolume(*mesh, gname);
    }
}

// -----------------------------------------------------------------------
// Export all
// -----------------------------------------------------------------------
void InversionExporter::exportAll(const InversionResult& result,
                                   const GravityData& observed,
                                   double xmin, double xmax,
                                   double ymin, double ymax,
                                   double zmin, double zmax,
                                   double cellSize)
{
    if (!result.finalModel) {
        std::cerr << "exporter: No model in result, skipping exportAll" << std::endl;
        return;
    }

    // Ensure output subdirectory exists (all exports below need it)
    if (!m_subfolder.empty()) {
        std::string dir = m_outDir + "/" + m_subfolder;
        std::error_code ec;
        std::filesystem::create_directories(dir, ec);
    }

    int nGroups = result.finalModel->groupCount();
    bool useNaming = !m_groupExportNames.empty()
                     && static_cast<int>(m_groupExportNames.size()) == nGroups;

    // Export each group mesh — interior (survey-area) and full-model (debug)
    for (int i = 0; i < result.finalModel->groupMeshCount(); ++i) {
        const SurfaceMesh* surf = result.finalModel->groupMesh(i);
        std::string surfName;
        if (useNaming && i < nGroups) {
            surfName = m_groupExportNames[i];
        } else {
            surfName = surf->name().empty()
                ? ("group_" + std::to_string(i))
                : surf->name();
        }

        // Survey-area exports (interior only, production)
        if (surf->paddingRings() > 0) {
            exportInteriorTS(*surf, surfName);
            exportInteriorOBJ(*surf, surfName);
            // Full-model debug exports with _full suffix
            exportTS(*surf, surfName + "_full");
            exportOBJ(*surf, surfName + "_full");
        } else {
            // No padding — same as before
            exportTS(*surf, surfName);
            exportOBJ(*surf, surfName);
        }

        exportRawVertices(*surf, surfName + "_vertices");
        exportRawTriangles(*surf, surfName + "_triangles");
    }

    exportUBCGIF(*result.finalModel, xmin, xmax, ymin, ymax, zmin, zmax, cellSize);
    exportLithoCSV(*result.finalModel, xmin, xmax, ymin, ymax, zmin, zmax, cellSize);
    exportPredictedCSV(observed, result.predictedData);
    exportLog(result);

    // Export closed volume per litho group — each group mesh is already a closed
    // triangulated boundary from the pipeline, no resampling/wall-stitching needed.
    for (int g = 0; g < nGroups; ++g) {
        const SurfaceMesh* mesh = result.finalModel->groupMesh(g);
        if (!mesh || mesh->vertexCount() == 0) continue;

        std::string gname;
        if (useNaming) {
            gname = m_groupExportNames[g];
        } else {
            gname = "group_" + std::to_string(g) + "_"
                   + result.finalModel->group(g).name;
            for (auto& ch : gname) {
                if (ch == ' ' || ch == '-' || ch == '.' || ch == '(' || ch == ')')
                    ch = '_';
            }
        }
        exportClosedVolume(*mesh, gname);
    }
}

} // namespace litho_invert
