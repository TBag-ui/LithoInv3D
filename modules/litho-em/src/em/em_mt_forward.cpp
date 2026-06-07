#include <litho_invert/em/em_mt_forward.h>
#include <litho_invert/litho/lithology_model.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace litho_invert {

EMMTForward::EMMTForward(std::shared_ptr<LithologyModel> model,
                         const std::vector<MTStation>& stations,
                         const MTData& data,
                         const EMConfig& config)
    : EMForwardModel(model, config)
    , m_stations(stations)
    , m_data(data)
{
    if (!m_solver) {
        m_solver = createEMSolver(config.solverMethod);
    }
}

size_t EMMTForward::dataCount() const {
    return m_data.size();
}

size_t EMMTForward::parameterCount() const {
    if (!m_model) return 0;
    return static_cast<size_t>(m_model->totalDofCount());
}

VectorXd EMMTForward::compute(const VectorXd& params) {
    if (!m_model) {
        return VectorXd::Zero(static_cast<Index>(m_data.size()));
    }

    // Apply parameters to model geometry
    m_model->applyParameterVector(params);

    // Collect unique station positions and frequencies for solver calls.
    // Group data indices by (stationIndex, frequencyIndex).
    std::map<std::pair<int, int>, std::vector<size_t>> sfGroups;
    for (size_t i = 0; i < m_data.size(); ++i) {
        sfGroups[{m_data[i].stationIndex, m_data[i].frequencyIndex}].push_back(i);
    }

    VectorXd result(static_cast<Index>(m_data.size()));

    // For MT there is no subsetting — the whole mesh contributes to every
    // station.  Plane waves illuminate the entire subsurface.
    if (m_solver && m_solver->supportsMT()) {
        // Collect all unique station positions
        std::vector<Vector3d> stationPositions;
        std::vector<double> frequencies_Hz;
        for (const auto& stn : m_stations) {
            stationPositions.push_back(stn.position);
        }

        // Collect all unique frequencies
        std::set<double> freqSet;
        for (size_t i = 0; i < m_data.size(); ++i) {
            int fi = m_data[i].frequencyIndex;
            for (const auto& stn : m_stations) {
                if (fi >= 0 && fi < static_cast<int>(stn.frequencies_Hz.size())) {
                    freqSet.insert(stn.frequencies_Hz[fi]);
                }
            }
        }
        frequencies_Hz.assign(freqSet.begin(), freqSet.end());

        if (!frequencies_Hz.empty()) {
            VectorXd solverResult = m_solver->solveMT(
                stationPositions, frequencies_Hz, *m_model, m_config);

            // solverResult layout: 8 values per station per frequency
            // Map back to individual data rows
            int nFreq = static_cast<int>(frequencies_Hz.size());
            for (const auto& kvSolver : sfGroups) {
                int stnIdx = kvSolver.first.first;
                int freqIdx = kvSolver.first.second;
                const std::vector<size_t>& dataIndices = kvSolver.second;

                // Find the frequency in our sorted list
                double targetFreq = 0.0;
                if (stnIdx >= 0 && stnIdx < static_cast<int>(m_stations.size()) &&
                    freqIdx >= 0 && freqIdx < static_cast<int>(m_stations[stnIdx].frequencies_Hz.size())) {
                    targetFreq = m_stations[stnIdx].frequencies_Hz[freqIdx];
                }

                int freqPos = -1;
                for (int k = 0; k < nFreq; ++k) {
                    if (std::abs(frequencies_Hz[k] - targetFreq) < 1e-6) {
                        freqPos = k;
                        break;
                    }
                }

                if (freqPos >= 0) {
                    int baseIdx = (stnIdx * nFreq + freqPos) * IMPEDANCE_COMPONENTS;
                    for (size_t di : dataIndices) {
                        int iComp = m_data[di].iComp;  // 0=x, 1=y
                        int jComp = m_data[di].jComp;
                        int compIdx = baseIdx + (iComp * 2 + jComp) * 2;  // real part
                        if (compIdx + 1 < solverResult.size()) {
                            result[static_cast<Index>(di)] = solverResult[compIdx];  // Re(Z)
                        }
                    }
                }
            }
            return result;
        }
    }

    // --- Fallback: simple half-space MT response ---
    // When no solver is available, compute the analytic 1D half-space
    // impedance: Z = sqrt(i·ω·μ₀/σ) for each station/frequency.
    // Zxy = -Zyx = Z (in Ω), Zxx = Zyy = 0.
    //
    // This is dimensionally correct (a homogeneous Earth) but ignores
    // all 3D structure.  The real solver replaces this.
    const double mu0 = 4.0 * M_PI * 1e-7;  // H/m

    for (const auto& kv : sfGroups) {
        int stnIdx = kv.first.first;
        int freqIdx = kv.first.second;
        const std::vector<size_t>& dataIndices = kv.second;

        double freq_Hz = 1.0;  // default
        if (stnIdx >= 0 && stnIdx < static_cast<int>(m_stations.size()) &&
            freqIdx >= 0 && freqIdx < static_cast<int>(m_stations[stnIdx].frequencies_Hz.size())) {
            freq_Hz = m_stations[stnIdx].frequencies_Hz[freqIdx];
        }

        // Use the average crustal conductivity for half-space estimate
        double sigma = 1e-4;  // S/m default
        if (m_model && m_model->groupCount() > 0) {
            sigma = 0.0;
            for (int g = 0; g < m_model->groupCount(); ++g) {
                sigma += m_model->group(g).conductivity;
            }
            sigma /= m_model->groupCount();
            if (sigma <= 0.0) sigma = 1e-4;
        }

        // Z = (1 + i) * sqrt(ω·μ₀ / (2σ))
        double omega = 2.0 * M_PI * freq_Hz;
        double zMag = std::sqrt(omega * mu0 / (2.0 * sigma));  // |Z| in Ω

        for (size_t di : dataIndices) {
            int iComp = m_data[di].iComp;
            int jComp = m_data[di].jComp;

            // Diagonal: Zxx = Zyy = 0 for half-space
            // Off-diagonal: Zxy = -Zyx = Z
            if (iComp == jComp) {
                result[static_cast<Index>(di)] = 0.0;  // real part is 0
            } else if (iComp == 0 && jComp == 1) {
                // Zxy = +Z (real part)
                result[static_cast<Index>(di)] = zMag;
            } else {
                // Zyx = -Z (real part)
                result[static_cast<Index>(di)] = -zMag;
            }
        }
    }

    return result;
}

VectorXd EMMTForward::computeGroupUnitResponse(int groupIndex) {
    if (!m_model || groupIndex < 0 || groupIndex >= m_model->groupCount()) {
        return VectorXd::Zero(static_cast<Index>(m_data.size()));
    }

    // Store original conductivities, set target group to 1.0 S/m
    std::vector<double> originalConductivities(m_model->groupCount());
    for (int g = 0; g < m_model->groupCount(); ++g) {
        originalConductivities[g] = m_model->group(g).conductivity;
        m_model->group(g).conductivity = 0.0;
    }
    m_model->group(groupIndex).conductivity = 1.0;

    VectorXd unitResponse = compute(m_model->assembleParameterVector());

    // Restore
    for (int g = 0; g < m_model->groupCount(); ++g) {
        m_model->group(g).conductivity = originalConductivities[g];
    }

    return unitResponse;
}

VectorXd EMMTForward::computePaddingUnitResponse() const {
    return VectorXd::Zero(static_cast<Index>(m_data.size()));
}

void EMMTForward::enablePadding(bool enabled, double paddingDepth) {
    m_paddingEnabled = enabled;
    m_paddingDepth = paddingDepth;
    if (enabled) {
        buildFlatDeepBottom();
    }
}

void EMMTForward::buildFlatDeepBottom() {
    if (!m_model || m_model->groupMeshCount() == 0) return;
    const SurfaceMesh* deepestSurf = m_model->groupMesh(m_model->groupMeshCount() - 1);
    m_flatDeepBottom = std::make_shared<SurfaceMesh>();
    m_flatDeepBottom->setName("mt_padding_bottom");
    for (uint32_t v = 0; v < deepestSurf->vertexCount(); ++v) {
        Vector3d pos = deepestSurf->vertex(v).position;
        pos.z() = m_paddingDepth;
        m_flatDeepBottom->addVertex(pos);
    }
    for (uint32_t t = 0; t < deepestSurf->triangleCount(); ++t) {
        const Triangle& tri = deepestSurf->triangle(t);
        m_flatDeepBottom->addTriangle(tri.v0, tri.v1, tri.v2);
    }
    m_flatDeepBottom->buildNeighbors();
}

} // namespace litho_invert
