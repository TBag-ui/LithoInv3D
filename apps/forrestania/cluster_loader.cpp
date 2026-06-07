#include "cluster_loader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <unordered_map>
#include <set>
#include <limits>

namespace litho_invert {

// Helper: trim whitespace and quotes
static std::string trimToken(const std::string& raw) {
    size_t start = raw.find_first_not_of(" \t\r\"");
    size_t end = raw.find_last_not_of(" \t\r\"");
    if (start == std::string::npos) return "";
    return raw.substr(start, end - start + 1);
}

// Parse a token as double, returning NaN for "n/a" or empty.
static double parseDoubleOrNaN(const std::string& token) {
    std::string t = trimToken(token);
    if (t.empty() || t == "n/a" || t == "NA" || t == "N/A") {
        return std::numeric_limits<double>::quiet_NaN();
    }
    try {
        return std::stod(t);
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

// Parse a token as yes/no bool.
static bool parseYesNo(const std::string& token) {
    std::string t = trimToken(token);
    return t == "yes" || t == "Yes" || t == "YES" || t == "true" || t == "True" || t == "1";
}

std::vector<ClusterProperties> loadClusterProperties(const std::string& filepath) {
    std::vector<ClusterProperties> clusters;

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "loadClusterProperties: Could not open " << filepath << std::endl;
        return clusters;
    }

    std::string line;
    int lineNum = 0;
    while (std::getline(file, line)) {
        ++lineNum;

        if (line.empty() || line[0] == '#') continue;

        // Split line into tokens
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;
        while (std::getline(ss, token, ',')) {
            tokens.push_back(trimToken(token));
        }

        // Skip header
        if (lineNum == 1 && !tokens.empty() && tokens[0] == "cluster_id") {
            continue;
        }

        if (tokens.size() < 11) {
            std::cerr << "loadClusterProperties: Skipping line " << lineNum
                      << ": expected at least 11 columns, got " << tokens.size() << std::endl;
            continue;
        }

        ClusterProperties cp;
        try {
            cp.cluster_id = std::stoi(tokens[0]);
        } catch (...) {
            std::cerr << "loadClusterProperties: Skipping line " << lineNum
                      << ": bad cluster_id" << std::endl;
            continue;
        }

        cp.working_name = tokens[1];
        cp.sample_count = std::stoi(tokens[2]);

        cp.density_median = std::stod(tokens[3]);
        cp.density_p10   = std::stod(tokens[4]);
        cp.density_p90   = std::stod(tokens[5]);

        cp.susceptibility_median = std::stod(tokens[6]);
        cp.susceptibility_p10   = std::stod(tokens[7]);
        cp.susceptibility_p90   = std::stod(tokens[8]);

        cp.has_measured_density        = parseYesNo(tokens[9]);
        cp.has_measured_susceptibility = parseYesNo(tokens[10]);

        if (tokens.size() >= 12)
            cp.remanence_magnitude    = parseDoubleOrNaN(tokens[11]);
        if (tokens.size() >= 13)
            cp.remanence_inclination  = parseDoubleOrNaN(tokens[12]);
        if (tokens.size() >= 14)
            cp.remanence_declination  = parseDoubleOrNaN(tokens[13]);

        if (tokens.size() >= 15) {
            cp.has_measured_remanence = parseYesNo(tokens[14]);
        }

        clusters.push_back(cp);
    }

    std::cout << "loadClusterProperties: loaded " << clusters.size()
              << " clusters from " << filepath << std::endl;
    return clusters;
}

// ---------------------------------------------------------------------------
// Build LithologyModel from cluster properties + cluster→group mapping
// ---------------------------------------------------------------------------
std::shared_ptr<LithologyModel> buildModelFromClusters(
    const std::vector<ClusterProperties>& clusters,
    const std::vector<ClusterGroupMapping>& mapping,
    double bottomDepth) {

    // Build reverse map: litho_group_id → vector of cluster indices
    std::unordered_map<int, std::vector<size_t>> groupToClusters;
    for (const auto& m : mapping) {
        // Find the cluster index in clusters vector
        for (size_t i = 0; i < clusters.size(); ++i) {
            if (clusters[i].cluster_id == m.cluster_id) {
                groupToClusters[m.litho_group_id].push_back(i);
                break;
            }
        }
    }

    // Determine max group id
    int maxGroup = -1;
    for (const auto& m : mapping) {
        if (m.litho_group_id > maxGroup) maxGroup = m.litho_group_id;
    }

    auto model = std::make_shared<LithologyModel>();
    model->setBottomDepth(bottomDepth);

    // Helper: average properties for a set of cluster indices
    auto avgProps = [&](const std::vector<size_t>& cis) -> std::tuple<double,double,double,double,double,std::string> {
        double r=0,c=0,mr=0,inc=0,dec=0; int nr=0,nc=0,nmr=0,ninc=0,ndec=0;
        std::string name;
        for (size_t ci : cis) {
            const auto& cp = clusters[ci];
            if (name.empty()) name = cp.working_name;
            r += cp.density_median; nr++;
            c += cp.susceptibility_median; nc++;
            if (!std::isnan(cp.remanence_magnitude))    { mr += cp.remanence_magnitude; nmr++; }
            if (!std::isnan(cp.remanence_inclination))  { inc += cp.remanence_inclination; ninc++; }
            if (!std::isnan(cp.remanence_declination))  { dec += cp.remanence_declination; ndec++; }
        }
        return {nr>0?r/nr:2.67, nc>0?c/nc:0.0, nmr>0?mr/nmr:std::numeric_limits<double>::quiet_NaN(),
                ninc>0?inc/ninc:std::numeric_limits<double>::quiet_NaN(),
                ndec>0?dec/ndec:std::numeric_limits<double>::quiet_NaN(), name};
    };

    for (int gid = 0; gid <= maxGroup; ++gid) {
        auto it = groupToClusters.find(gid);
        auto [rho, chi, mrem, inc, dec, name] =
            (it != groupToClusters.end() && !it->second.empty())
            ? avgProps(it->second)
            : [&]() -> std::tuple<double,double,double,double,double,std::string> {
                // Interpolate from nearest neighboring groups
                double rAbove=0,cAbove=0; bool foundAbove=false;
                double rBelow=0,cBelow=0; bool foundBelow=false;
                for (int ga=gid-1; ga>=0; --ga) {
                    auto ita=groupToClusters.find(ga);
                    if (ita!=groupToClusters.end()&&!ita->second.empty()) {
                        auto [ra,ca,_,__,___,____]=avgProps(ita->second);
                        rAbove=ra; cAbove=ca; foundAbove=true; break;
                    }
                }
                for (int gb=gid+1; gb<=maxGroup; ++gb) {
                    auto itb=groupToClusters.find(gb);
                    if (itb!=groupToClusters.end()&&!itb->second.empty()) {
                        auto [rb,cb,_,__,___,____]=avgProps(itb->second);
                        rBelow=rb; cBelow=cb; foundBelow=true; break;
                    }
                }
                double ri, ci;
                if (foundAbove&&foundBelow)         { ri=(rAbove+rBelow)/2.0; ci=(cAbove+cBelow)/2.0; }
                else if (foundAbove)                { ri=rAbove*1.05; ci=cAbove*1.5; }
                else if (foundBelow)                { ri=rBelow*0.95; ci=cBelow*0.5; }
                else                                { ri=2.67; ci=0.0; }
                std::cerr << "buildModelFromClusters: no clusters map to group " << gid
                          << " — interpolated rho=" << ri << " chi=" << ci << std::endl;
                return {ri, ci, std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN(),
                        std::numeric_limits<double>::quiet_NaN(),
                        std::string("Interpolated_") + std::to_string(gid)};
            }();

        LithoGroup group(gid, name, rho, chi);
        group.remanence_magnitude    = mrem;
        group.remanence_inclination  = inc;
        group.remanence_declination  = dec;
        model->addGroup(group);
    }

    std::cout << "buildModelFromClusters: " << model->groupCount()
              << " groups from " << clusters.size() << " clusters" << std::endl;
    for (int g = 0; g < model->groupCount(); ++g) {
        const auto& grp = model->group(g);
        std::cout << "  [" << grp.id << "] " << grp.name
                  << ": rho=" << grp.density
                  << "  chi=" << grp.susceptibility
                  << "  M_rem=" << grp.remanence_magnitude
                  << " A/m  I=" << grp.remanence_inclination
                  << "  D=" << grp.remanence_declination << std::endl;
    }

    return model;
}

// ---------------------------------------------------------------------------
// Configure remanence mode from cluster properties
// ---------------------------------------------------------------------------
bool configureRemanenceFromClusters(
    std::shared_ptr<LithologyModel> model,
    const std::vector<ClusterProperties>& clusters,
    const std::vector<ClusterGroupMapping>& mapping,
    RemanentMagnetizationMode& outMode) {

    constexpr double FIELD_NT = 55000.0;
    constexpr double FIELD_INC = 75.0;
    constexpr double FIELD_DEC = -20.0;

    // Determine per-group: any known / any unknown remanence
    int nGroups = model->groupCount();
    std::vector<bool> groupHasKnown(nGroups, false);
    std::vector<bool> groupHasUnknown(nGroups, false);

    for (const auto& m : mapping) {
        if (m.litho_group_id < 0 || m.litho_group_id >= nGroups) continue;
        const ClusterProperties* cp = nullptr;
        for (const auto& c : clusters) {
            if (c.cluster_id == m.cluster_id) { cp = &c; break; }
        }
        if (!cp) continue;

        if (cp->has_measured_remanence) {
            groupHasKnown[m.litho_group_id] = true;
        } else {
            groupHasUnknown[m.litho_group_id] = true;
        }
    }

    bool anyUnknownRemanence = false;
    bool anyKnownRemanence = false;
    for (int g = 0; g < nGroups; ++g) {
        if (groupHasKnown[g]) anyKnownRemanence = true;
        if (groupHasUnknown[g]) anyUnknownRemanence = true;
    }

    if (anyUnknownRemanence) {
        outMode = RemanentMagnetizationMode::FixedVectorPerGroup;

        // Set initial M_rem guess for groups where remanence is unknown.
        // Use group's averaged susceptibility (already set by buildModelFromClusters)
        // with Q≈1 as a cautious initial guess.
        for (int g = 0; g < nGroups; ++g) {
            auto& grp = model->group(g);
            if (groupHasUnknown[g] && grp.susceptibility > 0.0) {
                // M_equiv (A/m) = chi (SI) * F (nT) / 100
                double m_ind = grp.susceptibility * FIELD_NT / 100.0;
                grp.remanence_magnitude = m_ind; // Q≈1 initial guess

                if (std::isnan(grp.remanence_inclination)) {
                    grp.remanence_inclination = FIELD_INC;
                }
                if (std::isnan(grp.remanence_declination)) {
                    grp.remanence_declination = FIELD_DEC;
                }
            } else if (groupHasUnknown[g] && grp.susceptibility == 0.0) {
                // Non-magnetic group — remanence is zero
                grp.remanence_magnitude = 0.0;
                grp.remanence_inclination = FIELD_INC;
                grp.remanence_declination = FIELD_DEC;
            }
        }

        std::cout << "configureRemanence: FixedVectorPerGroup mode "
                  << "(known=" << (anyKnownRemanence ? "mixed" : "none")
                  << " unknown=" << (anyUnknownRemanence ? "mixed" : "none") << ")"
                  << std::endl;
    } else if (anyKnownRemanence) {
        outMode = RemanentMagnetizationMode::FixedVectorPerGroup;
        std::cout << "configureRemanence: FixedVectorPerGroup (all prescribed)" << std::endl;
    } else {
        outMode = RemanentMagnetizationMode::EffectiveSusceptibility;
        std::cout << "configureRemanence: EffectiveSusceptibility (no remanence data)" << std::endl;
    }

    // Fill any remaining NaN in remanence parameters with sensible defaults.
    // This handles groups where has_measured_remanence=yes but inc/dec are n/a
    // (e.g. M_rem=0 with unmeasured direction).
    if (outMode != RemanentMagnetizationMode::EffectiveSusceptibility) {
        for (int g = 0; g < nGroups; ++g) {
            auto& grp = model->group(g);
            if (std::isnan(grp.remanence_magnitude))
                grp.remanence_magnitude = 0.0;
            if (std::isnan(grp.remanence_inclination))
                grp.remanence_inclination = FIELD_INC;
            if (std::isnan(grp.remanence_declination))
                grp.remanence_declination = FIELD_DEC;
        }
    }

    return anyUnknownRemanence;
}

// ---------------------------------------------------------------------------
// Build LithologyModel directly from clusters — one group per inv cluster_id.
// ordering: geo_group_index → cluster_id. If empty, uses cluster_id ascending.
// ---------------------------------------------------------------------------
std::shared_ptr<LithologyModel> buildModelFromClustersDirect(
    const std::vector<ClusterProperties>& clusters,
    double bottomDepth,
    const std::vector<int>& ordering) {

    std::vector<int> kGeoToInv = ordering;
    if (kGeoToInv.empty()) {
        std::set<int> uniqueIds;
        for (const auto& c : clusters) uniqueIds.insert(c.cluster_id);
        kGeoToInv.assign(uniqueIds.begin(), uniqueIds.end());
    }

    auto findCluster = [&](int cid) -> const ClusterProperties* {
        for (const auto& c : clusters) {
            if (c.cluster_id == cid) return &c;
        }
        return nullptr;
    };

    auto model = std::make_shared<LithologyModel>();
    model->setBottomDepth(bottomDepth);

    int nGroups = static_cast<int>(kGeoToInv.size());
    for (int gi = 0; gi < nGroups; ++gi) {
        int invCid = kGeoToInv[gi];
        const ClusterProperties* cp = findCluster(invCid);
        if (!cp) {
            std::cerr << "buildModelFromClustersDirect: inv cluster " << invCid
                      << " not found in cluster properties" << std::endl;
            continue;
        }

        std::string name = "Cluster_ID_" + std::to_string(invCid);
        LithoGroup group(gi, name, cp->density_median, cp->susceptibility_median);
        group.remanence_magnitude    = cp->remanence_magnitude;
        group.remanence_inclination  = cp->remanence_inclination;
        group.remanence_declination  = cp->remanence_declination;
        model->addGroup(group);
    }

    std::cout << "buildModelFromClustersDirect: " << model->groupCount()
              << " groups (";
    if (!ordering.empty()) std::cout << "borehole-derived:";
    else std::cout << "cluster_id ascending:";
    for (int cid : kGeoToInv) std::cout << " c" << cid;
    std::cout << ")" << std::endl;
    for (int g = 0; g < model->groupCount(); ++g) {
        const auto& grp = model->group(g);
        std::cout << "  [" << grp.id << "] " << grp.name
                  << ": rho=" << grp.density
                  << "  chi=" << grp.susceptibility
                  << "  M_rem=" << grp.remanence_magnitude
                  << " A/m  I=" << grp.remanence_inclination
                  << "  D=" << grp.remanence_declination << std::endl;
    }

    return model;
}

// ---------------------------------------------------------------------------
// Configure remanence for direct mode.
// ordering: geo_group_index → cluster_id. If empty, uses cluster_id ascending.
// ---------------------------------------------------------------------------
bool configureRemanenceFromClustersDirect(
    std::shared_ptr<LithologyModel> model,
    const std::vector<ClusterProperties>& clusters,
    RemanentMagnetizationMode& outMode,
    const std::vector<int>& ordering) {

    constexpr double FIELD_NT = 55000.0;
    constexpr double FIELD_INC = 75.0;
    constexpr double FIELD_DEC = -20.0;

    std::vector<int> kGeoToInv = ordering;
    if (kGeoToInv.empty()) {
        std::set<int> uniqueIds;
        for (const auto& c : clusters) uniqueIds.insert(c.cluster_id);
        kGeoToInv.assign(uniqueIds.begin(), uniqueIds.end());
    }

    int nGroups = model->groupCount();
    std::vector<bool> groupHasKnown(nGroups, false);
    std::vector<bool> groupHasUnknown(nGroups, false);

    auto findCluster = [&](int cid) -> const ClusterProperties* {
        for (const auto& c : clusters) {
            if (c.cluster_id == cid) return &c;
        }
        return nullptr;
    };

    for (int gi = 0; gi < nGroups && gi < static_cast<int>(kGeoToInv.size()); ++gi) {
        int invCid = kGeoToInv[gi];
        const ClusterProperties* cp = findCluster(invCid);
        if (!cp) continue;

        if (cp->has_measured_remanence) {
            groupHasKnown[gi] = true;
        } else {
            groupHasUnknown[gi] = true;
        }
    }

    bool anyUnknownRemanence = false;
    bool anyKnownRemanence = false;
    for (int g = 0; g < nGroups; ++g) {
        if (groupHasKnown[g]) anyKnownRemanence = true;
        if (groupHasUnknown[g]) anyUnknownRemanence = true;
    }

    if (anyUnknownRemanence) {
        outMode = RemanentMagnetizationMode::FixedVectorPerGroup;

        for (int g = 0; g < nGroups; ++g) {
            auto& grp = model->group(g);
            if (groupHasUnknown[g] && grp.susceptibility > 0.0) {
                double m_ind = grp.susceptibility * FIELD_NT / 100.0;
                grp.remanence_magnitude = m_ind;
                if (std::isnan(grp.remanence_inclination))
                    grp.remanence_inclination = FIELD_INC;
                if (std::isnan(grp.remanence_declination))
                    grp.remanence_declination = FIELD_DEC;
            } else if (groupHasUnknown[g] && grp.susceptibility == 0.0) {
                grp.remanence_magnitude = 0.0;
                grp.remanence_inclination = FIELD_INC;
                grp.remanence_declination = FIELD_DEC;
            }
        }

        std::cout << "configureRemanenceDirect: FixedVectorPerGroup mode "
                  << "(known=" << (anyKnownRemanence ? "mixed" : "none")
                  << " unknown=" << (anyUnknownRemanence ? "mixed" : "none") << ")"
                  << std::endl;
    } else if (anyKnownRemanence) {
        outMode = RemanentMagnetizationMode::FixedVectorPerGroup;
        std::cout << "configureRemanenceDirect: FixedVectorPerGroup (all prescribed)" << std::endl;
    } else {
        outMode = RemanentMagnetizationMode::EffectiveSusceptibility;
        std::cout << "configureRemanenceDirect: EffectiveSusceptibility (no remanence data)" << std::endl;
    }

    if (outMode != RemanentMagnetizationMode::EffectiveSusceptibility) {
        for (int g = 0; g < nGroups; ++g) {
            auto& grp = model->group(g);
            if (std::isnan(grp.remanence_magnitude))
                grp.remanence_magnitude = 0.0;
            if (std::isnan(grp.remanence_inclination))
                grp.remanence_inclination = FIELD_INC;
            if (std::isnan(grp.remanence_declination))
                grp.remanence_declination = FIELD_DEC;
        }
    }

    return anyUnknownRemanence;
}

} // namespace litho_invert

