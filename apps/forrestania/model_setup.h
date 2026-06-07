#pragma once
#include <litho_invert/litho/lithology_model.h>
#include <litho_invert/surface/surface_mesh.h>
#include <memory>
#include <algorithm>

namespace litho_invert {

inline void finalizeModelSetup(std::shared_ptr<LithologyModel> model,
                                VertexFreedom freedom = VertexFreedom::XYZ_FREE,
                                double marginMultiplier = 1.0) {
    const int nGroups = model->groupMeshCount();
    if (nGroups < 1) return;

    // 1. Per-mesh wide depth bounds + vertex freedom reset
    for (int g = 0; g < nGroups; ++g) {
        auto* mesh = model->groupMesh(g);
        if (!mesh) continue;
        mesh->setBounds(-10000.0, 100.0);
        mesh->setAllFreedom(freedom);
    }

    // 2. Exterior face classification
    model->fixExteriorFaces();

    // 3. Model bounds from actual AABB
    model->applyModelBounds(marginMultiplier);

    // 4. Bottom depth from deepest vertex
    double minZ = 1e30;
    for (int g = 0; g < nGroups; ++g) {
        const auto* m = model->groupMesh(g);
        if (!m) continue;
        for (uint32_t vi = 0; vi < m->vertexCount(); ++vi) {
            minZ = std::min(minZ, m->vertex(vi).position.z());
        }
    }
    if (minZ < 1e29) {
        model->setBottomDepth(minZ - 100.0);
    }
}

} // namespace litho_invert

