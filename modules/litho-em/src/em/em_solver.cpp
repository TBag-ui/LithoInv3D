#include <litho_invert/em/em_solver.h>
#include <litho_invert/em/em_active_forward.h>
#include <stdexcept>

namespace litho_invert {

std::shared_ptr<EMSolver> createEMSolver(const std::string& method) {
    // =========================================================================
    // EXTENSION POINT: Add new solver types here.
    //
    // When implementing a new EMSolver subclass:
    //   1. Add an #include for the header above
    //   2. Add an else-if branch below
    //   3. The config string (e.g. "fdfd") is set by the user in EMConfig
    //
    // Currently only "ie" (integral equation on polyhedra) is stubbed.
    // The actual IE solver will be implemented in a separate file
    // (src/em/ie_solver.cpp) that depends on the geometry routines
    // already available (solidAngle, lineIntegralTerm in geometry.h).
    // =========================================================================

    if (method == "ie") {
        // Placeholder: return nullptr for now — the forward models will
        // fall back to a built-in approximate computation until the full
        // IE solver is implemented.
        return nullptr;
    }

    throw std::runtime_error("Unknown EMSolver method: " + method);
}

} // namespace litho_invert

