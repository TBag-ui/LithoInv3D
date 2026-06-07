#pragma once
#include <litho_invert/core/common.h>

namespace litho_invert {

// Triple product: a · (b × c)
inline double tripleProduct(const Vector3d& a, const Vector3d& b, const Vector3d& c) {
    return a.dot(b.cross(c));
}

// Signed volume of tetrahedron (a,b,c,d)
inline double tetraVolume(const Vector3d& a, const Vector3d& b,
                           const Vector3d& c, const Vector3d& d) {
    return tripleProduct(b - a, c - a, d - a) / 6.0;
}

// Solid angle subtended by triangle (a,b,c) at observation point o
// Uses the formula from van Oosterom & Strackee (1983), also used in Nagy (2000)
double solidAngle(const Vector3d& a, const Vector3d& b,
                  const Vector3d& c, const Vector3d& o);

// Line integral term L_i from Nagy (2000), Eq. 6
// Integral along edge from a to b at observation point o
double lineIntegralTerm(const Vector3d& a, const Vector3d& b,
                        const Vector3d& o);

// Surface integral term S_f from Nagy (2000), Eq. 7
// Integral over triangle (a,b,c) at observation point o
double surfaceIntegralTerm(const Vector3d& a, const Vector3d& b,
                           const Vector3d& c, const Vector3d& o);

} // namespace litho_invert

