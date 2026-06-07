#include <litho_invert/core/geometry.h>
#include <cmath>

namespace {
    constexpr double kEpsilon = 1e-15;
}

namespace litho_invert {

double solidAngle(const Vector3d& a, const Vector3d& b,
                  const Vector3d& c, const Vector3d& o) {
    Vector3d ra = a - o;
    Vector3d rb = b - o;
    Vector3d rc = c - o;

    double ra_len = ra.norm();
    double rb_len = rb.norm();
    double rc_len = rc.norm();

    // If any vector is zero, the point is on a vertex => solid angle is 0
    if (ra_len < kEpsilon || rb_len < kEpsilon || rc_len < kEpsilon) return 0.0;

    double numerator = tripleProduct(ra, rb, rc);
    double denominator = ra_len * rb_len * rc_len
                       + ra.dot(rb) * rc_len
                       + ra.dot(rc) * rb_len
                       + rb.dot(rc) * ra_len;

    return 2.0 * std::atan2(numerator, denominator);
}

double lineIntegralTerm(const Vector3d& a, const Vector3d& b,
                        const Vector3d& o) {
    Vector3d ra = a - o;
    Vector3d rb = b - o;
    Vector3d e = b - a;
    double L = e.norm();

    if (L < kEpsilon) return 0.0;  // degenerate edge

    double s1 = ra.dot(e) / L;
    double s2 = rb.dot(e) / L;
    double d1 = ra.norm();
    double d2 = rb.norm();

    double num = s2 + d2;
    double denom = s1 + d1;
    if (std::abs(denom) < kEpsilon && std::abs(num) < kEpsilon) return 0.0;
    if (std::abs(denom) < kEpsilon) return 0.0;
    if (num <= 0.0 || denom <= 0.0) return 0.0;
    return std::log(num / denom);
}

double surfaceIntegralTerm(const Vector3d& a, const Vector3d& b,
                           const Vector3d& c, const Vector3d& o) {
    double omega = solidAngle(a, b, c, o);

    // Face normal
    Vector3d n = (b - a).cross(c - a);
    double n_norm = n.norm();
    if (n_norm < kEpsilon) return 0.0;

    // Perpendicular distance from observation point to plane
    Vector3d unit_n = n / n_norm;
    double h = (a - o).dot(unit_n);

    return omega * h;
}

} // namespace litho_invert

