#pragma once

#include <algorithm>
#include <cmath>

namespace fastnsa {

// Geometry libraries can differ by a few ULPs when summing segment lengths.
// Canonicalize endpoint round-off, while rejecting genuine mapping errors.
inline bool canonical_edge_offset(double value, double length, double& canonical) noexcept {
    if (!std::isfinite(value) || !std::isfinite(length) || length < 0.0)
        return false;
    const double tolerance = std::max(1e-9, 1e-12 * length);
    if (value < -tolerance || value > length + tolerance)
        return false;
    canonical = std::clamp(value, 0.0, length);
    return true;
}

}  // namespace fastnsa
