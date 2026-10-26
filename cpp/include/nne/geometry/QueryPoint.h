#pragma once

#include "nne/geometry/Polyline.h"

namespace fastnsa::nne {

struct QueryPoint : Point2D {
    QueryPoint(double x_value, double y_value) : Point2D{x_value, y_value} {}
};

}  // namespace fastnsa::nne
